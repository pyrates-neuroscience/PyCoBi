import os
import pickle
import warnings
from dataclasses import dataclass, field

import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame, MultiIndex, Series
from mpl_toolkits.mplot3d import Axes3D
from pyrates import CircuitTemplate, clear
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from typing import Union, Any, Optional, List

from .utility import get_solution_keys, get_branch_info, get_solution_variables, \
    get_solution_params, get_lyapunov_exponents, parse_point_diagnostics


@dataclass
class Continuation:
    """One continuation tracked by `ODESystem` — bundles together what was
    previously scattered across the `auto_solutions`, `results`, `_results_map`,
    and `_branches` dicts.

    Fields
    ------
    key
        The pycobi key (an int) under which this continuation is stored.
    name
        User-supplied name passed to `run(..., name=...)`, or `None`.
    branch_id
        auto-07p's BR identifier for the branch this continuation lives on.
        Multiple continuations may share a branch_id (e.g. when one extends
        another or when bidirectional merges into the same branch).
    icps
        Continuation parameter index tuples used to grow this branch. A
        bidirectional run that goes both ways in PAR(4) records `[(4,)]`
        once; a branch that was extended a second time in PAR(5) appends
        `(5,)`.
    auto_solution
        The auto-07p ``bifDiag`` object holding the actual solution data.
        Not pickle-safe — excluded from `to_file`; on `from_file` this
        will be `None` for loaded instances.
    summary
        PyCoBi's parsed `DataFrame` summary of the continuation. Populated
        by `_create_summary` after `auto.run` returns.
    """
    key: int
    name: Optional[str]
    branch_id: int
    icps: List[tuple] = field(default_factory=list)
    auto_solution: Any = None
    summary: Optional[DataFrame] = None


# auto-07p reserves PAR(11)..PAR(14) for internal use (period, time, etc.), so
# user parameters must skip that slot range. Read the canonical value from
# PyRates' FortranBackend (the layer that actually allocates PAR slots), with a
# hard-coded fallback for older PyRates that doesn't expose it. Keeping a
# single source of truth so PyCoBi's _var_map and PyRates' parnames stay in
# lockstep.
try:
    from pyrates.backend.fortran.fortran_backend import FortranBackend as _FortranBackend
    _AUTO_BLOCKED_PAR_RANGE = _FortranBackend._AUTO_BLOCKED_PAR_RANGE
except (ImportError, AttributeError):
    _AUTO_BLOCKED_PAR_RANGE = (10, 15)


class ODESystem:

    __slots__ = ["auto_solutions", "results", "_orig_dir", "dir", "_auto", "_last_cont", "_cont_num", "_results_map",
                 "_branches", "_bifurcation_styles", "_temp", "additional_attributes", "_eq", "_var_map",
                 "_var_map_inv", "continuations"]

    blocked_indices = _AUTO_BLOCKED_PAR_RANGE

    def __init__(self, eq_file: str, working_dir: str = None, auto_dir: str = None, init_cont: bool = False,
                 params: list = None, state_vars: list = None, **kwargs) -> None:
        """

        Parameters
        ----------
        eq_file
            Equation file that this instance of PyCoBi will use for all calls to `PyCoBi.run`
        working_dir
            Directory in which all the fortran equation and auto-07p constant files are saved.
        auto_dir
            Installation directory of auto-07p.
        init_cont
            If true, an initial-value integration with respect to time is performed at instantiation, using the
            equation file provided via the keyword argument `e=<fname>` (a file named `<fname>.f90` should exist in
            `working_dir`) and the auto constants provided via the keyword argument `c=<fname>` (a file named
            `c.<fname>` should exist in `working_dir`). Defaults to false — many use cases start from a pre-converged
            steady state and don't need an IVP, and time integration can be slow or fail to converge for stiff
            systems. Set to true to opt into the legacy behaviour.
        params
            Optional ordered list with names of all parameters in the model equations. Can be used to refer to model
            parameters.
        state_vars
            Optional ordered list that provides a name for each entry in the state vector of the model equations.
        kwargs
            Additional keyword arguments. When `init_cont=True`, these are forwarded to `ODESystem.run` for the
            initial time-integration call. Ignored otherwise.
        """
        
        # make sure that auto-07p environment variables are set
        if 'AUTO_DIR' not in os.environ:
            if auto_dir is None:
                raise ValueError('Auto-07p directory has not been set as environment variable. '
                                 'Please provide path to cmds/auto.env.sh or set environment variable yourself.')
            else:
                auto_dir = auto_dir.replace('$HOME', '~')
                auto_dir = os.path.expanduser(auto_dir)
                os.environ['AUTO_DIR'] = auto_dir
                path = f"{auto_dir}/cmds:{auto_dir}/bin:{os.environ['PATH']}"
                os.environ['PATH'] = path

        import auto as a

        # open attributes
        self.auto_solutions = {}
        self.results = {}
        # `continuations` is the canonical store; `auto_solutions`, `results`,
        # `_results_map`, and `_branches` are kept as mirrors of its fields for
        # backward compatibility with external code that reads them directly.
        # Mirrors are populated by `_register_continuation` / `_record_summary`;
        # direct mirror writes by external code don't propagate back here.
        self.continuations = {}
        self._orig_dir = os.getcwd()
        if working_dir:
            try:
                os.chdir(working_dir)
            except FileNotFoundError:
                os.chdir(f"{os.getcwd()}/{working_dir}")
        self.dir = os.getcwd()
        self.additional_attributes = {}

        # private attributes
        self._auto = a
        self._eq = eq_file
        self._last_cont = 0
        self._cont_num = 0
        self._results_map = {}
        self._branches = {}
        self._bifurcation_styles = {'LP': {'marker': 'v', 'color' : '#5D6D7E'},
                                    'HB': {'marker': 'o', 'color': '#148F77'},
                                    'CP': {'marker': 'd', 'color': '#5D6D7E'},
                                    'PD': {'marker': 'h', 'color': '#5D6D7E'},
                                    'BT': {'marker': 's', 'color': 'k'},
                                    'GH': {'marker': 'o', 'color': '#148F77'}
                                    }
        self._temp = kwargs.pop("template", None)

        # Build name <-> PAR/U index maps. Two paths feed into this code:
        #
        #   (1) Hand-written .f90 + c.* without parnames/unames. The user
        #       passes `params=[...]`, `state_vars=[...]` so PyCoBi can
        #       translate user-facing names ("eta", "r") to auto-07p's
        #       internal "PAR(i)" / "U(i)" form for both inputs (ICP, UZR,
        #       UZSTOP keys) and outputs (DataFrame column relabelling).
        #
        #   (2) PyRates-generated .f90 + c.* with parnames/unames. The user
        #       passes namespaced names ("p/qif_op/eta") via `from_template`,
        #       but auto-07p's solution exposes the bare local names ("eta")
        #       through `solution.coordnames`. Since the namespaced and bare
        #       names don't collide, every lookup misses and the input/output
        #       strings pass through unchanged — auto-07p resolves them via
        #       parnames/unames internally. `_var_map` is dead code on this
        #       path but doesn't get in the way.
        #
        # Per-entry storage is a (kind, idx) tuple — kind 'P' for parameters
        # (mapped to "PAR(idx)") or 'U' for state variables (mapped to
        # "U(idx)"). The "plot" string form is derived on demand by `_map_var`
        # rather than stored alongside the int. PAR(14) is auto-07p's reserved
        # time slot, always present.
        self._var_map = {"t": ("P", 14)}
        self._var_map_inv = {}
        if params:
            increment = 1
            for i, key in enumerate(params):
                idx = i + increment
                if self.blocked_indices[0] <= idx <= self.blocked_indices[1]:
                    idx -= increment
                    increment += self.blocked_indices[1] - self.blocked_indices[0]
                    idx += increment
                self._var_map[key] = ("P", idx)
        if state_vars:
            for i, key in enumerate(state_vars):
                self._var_map[key] = ("U", i + 1)
        # Only the "plot" string -> user-facing-name direction is ever read
        # (by `_create_summary`'s column remapping and `extract`). Derive each
        # plot string from the (kind, idx) tuple — `_map_var(name, "plot")`
        # produces the same form.
        for name in self._var_map:
            self._var_map_inv[self._map_var(name, "plot")] = name

        # perform initial continuation in time to ensure convergence to steady-state solution
        if init_cont:
            _ = self.run(ICP=[14], **kwargs)

    def __getitem__(self, item):
        # Direct int-key lookup, then a named-continuation lookup via
        # _results_map. On a double-miss, raise a single KeyError listing every
        # name and int key currently registered — much friendlier than the
        # opaque `KeyError: <item>` you'd get from the inner dict.
        try:
            return self.results[item]
        except KeyError:
            pass
        try:
            return self.results[self._results_map[item]]
        except KeyError:
            known_names = sorted(self._results_map.keys())
            known_keys = sorted(self.results.keys())
            raise KeyError(
                f"{item!r} is neither a registered continuation name nor a stored "
                f"pyauto key. Known names: {known_names}. Known keys: {known_keys}."
            )

    @property
    def pyrates_template(self):
        return self._temp

    def close_session(self, clear_files: bool = False, **kwargs):
        if clear_files:
            clear(self._temp, **kwargs)
        os.chdir(self._orig_dir)

    @staticmethod
    def reset_auto_state() -> None:
        """Clear auto-07p's per-process cross-run state (``parnames``, ``unames``).

        auto-07p's Python wrapper deliberately leaves the ``parnames`` /
        ``unames`` entries on its global runner intact between successive
        ``run()`` calls — see ``auto/runAUTO.py``::

            # do not completely replace existing constants data but
            # leave the special keys such as unames, parnames, etc, intact

        That's helpful when iterating on a single model, but it leaks across
        unrelated model loads. Concretely: an `ODESystem` for model A whose
        generated c.* declares ``unames = {1: 'r', 2: 'v'}`` populates the
        global runner; instantiating model B whose c.* declares no unames will
        inherit ``{1: 'r', 2: 'v'}`` and silently relabel B's DataFrame
        columns with A's state-variable names. The same applies to
        ``parnames``.

        Call this between unrelated model loads (typically in test teardown,
        or right before constructing a fresh `ODESystem` from a different
        model). No-ops if `auto` was never imported, or if its internal layout
        differs from what we expect.
        """
        runner = ODESystem._get_auto_runner()
        if runner is None:
            return
        constants = runner.options.get('constants')
        if constants is None:
            return
        for key in ('parnames', 'unames'):
            constants[key] = None

    @staticmethod
    def _get_auto_runner():
        """Locate auto-07p's global ``runAUTO`` instance via the ``withrunner``
        closure that ``AUTOSimpleFunctions`` binds into every command.

        Auto's package-level ``run`` / ``load`` / etc. are FunctionType copies
        whose globals carry a ``withrunner`` closure pointing at the
        AUTOSimpleFunctions singleton; the singleton's ``_runner`` is the
        process-wide runAUTO instance that holds the persisted ``constants``
        dict. Returns ``None`` if auto isn't imported or the layout changes
        upstream — callers must tolerate that.
        """
        try:
            import auto as a
        except ImportError:
            return None
        run_fn = getattr(a, 'run', None)
        if run_fn is None:
            return None
        withrunner = run_fn.__globals__.get('withrunner')
        if withrunner is None or withrunner.__closure__ is None:
            return None
        for cell in withrunner.__closure__:
            simple_funcs = cell.cell_contents
            if hasattr(simple_funcs, '_runner'):
                return simple_funcs._runner
        return None

    @classmethod
    def from_yaml(cls, path: str, working_dir: str = None, auto_dir: str = None, init_cont: bool = False,
                  init_kwargs: dict = None, analytical_jacobian: bool = True,
                  auto_constants: Union[str, tuple, list] = ('ivp',), **kwargs):
        """Instantiates `ODESystem` from a YAML definition file.

        Parameters
        ----------
        path
            Full path to a YAML model definition file for a `pyrates.CircuitTemplate`.
        working_dir
            Directory in which all the fortran equation and auto-07p constant files are saved.
        auto_dir
            Installation directory of auto-07p.
        init_cont
            If true, an IVP time integration is run at instantiation against the freshly generated c.ivp file.
            Defaults to false — set to true to opt into the legacy behaviour. See `__init__` for details.
        init_kwargs
            Additional keyword arguments that will be provided to the `ODESystem.run` method for performing the time
            integration.
        analytical_jacobian
            If true (default), instruct PyRates to symbolically differentiate the vector field and emit DFDU/DFDP
            inside the generated `func` subroutine; the generated `c.*` file will set `JAC=1` so auto-07p uses the
            analytical Jacobian. Set to false to fall back to auto-07p's finite-difference Jacobian (useful when
            symbolic differentiation is slow or produces unwieldy expressions for the model at hand).
        auto_constants
            Name (or iterable of names) of auto-07p continuation scenarios to generate `c.<name>` files for. See
            `from_template` for the recognised scenarios and their default constants. Defaults to `('ivp',)` for
            backward compatibility.
        kwargs
            Additional keyword arguments provided to the `pyrates.CircuitTemplate.get_run_func` method that is used to
            generate the fortran equation file and the auto constants file that will be used to initialize `ODESystem`.

        Returns
        -------
        ODESystem
            `ODESystem` instance.
        """

        return cls.from_template(CircuitTemplate.from_yaml(path), working_dir=working_dir, auto_dir=auto_dir,
                                 init_cont=init_cont, init_kwargs=init_kwargs,
                                 analytical_jacobian=analytical_jacobian,
                                 auto_constants=auto_constants, **kwargs)

    @classmethod
    def from_template(cls, template: CircuitTemplate, working_dir: str = None, auto_dir: str = None,
                      init_cont: bool = False, init_kwargs: dict = None, analytical_jacobian: bool = True,
                      auto_constants: Union[str, tuple, list] = ('ivp',), **kwargs):
        """Instantiates `ODESystem` from a `pyrates.CircuitTemplate`.

        Parameters
        ----------
        template
            Instance of the class `pyrates.CircuitTemplate`.
        working_dir
            Directory in which all the fortran equation and auto-07p constant files are saved.
        auto_dir
            Installation directory of auto-07p.
        init_cont
            If true, an IVP time integration is run at instantiation against the freshly generated c.ivp file.
            Defaults to false — set to true to opt into the legacy behaviour. See `__init__` for details.
        init_kwargs
            Additional keyword arguments that will be provided to the `ODESystem.run` method for performing the time
            integration.
        analytical_jacobian
            If true (default), instruct PyRates to symbolically differentiate the vector field and emit DFDU/DFDP
            inside the generated `func` subroutine; the generated `c.*` file will set `JAC=1` so auto-07p uses the
            analytical Jacobian. Set to false to fall back to auto-07p's finite-difference Jacobian (useful when
            symbolic differentiation is slow or produces unwieldy expressions for the model at hand). Can be
            overridden on a per-continuation basis by passing `JAC=0` or `JAC=1` to `ODESystem.run`.
        auto_constants
            Name (or iterable of names) of auto-07p continuation scenarios PyRates should emit `c.<name>` files for.
            One file is written per requested scenario, each pre-configured with auto-07p constants appropriate for
            that mode. Recognised scenarios:

              * ``'ivp'`` — initial-value problem / time integration (``IPS=-2``). Required when ``init_cont=True``.
              * ``'eq'``  — equilibrium continuation in one parameter (``IPS=1``).
              * ``'lc'``  — limit-cycle continuation in one parameter, with PAR(11) as the period (``IPS=2``).
              * ``'bvp'`` — boundary-value problem (``IPS=4``).

            Pass e.g. ``auto_constants=('ivp', 'eq', 'lc')`` to set up all three at once and then switch scenarios on
            a per-call basis via ``ode.run(c='eq', ...)``. Auto-07p constants passed as kwargs (``NMX``, ``DSMAX``,
            ``UZSTOP``, ...) apply to every requested scenario; per-scenario overrides can be applied at run-time on
            the corresponding ``ODESystem.run`` call. Defaults to ``('ivp',)`` for backward compatibility.
        kwargs
            Additional keyword arguments provided to the `pyrates.CircuitTemplate.get_run_func` method that is used to
            generate the fortran equation file and the auto constants file that will be used to initialize `ODESystem`.

        Returns
        -------
        ODESystem
            `ODESystem` instance.
        """

        # normalise & validate auto_constants (before any I/O so an obviously
        # inconsistent combo errors out cleanly rather than as a stale
        # working-dir FileNotFoundError from chdir downstream).
        scenarios = (auto_constants,) if isinstance(auto_constants, str) else tuple(auto_constants)
        if init_cont and 'ivp' not in scenarios:
            raise ValueError(
                f"init_cont=True performs an IVP integration against c.ivp, but 'ivp' is missing from "
                f"auto_constants={scenarios!r}. Either include 'ivp' (e.g. auto_constants=('ivp', 'eq')) or "
                f"set init_cont=False."
            )

        # change working directory
        if working_dir:
            try:
                os.chdir(working_dir)
            except FileNotFoundError:
                os.chdir(f"{os.getcwd()}/{working_dir}")

        # preparations
        func_name = kwargs.pop("func_name", "vector_field")
        file_name = kwargs.pop("file_name", "system_equations")
        dt = kwargs.pop("step_size", 1e-3)
        solver = kwargs.pop("solver", "scipy")
        if init_kwargs is None:
            init_kwargs = {}

        # update circuit template variables
        if "node_vars" in kwargs:
            template.update_var(node_vars=kwargs.pop("node_vars"))
        if "edge_vars" in kwargs:
            template.update_var(edge_vars=kwargs.pop("edge_vars"))

        # generate fortran files
        prec = kwargs.pop("float_precision", "float64")
        _, _, params, state_vars = template.get_run_func(func_name, dt, file_name=file_name, backend="fortran",
                                                         float_precision=prec, auto=True, auto_jac=analytical_jacobian,
                                                         auto_constants=scenarios,
                                                         vectorize=False, solver=solver, **kwargs)

        # PyRates returns the full positional argument list for the run function,
        # which prepends some non-parameter args (state vector ``y``, derivative
        # ``dy``, time ``t``, optional history function ``hist`` for DDEs) before
        # the actual model parameters. Filter by name rather than slicing by
        # position — survives DDE models (where ``hist`` shifts the offset) and
        # is robust to upstream signature reordering.
        non_param_args = {'t', 'y', 'dy', 'hist'}
        param_names = tuple(p for p in params if p not in non_param_args)

        # initialize ODESystem
        return cls(auto_dir=auto_dir, init_cont=init_cont, c="ivp", eq_file=file_name, template=template,
                   params=param_names, state_vars=list(state_vars), **init_kwargs)

    @classmethod
    def from_file(cls, filename: str, auto_dir: str = None):
        """Load `ODESystem` from a pickle file written by `to_file`.

        Parameters
        ----------
        filename
            Path to the pickle file written by `ODESystem.to_file`.
        auto_dir
            Installation directory of `auto-07p`.

        Returns
        -------
        ODESystem
            ODESystem instance with state restored from the file. Slots not stored on disk
            (the live `auto` module, the PyRates `CircuitTemplate`, the cwd snapshots) are
            re-initialised by `__init__` rather than from the file.
        """
        pyauto_instance = cls('', auto_dir=auto_dir, init_cont=False)
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        for key, val in data.items():
            attr = getattr(pyauto_instance, key, None)
            # Merge dicts in place to preserve any fresh init values; otherwise
            # just bind the loaded value. This makes from_file work for the
            # whole-state pickle as well as the results-only one.
            if isinstance(attr, dict) and isinstance(val, dict):
                attr.update(val)
            else:
                setattr(pyauto_instance, key, val)
        # `continuations` is in _PICKLE_EXCLUDE — rebuild it from the mirror
        # dicts we just restored so `get_continuation(...)` works after load.
        pyauto_instance._rebuild_continuations_from_mirrors()
        return pyauto_instance

    # Slots that to_file deliberately omits. ``dir`` / ``_orig_dir`` are
    # session-local cwd snapshots that should be re-derived on load.
    # ``_auto`` is a live Python module (the auto-07p package) — pickle can't
    # serialise modules, and __init__ re-imports it anyway. ``_temp`` is a
    # PyRates CircuitTemplate that holds lambdified sympy functions which
    # don't pickle; users that want the template on disk should pickle it
    # separately or save the YAML path alongside. ``auto_solutions`` contains
    # auto-07p bifDiag objects that hold open BufferedReader handles on the
    # fort.* files and refuse to pickle. ``continuations`` mirrors hold the
    # same bifDiag objects on their `.auto_solution` field — also unpicklable;
    # rebuilt by `_rebuild_continuations_from_mirrors` on load. ``_last_cont``
    # is normally an int but `merge()` rebinds it to a solution object — same
    # problem; skip it and rely on `_cont_num` to track the count.
    _PICKLE_EXCLUDE = frozenset({
        "dir", "_orig_dir", "_auto", "_temp",
        "auto_solutions", "_last_cont", "continuations",
    })

    def to_file(self, filename: str, results_only: bool = True, **kwargs) -> None:
        """Save the instance state on disc via pickle.

        Parameters
        ----------
        filename
            Path to write the pickle file to. If a file already exists at this path it
            will be overwritten.
        results_only
            When true (default), only the PyCoBi-side bookkeeping (`results`, `_branches`,
            `_results_map`) is saved — enough to reproduce DataFrames and plots without
            rerunning auto-07p. When false, all pickle-safe slots are saved (see
            `_PICKLE_EXCLUDE` for the slots intentionally omitted because they can't
            round-trip).
        kwargs
            Extra metadata to attach to the dump. Restored as `additional_attributes`.

        Returns
        -------
        None
        """

        if results_only:
            data = {'results': self.results, '_branches': self._branches, '_results_map': self._results_map}
        else:
            data = {key: getattr(self, key)
                    for key in self.__slots__ if key not in self._PICKLE_EXCLUDE}
        data.update({'additional_attributes': kwargs})

        try:
            with open(filename, 'xb') as f:
                pickle.dump(data, f)
        except FileExistsError:
            with open(filename, 'wb') as f:
                pickle.dump(data, f)

    def run(self, origin: Union[int, str, object] = None, starting_point: Union[str, int] = None, variables: list = None,
            params: list = None, get_stability: bool = True, get_period: bool = False, get_timeseries: bool = False,
            get_eigenvals: bool = False, get_lyapunov_exp: bool = False, reduce_limit_cycle: bool = True,
            bidirectional: bool = False, name: str = None, _reverse_direction: bool = False, **auto_kwargs) -> tuple:
        """
        Wraps auto-07p command `run` and stores requested solution details on instance.

        Parameters
        ----------
        origin
            Key of the solution branch that contains the solution `starting_point`, from which the new continuation will
            be started.
        starting_point
            Solution on the origin branch to start the new continuation from. Accepted forms:

              * Auto-07p label string — ``'EP'``, ``'LP1'``, ``'HB2'`` etc. The first two characters
                are the bifurcation type; the optional trailing integer disambiguates when the branch
                carries several solutions of the same type (1-based, defaults to 1). The IVP that
                `init_cont=True` runs produces ``'EP1'`` for the initial state and ``'EP2'`` for the
                converged steady state — use ``'EP2'`` when starting an equilibrium continuation from
                the IVP's terminal state.
              * Bare bifurcation type (no number) — ``'EP'`` is equivalent to ``'EP1'``.
              * Integer point index — a 1-based auto-07p point number on the branch (the ``PT``
                column of the printed table). Useful when auto produced unlabeled regular points
                that you want to continue from.
              * ``None`` — only valid on the very first call against a fresh `ODESystem` (when
                no prior continuation exists to extend); subsequent calls require an explicit
                starting point so PyCoBi knows which branch to extend.
        variables
            Keys of the state variables that should be recorded for each continuation recording step.
        params
            Keys of the parameters that should be recorded for each continuation recording step.
        get_stability
            If true, the stability of each solution will be stored in the results under the key 'stability'.
        get_period
            If true, the period of periodic solutions will be stored in the results under the key 'period'.
        get_timeseries
            If true, the time vector associated with the state variables of a periodic solution will be stored under the
            key 'time'.
        get_eigenvals
            If true, the eigenvalues (floquet multipliers) or steady-state (periodic) solutions will be stored under the
            key 'eigenvalues'.
        get_lyapunov_exp
            If true, the local lyapunov exponents of solutions will be stored under the key 'lyapunov'.
        reduce_limit_cycle
            If true, the values of each state variable will be reduced to the minimum and maximum for limit cycle 
            solutions. Else, the state variable values will be stored for multiple discretized points along the limit
            cycle solution (number depends on the arguments passed to Auto).
        bidirectional
            If true, parameter continuation will be performed into both directions for a given continuation parameter.
        name
            Name, under which the resulting solution branch will be accessible for future continuations.
        _reverse_direction
            Private flag set internally by the recursive call that `bidirectional=True` makes. Tells `run()` that
            this invocation is the reverse-direction half of a bidirectional continuation, so the result is merged
            into the forward branch rather than registered as a fresh continuation. Don't pass manually.
        auto_kwargs
            Additional keyword arguments to be passed to the auto command `run`. All auto-07p constants can be
            overridden here (e.g. `NMX`, `DSMAX`, `ICP`, ...). In particular, `JAC=0`/`JAC=1` overrides the
            Jacobian source for this single continuation: pass `JAC=0` to force finite-difference Jacobian even if
            `from_template` / `from_yaml` was instantiated with `analytical_jacobian=True`, and vice versa.

        Returns
        -------
        tuple
            DataFrame with the results, auto solution branch object.
        """

        # auto call
        ###########

        # extract starting point of continuation
        if self._last_cont == 0 and self._last_cont not in self.auto_solutions:
            auto_kwargs["e"] = self._eq
        if 'IRS' in auto_kwargs or 's' in auto_kwargs:
            raise ValueError('Usage of keyword arguments `IRS` and `s` is disabled in pycobi. To start from a previous'
                             'solution, use the `starting_point` keyword argument and provide a tuple of branch '
                             'number and point number as returned by the `run` method.')
        if not starting_point and self._last_cont > 0:
            raise ValueError('A starting point is required for further continuation. Either provide a solution to start'
                             ' from via the `starting_point` keyword argument or create a fresh `ODESystem` instance.')
        if origin is None:
            origin = self._last_cont
        elif type(origin) is str:
            origin = self._results_map[origin]
        elif type(origin) is not int:
            origin = origin.pycobi_key

        # call to auto
        auto_kwargs = self._map_auto_kwargs(auto_kwargs)
        solution = self._call_auto(starting_point, origin, **auto_kwargs)

        # extract information from auto solution
        ########################################

        # extract branch and solution info
        new_branch, new_icp = get_branch_info(solution)
        new_points = get_solution_keys(solution)

        # get all passed variables and params
        solution_tmp, *_ = self.get_solution(point=new_points[0], cont=solution)
        if variables is None:
            variables = self._get_all_var_keys(solution_tmp)
        variables = [self._map_var(v, mode="plot") for v in variables]
        if params is None:
            try:
                params = self._get_all_param_keys(solution_tmp)
            except KeyError:
                n_params = auto_kwargs['NPAR']
                params = [f"PAR({i})" for i in range(1, n_params+1)]
        params = [self._map_var(p, mode="plot") for p in params]

        # store solution and extracted information in pycobi
        ####################################################

        # Decide whether this continuation extends an existing branch (merge
        # path) or starts a fresh one. Three cases:
        #   1. Same (branch, origin, icp) tuple seen before — auto extended
        #      the existing branch; merge into the previous result.
        #   2. The reverse-direction half of a bidirectional run — merge into
        #      the forward branch identified by `_last_cont`.
        #   3. Otherwise — allocate a fresh pyauto key.
        if new_branch in self._branches and origin in self._branches[new_branch] \
                and new_icp in self._branches[new_branch][origin]:
            solution_old, *_ = self.get_solution(origin)
            pyauto_key = solution_old.pycobi_key
            solution, new_points = self.merge(pyauto_key, solution, new_icp)
        elif _reverse_direction and 'DS' in auto_kwargs and auto_kwargs['DS'] == '-':
            solution_old = self.auto_solutions[self._last_cont]
            pyauto_key = solution_old.pycobi_key
            solution, new_points = self.merge(pyauto_key, solution, new_icp)
        else:
            pyauto_key = self._cont_num + 1 if self._cont_num in self.auto_solutions else self._cont_num
            solution.pycobi_key = pyauto_key

        # The reverse-direction half of a bidirectional run doesn't register
        # a fresh name — it merges into the forward branch which is already
        # in _results_map under the user's name.
        registered_name = name if (name and not _reverse_direction) else None
        self._register_continuation(
            key=pyauto_key, name=registered_name, branch_id=new_branch,
            icp=new_icp, auto_solution=solution,
        )

        # if continuation should be bidirectional, call this method again with reversed continuation direction
        ######################################################################################################

        if bidirectional:

            # perform continuation in opposite direction; the recursive call's
            # _reverse_direction=True tells it to merge into this branch rather
            # than register itself as a fresh continuation.
            ds = auto_kwargs.pop('DS', None)
            _, solution = self.run(origin, starting_point, variables=variables, params=params,
                                   get_stability=get_stability, get_period=get_period, get_timeseries=get_timeseries,
                                   get_eigenvals=get_eigenvals, get_lyapunov_exp=get_lyapunov_exp, bidirectional=False,
                                   _reverse_direction=True, DS=1e-3 if ds == '-' else '-', **auto_kwargs)

        else:

            # store summary of continuation results
            if new_icp[0] == 14:
                get_stability = False
            summary = self._create_summary(solution=solution, points=new_points, variables=variables,
                                           params=params, timeseries=get_timeseries, stability=get_stability,
                                           period=get_period, eigenvals=get_eigenvals, lyapunov_exp=get_lyapunov_exp,
                                           reduce_limit_cycle=reduce_limit_cycle)
            self._record_summary(pyauto_key, summary)

        return self.results[pyauto_key], solution

    def merge(self, key: int, cont, icp: tuple):
        """Merges two solutions from two separate auto continuations.

        Parameters
        ----------
        key
            PyCoBi identifier under which the merged solution should be stored. Must be equal to identifier of first
            continuation.
        cont
            auto continuation object that should be merged with the continuation object under `key`.
        icp
            Continuation parameter that was used in both continuations that are to be merged.
        """

        # call merge in auto
        solution = self._auto.merge(self.auto_solutions[key] + cont)
        solution.pycobi_key = key

        # mirror updates (idempotent — `run()` re-syncs through
        # `_register_continuation`, but keeping them here lets external
        # callers use `merge` without the surrounding bookkeeping)
        self.auto_solutions[key] = solution
        self._last_cont = key

        # also reflect the merged solution on the canonical Continuation
        # if one is registered for this key
        if key in self.continuations:
            self.continuations[key].auto_solution = solution

        # extract solution points
        points = list(solution.data[0].labels.by_index.keys())

        return solution, points

    # ------------------------------------------------------------------
    # Centralised continuation bookkeeping (replaces the scattered writes
    # to auto_solutions / results / _results_map / _branches that used to
    # live inline in `run`).
    # ------------------------------------------------------------------

    def _register_continuation(self, key: int, name: Optional[str], branch_id: int,
                               icp: tuple, auto_solution: Any) -> "Continuation":
        """Add a new `Continuation` or update an existing one, syncing all
        four legacy mirror dicts in the process.

        Idempotent on `key`: when the entry already exists (typical merge /
        bidirectional-reverse path) the auto_solution is replaced, the icp
        appended to both the dataclass and the `_branches` mirror, and any
        non-None `name` is set if not already present.
        """
        existing = self.continuations.get(key)
        if existing is None:
            cont = Continuation(
                key=key, name=name, branch_id=branch_id,
                icps=[icp], auto_solution=auto_solution, summary=None,
            )
            self.continuations[key] = cont
        else:
            cont = existing
            cont.auto_solution = auto_solution
            if icp not in cont.icps:
                cont.icps.append(icp)
            if name and not cont.name:
                cont.name = name

        # ---- mirror sync ----
        self.auto_solutions[key] = auto_solution
        self._last_cont = key
        if name:
            self._results_map[name] = key
        # `_branches` carries an icp list per (branch_id, key) — kept
        # append-only (with duplicates allowed) so the merge-detection
        # condition in `run()` keeps matching exactly as before.
        if branch_id not in self._branches:
            self._branches[branch_id] = {key: []}
        elif key not in self._branches[branch_id]:
            self._branches[branch_id][key] = []
        self._branches[branch_id][key].append(icp)
        self._cont_num = len(self.auto_solutions)
        return cont

    def _record_summary(self, key: int, summary: DataFrame) -> None:
        """Attach a parsed summary DataFrame to a Continuation and its
        `results` mirror."""
        if key in self.continuations:
            self.continuations[key].summary = summary
        self.results[key] = summary

    def _rebuild_continuations_from_mirrors(self) -> None:
        """Reconstruct `self.continuations` from the legacy mirror dicts.
        Called by `from_file` after the mirrors have been restored from
        disk — `continuations` is in `_PICKLE_EXCLUDE` (because each entry
        carries an unpicklable auto_solution), so it has to be rebuilt.
        Loaded continuations have `auto_solution=None`; callers wanting
        to drive auto from a loaded instance need to re-run the model.
        """
        self.continuations.clear()
        key_to_name = {key: name for name, key in self._results_map.items()}
        key_to_branch: dict = {}
        key_to_icps: dict = {}
        for branch_id, by_key in self._branches.items():
            for key, icps in by_key.items():
                key_to_branch[key] = branch_id
                key_to_icps.setdefault(key, []).extend(icps)

        all_keys = set(self.results) | set(self.auto_solutions) | set(key_to_branch)
        for key in all_keys:
            # Dedupe icps — `_branches[branch_id][key]` keeps duplicates
            # (append-only by design; the merge-detection condition in `run`
            # uses `in` on it), but the dataclass surface dedupes so a
            # bidirectional run's `icps` stays `[(4,)]` rather than `[(4,), (4,)]`
            # both pre- and post-pickle.
            seen = []
            for icp in key_to_icps.get(key, []):
                if icp not in seen:
                    seen.append(icp)
            self.continuations[key] = Continuation(
                key=key,
                name=key_to_name.get(key),
                branch_id=key_to_branch.get(key, 0),
                icps=seen,
                auto_solution=self.auto_solutions.get(key),
                summary=self.results.get(key),
            )

    def get_continuation(self, key_or_name: Union[int, str]) -> "Continuation":
        """Return the `Continuation` dataclass for a stored continuation,
        looked up by user-supplied name or by pyauto-key int.

        Examples
        --------
        >>> sols, _ = ode.run(starting_point='EP2', name='eta_branch', ICP='eta', ...)
        >>> cont = ode.get_continuation('eta_branch')
        >>> cont.branch_id, cont.icps, len(cont.summary)
        (1, [(4,)], 30)
        """
        if isinstance(key_or_name, str):
            try:
                key = self._results_map[key_or_name]
            except KeyError:
                raise KeyError(
                    f"No continuation named {key_or_name!r}; "
                    f"known names: {sorted(self._results_map)}"
                )
        else:
            key = key_or_name
        try:
            return self.continuations[key]
        except KeyError:
            raise KeyError(
                f"No continuation with key {key!r}; "
                f"known keys: {sorted(self.continuations)}"
            )

    def get_summary(self, cont: Optional[Union[Any, str, int]] = None, point=None) -> DataFrame:
        """Extract summary of continuation from PyCoBi.

        Parameters
        ----------
        cont
            Key of the solution branch.
        point
            Key of the solution on the branch.

        Returns
        -------
        DataFrame
            All recorded state variables, parameters, etc. for the solution/solution branch.
        """

        # get continuation summary
        if type(cont) is int:
            summary = self.results[cont]
        elif type(cont) is str:
            summary = self.results[self._results_map[cont]]
        elif cont is None:
            summary = self.results[self._last_cont]
        else:
            summary = self.results[cont.pycobi_key]

        # return continuation or point summary
        if not point:
            return summary
        elif type(point) is str:
            n = int(point[2:]) if len(point) > 2 else 1
            i = 1
            for p in summary.index:
                if point[:2] == summary.loc[p, 'bifurcation']:
                    if i == n:
                        return summary.loc[p, :]
                    i += 1
            else:
                raise KeyError(f'Invalid point: {point} was not found on continuation {cont}.')

        return summary.loc[point, :]

    def get_solution(self, cont: Union[Any, str, int], point: Union[str, int] = None) -> Union[Any, tuple]:
        """Extract auto solution object of a given solution/solution branch.

        Parameters
        ----------
        cont
            Key of the solution branch.
        point
            Key of the solution on the branch.

        Returns
        -------
        Union[Any, tuple]
            Solution type (only if `point` is provided), auto solution object.
        """

        # extract continuation object
        if type(cont) is int:
            cont = self.auto_solutions[cont]
        elif type(cont) is str:
            cont = self.auto_solutions[self._results_map[cont]]

        if point is None:
            return cont, None, None

        # extract solution point from continuation object and its solution type
        try:

            # extract solution point via string label
            s = cont(point)
            solution_name, solution_idx = point[:2], point[2:]
            solution_idx = int(solution_idx) if len(solution_idx) > 0 else 0

        except (AttributeError, KeyError, TypeError):

            # extract solution point via integer index — iterate the branch's
            # data entries until one of them holds the requested label.
            for bd in cont.data:
                try:
                    if type(point) is int:
                        s = bd.labels.by_index[point]
                        solution_name = list(s.keys())[0]
                        idx = np.argwhere([p == point for p in bd.labels.by_label[solution_name]]).squeeze()
                        solution_idx = int(idx + 1)
                        break
                    else:
                        s = bd.labels.by_label[point]
                        solution_name, solution_idx = point[:2], point[2:]
                        break
                except (KeyError, IndexError):
                    continue
            else:
                s = None
                solution_name = 'No Label'
                solution_idx = 0

            # make sure a proper solution was extracted, else return an unlabeled solution
            if solution_name != 'No Label':
                try:
                    s = s[solution_name]['solution']
                except KeyError:
                    solution_name = 'No Label'

        return s, solution_name, solution_idx

    def extract(self, keys: list, cont: Union[Any, str, int], point: Union[str, int] = None) -> tuple:
        """Extract properties from a solution.

        Parameters
        ----------
        keys
            Keys of the properties (e.g. state variable names, parameter names, ...).
        cont
            Key of the solution branch.
        point
            Key of the solution on the branch. When omitted (default), all rows of the
            summary are returned; when provided, only the row at that label.

        Returns
        -------
        tuple
            Tuple with 2 entries: (1) the requested slice of the continuation summary,
            either a `DataFrame` (when `point` is None — multiple rows × `len(keys)` cols)
            or a `Series` (when `point` is provided — a single labeled row indexed by
            the requested keys); (2) a `dict` mapping each input key to the column name
            that was actually used to look it up in the summary (these can differ when
            PyCoBi's `_var_map_inv` translates `'PAR(i)'` / `'U(i)'` back to user-facing
            names).
        """
        summary = self.get_summary(cont)
        columns = [k for k, _ in list(summary.keys())]
        keys_new = [self._resolve_summary_key(key, columns) for key in keys]
        key_map = {key_old: key_new for key_old, key_new in zip(keys, keys_new)}
        if point:
            return summary.loc[point, keys_new], key_map
        return summary.loc[:, keys_new], key_map

    def _resolve_summary_key(self, key: str, columns: list) -> str:
        """Resolve a user-supplied key to a column name actually present in the summary.

        Resolution order (each step falls through to the next on miss):

          1. ``key`` itself is a column — common when the user already passes a
             summary-native name (e.g. ``'eta'``, ``'r'``).
          2. ``key`` lives in ``_var_map_inv`` — typical when PyCoBi was set up
             with explicit ``params=[...]`` / ``state_vars=[...]`` and the user
             passes the auto-07p-native form (``'PAR(4)'``, ``'U(1)'``).
          3. Strip the namespace prefix from step 2's result and retry — needed
             when PyRates' ``parnames`` / ``unames`` emit the bare local name
             (``'eta'``) while PyCoBi's ``_var_map_inv`` carries the namespaced
             form (``'p/qif_op/eta'``).

        Raises a ``KeyError`` that lists what was tried if nothing matches.
        """
        if key in columns:
            return key
        mapped = self._var_map_inv.get(key)
        if mapped is not None and mapped in columns:
            return mapped
        bare = mapped.rsplit('/', 1)[-1] if isinstance(mapped, str) else None
        if bare is not None and bare in columns:
            return bare
        tried = [key]
        if mapped is not None:
            tried.append(mapped)
        if bare is not None and bare not in tried:
            tried.append(bare)
        raise KeyError(
            f"{key!r} is not a recognised summary column. Tried: {tried}. "
            f"Available columns: {columns}."
        )

    def plot_continuation(self, x: str, y: str, cont: Union[Any, str, int], ax: plt.Axes = None,
                          force_axis_lim_update: bool = False, bifurcation_legend: bool = True,
                          get_stability: bool = True, **kwargs) -> LineCollection:
        """Line plot of 1D/2D parameter continuations and the respective codimension 1/2 bifurcations.

        Parameters
        ----------
        x
            Key of the parameter/variable plotted on the x-axis.
        y
            Key of the variable/parameter plotted on the y-axis.
        cont
            Key of the solution branch to be plotted.
        ax
            Axis in which to plot the data. If not provided, a new figure will be created.
        force_axis_lim_update
            If true, the axis limits of x and y axis will be updated after creating the line plots.
        bifurcation_legend
            If true, a legend will be plotted that lists the type of all special solutions on a continuation curve.
        get_stability
            If true, the stability of the solutions will be indicated via different line styles.
        kwargs
            Additional keyword arguments that allow to control the appearance of the line plot.

        Returns
        -------
        LineCollection
            Line object that was created.
        """

        if ax is None:
            fig, ax = plt.subplots()
        label_pad = kwargs.pop('labelpad', 5)
        tick_pad = kwargs.pop('tickpad', 5)
        axislim_pad = kwargs.pop('axislimpad', 0)

        # extract information from branch solutions
        if x in ["PAR(14)", "t"]:
            x = "t"
            results, vmap = self.extract([x, y], cont=cont)
            results['stability'] = np.asarray([True] * len(results[x]))
            results['bifurcation'] = np.asarray(['RG'] * len(results[x]))
        elif get_stability:
            results, vmap = self.extract([x, y, 'stability', 'bifurcation'], cont=cont)
        else:
            results, vmap = self.extract([x, y, 'bifurcation'], cont=cont)
            results['stability'] = np.asarray([True] * len(results[vmap[x]]))
        x, y = vmap[x], vmap[y]

        # plot bifurcation points
        bifurcation_point_kwargs = ['default_color', 'default_marker', 'default_size', 'custom_bf_styles',
                                    'ignore']
        kwargs_tmp = {key: kwargs.pop(key) for key in bifurcation_point_kwargs if key in kwargs}
        self.plot_bifurcation_points(solution_types=results['bifurcation'], x_vals=results[x],
                                     y_vals=results[y], ax=ax, **kwargs_tmp)

        # set title variable if passed
        tvar = kwargs.pop('title_var', None)
        if tvar:
            tvar_results, tmap = self.extract([tvar], cont=cont)
            tval = tvar_results[tmap[tvar]][0]
            ax.set_title(f"{tvar} = {tval}")

        # plot main continuation
        x_data, y_data = results[x], results[y]
        line_col = self._get_line_collection(x=x_data.values, y=y_data.values, stability=results['stability'], **kwargs)
        ax.add_collection(line_col)
        ax.autoscale()

        # cosmetics
        ax.tick_params(axis='both', which='major', pad=tick_pad)
        ax.set_xlabel(x, labelpad=label_pad)
        ax.set_ylabel(y, labelpad=label_pad)
        self._update_axis_lims(ax, ax_data=[x_data, y_data], padding=axislim_pad, force_update=force_axis_lim_update)
        # Skip `ax.legend()` when no labeled artists exist — otherwise
        # matplotlib emits "No artists with labels found to put in legend"
        # for every continuation that recorded zero bifurcation points.
        if bifurcation_legend and ax.get_legend_handles_labels()[1]:
            ax.legend()

        return line_col

    def plot_trajectory(self, variables: Union[list, tuple], cont: Union[Any, str, int], point: Union[str, int] = None,
                        ax: plt.Axes = None, force_axis_lim_update: bool = False, cutoff: float = None,
                        colorbar: bool = False, colorbar_label: str = None, **kwargs
                        ) -> LineCollection:
        """Plot trajectory of state variables through phase space over time.

        Parameters
        ----------
        variables
            State variables for which to create the trajectory. If 2, a 2D plot will be created, if 3, a 3D plot.
        cont
            Key of the solution branch to be plotted.
        point
            Key of the solution on the solution branch for which to plot the trajectories.
        ax
            Axis in which to plot the data. If not provided, a new figure will be created.
        force_axis_lim_update
            If true, the axis limits of x and y-axis will be updated after creating the line plots.
        cutoff
            Initial time to be disregarded for plotting.
        colorbar
            For 3D plots only: if true, attach a colorbar to the figure showing the
            scalar that ``_get_3d_line_collection`` mapped onto the LineCollection's
            color (default: the projected x-axis variable). Useful when the
            colour gradient already encodes time or a state variable. Ignored
            for 2D plots.
        colorbar_label
            Label for the colorbar. Defaults to the array key
            (``'x'`` / ``'y'`` / ``'z'`` — whichever the LineCollection's
            ``array=`` kwarg points at).
        kwargs
            Additional keyword arguments that allow to control the appearance of the line plot.

        Returns
        -------
        LineCollection
            Line object that was created.
        """

        # extract information from branch solutions
        try:
            results, vmap = self.extract(list(variables) + ['stability'], cont=cont, point=point)
        except KeyError:
            results, vmap = self.extract(list(variables), cont=cont, point=point)
            results['stability'] = None
        variables = [vmap[v] for v in variables]

        # apply cutoff, if passed
        if cutoff:
            try:
                time, _ = self.extract(['t'], cont=cont, point=point)
                time = time['t']
            except KeyError:
                try:
                    time, _ = self.extract(['time'], cont=cont, point=point)
                    time = time['time']
                except KeyError:
                    raise ValueError("Could not find time variable on solution to apply cutoff to. Please consider "
                                     "adding the keyword argument `get_timeseries` to the `PyCoBi.run()` call for which"
                                     "the phase space trajectory should be plotted.")
            idx = np.where(time > cutoff)
            for key, val in results.items():
                if hasattr(val, 'shape') and val.shape:
                    results[key] = val[idx]

        if len(variables) == 2:

            # create 2D plot
            if ax is None:
                fig, ax = plt.subplots()

            # plot phase trajectory
            line_col = self._get_line_collection(x=results[variables[0]], y=results[variables[1]],
                                                 stability=results['stability'], **kwargs)
            ax.add_collection(line_col)
            ax.autoscale()

            # cosmetics
            ax.set_xlabel(variables[0])
            ax.set_ylabel(variables[1])

        elif len(variables) == 3:

            # create 3D plot
            if ax is None:
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
            label_pad = kwargs.pop('labelpad', 30)
            tick_pad = kwargs.pop('tickpad', 20)
            axislim_pad = kwargs.pop('axislimpad', 0.1)

            # plot phase trajectory
            x, y, z = results[variables[0]], results[variables[1]], results[variables[2]]
            # `array=` controls which projected coordinate the LineCollection's
            # color map encodes; capture it before _get_3d_line_collection pops
            # it from kwargs so we can label the colorbar correctly.
            array_key = kwargs.get('array', 'x')
            line_col = self._get_3d_line_collection(x=x, y=y, z=z, stability=results['stability'], **kwargs)
            ax.add_collection3d(line_col)
            ax.autoscale()

            # cosmetics
            ax.tick_params(axis='both', which='major', pad=tick_pad)
            ax.set_xlabel(variables[0], labelpad=label_pad)
            ax.set_ylabel(variables[1], labelpad=label_pad)
            ax.set_zlabel(variables[2], labelpad=label_pad)
            self._update_axis_lims(ax, [x, y, z], padding=axislim_pad, force_update=force_axis_lim_update)

            if colorbar:
                # Resolve the array-key string to the user-facing variable
                # name for a sensible default colorbar label.
                array_to_label = {
                    'x': variables[0], 'y': variables[1], 'z': variables[2],
                }
                label = colorbar_label or array_to_label.get(array_key, array_key)
                ax.figure.colorbar(line_col, ax=ax, label=label, shrink=0.7)

        else:

            raise ValueError('Invalid number of state variables to plot. First argument can only take 2 or 3 state'
                             'variable names as input.')

        return line_col

    def plot_timeseries(self, var: str, cont: Union[Any, str, int], points: list = None, ax: plt.Axes = None,
                        linespecs: list = None, **kwargs) -> plt.Axes:
        """Plot state variable of a periodic solution over time.

        Parameters
        ----------
        var
            Key of the state variable.
        cont
            Key of the solution branch.
        points
            List with keys of the solutions for which to create time series plots.
            When ``None`` (default), every labelled point on the continuation
            is plotted as a separate trace.
        ax
            Axis in which to plot the data. If not provided, a new figure will be created.
        linespecs
            Per-trace keyword overrides; ``linespecs[i]`` is merged into
            ``kwargs`` for the i-th point.
        kwargs
            Additional keyword arguments that control the appearance of the plot.

        Returns
        -------
        plt.Axes
            Axis object that contains the plotted timeseries.
        """

        # Resolve `points`. The pre-1.0 `points=None` branch built an N-row
        # results dict but then iterated `range(len(['RG']))`, silently
        # dropping all but the first trace; worse, an unpacking-only-works-
        # for-2-elements bug broke the path entirely. Replaced with an
        # explicit enumeration of every labelled point on the continuation.
        if not points:
            cont_key = self._results_map[cont] if isinstance(cont, str) else cont
            try:
                stored = self.results[cont_key]
            except KeyError as exc:
                raise KeyError(
                    f"plot_timeseries: continuation {cont!r} not found "
                    f"in self.results"
                ) from exc
            points = list(stored.keys())
            if not points:
                raise ValueError(
                    f"plot_timeseries: continuation {cont!r} has no recorded "
                    f"points to plot. Pass `points=[...]` explicitly or rerun "
                    f"the continuation with NPR set so points are labelled."
                )

        # extract information from branch solutions
        results = []
        vmap: dict = {}
        for p in points:
            r, vmap = self.extract([var, 'time'], cont=cont, point=p)
            results.append(r)
        var_col = vmap[var]
        time_col = vmap['time']

        # create plot
        if ax is None:
            _, ax = plt.subplots()

        # plot phase trajectory
        if not linespecs:
            linespecs = [dict() for _ in range(len(points))]
        for i in range(len(points)):
            time = np.atleast_1d(results[i][time_col]).squeeze()
            y = np.atleast_1d(results[i][var_col]).squeeze()
            kwargs_tmp = dict(kwargs)
            kwargs_tmp.update(linespecs[i])
            line_col = self._get_line_collection(x=time, y=y, **kwargs_tmp)
            ax.add_collection(line_col)
        ax.autoscale()
        ax.legend([str(p) for p in points])

        return ax

    def plot_bifurcation_points(self, solution_types: DataFrame, x_vals: DataFrame, y_vals: DataFrame, ax: plt.Axes,
                                default_color: str = 'k', default_marker: str = '*', default_size: float = 10,
                                ignore: list = None, custom_bf_styles: dict = None) -> tuple:
        """Plot markers for special solutions at coordinates in 2D space.

        Parameters
        ----------
        solution_types
            Type of each solution, entries of DataFrame should be strings.
        x_vals
            X-coordinates of each solution.
        y_vals
            Y-coordinates of each special solution
        ax
            Axis in which to plot the data. If not provided, a new figure will be created.
        default_color
            Default color to be used if bifurcation style is not known.
        default_marker
            Default marker style to be used if bifurcation style is not known.
        default_size
            Default marker size.
        ignore
            List of solution types that should not be displayed.
        custom_bf_styles
            Dictionary containing adjustments to the default bifurcation markers and colors.

        Returns
        -------
        tuple
            A 2-entry tuple of (1) a list of PathCollections that correspond to bifurcation points, and (2) a list of
            corresponding bifurcation types.
        """

        if not ignore:
            ignore = []

        # set bifurcation styles
        if custom_bf_styles:
            for key, args in custom_bf_styles.items():
                self.update_bifurcation_style(key, **args)
        bf_styles = self._bifurcation_styles.copy()

        # draw bifurcation points. Pre-1.0 used `plt.sca(ax) + plt.plot`
        # which silently couples to matplotlib's global "current axes" state
        # — calling this from a function that itself activates a different
        # axes would draw on the wrong figure. Routed through `ax.plot`
        # directly.
        points, labels = ax.get_legend_handles_labels()
        for bf, x, y in zip(solution_types.values, x_vals.values, y_vals.values):
            if bf not in "EPMXRG" and bf not in ignore:
                if bf in bf_styles:
                    m = bf_styles[bf]['marker']
                    c = bf_styles[bf]['color']
                else:
                    m = default_marker
                    c = default_color
                # Limit-cycle case: y is a (min, max) pair — draw a marker
                # on each envelope. Equilibrium case: y is scalar.
                if y.shape and np.sum(y.shape) > 1:
                    if bf not in labels:
                        line = ax.plot(x, y.max(), markersize=default_size, marker=m, c=c, label=bf)
                        points.append(line[0])
                        labels.append(bf)
                    else:
                        ax.plot(x, y.max(), markersize=default_size, marker=m, c=c)
                    ax.plot(x, y.min(), markersize=default_size, marker=m, c=c)
                else:
                    if bf not in labels:
                        line = ax.plot(x, y, markersize=default_size, marker=m, c=c, label=bf)
                        points.append(line[0])
                        labels.append(bf)
                    else:
                        ax.plot(x, y, markersize=default_size, marker=m, c=c)
        return points, labels

    def plot_continuation_grid(self, plots: list, ncols: int = 2, figsize: tuple = None,
                                sharex: bool = False, sharey: bool = False,
                                **shared_kwargs) -> tuple:
        """Lay out multiple 1D/2D continuations as a grid of subplots.

        Convenience helper to compare continuations side-by-side (e.g. a
        codim-1 scan in eta next to its codim-2 fold curve in (eta, Delta),
        or several "same x/y, different parameter setting" diagrams).

        Parameters
        ----------
        plots
            List of dicts, one per subplot. Each dict must contain
            ``'x'``, ``'y'``, and ``'cont'`` (forwarded to
            :meth:`plot_continuation` as the corresponding positional
            args). Any additional keys override ``shared_kwargs`` for that
            specific subplot — except ``'title'``, which is set on the
            subplot's axes via ``ax.set_title``. The optional
            ``'panel_label'`` key is drawn in the upper-left corner of
            the subplot (useful for figure-quality (a), (b), (c)
            annotations).
        ncols
            Number of columns in the grid. Rows are derived from
            ``ceil(len(plots) / ncols)``. Default 2.
        figsize
            Forwarded to :func:`matplotlib.pyplot.subplots`. Defaults to
            ``(5 * ncols, 4 * nrows)`` — i.e. each subplot gets ~5x4 inches.
        sharex, sharey
            Forwarded to :func:`matplotlib.pyplot.subplots`. Useful when
            all panels live on the same parameter range.
        shared_kwargs
            Keyword arguments applied to every subplot (e.g.
            ``bifurcation_legend=False`` to suppress the per-panel
            legend). Per-plot keys override these.

        Returns
        -------
        tuple
            ``(fig, axes, line_cols)``. ``axes`` is the flat list of
            ``Axes`` (length ``len(plots)``; trailing positions in the
            grid that have no plot are deleted). ``line_cols`` is the
            list of LineCollections returned by each ``plot_continuation``
            call, in the same order as ``plots``.
        """
        if not plots:
            raise ValueError("plot_continuation_grid requires at least one plot spec")

        n = len(plots)
        nrows = (n + ncols - 1) // ncols  # ceil(n / ncols)
        if figsize is None:
            figsize = (5 * ncols, 4 * nrows)

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize,
                                 sharex=sharex, sharey=sharey, squeeze=False)
        axes_flat = list(axes.flatten())

        line_cols = []
        for i, spec in enumerate(plots):
            spec = dict(spec)
            title = spec.pop('title', None)
            panel_label = spec.pop('panel_label', None)
            x = spec.pop('x')
            y = spec.pop('y')
            cont = spec.pop('cont')

            # per-plot spec overrides the shared defaults
            kwargs_tmp = dict(shared_kwargs)
            kwargs_tmp.update(spec)

            ax = axes_flat[i]
            line_col = self.plot_continuation(x=x, y=y, cont=cont, ax=ax, **kwargs_tmp)
            line_cols.append(line_col)

            if title is not None:
                ax.set_title(title)
            if panel_label is not None:
                # axes-fraction text in the upper-left, bold, slightly
                # larger than the tick labels — matches the convention
                # used in most published bifurcation figures.
                ax.text(0.02, 0.96, panel_label, transform=ax.transAxes,
                        ha='left', va='top', fontweight='bold', fontsize='large')

        # Hide any trailing subplot positions we didn't use.
        for j in range(n, len(axes_flat)):
            fig.delaxes(axes_flat[j])

        fig.tight_layout()
        return fig, axes_flat[:n], line_cols

    def update_bifurcation_style(self, bf_type: str, marker: str = None, color: str = None) -> None:
        """Update the default marker and color of a given special solution type.

        Parameters
        ----------
        bf_type
            Type of the special solution.
        marker
            New marker type.
        color
            New color.

        Returns
        -------
        None
        """

        if bf_type in self._bifurcation_styles:
            if marker:
                self._bifurcation_styles[bf_type]['marker'] = marker
            if color:
                self._bifurcation_styles[bf_type]['color'] = color
        else:
            if marker is None:
                marker = 'o'
            if color is None:
                color = 'k'
            self._bifurcation_styles.update({bf_type: {'marker': marker, 'color': color}})

    def _create_summary(self, solution: Union[Any, dict], points: list, variables: list, params: list,
                        timeseries: bool, stability: bool, period: bool, eigenvals: bool, lyapunov_exp: bool,
                        reduce_limit_cycle: bool
                        ) -> DataFrame:
        """Creates summary of auto continuation and stores it in dictionary.

        Builds a single dict-of-lists keyed by `(column_name, sub_index)` tuples
        (scalar columns use `''` for the sub-index, matching the legacy shape
        produced by merging a flat-column DataFrame into a MultiIndex one).
        One DataFrame construction at the end replaces the previous two-build
        + column-by-column merge dance.

        Parameters
        ----------
        solution
        points
        variables
        params
        timeseries
        stability
        period
        eigenvals
        lyapunov_exp
        reduce_limit_cycle

        Returns
        -------
        DataFrame
            Continuation summary with a `MultiIndex` columns axis; vector-valued
            quantities (state-var min/max, eigenvalues, lyapunov exponents) use
            integer sub-indices, scalar quantities use empty-string sub-indices.
        """

        # ``col_values`` is the single mutable structure built in the loop.
        # Keys are the (name, sub_index) tuples that go straight into the
        # final MultiIndex; values are per-row lists. dict insertion order
        # determines the final column order.
        col_values: dict = {}
        indices: list = []

        for point in points:
            s, solution_type, solution_idx = self.get_solution(cont=solution, point=point)
            if solution_type == 'No Label' or solution_type == 'MX':
                continue

            indices.append(point)
            var_vals = get_solution_variables(s, variables, timeseries)
            param_vals = get_solution_params(s, params)
            # Parse the diagnostic block at most once per point; reused for
            # stability + eigenvalues + lyapunov below.
            diag = parse_point_diagnostics(s) if (stability or eigenvals or lyapunov_exp) else None
            period_val = (get_solution_params(s, ['PAR(11)'])[0]
                          if (period or lyapunov_exp or eigenvals) else None)

            # Column insertion order is chosen to match the legacy output
            # column order produced by the old data_2d-first / data_1d-then-
            # merged build, so existing user scripts that rely on positional
            # access keep working.

            # --- state-variable values (vector for limit cycles) ---
            for var, val in zip(variables, var_vals):
                if len(val) > 1 and reduce_limit_cycle:
                    col_values.setdefault((var, 0), []).append(np.min(val))
                    col_values.setdefault((var, 1), []).append(np.max(val))
                else:
                    for i, v in enumerate(val):
                        col_values.setdefault((var, i), []).append(v)

            # --- eigenvalues / Floquet multipliers ---
            if eigenvals:
                for i, v in enumerate(diag['eigenvalues']):
                    col_values.setdefault(('eigenvalues', i), []).append(v)

            # --- lyapunov exponents ---
            if lyapunov_exp:
                for i, lyap in enumerate(get_lyapunov_exponents(diag['eigenvalues'], period_val)):
                    col_values.setdefault(('lyapunov_exponents', i), []).append(lyap)

            # --- bifurcation type / index ---
            col_values.setdefault(('bifurcation', ''), []).append(solution_type)
            col_values.setdefault(('bifurcation_index', ''), []).append(solution_idx)

            # --- parameter values ---
            for param, val in zip(params, param_vals):
                col_values.setdefault((param, ''), []).append(val)

            # --- time vector (when get_timeseries=True) ---
            if len(var_vals) > len(variables) and timeseries:
                col_values.setdefault(('time', ''), []).append(var_vals[-1])

            if stability:
                col_values.setdefault(('stability', ''), []).append(bool(diag['stable']))

            if period:
                col_values.setdefault(('period', ''), []).append(period_val)

        # ---- single DataFrame construction ----
        if not col_values:
            return DataFrame(index=indices)

        # Apply _var_map_inv remapping to translate "PAR(i)" / "U(i)" column
        # names back to user-facing names where one is registered. Done once,
        # over the discovered columns, rather than per-point.
        remapped = [
            (self._var_map_inv[name] if name in self._var_map_inv else name, sub)
            for (name, sub) in col_values
        ]
        columns = MultiIndex.from_tuples(remapped)

        # `index` may have one extra element when the last `points` entry was
        # dropped mid-row (legacy `_to_dataframe` fallback handled this by
        # trimming `data[:-1]`). Here every dropped point is filtered out at
        # the top of the loop so `indices` and the per-column lists agree by
        # construction.
        return DataFrame(dict(zip(columns, col_values.values())), index=indices)

    def _call_auto(self, starting_point: Union[str, int], origin: Union[Any, dict], **auto_kwargs) -> Any:
        if starting_point:
            s, solution_name, _ = self.get_solution(point=starting_point, cont=origin)
            if solution_name == "No Label":
                raise KeyError(f"Starting point {starting_point} could not be found on the provided origin branch.")
            solution = self._auto.run(s, **auto_kwargs)
        else:
            solution = self._auto.run(**auto_kwargs)
        return self._start_from_solution(solution)

    def _update_axis_lims(self, ax: Union[plt.Axes, Axes3D], ax_data: list, padding: float = 0.,
                          force_update: bool = False) -> None:
        ax_names = ['x', 'y', 'z']
        for i, data in enumerate(ax_data):
            axis_limits = self._get_axis_lims(np.asarray(data), padding=padding)
            if force_update:
                min_val, max_val = axis_limits
            else:
                min_val, max_val = eval(f"ax.get_{ax_names[i]}lim()")
                min_val, max_val = np.min([min_val, axis_limits[0]]), np.max([max_val, axis_limits[1]])
            eval(f"ax.set_{ax_names[i]}lim(min_val, max_val)")

    def _map_auto_kwargs(self, kwargs: dict) -> dict:

        # handle the continuation parameter
        if "ICP" in kwargs:
            val = kwargs.pop("ICP")
            if type(val) is str:
                kwargs["ICP"] = self._map_var(val)
            elif type(val) in [list, tuple]:
                kwargs["ICP"] = [self._map_var(v) if type(v) is str else v for v in val]
            else:
                kwargs["ICP"] = val

        # handle PAR-keyed dict constants (named -> integer index). On the
        # PyRates-generated path auto-07p resolves names itself via parnames,
        # so unmapped strings pass through harmlessly; on the hand-written
        # path the explicit translation here is what makes string keys work.
        for key in ("UZR", "UZSTOP", "THL", "THU"):
            if key in kwargs:
                d = kwargs.pop(key)
                kwargs[key] = {self._map_var(k) if type(k) is str else k: v for k, v in d.items()}

        return kwargs

    def _map_var(self, var: str, mode: str = "cont"):
        """Translate a user-facing var name to auto-07p's internal form.

        With ``mode="cont"`` returns the integer index (PAR slot for
        parameters, U slot for state variables) for use in ICP / UZR / UZSTOP
        / THL / THU. With ``mode="plot"`` returns the ``"PAR(i)"`` or
        ``"U(i)"`` string used as a DataFrame column key. Unknown names pass
        through unchanged so non-PyCoBi-managed keys (raw ``PAR(i)`` strings,
        bare ints, etc.) still work.
        """
        entry = self._var_map.get(var)
        if entry is None:
            return var
        kind, idx = entry
        if mode == "cont":
            return idx
        # "plot" mode (the only other mode currently used)
        return f"PAR({idx})" if kind == "P" else f"U({idx})"

    @staticmethod
    def _get_all_var_keys(solution):
        # Prefer the solution's own coord names — when PyRates emits `unames`
        # into the auto-07p c.* file (post-pyrates>=1.1), auto exposes state
        # variables under their user-facing name rather than as ``U(i)``.
        # Fall back to the historical ``U(i)`` form for hand-written fortran
        # systems or older PyRates without unames.
        coords = getattr(solution, 'coordnames', None)
        if coords:
            return list(coords)
        return [f'U({i+1})' for i in range(solution['NDIM'])]

    @staticmethod
    def _get_all_param_keys(solution):
        return solution.PAR.coordnames

    def _start_from_solution(self, solution: Any) -> Any:
        """Auto-retry hook for runs that produced only a starting direction.

        When auto-07p's first run returns just a starting-direction diagnostic
        with a single ``EP`` label (no continuation steps taken), this method
        re-invokes ``auto.run`` from that EP to actually kick off the
        continuation. Emits a UserWarning so callers know a second auto.run
        fired — the second call uses no kwargs from the original, so a user
        relying on specific NMX / DS / DSMAX overrides may want to know.

        The retry also pops the EP from the original solution's labels dict
        (auto's API quirk: that popitem is how we extract the seed solution),
        leaving the original solution object structurally modified.
        """
        diag = str(solution[0].diagnostics)
        sol_keys = get_solution_keys(solution)
        if 'Starting direction of the free parameter(s)' in diag and len(sol_keys) == 1 and \
                "EP" in list(solution[0].labels.by_index[sol_keys[0]])[0]:
            warnings.warn(
                "auto-07p's first run took no continuation steps (only a starting-direction "
                "diagnostic was produced); restarting auto.run from the single EP label without "
                "the original run's keyword arguments. If this is unexpected, check that your "
                "auto constants (DS, DSMAX, ICP, ...) are appropriate for the model.",
                UserWarning,
                stacklevel=3,
            )
            _, s = solution[0].labels.by_index.popitem()
            solution = self._auto.run(s['EP']['solution'])
        return solution

    @staticmethod
    def _get_line_collection(x, y, stability=None, line_style_stable='solid', line_style_unstable='dotted',
                             line_color_stable='k', line_color_unstable='gray', **kwargs) -> LineCollection:
        """

        Parameters
        ----------
        x
        y
        stability
        line_style_stable
        line_style_unstable
        line_color_stable
        line_color_unstable
        kwargs

        Returns
        -------
        LineCollection
        """

        # Combine x and y into segment-friendly (N, 2) arrays. The squeeze-
        # then-take-shape[0] pattern raises IndexError on degenerate length-1
        # inputs (the squeezed array is 0-D); a `reshape(-1, 1)` fallback
        # handles both 1-D-length-1 and N-D-length-1 cases without crashing
        # the line-collection construction. The pre-1.0 code had this guard
        # on `x` but not on `y`.
        x = np.asarray(x).reshape(-1, 1)
        if hasattr(y[0], "shape") and sum(y[0].shape) > 1:
            y = np.asarray([y[i] for i in range(y.shape[0])])
            y_max = np.reshape(y.max(axis=1), (y.shape[0], 1))
            y_min = np.reshape(y.min(axis=1), (y.shape[0], 1))
            y_min = np.append(x, y_min, axis=1)
            y = y_max
            add_min = True
        else:
            y = np.asarray(y).reshape(-1, 1)
            add_min = False
        y = np.append(x, y, axis=1)

        # if stability was passed, collect indices for stable line segments
        ###################################################################

        # The size>1 guard skips the segmentation logic for degenerate
        # 1-element stability arrays (which would yield a single empty
        # segment via the diff). The pre-1.0 form used
        # `np.sum(stability.shape) > 1` — semantically equivalent for 1-D
        # arrays but cryptic; switched to the explicit `.size > 1`.
        if stability is not None and np.asarray(stability).size > 1:

            # collect indices
            stability = np.asarray(stability, dtype='int')
            stability_changes = np.concatenate([np.zeros((1,)), np.diff(stability)])
            idx_changes = np.sort(np.argwhere(stability_changes != 0))
            idx_changes = np.append(idx_changes, len(stability_changes))

            # create line segments
            lines, styles, colors = [], [], []
            idx_old = 1
            for idx in idx_changes:
                lines.append(y[idx_old-1:idx, :])
                styles.append(line_style_stable if stability[idx_old] else line_style_unstable)
                colors.append(line_color_stable if stability[idx_old] else line_color_unstable)
                if add_min:
                    lines.append(y_min[idx_old - 1:idx, :])
                    styles.append(line_style_stable if stability[idx_old] else line_style_unstable)
                    colors.append(line_color_stable if stability[idx_old] else line_color_unstable)
                idx_old = idx

        else:

            lines = [y, y_min] if add_min else [y]
            styles = [line_style_stable, line_style_stable] if add_min else [line_style_stable]
            colors = [line_color_stable, line_color_stable] if add_min else [line_color_stable]

        # Pop the two LineCollection kwargs we set explicitly so a user
        # passing `colors=` or `linestyles=` through one of the plot_* helpers
        # overrides the computed values rather than colliding with them
        # (passing both via `**kwargs` raises TypeError).
        colors = kwargs.pop('colors', colors)
        styles = kwargs.pop('linestyles', styles)
        return LineCollection(segments=lines, linestyles=styles, colors=colors, **kwargs)

    @staticmethod
    def _get_3d_line_collection(x, y, z, stability=None, line_style_stable='solid', line_style_unstable='dotted',
                                **kwargs) -> Line3DCollection:
        """

        Parameters
        ----------
        x
        y
        z
        stability
        line_style_stable
        line_style_unstable
        kwargs

        Returns
        -------
        Line3DCollection
        """

        # combine y and param vals
        x = np.reshape(x, (x.squeeze().shape[0], 1))
        y = np.reshape(y, (y.squeeze().shape[0], 1))
        z = np.reshape(z, (z.squeeze().shape[0], 1))
        y = np.append(x, y, axis=1)
        y = np.append(y, z, axis=1)

        # if stability was passed, collect indices for stable line segments
        ###################################################################

        if stability is not None and np.asarray(stability).size > 1:

            # collect indices
            stability = np.asarray(stability, dtype='int')
            stability_changes = np.concatenate([np.zeros((1,)), np.diff(stability)])
            idx_changes = np.sort(np.argwhere(stability_changes != 0))
            idx_changes = np.append(idx_changes, len(stability_changes))

            # create line segments
            lines, styles = [], []
            idx_old = 1
            for idx in idx_changes:
                lines.append(y[idx_old - 1:idx, :])
                styles.append(line_style_stable if stability[idx_old] else line_style_unstable)
                idx_old = idx

        else:

            lines = [y]
            styles = [line_style_stable]

        # create line collection
        array = kwargs.pop('array', 'x')
        # Same as `_get_line_collection`: pop `linestyles` so a user-supplied
        # value overrides the per-stability-block styles rather than colliding
        # with the explicit kwarg below.
        styles = kwargs.pop('linestyles', styles)
        # NOTE: Line3DCollection takes `lines` as a positional argument (not
        # `segments=` like the 2D LineCollection). Passing `segments=lines`
        # used to silently fail on older matplotlib and now raises outright.
        line_col = Line3DCollection(lines, linestyles=styles, **kwargs)

        # post-processing
        if array == 'x':
            array = x.squeeze()
        elif array == 'y':
            array = y[:, 1].squeeze()
        elif array == 'z':
            array = z.squeeze()
        line_col.set_array(array)

        return line_col

    @staticmethod
    def _get_axis_lims(x: np.array, padding: float = 0.) -> tuple:
        x_min, x_max = x.min(), x.max()
        x_pad = (x_max - x_min) * padding
        return x_min - x_pad, x_max + x_pad


def _extract_merge_point(p: int, df: DataFrame) -> Series:
    p_tmp = df.loc[p, :]
    if len(p_tmp.shape) > 1 and p_tmp.shape[0] > 1:
        return p_tmp.iloc[0, :]
    return p_tmp
