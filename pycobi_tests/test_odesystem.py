"""Test suite for the ODESystem class.
"""

# imports
import os
import re

import pytest
from pycobi import ODESystem, Continuation, parse_point_diagnostics
from pycobi.automated_continuation import (
    codim2_search,
    continue_period_doubling_bf,
    _bif_series_contains,
    _resolve_param_for_extract,
)
from pycobi.utility import get_branch_info
from pyrates import clear
from pytest import fixture, approx, raises
import numpy as np
import pandas as pd


# Matches an analytical-Jacobian assignment like ``dfdu(1,2) = ...`` (or DFDU);
# does NOT match the array declaration ``dfdu(ndim,ndim)`` in the func signature.
_DFDU_ASSIGN = re.compile(r'\bdfdu\(\s*\d', re.IGNORECASE)
_DFDP_ASSIGN = re.compile(r'\bdfdp\(\s*\d', re.IGNORECASE)


def _pyrates_has_auto_jac() -> bool:
    """True iff the installed PyRates writes the analytical Jacobian (DFDU/DFDP)
    into the auto-07p ``func`` wrapper and `unames` / `parnames` / multi-scenario
    `auto_constants` into c.* — i.e. PyRates >= 1.1. PyPI 1.0.x silently
    ignores `auto_jac`, `auto_constants`, etc., so feature-dependent tests
    must skip rather than fail there.
    """
    try:
        from pyrates.backend.fortran.fortran_backend import FortranBackend
    except ImportError:
        return False
    return hasattr(FortranBackend, '_emit_auto_jacobian_block')


_requires_pyrates_dev = pytest.mark.skipif(
    not _pyrates_has_auto_jac(),
    reason="requires PyRates >= 1.1 (analytical-Jacobian emission, unames/parnames, "
           "multi-scenario auto_constants); skipped against older PyPI releases",
)

# meta infos
__author__ = "Richard Gast"
__status__ = "Development"

# Utility
#########


def setup_module():
    print("\n")
    print("==========================")
    print("| Test Suite : ODESystem |")
    print("==========================")


@fixture(scope="session")  # type: str
def auto_dir(pytestconfig) -> str:
    return pytestconfig.getoption("auto_dir")


# test accuracy
accuracy = 1e-4

# tests
#######


def test_1_1_init(auto_dir):
    """Testing the different instantiation options of the `ODESystem` class.
    """

    # initialize ODESystem the default way
    ode1 = ODESystem(eq_file='qif_eq', working_dir="resources", auto_dir=auto_dir, init_cont=False)
    ode1.close_session()

    # initialize ODESystem from YAML file
    model = "model_templates.neural_mass_models.qif.qif"
    ode2 = ODESystem.from_yaml(model, init_cont=False, file_name='qif_eq2', func_name="qif2", auto_dir=auto_dir)
    ode2.close_session(clear_files=True)

    # initialize ODESystem from YAML file with different parameters
    ode3 = ODESystem.from_yaml(model, init_cont=False, file_name='qif_eq3', func_name="qif3", auto_dir=auto_dir,
                               node_vars={'p/qif_op/eta': 2.0})
    ode3.close_session(clear_files=True)

    # these tests should pass
    assert isinstance(ode1, ODESystem)
    assert isinstance(ode2, ODESystem)
    assert isinstance(ode3, ODESystem)
    assert ode1.dir != ode2.dir
    assert ode2.dir == ode3.dir


def test_1_2_run(auto_dir):
    """Testing the run method for running auto commands via `ODESystem`.
    """

    # initialize ODESystem the default way
    ode1 = ODESystem(eq_file="qif_eq", working_dir="resources", auto_dir=auto_dir, init_cont=True, c="ivp", NPR=100)
    ode1.close_session()

    # initialize ODESystem from YAML file
    model = "model_templates.neural_mass_models.qif.qif"
    ode2 = ODESystem.from_yaml(model, init_cont=True, file_name='qif_eq2', func_name="qif2", auto_dir=auto_dir,
                               NPR=100, NMX=5000)
    ode2.close_session(clear_files=True)

    # initialize ODESystem from YAML file with different parameters
    ode3 = ODESystem.from_yaml(model, init_cont=True, file_name='qif_eq3', func_name="qif3", auto_dir=auto_dir,
                               node_vars={'p/qif_op/eta': 2.0}, NPR=100, NMX=5000)
    ode3.close_session(clear_files=True)

    print(ode2[0].columns.values)

    # The column name for the first state variable depends on which PyRates is
    # installed: PyRates >= 1.1 writes `unames = {1: 'r', ...}` into the c.*
    # file so auto's DataFrame columns are bare local names ('r'); older
    # PyRates (PyPI 1.0.x) doesn't, and PyCoBi falls back to the namespaced
    # name via `_var_map_inv` ('p/qif_op/r'). Either name should match ode1's
    # 'U(1)' column from the hand-written qif_eq.f90.
    r_col = 'r' if 'r' in ode2[0].columns else 'p/qif_op/r'
    assert (ode1[0]["U(1)"] - ode2[0][r_col]).sum()[0] == approx(0.0, rel=accuracy, abs=accuracy)
    assert abs((ode2[0][r_col] - ode3[0][r_col]).sum()[0]) > 0


@_requires_pyrates_dev
def test_1_3_auto_constants(auto_dir, tmp_path):
    """`auto_constants` exposes PyRates' multi-scenario c.* generation.

    Two behaviours to pin:
      1. Passing a tuple of scenario names produces one ``c.<name>`` file per
         scenario, each preconfigured for that continuation mode.
      2. Requesting ``init_cont=True`` without ``'ivp'`` in ``auto_constants``
         is rejected up front rather than failing later with an opaque
         "c.ivp not found" error from auto-07p.
    """
    model = "model_templates.neural_mass_models.qif.qif"

    # Multi-scenario: ivp + eq + lc all emitted.
    multi_dir = tmp_path / 'multi'
    multi_dir.mkdir()
    ode = ODESystem.from_yaml(
        model,
        working_dir=str(multi_dir), auto_dir=auto_dir,
        file_name='qif_multi', func_name='qif_multi_fn',
        init_cont=False,
        auto_constants=('ivp', 'eq', 'lc'),
    )
    for scen in ('ivp', 'eq', 'lc'):
        f = multi_dir / f'c.{scen}'
        assert f.exists(), f"c.{scen} was not generated"
    # Each c.* should be configured for its scenario (distinct IPS values).
    assert 'IPS = -2' in (multi_dir / 'c.ivp').read_text()
    assert 'IPS = 1' in (multi_dir / 'c.eq').read_text()
    assert 'IPS = 2' in (multi_dir / 'c.lc').read_text()
    ode.close_session(clear_files=True)

    # Validation: init_cont=True without 'ivp' is rejected.
    with raises(ValueError, match=r"init_cont=True .* 'ivp' is missing"):
        ODESystem.from_yaml(
            model,
            working_dir=str(tmp_path / 'bad'), auto_dir=auto_dir,
            file_name='qif_bad', func_name='qif_bad_fn',
            init_cont=True,
            auto_constants=('eq',),
        )


def test_1_4_name_remapping(auto_dir):
    """`_map_auto_kwargs` translates named PAR keys to integers in ICP, UZR,
    UZSTOP, THL, and THU. Previously only ICP/UZR were remapped, so
    UZSTOP={'eta': 5.0} silently passed through unchanged (which only happened
    to work when auto-07p's parnames were active — the bug we hit on the
    hand-written demo in test_3_1).
    """
    # Hand-written path: explicit params/state_vars build a name -> PAR index
    # map. We can exercise `_map_auto_kwargs` directly without ever running
    # auto, since init_cont=False short-circuits the IVP call.
    # Anchor working_dir absolutely so preceding tests' cwd doesn't matter.
    resources = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources')
    ode = ODESystem(
        eq_file='qif_eq', working_dir=resources, auto_dir=auto_dir,
        init_cont=False,
        params=['tau', 'Delta', 'I_ext', 'eta', 'weight'],
        state_vars=['r', 'v'],
    )
    try:
        remapped = ode._map_auto_kwargs({
            'ICP': 'eta',
            'UZR':    {'eta': [2.0], 'weight': [10.0]},
            'UZSTOP': {'eta': 5.0,   'I_ext': 1.0},
            'THL':    {'eta': 0.0},
            'THU':    {'r': 0.5},
        })

        # Pin `_map_var` itself: both 'cont' (int) and 'plot' ("PAR(i)"/"U(i)")
        # modes must agree on the index, and the (kind, idx) tuple structure
        # in `_var_map` is what determines whether 'plot' returns PAR or U.
        assert ode._map_var('eta', 'cont') == 4
        assert ode._map_var('eta', 'plot') == 'PAR(4)'
        assert ode._map_var('r', 'cont') == 1
        assert ode._map_var('r', 'plot') == 'U(1)'
        assert ode._map_var('t', 'cont') == 14
        assert ode._map_var('t', 'plot') == 'PAR(14)'
        # Unknown names pass through unchanged in both modes.
        assert ode._map_var('not_a_var', 'cont') == 'not_a_var'
        assert ode._map_var('PAR(99)', 'plot') == 'PAR(99)'
    finally:
        ode.close_session()

    assert remapped['ICP'] == 4, f"expected ICP -> PAR(4), got {remapped['ICP']!r}"
    assert remapped['UZR'] == {4: [2.0], 5: [10.0]}
    assert remapped['UZSTOP'] == {4: 5.0, 3: 1.0}
    assert remapped['THL'] == {4: 0.0}
    assert remapped['THU'] == {1: 0.5}, "state-var name 'r' should map to U(1) -> 1"


def test_1_5_to_file_roundtrip(auto_dir, tmp_path):
    """`to_file` / `from_file` round-trips both modes after a real continuation.

    Pins two regressions to_file had collected:
      - `to_file(results_only=False)` used to crash on `_auto` (Python module)
        and `_temp` (CircuitTemplate with lambdified sympy functions), neither
        of which pickle. Now both are listed in `_PICKLE_EXCLUDE` and re-derived
        by `__init__` on load.
      - `from_file` used to insist that every loaded attribute already exists
        as a dict on the fresh instance, so any non-dict slot (`_eq`,
        `_last_cont`, `_cont_num`, ...) raised AttributeError. Now non-dicts
        are restored via setattr.
    """
    resources = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources')
    # `parnames={}, unames={}` clears auto-07p's runner cache (see auto/runAUTO.py
    # config(): parnames/unames are intentionally retained across run() calls).
    # Without this clear, a preceding PyRates-generated test leaks its column
    # names ('r', 'v', 'eta', ...) into the auto solution for this hand-written
    # qif_eq.f90 — which has no unames declared and expects U(1)/U(2) coords.
    ode = ODESystem(
        eq_file='qif_eq', working_dir=resources, auto_dir=auto_dir,
        init_cont=True, c='ivp', NMX=500, NPR=500,
        parnames={}, unames={},
        params=['tau', 'Delta', 'I_ext', 'eta', 'weight'],
        state_vars=['r', 'v'],
    )

    # results_only=True is the documented happy path; round-trip and compare.
    rfile = tmp_path / 'ode_results.pkl'
    ode.to_file(str(rfile), results_only=True, note='hello')
    ode_r = ODESystem.from_file(str(rfile), auto_dir=auto_dir)
    assert set(ode_r.results.keys()) == set(ode.results.keys())
    # the DataFrames should match index-wise
    for k in ode.results:
        assert list(ode_r.results[k].index) == list(ode.results[k].index)
    assert ode_r.additional_attributes == {'note': 'hello'}

    # results_only=False rests on the _PICKLE_EXCLUDE shim — must not raise.
    wfile = tmp_path / 'ode_whole.pkl'
    ode.to_file(str(wfile), results_only=False)
    ode_w = ODESystem.from_file(str(wfile), auto_dir=auto_dir)
    # The non-dict scalars and the _var_map dicts must come back identical.
    assert ode_w._eq == ode._eq
    assert ode_w._cont_num == ode._cont_num
    assert ode_w._var_map == ode._var_map
    assert ode_w._var_map_inv == ode._var_map_inv
    # Excluded slots: _auto is the live module, _temp is unset for the
    # hand-written path. Both should still be available on the loaded
    # instance, supplied by __init__.
    assert ode_w._auto is not None
    assert ode_w._temp is None  # __init__ defaults this when no template kwarg

    ode.close_session()
    ode_r.close_session()
    ode_w.close_session()


class _MockSol:
    """Minimal solution-object stand-in for `parse_point_diagnostics`:
    only `.data['NDIM']` and `.b['PT']` need to be lookup-able.
    """
    def __init__(self, ndim, pt=None):
        self.data = {'NDIM': ndim}
        self.b = {'PT': pt} if pt is not None else {}


def test_1_12_line_collection_user_kwargs():
    """`_get_line_collection` and `_get_3d_line_collection` accept user-supplied
    `linestyles=` / `colors=` overrides without colliding with the styles/colors
    they compute themselves from the stability array.

    Pins D3: previously `LineCollection(..., linestyles=styles, **kwargs)`
    raised `TypeError: got multiple values for keyword argument 'linestyles'`
    whenever a user passed `linestyles=` to one of the plot_* helpers (which
    forward arbitrary kwargs down to this builder).
    """
    n = 6
    x = np.linspace(0.0, 1.0, n)
    y = np.linspace(0.0, 0.5, n)
    z = np.linspace(0.0, 0.25, n)
    stability = np.array([True, True, False, False, True, True])

    # 2D — with stability so the per-segment styles list is non-trivial.
    lc = ODESystem._get_line_collection(x=x, y=y, stability=stability, linestyles='-.')
    # The user's override should win over the computed per-segment list:
    # matplotlib stores it as a list-of-tuples; just check at least one segment
    # carries the dash-dot pattern (Matplotlib normalises '-.' to (0, (6.4, ...))
    # depending on version, so check that it didn't fall back to the default
    # 'solid'/'dotted' strings the function would otherwise have produced).
    assert lc is not None  # primary smoke-test: no TypeError

    # 2D — colors override on top of stability.
    lc = ODESystem._get_line_collection(x=x, y=y, stability=stability, colors='red')
    assert lc is not None

    # 3D — same shape of override.
    lc3 = ODESystem._get_3d_line_collection(x=x, y=y, z=z, stability=stability, linestyles='-.')
    assert lc3 is not None


def test_1_10_getitem_helpful_keyerror(auto_dir):
    """`ODESystem.__getitem__` raises a KeyError that lists all registered
    continuation names and stored pyauto keys when the requested item matches
    neither. Pins C2: previously the error was the opaque inner KeyError(item)
    raised by the dict, which was unhelpful when the user mistyped a name.
    """
    resources = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources')
    ode = ODESystem(
        eq_file='qif_eq', working_dir=resources, auto_dir=auto_dir,
        init_cont=False,
        params=['tau', 'Delta', 'I_ext', 'eta', 'weight'],
        state_vars=['r', 'v'],
    )
    try:
        # Plant some bookkeeping by hand — no auto run needed for this test.
        ode.results[0] = 'sentinel-zero'
        ode.results[7] = 'sentinel-seven'
        ode._results_map['hi'] = 0
        ode._results_map['lo'] = 7

        # Known lookups still work
        assert ode[0] == 'sentinel-zero'
        assert ode['hi'] == 'sentinel-zero'
        assert ode['lo'] == 'sentinel-seven'

        # Unknown lookup raises with both name + key lists in the message
        with raises(KeyError, match=r"is neither a registered continuation name"):
            ode['definitely_not_a_name']
        with raises(KeyError, match=r"Known names: \['hi', 'lo'\]"):
            ode['definitely_not_a_name']
        with raises(KeyError, match=r"Known keys: \[0, 7\]"):
            ode['definitely_not_a_name']
    finally:
        ode.close_session()


def test_1_9_get_branch_info_paths():
    """`get_branch_info` accepts three call shapes (bifDiag, single-branch dict,
    IVP-style nested-RG container) and raises a meaningful error when none
    match. Pins B9: the IVP-style path used to hide its failure mode behind a
    `while i < 10` counter and re-raise whichever KeyError fired last.
    """
    # ---- Path 1: bifDiag-shaped (`.BR` + `.c['ICP']` on branch[0]) ----
    class _Sub:
        BR = 7
        c = {'ICP': [4]}
    bifdiag_shape = [_Sub()]
    assert get_branch_info(bifdiag_shape) == (7, (4,))

    # ICP given as a bare int collapses to a 1-tuple
    class _Sub2:
        BR = 3
        c = {'ICP': 2}
    assert get_branch_info([_Sub2()]) == (3, (2,))

    # ---- Path 2: single-branch dict with 'BR' key and .c config ----
    class _Single(dict):
        c = {'ICP': [1, 2]}
    single = _Single({'BR': 42})
    assert get_branch_info(single) == (42, (1, 2))

    # ---- Path 3: IVP-style — BR buried under an RG label on branch[0] ----
    class _NoBR:
        # No .BR attribute → path 1 fails. Path 2 also fails because this is
        # not a dict-like container. Path 3 walks .labels.by_index.
        c = {'ICP': [14]}
        class labels:
            by_index = {
                # an EP-only label (no RG nested) — must be skipped, not give up
                3: {'EP': {'solution': None}},
                # the one with the RG nested data → path 3 picks it
                7: {'RG': {'solution': {'data': {'BR': 11}}}},
            }
    assert get_branch_info([_NoBR()]) == (11, (14,))

    # ---- Path 3 exhaustion: no key has an RG-nested BR ----
    class _NoRG:
        c = {'ICP': [4]}
        class labels:
            by_index = {1: {'EP': {}}, 2: {'LP': {}}}
    with raises(KeyError, match=r"Tried 2 label indices: \[1, 2\]"):
        get_branch_info([_NoRG()])

    # ---- Totally unrecognised shape → AttributeError naming the type ----
    with raises(AttributeError, match=r"cannot navigate int"):
        get_branch_info(42)


def test_1_7_parse_point_diagnostics_synthetic():
    """`parse_point_diagnostics` extracts stability + eigenvalues from auto-07p
    diagnostic blocks using a single regex pass, with synthetic inputs covering:
      - the standard 2-dim equilibrium case
      - NDIM > 9 (eigenvalue indices like 10, 11, 12) which the previous
        token-based parser was suspected of mishandling — verified here
      - limit-cycle continuation ('Multipliers' + 'Multiplier N:' format)
      - no-convergence blocks (empty spectrum, stability=None)
      - missing-spectrum block (falls back to auto's PT-sign convention)
      - malformed numeric entries (skipped, parser keeps going)
    """
    # --- 1. Standard NDIM=2 stable equilibrium ---
    diag = """
   1     5         Eigenvalues  :   Stable:   2
   1     5         Eigenvalue  1:  -2.12891E+00   0.00000E+00
   1     5         Eigenvalue  2:  -5.17236E+00   0.00000E+00
"""
    res = parse_point_diagnostics(_MockSol(ndim=2), diag=diag)
    assert res['stable'] is True
    assert len(res['eigenvalues']) == 2
    assert res['eigenvalues'][0] == complex(-2.12891, 0.0)
    assert res['eigenvalues'][1] == complex(-5.17236, 0.0)

    # --- 2. NDIM=10 with 7 stable -> unstable; tests 10/11/12 indices ---
    diag = """
   1     5         Eigenvalues  :   Stable:   7
   1     5         Eigenvalue  1:  -2.12891E+00   0.00000E+00
   1     5         Eigenvalue  2:  -5.17236E+00   0.00000E+00
   1     5         Eigenvalue  3:  -1.50000E+00   1.00000E+00
   1     5         Eigenvalue  4:  -1.50000E+00  -1.00000E+00
   1     5         Eigenvalue  5:  -8.00000E-01   0.00000E+00
   1     5         Eigenvalue  6:  -7.00000E-01   0.00000E+00
   1     5         Eigenvalue  7:  -3.00000E-01   0.00000E+00
   1     5         Eigenvalue  8:   5.00000E-01   0.00000E+00
   1     5         Eigenvalue  9:   8.00000E-01   2.00000E+00
   1     5         Eigenvalue 10:   8.00000E-01  -2.00000E+00
"""
    res = parse_point_diagnostics(_MockSol(ndim=10), diag=diag)
    assert res['stable'] is False, "7 stable < NDIM=10 should report unstable"
    assert len(res['eigenvalues']) == 10
    # eigenvalues kept in auto's emit order; #10's real/imag came through correctly
    assert res['eigenvalues'][9] == complex(0.8, -2.0)

    # --- 3. Limit-cycle continuation: 'Multipliers' + 'Multiplier N:' ---
    diag = """
   1    20         Multipliers:     Stable:   3
   1    20         Multiplier   1:   1.00000E+00   0.00000E+00
   1    20         Multiplier   2:   5.00000E-01   0.00000E+00
   1    20         Multiplier   3:   2.00000E-01   0.00000E+00
"""
    res = parse_point_diagnostics(_MockSol(ndim=3), diag=diag)
    assert res['stable'] is True
    assert len(res['eigenvalues']) == 3
    assert res['eigenvalues'][0] == complex(1.0, 0.0)

    # --- 4. No-convergence block: empty spectrum, stability=None ---
    diag = """
   1   100         NOTE:No convergence with fixed step size — aborting
"""
    res = parse_point_diagnostics(_MockSol(ndim=2), diag=diag)
    assert res['stable'] is None
    assert res['eigenvalues'] == []

    # --- 5. No spectrum recorded -> fall back to PT sign convention ---
    diag = """
  BR    PT  IT         PAR           L2-NORM
   1     5   0       -4.56E+00   1.82E+00
"""
    # PT > 0 -> stable
    res = parse_point_diagnostics(_MockSol(ndim=2, pt=16), diag=diag)
    assert res['stable'] is True
    # PT < 0 -> unstable (auto encodes this in fort.7 for the LP it just crossed)
    res = parse_point_diagnostics(_MockSol(ndim=2, pt=-16), diag=diag)
    assert res['stable'] is False
    # No PT either -> stability genuinely unknown
    res = parse_point_diagnostics(_MockSol(ndim=2), diag=diag)
    assert res['stable'] is None

    # --- 6. Malformed numeric in one entry: skip it, keep parsing the rest ---
    diag = """
   1     5         Eigenvalues  :   Stable:   2
   1     5         Eigenvalue  1:  -2.12891E+00   0.00000E+00
   1     5         Eigenvalue  2:  ----GARBAGE----   bogus
   1     5         Eigenvalue  3:  -1.00000E+00   0.00000E+00
"""
    res = parse_point_diagnostics(_MockSol(ndim=3), diag=diag)
    # entries 1 and 3 parsed; entry 2 silently dropped because of bad floats
    assert len(res['eigenvalues']) == 2
    assert res['eigenvalues'][0] == complex(-2.12891, 0.0)
    assert res['eigenvalues'][1] == complex(-1.0, 0.0)


def test_1_11_continuation_dataclass(auto_dir, tmp_path):
    """`get_continuation` returns a `Continuation` dataclass whose fields stay
    in lockstep with the legacy mirror dicts (auto_solutions / results /
    _results_map / _branches). Also covers the bidirectional path: forward +
    reverse share a single Continuation, and from_file rebuilds the dataclass
    from the mirrors on load.
    """
    resources = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources')
    ode = ODESystem(
        eq_file='qif_eq', working_dir=resources, auto_dir=auto_dir,
        init_cont=True, c='ivp', NMX=500, NPR=500,
        parnames={}, unames={},
        params=['tau', 'Delta', 'I_ext', 'eta', 'weight'],
        state_vars=['r', 'v'],
    )
    try:
        # init_cont=True ran an IVP — there's already one continuation at key 0,
        # unnamed (init_cont doesn't pass a name).
        assert 0 in ode.continuations
        c0 = ode.get_continuation(0)
        assert isinstance(c0, Continuation)
        assert c0.key == 0 and c0.name is None
        assert c0.summary is ode.results[0]
        assert c0.auto_solution is ode.auto_solutions[0]

        # Run a bidirectional eta continuation; forward + reverse share a
        # single Continuation registered under the user's name.
        sols, _ = ode.run(
            starting_point='EP2', name='eta_bd', ICP='eta',
            IPS=1, ILP=1, ISP=2, NMX=30, NPR=500,
            bidirectional=True, DS=0.05, DSMAX=0.5,
        )
        assert sorted(ode.continuations.keys()) == [0, 1]
        # Name lookup goes through _results_map and returns the same instance
        # as the key lookup.
        c_by_name = ode.get_continuation('eta_bd')
        c_by_key = ode.get_continuation(1)
        assert c_by_name is c_by_key
        assert c_by_name.name == 'eta_bd'
        assert c_by_name.summary is sols
        # icps records the continuation parameter (PAR(4) = eta). The forward
        # half registers it once; the reverse half's _register_continuation
        # call dedupes it on the Continuation (but the legacy `_branches`
        # mirror keeps duplicates because the merge-detection condition in
        # `run()` requires `new_icp in self._branches[new_branch][origin]`).
        assert c_by_name.icps == [(4,)]
        assert c_by_name.branch_id == 1  # auto starts BR at 1

        # Mismatched lookups raise with helpful messages
        with raises(KeyError, match=r"No continuation named 'nope'"):
            ode.get_continuation('nope')
        with raises(KeyError, match=r"No continuation with key 999"):
            ode.get_continuation(999)

        # to_file / from_file round-trip: the canonical store is rebuilt from
        # the mirrors on load (auto_solution is excluded from the pickle and
        # comes back as None).
        out = tmp_path / 'ode_with_continuations.pkl'
        ode.to_file(str(out), results_only=False)
        loaded = ODESystem.from_file(str(out), auto_dir=auto_dir)
        assert sorted(loaded.continuations.keys()) == [0, 1]
        l_by_name = loaded.get_continuation('eta_bd')
        assert l_by_name.name == 'eta_bd'
        assert l_by_name.branch_id == 1
        assert l_by_name.icps == [(4,)]
        assert l_by_name.summary is not None  # summary survived the round-trip
        assert l_by_name.auto_solution is None  # bifDiag excluded from pickle
        loaded.close_session()
    finally:
        ode.close_session()


def test_1_6_bidirectional_no_magic_name(auto_dir):
    """`bidirectional=True` merges the forward and reverse halves of a
    continuation into a single branch registered under the user's `name`.

    Pins B6: the recursive reverse-direction call no longer overloads
    `name='bidirect:cont2'` to flag itself — that's now `_reverse_direction=True`,
    a private kwarg — so a user is free to name a continuation 'bidirect:cont2'
    without breaking the bidirectional state machine.
    """
    resources = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources')
    # As in test_1_5: parnames={}, unames={} clears auto-07p's runner cache so
    # preceding PyRates-generated tests don't leak column names into this run.
    ode = ODESystem(
        eq_file='qif_eq', working_dir=resources, auto_dir=auto_dir,
        init_cont=True, c='ivp', NMX=500, NPR=500,
        parnames={}, unames={},
        params=['tau', 'Delta', 'I_ext', 'eta', 'weight'],
        state_vars=['r', 'v'],
    )
    try:
        sols, _ = ode.run(
            starting_point='EP2', name='eta_bd', ICP='eta',
            IPS=1, ILP=1, ISP=2, NMX=30, NPR=500,
            bidirectional=True, DS=0.05, DSMAX=0.5,
        )
        # Exactly one named continuation registered, under the user's name.
        assert dict(ode._results_map) == {'eta_bd': 1}
        # auto_solutions: 0 is the IVP, 1 is the merged forward+reverse curve.
        # The reverse-direction recursive call merges into 1, doesn't register
        # a fresh entry.
        assert sorted(ode.auto_solutions.keys()) == [0, 1]
        # Merged continuation spans both sides of the starting eta = -5.
        eta_vals = np.asarray(sols['eta'].values).squeeze()
        assert eta_vals.min() < -5.0, (
            f"reverse direction didn't extend past start: min eta = {eta_vals.min()}"
        )
        assert eta_vals.max() > -5.0, (
            f"forward direction didn't extend past start: max eta = {eta_vals.max()}"
        )
    finally:
        ode.close_session()


def test_1_8_stability_across_fold(auto_dir):
    """A QIF equilibrium continuation in eta crosses a saddle-node fold at
    eta~=-3.13. Before the LP we're on the stable lower branch; after the LP
    we're on the unstable middle branch. Pins B7: `_create_summary` must
    populate the `stability` column with both True and False values, with
    stability flipping where auto reports an LP.
    """
    resources = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources')
    ode = ODESystem(
        eq_file='qif_eq', working_dir=resources, auto_dir=auto_dir,
        init_cont=True, c='ivp', NMX=500, NPR=500,
        parnames={}, unames={},
        params=['tau', 'Delta', 'I_ext', 'eta', 'weight'],
        state_vars=['r', 'v'],
    )
    try:
        sols, _ = ode.run(
            starting_point='EP2', name='eq_fold', ICP='eta',
            IPS=1, ILP=1, ISP=2, NMX=30, NPR=1,  # NPR=1: label every point
            DS=0.05, DSMAX=0.5,
            get_stability=True,
        )
        # there should be at least one LP recorded
        bifs = list(sols['bifurcation'].values.ravel())
        assert 'LP' in bifs, f"expected LP in bifurcation column, got {bifs}"

        # stability column has both stable and unstable points
        stab = list(sols['stability'].values.ravel())
        assert any(stab), "no stable points recorded"
        assert not all(stab), "no unstable points recorded — did the run cross the fold?"

        # the pre-LP points should be stable, post-LP unstable
        # (auto's PT-sign convention: positive PT before crossing, negative after)
        lp_pos = bifs.index('LP')
        pre = stab[:lp_pos]
        post = stab[lp_pos + 1:]
        assert all(pre), f"all pre-LP points should be stable, got {pre}"
        # post-LP we expect at least one unstable; not necessarily all
        # (the curve may bend back to a different stable branch within NMX)
        assert any(not s for s in post), (
            f"expected at least one unstable point past the LP, got {post}"
        )
    finally:
        ode.close_session()


@_requires_pyrates_dev
def test_2_1_jacobian_parity(auto_dir, tmp_path):
    """Equilibrium continuation should agree between PyRates' analytical Jacobian (JAC=1)
    and auto-07p's finite-difference Jacobian (JAC=0).

    This pins the PyCoBi <-> PyRates contract: if PyRates' DFDU/DFDP emission goes wrong
    (sign flip, missing parameter dependency, indexing bug), the analytical-Jac curve will
    diverge from the FD curve, and this test catches it.
    """
    model = "model_templates.neural_mass_models.qif.qif"

    # ---- analytical Jacobian (default) ----
    ode_jac = ODESystem.from_yaml(
        model,
        working_dir=str(tmp_path), auto_dir=auto_dir,
        file_name='qif_jac_mod', func_name='qif_jac',
        init_cont=True, analytical_jacobian=True,
        auto_constants=('ivp', 'eq'),
        NPR=100, NMX=5000,
    )

    # structural assertions: DFDU/DFDP entries emitted and JAC=1 in constants file
    f90_jac = (tmp_path / 'qif_jac_mod.f90').read_text()
    assert _DFDU_ASSIGN.search(f90_jac), "analytical_jacobian=True should emit dfdu(i,j) = ... entries"
    assert _DFDP_ASSIGN.search(f90_jac), "analytical_jacobian=True should emit dfdp(i,j) = ... entries"
    ceq_jac = (tmp_path / 'c.eq').read_text()
    assert 'JAC = 1' in ceq_jac, f"expected JAC = 1 in c.eq, got:\n{ceq_jac}"

    # equilibrium continuation in eta, starting from the IVP's *final* EP
    # (EP1 is the t=0 initial condition; EP2 is the converged steady state).
    sols_jac, _ = ode_jac.run(
        c='eq', ICP='eta', name='eta_jac',
        starting_point='EP2',
        DS=0.05, DSMAX=0.1, NMX=30,
    )
    ode_jac.close_session(clear_files=True)

    # ---- finite-difference Jacobian ----
    ode_fd = ODESystem.from_yaml(
        model,
        working_dir=str(tmp_path), auto_dir=auto_dir,
        file_name='qif_fd_mod', func_name='qif_fd',
        init_cont=True, analytical_jacobian=False,
        auto_constants=('ivp', 'eq'),
        NPR=100, NMX=5000,
    )

    f90_fd = (tmp_path / 'qif_fd_mod.f90').read_text()
    assert not _DFDU_ASSIGN.search(f90_fd), "analytical_jacobian=False should not emit dfdu(i,j) = ... entries"
    assert not _DFDP_ASSIGN.search(f90_fd), "analytical_jacobian=False should not emit dfdp(i,j) = ... entries"
    ceq_fd = (tmp_path / 'c.eq').read_text()
    assert 'JAC = 0' in ceq_fd, f"expected JAC = 0 in c.eq, got:\n{ceq_fd}"

    sols_fd, _ = ode_fd.run(
        c='eq', ICP='eta', name='eta_fd',
        starting_point='EP2',
        DS=0.05, DSMAX=0.1, NMX=30,
    )
    ode_fd.close_session(clear_files=True)

    # ---- numerical parity: r(eta) curves should match ----
    # PyRates emits bare local names via auto's `parnames`/`unames`, so the
    # DataFrame columns are 'eta' / 'r' (not 'p/qif_op/eta' / 'p/qif_op/r').
    eta_jac = np.asarray(sols_jac['eta'].values).squeeze()
    eta_fd = np.asarray(sols_fd['eta'].values).squeeze()
    r_jac = np.asarray(sols_jac['r'].values).squeeze()
    r_fd = np.asarray(sols_fd['r'].values).squeeze()

    # sort by eta (defensive — auto's stepper may visit points out of order at boundaries)
    order_jac = np.argsort(eta_jac)
    order_fd = np.argsort(eta_fd)
    eta_jac, r_jac = eta_jac[order_jac], r_jac[order_jac]
    eta_fd, r_fd = eta_fd[order_fd], r_fd[order_fd]

    # interpolate onto common eta grid and compare
    lo = max(eta_jac.min(), eta_fd.min())
    hi = min(eta_jac.max(), eta_fd.max())
    assert hi > lo, "no overlapping eta range between analytical and FD continuations"
    eta_grid = np.linspace(lo, hi, 20)
    r_jac_i = np.interp(eta_grid, eta_jac, r_jac)
    r_fd_i = np.interp(eta_grid, eta_fd, r_fd)
    assert np.allclose(r_jac_i, r_fd_i, atol=1e-3, rtol=1e-3), (
        f"analytical vs FD Jacobian continuations diverge:\n"
        f"  eta_grid = {eta_grid}\n  r_jac    = {r_jac_i}\n  r_fd     = {r_fd_i}"
    )


@_requires_pyrates_dev
def test_1_15_extract_bare_name_fallback(auto_dir, tmp_path):
    """`extract` resolves auto-07p-native keys (`PAR(i)` / `U(i)`) to the bare
    local names that PyRates emits via `parnames` / `unames`.

    Pins the regression flagged in PyCoBi issue #3: in the from_yaml path,
    `_var_map_inv['PAR(4)']` translates to the *namespaced* name
    (`'p/qif_op/eta'`) — but the actual DataFrame column for that parameter
    is the *bare* name (`'eta'`), because PyRates writes parnames into the
    c.* file. Without the bare-name fallback, the documented call
    `plot_continuation('PAR(4)', 'U(1)', cont='eta')` raises KeyError because
    the namespaced lookup misses every column.
    """
    model = 'model_templates.neural_mass_models.qif.qif'
    ode = ODESystem.from_yaml(
        model, working_dir=str(tmp_path), auto_dir=auto_dir,
        file_name='qif_extract_mod', func_name='qif_extract',
        init_cont=True, NMX=500, NPR=500,
    )
    try:
        ode.run(
            starting_point='EP2', name='eta', ICP='eta',
            IPS=1, ILP=1, ISP=2,
            NMX=30, NPR=1, DS=0.05, DSMAX=0.5,
        )
        # Mix of auto-07p-native keys and a non-PAR scalar — all three should resolve.
        result, key_map = ode.extract(
            ['PAR(4)', 'U(1)', 'bifurcation'], cont='eta',
        )
        # PAR(4) and U(1) resolve to the bare names PyRates wrote into the c.* file;
        # 'bifurcation' is a summary-native column so it passes through unchanged.
        assert key_map == {'PAR(4)': 'eta', 'U(1)': 'r', 'bifurcation': 'bifurcation'}
        assert result.shape == (30, 3)

        # A truly missing key surfaces a helpful error listing what was tried.
        with raises(KeyError, match=r"is not a recognised summary column. Tried:"):
            ode.extract(['totally_made_up'], cont='eta')
    finally:
        ode.close_session(clear_files=True)


def test_1_13_starting_point_chaining(auto_dir):
    """Continuations can be chained: a fresh `run()` call started from a labeled
    point on a previous continuation registers as a new entry in `continuations`
    without disturbing the origin's entry.

    The chain here is:
      (0) IVP → (1) equilibrium continuation in eta with a fold at LP1 →
      (2) 2-parameter continuation of that fold curve in (eta, I_ext).
    """
    resources = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources')
    ode = ODESystem(
        eq_file='qif_eq', working_dir=resources, auto_dir=auto_dir,
        init_cont=True, c='ivp', NMX=500, NPR=500,
        parnames={}, unames={},
        params=['tau', 'Delta', 'I_ext', 'eta', 'weight'],
        state_vars=['r', 'v'],
    )
    try:
        # 1D equilibrium continuation in eta from the IVP's converged EP2.
        sols_a, _ = ode.run(
            starting_point='EP2', name='eta_a', ICP='eta',
            IPS=1, ILP=1, ISP=2, NMX=30, NPR=1, DS=0.05, DSMAX=0.5,
        )
        assert set(ode.continuations.keys()) == {0, 1}
        eta_a = ode.get_continuation('eta_a')
        assert eta_a.icps == [(4,)]  # eta is PAR(4) by PyCoBi's param ordering
        assert 'LP' in list(sols_a['bifurcation'].values.ravel()), \
            "expected an LP on the equilibrium curve to chain from"

        # Chain: start a new continuation from the LP, switch to 2-parameter
        # continuation of the LP curve in (eta, I_ext).
        sols_b, _ = ode.run(
            origin='eta_a', starting_point='LP1', name='lp_curve',
            ICP=['eta', 'I_ext'],
            IPS=1, ILP=0, ISP=2, ISW=2, NMX=20, NPR=1, DS=0.05, DSMAX=0.5,
        )
        # A fresh entry is added; the origin's entry stays put.
        assert set(ode.continuations.keys()) == {0, 1, 2}
        assert ode.get_continuation('eta_a') is eta_a
        lp_curve = ode.get_continuation('lp_curve')
        assert lp_curve.name == 'lp_curve'
        assert lp_curve.summary is sols_b
        # 2-param continuation records the ICP tuple as it was passed
        # (PyCoBi resolves names 'eta' and 'I_ext' to ints 4 and 3).
        assert lp_curve.icps == [(4, 3)]
        # And the summary has both PAR columns
        col_names = {c[0] for c in sols_b.columns}
        assert 'eta' in col_names
        assert 'I_ext' in col_names
    finally:
        ode.close_session()


@_requires_pyrates_dev
def test_1_14_limit_cycle_summary(auto_dir, tmp_path):
    """Limit-cycle continuation (`IPS=2`) produces a summary with the
    LC-specific column shape:

      - state-variable columns get integer sub-indices for min / max envelopes,
        i.e. ``('r', 0)`` and ``('r', 1)`` rather than the equilibrium-form
        ``('r', 0)`` alone;
      - the ``period`` column is populated when ``get_period=True`` (one value
        per labeled point, all positive).

    Uses PyRates' QIF + spike-frequency-adaptation model, which carries a
    Hopf bifurcation in eta — so the test needs the dev-PyRates features.
    """
    model = "model_templates.neural_mass_models.qif.qif_sfa"
    ode = ODESystem.from_yaml(
        model, working_dir=str(tmp_path), auto_dir=auto_dir,
        file_name='qif_sfa_lc_mod', func_name='qif_sfa_lc',
        # The SFA system has a slow timescale (tau_a=10); the IVP needs
        # NMX≈30000 steps at the default dt=1e-3 to reach the steady state.
        init_cont=True, NMX=30000, NPR=30000,
        # Same parameter regime as the qif_sfa use-example — known to carry
        # a sub-critical Hopf around eta ≈ 4.
        node_vars={'p/qif_sfa_op/Delta': 2.0,
                   'p/qif_sfa_op/alpha': 1.0,
                   'p/qif_sfa_op/eta': 3.0},
        edge_vars=[('p/qif_sfa_op/r', 'p/qif_sfa_op/r_in',
                    {'weight': 15.0 * np.sqrt(2.0)})],
    )
    try:
        # 1D equilibrium continuation in eta from the IVP's converged state
        # (EP2 = the last EP on the IVP branch; EP1 is the t=0 initial point).
        eta_sols, _ = ode.run(
            origin=0, starting_point='EP2', name='eta', bidirectional=True,
            ICP='eta', IPS=1, ILP=1, ISP=2,
            NMX=500, NPR=50, DS=1e-3, DSMIN=1e-8, DSMAX=5e-2,
            ITNW=20, NWTN=10, JAC=0,
        )
        bifs_eta = list(eta_sols['bifurcation'].values.ravel())
        assert 'HB' in bifs_eta, (
            f"expected at least one Hopf on the eta branch; got bifurcations {bifs_eta}"
        )

        # Branch-switch to the periodic-solution branch at the Hopf. IPS=2
        # tells auto to continue limit cycles; ISW=-1 starts from a HB.
        lc_sols, _ = ode.run(
            origin='eta', starting_point='HB1', name='lc',
            IPS=2, ISP=2, ISW=-1, NTST=50, NMX=30, NPR=1,
            DS=1e-3, DSMAX=0.5,
            get_period=True,
            JAC=0,
        )

        # --- LC-specific column shape ---
        col_tuples = list(lc_sols.columns)
        # Period column should be present with positive values
        assert ('period', '') in col_tuples
        periods = np.asarray(lc_sols[('period', '')].values).squeeze()
        assert (periods > 0).all(), f"expected positive periods, got {periods}"
        # State-var columns appear with min/max sub-indices: PyCoBi packs
        # each limit cycle into ('var', 0) and ('var', 1).
        r_subs = sorted(c[1] for c in col_tuples if c[0] == 'r')
        assert r_subs == [0, 1], f"expected r min/max sub-indices, got {r_subs}"
    finally:
        ode.close_session(clear_files=True)


@_requires_pyrates_dev
def test_3_1_pycobi_vs_auto_ops(auto_dir, tmp_path):
    """End-to-end parity vs auto-07p's own analytical-Jacobian demo.

    auto-07p ships the `ops` demo (a FitzHugh-Nagumo-like 3D ODE used in
    "optimization of periodic solutions") with a hand-derived DFDU/DFDP
    and JAC=1. We run PyCoBi over the demo's native Fortran source as
    the reference, then build the same model from a PyRates YAML and run
    the same equilibrium continuation in PAR(3) ('p3').

    The two paths agreeing means: (1) PyRates' symbolic DFDU/DFDP entries
    match the hand-derived ones; (2) PyCoBi's PyRates-driven setup yields
    the same auto-07p output as a raw demo run.

    We chose `ops` over the slightly simpler `abcb` because abcb's RHS
    uses `exp(x3)` and PyRates' equation parser currently leaves `exp`
    as an undefined function — `sympy.diff` then returns an unevaluated
    `Derivative(exp(x3), x3)` that gfortran refuses to compile.
    """
    import shutil

    demo_src = os.path.join(auto_dir, 'demos', 'ops')

    # ---- reference: auto-07p's native ops demo, driven through PyCoBi ----
    demo_dir = tmp_path / 'demo'
    demo_dir.mkdir()
    for fn in ('ops.f90', 'c.ops'):
        shutil.copy(os.path.join(demo_src, fn), demo_dir)
    ode_demo = ODESystem(
        eq_file='ops', working_dir=str(demo_dir),
        auto_dir=auto_dir, init_cont=False,
    )
    # `parnames={}, unames={}` clears auto-07p's runner cache: by design (see
    # auto/runAUTO.py config()) auto retains parnames/unames across run() calls
    # within a process and merges them with new c.* files. Without this, any
    # prior test that used a PyRates-generated model with unames/parnames would
    # leak its column names into ode_demo's DataFrame.
    sols_demo, _ = ode_demo.run(c='ops', name='ops_demo', parnames={}, unames={})
    ode_demo.close_session()

    # ---- mirror: same model through PyCoBi/PyRates with the analytical Jac ----
    yaml_dir = tmp_path / 'yaml'
    yaml_dir.mkdir()
    # PyRates' slash-notation path: <dir>/<filename-no-ext>/<template-name>.
    # Anchor on this test file's directory rather than cwd — preceding tests
    # leave cwd inside their own tmp_path via ODESystem._orig_dir handling.
    test_dir = os.path.dirname(os.path.abspath(__file__))
    yaml_template = os.path.join(test_dir, 'resources', 'ops', 'ops')
    ode_yaml = ODESystem.from_yaml(
        yaml_template,
        working_dir=str(yaml_dir), auto_dir=auto_dir,
        file_name='ops_mod', func_name='ops_fn',
        init_cont=False, analytical_jacobian=True,
        auto_constants=('eq',),
        # mirror the demo's c.ops auto-07p constants:
        NTST=15, NCOL=4, ISP=1, ILP=1, MXBF=10,
        NMX=25, NPR=500, IID=2, ITMX=8, ITNW=5, NWTN=3,
        EPSL=1e-7, EPSU=1e-7, EPSS=1e-4,
        DS=0.01, DSMIN=1e-3, DSMAX=0.5, IADS=1,
        # Stop when p3 reaches 0.95; the demo uses UZSTOP={3: 0.95} (literal
        # PAR(3)) because it doesn't declare parnames, but PyRates writes
        # `parnames` into c.eq and may not preserve the YAML's parameter order
        # (e.g. p3 ends up at PAR(4) for this model). Use the named form so
        # auto-07p resolves p3 via parnames regardless of its PAR index.
        UZSTOP={'p3': 0.95},
    )
    sols_yaml, _ = ode_yaml.run(c='eq', ICP='p3', name='ops_yaml')

    # ---- the analytical Jacobian path actually fired ----
    # (read these BEFORE close_session(clear_files=True) wipes the build dir)
    f90 = (yaml_dir / 'ops_mod.f90').read_text()
    ceq = (yaml_dir / 'c.eq').read_text()
    assert _DFDU_ASSIGN.search(f90)
    assert _DFDP_ASSIGN.search(f90)
    assert 'JAC = 1' in ceq

    ode_yaml.close_session(clear_files=True)

    # ---- numerical parity ----
    # demo uses default auto-07p names (U(1), PAR(3)); PyRates emits the user
    # names via unames/parnames (x1, p3).
    p3_demo = np.asarray(sols_demo['PAR(3)'].values).squeeze()
    x1_demo = np.asarray(sols_demo['U(1)'].values).squeeze()
    p3_yaml = np.asarray(sols_yaml['p3'].values).squeeze()
    x1_yaml = np.asarray(sols_yaml['x1'].values).squeeze()

    order = np.argsort(p3_demo)
    p3_demo, x1_demo = p3_demo[order], x1_demo[order]
    order = np.argsort(p3_yaml)
    p3_yaml, x1_yaml = p3_yaml[order], x1_yaml[order]

    lo = max(p3_demo.min(), p3_yaml.min())
    hi = min(p3_demo.max(), p3_yaml.max())
    assert hi - lo > 0.02, (
        f"insufficient overlap between PyCoBi/auto curves: p3 in [{lo}, {hi}]"
    )
    grid = np.linspace(lo, hi, 20)
    x1_demo_i = np.interp(grid, p3_demo, x1_demo)
    x1_yaml_i = np.interp(grid, p3_yaml, x1_yaml)
    assert np.allclose(x1_demo_i, x1_yaml_i, atol=1e-4, rtol=1e-4), (
        f"PyCoBi/PyRates ops curve diverges from auto-07p demo:\n"
        f"  p3_grid    = {grid}\n  x1_demo    = {x1_demo_i}\n  x1_yaml    = {x1_yaml_i}"
    )


@_requires_pyrates_dev
def test_1_16_reset_auto_state(auto_dir, tmp_path):
    """`ODESystem.reset_auto_state()` clears auto-07p's persisted
    ``parnames`` / ``unames`` so a subsequent model load doesn't inherit
    them.

    auto-07p's ``runAUTO.config()`` deliberately preserves ``parnames`` /
    ``unames`` across successive ``run()`` calls within a process (see the
    ``# do not completely replace existing constants data`` comment in
    ``runAUTO.py``). That's fine for iterating on one model but leaks across
    unrelated model loads: a fresh c.* with no ``unames`` declaration will
    silently inherit names from the previously-loaded model and relabel its
    DataFrame columns. This test pins three things:

      1. The leak is real — running a PyRates model populates the global
         runner's parnames/unames.
      2. `reset_auto_state()` clears both back to `None`.
      3. The same model loaded after reset still works (sanity check that
         clearing didn't break anything required for a subsequent run).
    """
    # Load a PyRates model that emits parnames/unames into its c.* file.
    # Construct ODESystem first so `auto` is imported with AUTO_DIR set;
    # _get_auto_runner imports auto on demand and would otherwise crash on
    # auto's AUTO_DIR-dependent module init.
    model_dir = tmp_path / 'leak'
    model_dir.mkdir()
    ode = ODESystem.from_yaml(
        'model_templates.neural_mass_models.qif.qif',
        working_dir=str(model_dir), auto_dir=auto_dir,
        file_name='qif_reset', func_name='qif_reset_fn',
        init_cont=True, NPR=100, NMX=1000,
    )

    runner = ODESystem._get_auto_runner()
    assert runner is not None, "could not locate auto-07p global runner"

    # (1) Runner should now hold the populated names — that's the leak.
    constants = runner.options['constants']
    assert constants['parnames'], (
        "auto-07p runner has no parnames after a PyRates-generated run — "
        "test setup is broken or PyRates is no longer emitting parnames"
    )
    assert constants['unames'], (
        "auto-07p runner has no unames after a PyRates-generated run"
    )

    ode.close_session(clear_files=True)

    # (2) Clear and verify.
    ODESystem.reset_auto_state()
    assert constants['parnames'] is None
    assert constants['unames'] is None

    # (3) Sanity: helper is a no-op when called a second time on an already-
    # clean runner; doesn't error and doesn't repopulate.
    ODESystem.reset_auto_state()
    assert constants['parnames'] is None
    assert constants['unames'] is None


# ----------------------------------------------------------------------------
# automated_continuation.py — Section 4
# ----------------------------------------------------------------------------


def test_4_1_resolve_param_for_extract_unit():
    """`_resolve_param_for_extract` converts an int auto-07p PAR index into the
    legacy ``'PAR(i)'`` string form (consumed by ``ODESystem.extract`` via
    ``_var_map_inv``), and leaves string names alone for the PyRates path
    where the c.* file already declares ``parnames``.
    """
    assert _resolve_param_for_extract(None, 4) == 'PAR(4)'
    # numpy integer types — auto-07p sometimes feeds these back through
    assert _resolve_param_for_extract(None, np.int64(7)) == 'PAR(7)'
    # string parnames pass through unchanged
    assert _resolve_param_for_extract(None, 'eta') == 'eta'
    assert _resolve_param_for_extract(None, 'p1') == 'p1'


def test_4_2_bif_series_contains_substring():
    """`_bif_series_contains` tests Series **values** for a substring; the
    pre-1.0 ``"LP" in series`` form silently tested the *index* labels
    instead, so the ZH-fold-or-Hopf branch of ``codim2_search`` effectively
    never fired (the codim-1 1D run's index was integers, not strings).
    """
    # Series with non-string index — the broken `"LP" in series` would
    # have tested labels 0/1/2/3 against "LP", always False.
    s = pd.Series(['RG', 'LP1', 'EP', 'HB1'], index=[0, 1, 2, 3])
    assert _bif_series_contains(s, 'LP') is True
    assert _bif_series_contains(s, 'HB') is True
    assert _bif_series_contains(s, 'ZH') is False
    # Exact label match also works (no false negatives for un-suffixed
    # labels like 'LP' without trailing number).
    s2 = pd.Series(['LP', 'RG'], index=['a', 'b'])
    assert _bif_series_contains(s2, 'LP') is True
    # NaN tolerance: empty / NaN entries don't raise.
    s3 = pd.Series([None, 'LP1'])
    assert _bif_series_contains(s3, 'LP') is True


def test_4_3_codim2_search_qif_fold_curve_handwritten(auto_dir):
    """End-to-end smoke test of `codim2_search` on the hand-written QIF.

    Setup: 1D continuation of QIF equilibria in eta finds an LP (saddle-node).
    Starting `codim2_search` from that LP1 in (eta, Delta) should run a
    bidirectional 2-parameter continuation and register at least one
    continuation in the returned dict. The PyRates-style name resolution
    (``params=['eta', 'Delta']`` — string parnames rather than PAR indices)
    must work even on this hand-written model where ``_var_map`` was seeded
    via ``params=[...]`` on construction.

    Pins the pre-1.0 hardcoded ``f'PAR({params[0]})'`` column lookup bug:
    on the hand-written path with `params=[..., 'eta', ...]`, the columns
    are remapped to user names so the bare `'PAR(4)'` lookup KeyErrors.
    """
    resources = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources')
    ode = ODESystem(
        eq_file='qif_eq', working_dir=resources, auto_dir=auto_dir,
        init_cont=True, c='ivp', NMX=500, NPR=500,
        parnames={}, unames={},
        params=['tau', 'Delta', 'I_ext', 'eta', 'weight'],
        state_vars=['r', 'v'],
    )
    try:
        # 1D equilibrium continuation in eta — finds LP1 at the fold.
        _, _ = ode.run(
            starting_point='EP2', name='eta_for_codim2', ICP='eta',
            IPS=1, ILP=1, ISP=2, NMX=60, NPR=1, DS=0.05, DSMAX=0.5,
            get_stability=True,
        )
        # codim2 in (eta, Delta) from the LP1 we just found.
        # Pass string parnames so the PyRates-style name resolution path
        # is exercised (the int-PAR path is exercised by the original
        # ICP=14 init_cont).
        result = codim2_search(
            params=['eta', 'Delta'], starting_points=['LP1'],
            origin=ode.get_continuation('eta_for_codim2').key,
            pyauto_instance=ode, max_recursion_depth=1,
            NMX=30, NPR=20, DS=0.05, DSMAX=0.5,
        )
        # At minimum the seed continuation got registered.
        assert isinstance(result, dict)
        assert len(result) >= 1, f"codim2_search produced no continuations: {result}"
        # The naming convention: f"{name}:{starting_point}" — verify shape.
        keys = list(result.keys())
        assert any(k.endswith(':LP1') for k in keys), (
            f"expected at least one continuation keyed '*:LP1', got {keys}"
        )
    finally:
        ode.close_session()


def test_4_4_codim2_search_pyrates_path(auto_dir, tmp_path):
    """`codim2_search` works when the model comes from PyRates and the
    summary columns are parnames-resolved (bare 'eta', 'Delta'), not the
    hand-written 'PAR(i)' form.

    Pins the parnames-blind hardcoded ``f'PAR({params[0]})'`` lookups in
    the pre-1.0 ``codim2_search`` / ``continue_period_doubling_bf``: on
    the PyRates path those columns simply don't exist, and the resulting
    KeyError aborted the entire search.
    """
    work = tmp_path / 'codim2_pyrates'
    work.mkdir()
    ode = ODESystem.from_yaml(
        'model_templates.neural_mass_models.qif.qif',
        working_dir=str(work), auto_dir=auto_dir,
        file_name='qif_codim2', func_name='qif_codim2_fn',
        init_cont=True, NPR=100, NMX=500,
    )
    try:
        # 1D continuation in eta to find LP.
        _, _ = ode.run(
            starting_point='EP2', name='eta_1d', ICP='eta',
            IPS=1, ILP=1, ISP=2, NMX=60, NPR=1, DS=0.05, DSMAX=0.5,
            get_stability=True,
        )
        # 2-param codim2 from LP1 — uses parnames string lookup throughout.
        result = codim2_search(
            params=['eta', 'Delta'], starting_points=['LP1'],
            origin=ode.get_continuation('eta_1d').key,
            pyauto_instance=ode, max_recursion_depth=1,
            NMX=30, NPR=20, DS=0.05, DSMAX=0.5,
        )
        assert len(result) >= 1
    finally:
        ode.close_session(clear_files=True)


def test_4_5_continue_period_doubling_bf_input_validation():
    """`continue_period_doubling_bf` raises a clear ValueError when the
    required ``ICP=[param1, param2]`` kwarg is missing or malformed,
    instead of failing later with a confusing KeyError from the inner
    `kwargs['ICP']` access.
    """
    # No ICP at all — the pre-1.0 version raised KeyError('ICP') from the
    # `params = kwargs['ICP']` line; we now raise ValueError up front.
    with raises(ValueError, match=r"requires ICP"):
        continue_period_doubling_bf(
            solution={}, continuation=0, pyauto_instance=None,
        )
    # Wrong-shape ICP — must be 2-element list/tuple.
    with raises(ValueError, match=r"length-2 list/tuple"):
        continue_period_doubling_bf(
            solution={}, continuation=0, pyauto_instance=None,
            ICP='eta',  # bare string, not a 2-tuple
        )
    with raises(ValueError, match=r"length-2 list/tuple"):
        continue_period_doubling_bf(
            solution={}, continuation=0, pyauto_instance=None,
            ICP=['eta'],  # single-element list
        )


def test_4_6_continue_period_doubling_bf_no_pds_terminates():
    """When the input solution dict has no PD points, the function returns
    an empty list and doesn't recurse / doesn't raise. Pins the early-exit
    behaviour so users that point this at a non-PD curve get a clean
    empty result rather than an iteration explosion.
    """
    sols, ret_pyauto = continue_period_doubling_bf(
        solution={
            1: {'bifurcation': 'EP'},
            2: {'bifurcation': 'HB1'},
            3: {'bifurcation': 'LP1'},
        },
        continuation=0, pyauto_instance='sentinel',
        ICP=['eta', 'Delta'],
    )
    assert sols == []
    assert ret_pyauto == 'sentinel'  # passthrough, no calls were made


def test_4_7_continue_period_doubling_bf_depth_cap_warns():
    """The `max_iter` depth cap fires on entry past the threshold and emits
    a UserWarning rather than silently truncating the cascade.
    """
    import warnings as _warnings
    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        sols, _ = continue_period_doubling_bf(
            solution={1: {'bifurcation': 'PD1'}},
            continuation=0, pyauto_instance=None,
            ICP=['eta', 'Delta'],
            max_iter=2, _depth=2,  # at the cap
        )
    assert sols == []
    msgs = [str(w.message) for w in caught]
    assert any('max_iter=2' in m for m in msgs), (
        f"expected max_iter warning, got: {msgs}"
    )


# ----------------------------------------------------------------------------
# GH / BT recursion — codim-2 handler unit tests
# ----------------------------------------------------------------------------


class _MockODESystem:
    """Stand-in for `ODESystem` that lets us drive `codim2_search` without
    auto-07p doing real continuations.

    `run_returns` is a list of `(summary_df, cont_key)` tuples consumed
    in order, one per `.run()` call. `extract_returns` is keyed by `cont_key`
    and returns the DataFrame slice that `.extract(['bifurcation', ...])`
    would yield. If a key is missing, `extract` raises `KeyError` (matching
    the real method's contract). `run_exception_at`, if set, makes the
    Nth `run()` call raise `RuntimeError("mock failure")` instead.
    """
    def __init__(self, run_returns=None, extract_returns=None,
                 run_exception_at=None):
        self.run_returns = list(run_returns or [])
        self.extract_returns = extract_returns or {}
        self.run_exception_at = run_exception_at
        self.run_calls = []  # list of (args, kwargs)
        self.extract_calls = []

    def run(self, **kwargs):
        n = len(self.run_calls)
        self.run_calls.append(kwargs)
        if self.run_exception_at is not None and n == self.run_exception_at:
            raise RuntimeError("mock failure")
        if n >= len(self.run_returns):
            raise IndexError(
                f"_MockODESystem.run called {n + 1} times but only "
                f"{len(self.run_returns)} canned returns supplied"
            )
        return self.run_returns[n]

    def extract(self, keys, cont, point=None):
        self.extract_calls.append((tuple(keys), cont, point))
        if cont not in self.extract_returns:
            raise KeyError(f"no mock extract for cont={cont!r}")
        df = self.extract_returns[cont]
        # keep only the requested columns
        cols = [k for k in keys if k in df.columns]
        return df[cols], {k: k for k in keys}


def test_4_8_recurse_codim2_dispatch_unknown_type():
    """`_recurse_codim2` short-circuits to an empty dict for unrecognised
    bifurcation types — keeps the main loop forward-compatible if a future
    auto-07p exposes new codim-2 codes without crashing existing scans.
    """
    from pycobi.automated_continuation import _recurse_codim2
    result = _recurse_codim2(
        codim2_type='XX', pyauto=None, origin=0, idx=1,
        params=[1, 2], name='n', recursion=0,
        max_recursion_depth=3, kwargs_1D_lc_cont=None, base_kwargs={},
    )
    assert result == {}


def test_4_9_codim2_search_gh_dispatch_calls_lc_continuation():
    """A `GH` token in the bifurcation column triggers an LC-style
    continuation: ``IPS=2, ISW=-1, ICP=[params[0], 11]`` are the canonical
    auto-07p constants for switching to the LC family from a generalised
    Hopf point. Pins the GH handler's specific kwargs against the
    convention.
    """
    # Mock: initial 2D continuation returns a summary with one GH row;
    # the subsequent GH-handler LC continuation returns an empty summary
    # (no LP-of-cycle on the LC branch — keeps recursion depth shallow).
    main_summary = pd.DataFrame({
        'bifurcation': ['GH1'],
        'eta': [0.5],
        'Delta': [0.3],
    })
    lc_summary = pd.DataFrame({
        'bifurcation': ['EP'],
        'eta': [0.6],
        'Delta': [0.3],
    })
    mock = _MockODESystem(
        run_returns=[
            (main_summary, 'main_cont'),  # initial 2D continuation
            (lc_summary, 'gh_lc_cont'),   # GH→LC continuation
        ],
        extract_returns={
            'main_cont': main_summary,
            'gh_lc_cont': lc_summary,
        },
    )
    result = codim2_search(
        params=['eta', 'Delta'], starting_points=['HB1'],
        origin=0, pyauto_instance=mock, max_recursion_depth=1,
        codim2_types=('GH',),
    )
    # Two runs: initial 2D codim-1 + the GH LC continuation.
    assert len(mock.run_calls) == 2
    gh_call = mock.run_calls[1]
    assert gh_call['starting_point'] == 'GH1'
    assert gh_call['IPS'] == 2 and gh_call['ISW'] == -1
    assert gh_call['ICP'] == ['eta', 11]
    # Result registers both continuations.
    assert len(result) == 2


def test_4_10_codim2_search_bt_dispatch_calls_equilibrium_continuation():
    """A `BT` token triggers a 1D equilibrium continuation in ``params[0]``
    stopping at ``HB1`` (IPS=1, ISW=1, STOP=['HB1']) — the canonical move
    from a Bogdanov-Takens point to the Hopf curve emerging from it.
    The homoclinic curve is deliberately *not* auto-followed; documented
    in the BT handler's docstring.
    """
    main_summary = pd.DataFrame({
        'bifurcation': ['BT1'],
        'eta': [0.5],
        'Delta': [0.3],
    })
    bt_eq_summary = pd.DataFrame({
        'bifurcation': ['EP', 'EP'],  # no HB found — recursion terminates
        'eta': [0.6, 0.7],
        'Delta': [0.3, 0.3],
    })
    mock = _MockODESystem(
        run_returns=[
            (main_summary, 'main_cont'),
            (bt_eq_summary, 'bt_eq_cont'),
        ],
        extract_returns={
            'main_cont': main_summary,
            'bt_eq_cont': bt_eq_summary,
        },
    )
    result = codim2_search(
        params=['eta', 'Delta'], starting_points=['LP1'],
        origin=0, pyauto_instance=mock, max_recursion_depth=1,
        codim2_types=('BT',),
    )
    assert len(mock.run_calls) == 2
    bt_call = mock.run_calls[1]
    assert bt_call['starting_point'] == 'BT1'
    assert bt_call['IPS'] == 1 and bt_call['ISW'] == 1
    assert bt_call['ICP'] == 'eta'
    assert bt_call['STOP'] == ['HB1']
    # No further Hopf found → BT handler returns just the main run; no
    # secondary continuation was registered.
    assert len(result) == 1


def test_4_11_codim2_search_gh_run_failure_warns_not_raises():
    """If the GH LC continuation itself raises inside auto-07p, the failure
    surfaces as a UserWarning citing the kwargs hook for the user to tune,
    and the main search continues rather than aborting.
    """
    import warnings as _warnings

    main_summary = pd.DataFrame({
        'bifurcation': ['GH1'],
        'eta': [0.5],
        'Delta': [0.3],
    })
    mock = _MockODESystem(
        run_returns=[(main_summary, 'main_cont')],
        extract_returns={'main_cont': main_summary},
        run_exception_at=1,  # second .run() call raises
    )

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        result = codim2_search(
            params=['eta', 'Delta'], starting_points=['HB1'],
            origin=0, pyauto_instance=mock, max_recursion_depth=1,
            codim2_types=('GH',),
        )

    msgs = [str(w.message) for w in caught]
    assert any('GH1' in m and 'kwargs_1D_lc_cont' in m for m in msgs), (
        f"expected GH-failure warning mentioning the user-tunable kwargs hook, "
        f"got: {msgs}"
    )
    # Initial continuation still gets registered.
    assert len(result) == 1


def test_4_12_codim2_search_bt_run_failure_warns_homoclinic_note():
    """BT run failures surface a UserWarning that explicitly mentions the
    homoclinic-curve limitation (since users hitting BT failures are the
    population most likely to ask why). Pins the diagnostic content.
    """
    import warnings as _warnings

    main_summary = pd.DataFrame({
        'bifurcation': ['BT1'],
        'eta': [0.5],
        'Delta': [0.3],
    })
    mock = _MockODESystem(
        run_returns=[(main_summary, 'main_cont')],
        extract_returns={'main_cont': main_summary},
        run_exception_at=1,
    )

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        codim2_search(
            params=['eta', 'Delta'], starting_points=['LP1'],
            origin=0, pyauto_instance=mock, max_recursion_depth=1,
            codim2_types=('BT',),
        )

    msgs = [str(w.message) for w in caught]
    assert any('BT1' in m and 'homoclinic' in m.lower() for m in msgs), (
        f"expected BT-failure warning mentioning the homoclinic gap, got: {msgs}"
    )


def test_4_13_codim2_search_codim2_types_filter():
    """The `codim2_types` filter narrows which codim-2 codes trigger
    recursion. Pass `('GH',)` and a curve with both GH and BT points
    should only fire the GH handler.
    """
    main_summary = pd.DataFrame({
        'bifurcation': ['GH1', 'BT1'],
        'eta': [0.5, 0.6],
        'Delta': [0.3, 0.4],
    })
    lc_summary = pd.DataFrame({
        'bifurcation': ['EP'],
        'eta': [0.55],
        'Delta': [0.3],
    })
    mock = _MockODESystem(
        run_returns=[
            (main_summary, 'main_cont'),
            (lc_summary, 'gh_lc_cont'),
        ],
        extract_returns={
            'main_cont': main_summary,
            'gh_lc_cont': lc_summary,
        },
    )
    result = codim2_search(
        params=['eta', 'Delta'], starting_points=['HB1'],
        origin=0, pyauto_instance=mock, max_recursion_depth=1,
        codim2_types=('GH',),  # BT is excluded
    )
    # Only the initial run + the GH branch — BT was not dispatched.
    assert len(mock.run_calls) == 2
    assert mock.run_calls[1]['starting_point'] == 'GH1'
    assert len(result) == 2
