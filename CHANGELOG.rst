Changelog
=========

1.0
---

1.0.1
~~~~~

Patch release fixing a silent name-resolution bug that affected multi-node
models with shared operators (the canonical Jansen-Rit case: three nodes
each carrying a ``v`` variable). PyRates disambiguates these collisions
with ``_v1`` / ``_v2`` / ``_v3`` suffixes in the c.* ``unames`` declaration
— so the four distinct ``v`` variables are stored as ``v``, ``v_v1``,
``v_v2``, ``v_v3`` in the DataFrame columns. Pre-1.0.1 the
``_resolve_summary_key`` fallback chain didn't know about PyRates'
uname/parname → column mapping; namespaced inputs like
``'pc/rpo_i/v'`` silently collapsed to the bare ``'v'`` column —
the same column for every node — without any error to flag the aliasing.

- :code:`_resolve_summary_key` now reads the active ``unames`` / ``parnames``
  list off auto-07p's process-global runner and uses it as a bridge between
  ``_var_map`` (namespaced name → auto-07p slot index) and the DataFrame
  columns (PyRates-emitted, possibly suffix-disambiguated). Each of the
  four Jansen-Rit ``v`` variables now resolves to its own distinct column.
  Falls through cleanly when the runner isn't yet populated (single-node
  models still work via the existing bare-strip fallback).
- :code:`documentation/use_examples/jansenrit_rg_updated.py` updated to
  use :code:`ODESystem.extract` for namespaced reads of DataFrame columns
  (direct ``t_sols["pc/rpo_e_in/v"]`` indexing bypasses the resolver and
  raises ``KeyError`` for namespaced keys on the post-PyRates-parnames
  schema). Inline comments document the underlying naming flow.

Pinned by :code:`test_5_13_resolve_summary_key_bridges_via_uname_parname`
(multi-node JRC layout with a stubbed auto-07p runner) and
:code:`test_5_14_resolve_summary_key_falls_through_when_no_runner`
(single-node fallback path).

1.0.0
~~~~~

First production-stable release. Consolidates the 0.10.0 maintenance /
API-tightening work with a deeper hardening pass on the automated
continuation helpers and visualization API, plus two new use examples
on top of a freshly built bi-exponential QIF-SFA model.

Headline changes:

- PyCoBi now surfaces PyRates' symbolic-Jacobian generation, multi-
  scenario ``c.*`` files, and unames/parnames so DFDU/DFDP entries land
  in the generated Fortran and auto-07p uses the analytical Jacobian by
  default.
- :func:`codim2_search` works on the PyRates-driven flow (parnames-named
  summary columns), reliably handles ZH / GH / BT codim-2 points, and
  surfaces auto-07p sub-run failures as :class:`UserWarning` rather than
  aborting the whole search.
- :meth:`ODESystem.plot_continuation_grid` lays out multiple
  continuations as a faceted grid in one call;
  :meth:`plot_trajectory(colorbar=True)` adds a colorbar for 3D phase
  trajectories; the default ``line_color_unstable`` is now ``'gray'`` so
  stable / unstable segments are visually distinguishable.

Behaviour changes (be aware on upgrade)
'''''''''''''''''''''''''''''''''''''''

- :code:`init_cont` now defaults to :code:`False` on :code:`ODESystem.__init__`,
  :code:`ODESystem.from_yaml`, and :code:`ODESystem.from_template`. Set
  :code:`init_cont=True` explicitly to opt into the legacy automatic-IVP
  behaviour.
- Default :code:`line_color_unstable` changed from :code:`'k'` to
  :code:`'gray'` (identical to stable previously; only the linestyle
  distinguished them). Pass :code:`line_color_unstable='k'` through any
  plot method to revert. The line-style distinction (``'solid'`` vs
  ``'dotted'``) is unchanged.
- :code:`continue_period_doubling_bf`'s :code:`max_iter` now bounds the
  *recursion depth* of the PD cascade — was previously a per-call
  iteration counter that only fired once on entry, which made the cap
  largely cosmetic on real cascades. Argument name preserved; semantics
  changed.

New features
''''''''''''

- :code:`analytical_jacobian: bool = True` named parameter on
  :code:`from_yaml` / :code:`from_template` forwards as PyRates' :code:`auto_jac`;
  the generated c.* file sets :code:`JAC=1` and DFDU/DFDP blocks are emitted in
  :code:`func`. Override per-call via :code:`JAC=0`/:code:`JAC=1` on :code:`run()`.
- :code:`auto_constants` (string or iterable, default :code:`('ivp',)`) on
  :code:`from_yaml` / :code:`from_template` exposes PyRates' multi-scenario c.*
  generation — request :code:`('ivp', 'eq', 'lc', 'bvp')` to set up all four
  scenarios at once and switch between them at run-time via :code:`run(c=...)`.
- :code:`Continuation` dataclass and :code:`ODESystem.get_continuation(key_or_name)`
  provide a single typed view over per-continuation state. The four legacy mirror
  dicts (:code:`auto_solutions`, :code:`results`, :code:`_results_map`,
  :code:`_branches`) keep working unchanged for backward-compatible reads.
- :code:`parse_point_diagnostics(s, diag=None)` exposes the per-point diagnostic
  parser (one-shot regex; returns :code:`{stable, eigenvalues, text}`).
- :code:`_map_auto_kwargs` now remaps named PAR keys to integers for
  :code:`UZSTOP`, :code:`THL`, and :code:`THU` in addition to :code:`ICP` and
  :code:`UZR`.
- :code:`ODESystem.reset_auto_state()` — staticmethod that clears auto-07p's
  per-process retained ``parnames`` / ``unames`` so a subsequent model load
  doesn't inherit column names from the previous one.
- :code:`codim2_search` recurses on ZH, GH (Bautin / generalised-Hopf), and
  BT (Bogdanov-Takens) codim-2 points (was: ZH only). Per-type recursion
  follows standard auto-07p conventions for the switch-and-continue (LC
  family from GH, equilibrium→Hopf from BT). Homoclinic continuation from
  BT is deliberately left as a documented manual step (requires HomCont /
  IPS=9). Recognised codim-2 types narrowable via :code:`codim2_types=`.
  New :code:`kwargs_1D_lc_cont=` hook tunes the GH LC continuation.
- :code:`ODESystem.plot_continuation_grid(plots, ncols=2, ...)` —
  faceted multi-subplot helper. Each plot spec is a dict carrying
  ``x``, ``y``, ``cont``, optional ``title`` / ``panel_label``, plus
  per-plot overrides for any :meth:`plot_continuation` keyword.
- :code:`ODESystem.plot_trajectory(colorbar=True, colorbar_label=...)` —
  for 3D phase trajectories, attaches a figure colorbar surfacing the
  coordinate that the Line3DCollection's colormap encodes.
- Two new gallery examples co-located with the new
  :code:`qif_biexp_sfa.yaml`: :code:`automated_continuation.py` (codim-2
  BT/GH scenario in the (eta, tau_r) plane plus a recipe for the PD
  cascade at small tau_r) and :code:`visualization_tools.py`
  (:code:`plot_continuation_grid` + :code:`plot_trajectory(colorbar=)`).

Bug fixes
'''''''''

- :code:`to_file` opened the pickle in text mode (:code:`'x'`) — first call
  always raised :code:`TypeError` and silently fell into the binary fallback.
  Now uses :code:`'xb'` / :code:`'wb'` via context managers, and
  :code:`from_file` reconstructs the new :code:`continuations` store from the
  mirror dicts on load (handles non-dict slots too).
- :code:`continue_period_doubling_bf` had a mutable default :code:`pds: list = []`
  that was shared across calls — replaced with :code:`pds=None` plus in-body
  initialisation.
- :code:`codim2_search` and :code:`continue_period_doubling_bf` iterated
  :code:`for bf, p1, p2 in <DataFrame>:` which actually iterates *column labels*
  — switched to :code:`.itertuples(index=False, name=None)` so the automated
  codim-2 / PD-cascade machinery actually runs over rows.
- :code:`codim2_search` and :code:`continue_period_doubling_bf` no longer
  hardcode :code:`f'PAR({params[0]})'` column lookups (they failed on the
  PyRates flow where columns are named via parnames, and on the hand-
  written flow when :code:`params=` was passed to ODESystem). Both now
  route column access through :meth:`ODESystem.extract` and accept int
  PAR indices or string names.
- :code:`codim2_search`'s ZH→fold/Hopf branch used :code:`"LP" in
  codim1_bifs` on a pandas Series — which silently tests *index labels*
  rather than values, so the recursion effectively never fired. Replaced
  with a substring match against the values.
- :code:`_create_summary` was appending the whole eigenvalue list per
  eigenvalue column rather than the scalar value; now stores scalars.
- :code:`Line3DCollection(segments=lines, ...)` was outright broken on current
  matplotlib (the 3D constructor takes :code:`lines` positionally). Fixed.
- :code:`_get_line_collection` and :code:`_get_3d_line_collection` now pop
  :code:`linestyles` from :code:`**kwargs` so user-supplied overrides don't
  collide with the per-segment styles the function computes. The size-1
  stability tolerance was tightened (uses :code:`.size > 1` rather than the
  cryptic :code:`np.sum(shape) > 1`). The :code:`dtype='int'` coerce is now
  in a try/except so all-None stability arrays (typical for IVP results)
  fall through to no-stability rendering rather than raising
  :code:`TypeError`.
- :code:`plot_continuation` no longer calls :code:`ax.legend()` when no
  labelled artists exist (was: matplotlib "No artists with labels"
  UserWarning per continuation that crossed zero special points).
- :code:`plot_timeseries(points=None)` now plots one trace per labelled
  point on the continuation — was previously silently dropping every
  trace past the first via :code:`range(len(['RG']))`. Empty
  continuations raise a clear ValueError pointing at NPR / explicit
  :code:`points=`.
- :code:`plot_bifurcation_points` no longer routes through
  :code:`plt.sca(ax) + plt.plot(...)` — uses :code:`ax.plot` directly so
  faceted layouts don't silently misroute to the current axes.
- :code:`plot_trajectory(cutoff=...)` no longer crashes with
  "key of type tuple not found and not a MultiIndex"
  (:code:`np.where` returns a 1-tuple of arrays; unpacked) and no longer
  fills the trimmed DataFrame with NaN via per-column slice-assign
  (now slices the whole container via :code:`.iloc[idx]`). Stability
  Series coerce in :code:`_get_line_collection` also tolerates label-
  reindexed inputs (uses :code:`.values` upfront).
- :code:`_resolve_summary_key` now strips the namespace prefix off the
  user's key directly (step 4 in its docstring) — covers the PyRates
  flow where the c.* file declares bare parnames but the user
  references the namespaced :code:`'p/qif_op/eta'` they declared the
  model with.
- :code:`_start_from_solution`'s silent second :code:`auto.run` retry now
  emits a :code:`UserWarning` so callers know the original :code:`run()`
  kwargs were dropped on the retry.

Internal refactors
''''''''''''''''''

- Centralised continuation bookkeeping in :code:`_register_continuation` and
  :code:`_record_summary`; the scattered writes that used to live in
  :code:`run()` / :code:`merge()` now happen in one place.
- Bidirectional continuation no longer uses the magic-string sentinel
  :code:`name == 'bidirect:cont2'` for re-entry detection — replaced with a
  private :code:`_reverse_direction` kwarg.
- Replaced the regex-fragile token-search in :code:`get_solution_eigenvalues`
  with a one-shot regex (:code:`parse_point_diagnostics`). The same parsed
  dict is reused for stability / eigenvalues / lyapunov within one summary
  build (was up to 3 redundant text walks per point).
- Reworked :code:`get_branch_info`'s magic-10 retry loop into explicit
  iteration with a clear error message naming every label index that was
  tried.
- :code:`_create_summary` no longer builds two parallel
  :code:`data_1d` / :code:`data_2d` structures and then merges them — a
  single :code:`dict[(name, sub_index), list]` produces the final
  :code:`MultiIndex` DataFrame in one shot. :code:`_to_dataframe` is removed.
- :code:`__getitem__` raises with a helpful :code:`KeyError` listing every
  registered name and stored key on a double-miss.
- :code:`ODESystem.blocked_indices` is now imported from PyRates'
  :code:`FortranBackend._AUTO_BLOCKED_PAR_RANGE` with a hard-coded
  :code:`(10, 15)` fallback — single source of truth.
- :code:`from_template` filters non-parameter args
  (:code:`{'t', 'y', 'dy', 'hist'}`) by name rather than slicing by position
  — robust against DDE models (where :code:`hist` shifts the offset).
- :code:`continue_period_doubling_bf` rewritten with explicit recursion-
  depth tracking (the old :code:`iteration` accumulator never accumulated
  across sibling branches). PD points dedupe via tuples of rounded
  coordinates; names follow a :code:`pd_d{depth}_n{i}` convention so
  parallel sub-cascades don't collide on shared labels.

0.9
---

0.9.1
~~~~~

- added python package `meson` to requirements to account for issue #3 (https://github.com/pyrates-neuroscience/PyCoBi/issues/3)
- updated the installation instructions in `README.md` and readthedocs, to account for issue #3 (https://github.com/pyrates-neuroscience/PyCoBi/issues/3)
- added a new use example for 1D and 2D parameter continuations of the IK neuron model
- updated the simple use example in the documentation

0.9.0
~~~~~

- cumulative update based on the user interface updates, complex variable support, and bug fixes since the introduction of `0.8.0`
- added a use example for the bifurcation analysis of the Jansen Rit model
- fixed a bug with `ODESystem.plot_continuation`: for `get_stability=False`, a `KeyError` was sometimes triggered by an incomplete mapping from the Auto-07p parameter names to the PyRates parameter names

0.8
---

0.8.10
~~~~~

0.8.10
~~~~~

- added requirement for numpy >= 2.0
- dropped support for python versions 3.8 and older due to new numpy requirements

0.8.9
~~~~~

- added the option to not extract stability information for plotting solutions to `ODESystem.plot_continuation`
- updated version dependency on numpy (enforced by issues between the newest version of numpy and Auto-07p)

0.8.8
~~~~~

- fixed bug where the plotting method `ODESystem.plot_continuation` would throw an error when `PAR(14)` was used as the time parameter for the x-axis
- fixed bug where the merging of branches caused the wrong diagnostic files to be pulled up for determining the solution stability / eigenvalue spectrum
- added a new utility function `get_point_diagnostics` in utility.py that extracts the correct diagnostics file from the Auto-07p output for a given solution
- changed `get_solution_stability` and `get_solution_eigenvalues` to use the new function in utility.py

0.8.7
~~~~~

- implemented support for complex variable types. Pass the keyword argument "float_precision=complex" to the `ODESystem.from_yaml` or `.from_template` commands in order to generate complex-valued model files.

0.8.6
~~~~~

- updated how summaries of solution branches are generated to save computation time (results are now merged first, such that a summary is created only once)
- fixed a minor bug with the extraction of solutions via their keys
- improved how the continuation direction for automated bidirectional continuations is handled

0.8.5
~~~~~

- added a method `ODESystem.from_template` that allows to instantiate `ODESystem` from a `pyrates.CircuitTemplate`
- method `ODESystem.from_yaml` now creates a `pyrates.CircuitTemplate` from a YAML file first, and then passes it to `ODESystem.from_template`
- updated CircleCI config: Dropped support for Python 3.6, added support for Python 3.10

0.8.4
~~~~~

- updated readthedocs configuration
- updates to `ODESystem.extract` and the plotting functions that use it: Users can now flexibly switch between using the parameter/variable indices when specifying the variables to plot, or whether they want to use the pyrates-style naming scheme
- added option to `ODESystem.run` to store only the minimum and maximum of each selected state variable for periodic solutions (set keyword argument `reduce_limit_cycle=True`)
- debugged scenario where initial condition was run, but not explicitly labeled, which caused issues with extracting starting points from the solution branch
- debugged issue with selecting specific variables and parameters that should be saved in the continuation summary
- added private method `ODESystem._map_var` that ensures stable mapping between user-provided variable names and intrinsic variable names

0.8.3
~~~~~

- Bug fix for the new parameter naming system: `ODESystem` now accounts for the blocked Auto-07p parameter vector entries

0.8.2
~~~~~

- updated plotting method `ODESystem.plot_continuation`: It can now automatically plot a legend for all bifurcation points in a bifurcation diagram
- Keyword argument added to `ODESystem.plot_continuation`: "bifurcation_legend" can be set to `True` or `False` to turn bifurcation type legends on/off
- Bifurcation markers are now plotted via the `matplotlib.pyplot.plot` function rather than the `matplotlib.pyplot.scatter` function

0.8.1
~~~~~

- updates of the readthedocs requirements
- minor bugfix in "pycobi.py", where the new variable/parameter naming changes from 0.8.0 were interfering with the old naming system

0.8.0
~~~~~

- added mapping functionalities that allow to use the pyrates-based names for variables and parameters in the model rather than the auto naming style
- implemented the changes with the `ODESystem.run` and `ODESystem.extract` methods
- dataframes return by `ODESystem.run` now contain the pyrates-like variable names in the column header
- `ODESystem.__init__` now takes a couple of new arguments: The "eq_file" is a positional argument that ties an `ODESystem` instance to a single fortran equation file. "params" and "state_vars" allow to provide the parameter and state variable names that can be used instead of the indexing style of auto.
- the parameter and state-variable names that can be provided to the `ODESystem.__init__` method can be directly obtained from the `CircuitTemplate.get_run_func` that is also used by the `ODESystem.from_yaml` method.

0.7
---

0.7.5
~~~~~

- added a bugfix that allows to merge solution branches that included multiple branches that might arise from a automated switching at a branch point
- fixed a bug with the readthedocs website not displaying math correctly
- updated version dependency on numpy (enforced by issues between the newest version of numpy and Auto-07p)

0.7.4
~~~~~

- debugged the QIF-SFA use example
- added zenodo doi to the readme and documentation
- removed bug with saving additional attributes on the ODESystem instance

0.7.3
~~~~~

- improved docstrings of all public methods of `ODESystem`
- moved most static methods of `ODESystem` to a separate `utility` package
- added API section to readthedocs documentation
- moved period doubling continuation and automated 2D bifurcation analysis to extra package `automated_continuation`

0.7.2
~~~~~

- added use example for the QIF-SFA model to the documentation
- improved support for `pandas.DataFrames` as the main results storage data type
- added the pyrates model template as an attribute to the `ODESystem`
- added the option of clearing all generated fortran/auto files via the `ODESystem.close_session()` method

0.7.1
~~~~~

- debugged circle CI config
- added readthedocs source files
- improved integration of pycobi and pandas

0.7.0
~~~~~

- summaries of parameter continuations are now stored and returned as `pandas.DataFrame` instances
- added `ODESystem.__getitem__` method that allows to directly access parameter continuation summaries via their keys
- added `ODESystem.close_session` method that changes the working directory of the system to the directory prior to initialization of `ODESystem`
- added possibility to call pytest scripts with auto-07p directory as additional command line argument (--auto_dir)
- updated gitignore
- added new tests

0.6
---

0.6.3
~~~~~

- added new badges to the readme
- added official support for Python 3.9
- added CircleCI config
- removed bug from the `ODESystem.from_yaml` method, where the run function was generated for a discrete time-step solver by default

0.6.2
~~~~~

debugged PyPI integation

0.6.1
~~~~~

changed package name from `PyAuto` to `PyCoBi`

0.6.0
~~~~~

- first official version
- wrapper to `Auto-07p`
- automated fortran file generation via `PyRates`
- simplified handling of auto environment variables
- simplified parameter continuation
- visualization functions
- save and load results of parameter continuation/bifurcation analysis
