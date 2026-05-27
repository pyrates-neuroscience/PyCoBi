"""Automated multi-step continuations and codimension-2 searches.

These helpers wrap :class:`ODESystem.run` with bookkeeping for two common
patterns:

* ``continue_period_doubling_bf``: chase a period-doubling cascade in two
  parameters, recursing on every new PD point encountered.
* ``codim2_search``: 2-parameter continuation of a codim-1 starting point
  with optional recursion into codim-2 points (ZH / GH / BT).

Stability note
--------------

These are convenience wrappers around the lower-level
:class:`ODESystem.run` API. Recursive codim-2 searches depend on the
underlying auto-07p run producing the expected branch structure; on
unfamiliar models you should validate the results against a manual run
before trusting them. The implementations follow standard auto-07p
conventions for switching between equilibrium / limit-cycle / fold
continuations, but they cannot guarantee numerical stability across all
model families. Sub-runs that fail with an exception inside
``ODESystem.run`` will be reported via :class:`UserWarning` rather than
abort the whole search.
"""

import warnings
from typing import Any, Union

import numpy as np

from .pycobi import ODESystem


def _resolve_param_for_extract(pyauto: ODESystem, param) -> str:
    """Convert a param identifier into a form ``ODESystem.extract`` can resolve.

    Accepts either:
      * an int auto-07p PAR index (e.g. ``4`` → ``'PAR(4)'``), or
      * a string parameter name (passed through unchanged).

    The string form covers the PyRates-driven flow where the generated
    ``c.*`` file declares ``parnames`` and the summary DataFrame columns
    carry the bare user names (``'eta'``, ``'p1'``, ...). The int form
    covers the hand-written Fortran flow where columns come out as
    ``'PAR(i)'`` and optionally get remapped to user names via
    ``ODESystem(params=[...])``.
    """
    if isinstance(param, (int, np.integer)):
        return f'PAR({int(param)})'
    return param


def _bif_series_contains(series, code: str) -> bool:
    """Return True iff any string in ``series`` contains ``code``.

    Replaces ``code in series`` (which silently tests the Series *index*,
    not the values) and ``code in series.values`` (which is exact equality,
    missing labels like ``'LP1'`` when probing for ``'LP'``).
    """
    return bool(series.astype(str).str.contains(code, na=False).any())


def continue_period_doubling_bf(solution: dict, continuation: Union[str, int, Any],
                                pyauto_instance: ODESystem,
                                max_iter: int = 1000, precision: int = 3,
                                _depth: int = 0, _pds: list = None,
                                **kwargs) -> tuple:
    """Recursively follow a period-doubling cascade in 2 parameters.

    For every ``PD`` (period-doubling) label on ``solution``, start a new
    continuation from it (forwarded to ``pyauto_instance.run`` with
    ``**kwargs``), inspect the resulting curve for further PD points, and
    recurse — bounded by ``max_iter`` recursion depth to prevent runaway
    cascades on chaotic models.

    Parameters
    ----------
    solution
        The previous continuation's solution dictionary (mapping point key
        to ``{'bifurcation': ..., ...}``).
    continuation
        Origin key of the parent continuation, forwarded to
        ``run(origin=)``.
    pyauto_instance
        The active ``ODESystem`` instance.
    max_iter
        Maximum recursion depth. Default 1000. Semantics changed in 1.0.0
        from per-call iteration count to recursion depth, which actually
        bounds runaway cascades.
    precision
        Decimal places used to dedupe already-visited PD points (so the
        same PD point at slightly different floating-point coordinates
        across sibling branches doesn't trigger infinite recursion).
    kwargs
        Remaining keyword arguments forwarded to ``pyauto_instance.run``.
        Must include ``ICP=[param1, param2]`` (two continuation parameters).
        Each entry may be an int auto-07p PAR index or a string name (the
        latter resolved via the active ``c.*`` ``parnames`` declaration).

    Returns
    -------
    tuple
        ``(continuation_names, pyauto_instance)``. ``continuation_names``
        lists the names registered on ``pyauto_instance`` for each PD
        continuation that was performed, in execution order, including
        recursive sub-cascades.
    """
    if _pds is None:
        _pds = []
    if _depth >= max_iter:
        warnings.warn(
            f"continue_period_doubling_bf: reached max_iter={max_iter} recursion "
            f"depth; stopping cascade. Increase max_iter if deeper cascades are "
            f"expected.",
            UserWarning, stacklevel=2,
        )
        return [], pyauto_instance

    if 'ICP' not in kwargs:
        raise ValueError(
            "continue_period_doubling_bf requires ICP=[param1, param2] in **kwargs"
        )
    params = kwargs['ICP']
    if not (isinstance(params, (list, tuple)) and len(params) == 2):
        raise ValueError(
            f"continue_period_doubling_bf requires ICP to be a length-2 list/tuple; "
            f"got {params!r}"
        )
    param_cols = [_resolve_param_for_extract(pyauto_instance, p) for p in params]

    solutions: list = []
    pd_count = 0  # auto-07p uses 1-indexed PD labels

    for _, point_info in solution.items():
        if 'PD' not in point_info.get('bifurcation', ''):
            continue
        pd_count += 1
        name = f'pd_d{_depth}_n{pd_count}'

        try:
            s_tmp, cont = pyauto_instance.run(
                starting_point=f'PD{pd_count}', name=name,
                origin=continuation, **kwargs,
            )
        except Exception as exc:
            warnings.warn(
                f"continue_period_doubling_bf: PD{pd_count} continuation failed "
                f"({type(exc).__name__}: {exc}); skipping this branch.",
                UserWarning, stacklevel=2,
            )
            continue
        solutions.append(name)

        try:
            bfs, _ = pyauto_instance.extract(['bifurcation'] + param_cols, cont=cont)
        except KeyError as exc:
            warnings.warn(
                f"continue_period_doubling_bf: could not extract bifurcation "
                f"summary for {name!r} ({exc}); skipping recursion on this branch.",
                UserWarning, stacklevel=2,
            )
            continue

        for bf, p1, p2 in bfs.itertuples(index=False, name=None):
            param_pos = (round(float(p1), precision), round(float(p2), precision))
            if 'PD' in str(bf) and param_pos not in _pds:
                _pds.append(param_pos)
                sub_sols, _ = continue_period_doubling_bf(
                    solution=s_tmp, continuation=cont,
                    pyauto_instance=pyauto_instance,
                    max_iter=max_iter, precision=precision,
                    _depth=_depth + 1, _pds=_pds, **kwargs,
                )
                solutions.extend(sub_sols)

    return solutions, pyauto_instance


def codim2_search(params: list, starting_points: list,
                  origin: Union[str, int, Any], pyauto_instance: ODESystem,
                  max_recursion_depth: int = 3, recursion: int = 0,
                  periodic: bool = False,
                  kwargs_2D_lc_cont: dict = None, kwargs_2D_cont: dict = None,
                  kwargs_1D_lc_cont: dict = None,
                  precision: int = 2,
                  codim2_types: tuple = ('ZH', 'GH', 'BT'),
                  **kwargs) -> dict:
    """Continue codim-1 bifurcation points in 2 parameters with recursive
    handling of codim-2 points along each curve.

    For every entry in ``starting_points`` (typically auto-07p labels like
    ``'LP1'``, ``'HB1'``, ``'PD1'``), run a 2-parameter continuation in
    ``params``. Then walk the resulting curve for codim-2 points and, for
    each new one, perform the appropriate switch-and-continue:

    * ``ZH`` (zero-Hopf): 1D equilibrium continuation in ``params[0]``
      stopping at the nearest fold or Hopf, then recursive 2D continuation
      of that codim-1 point.
    * ``GH`` (generalized Hopf / Bautin): switch to LC continuation in
      ``params[0]`` and ``PAR(11)`` (period), looking for a fold-of-cycles
      (LP) that, if found, is recursively continued in 2 parameters via
      the LC route.
    * ``BT`` (Bogdanov-Takens): 1D equilibrium continuation from BT
      stopping at the nearest Hopf, then 2D continuation of the resulting
      Hopf curve. *Does not* follow the homoclinic curve emerging from BT
      (that requires HomCont params with IPS=9; see auto-07p docs).

    The list of recognised codim-2 types can be narrowed via ``codim2_types``
    — e.g. pass ``('ZH',)`` to suppress GH/BT recursion entirely.

    Failures during sub-continuations are surfaced as ``UserWarning`` rather
    than propagated, on the assumption that a search across many points
    should tolerate individual misfires.

    Parameters
    ----------
    params
        Two-element list of continuation parameter identifiers (int PAR
        indices or string names; see ``_resolve_param_for_extract``).
    starting_points
        Codim-1 auto-07p labels (e.g. ``['LP1', 'HB1']``) to start from.
    origin
        Origin continuation key on which the starting points live.
    pyauto_instance
        The active ``ODESystem``.
    max_recursion_depth
        Maximum codim-2 recursion depth. Default 3.
    recursion
        Internal recursion counter; do not set directly.
    periodic
        If true, the initial continuation is run with ``IPS=2, ISW=2``
        (limit-cycle 2-parameter continuation); use this when the starting
        point sits on a limit-cycle branch (e.g. a PD or fold-of-cycle).
    kwargs_2D_lc_cont
        Extra kwargs merged into the 2D LC continuation (``periodic=True``).
    kwargs_2D_cont
        Extra kwargs merged into the 2D equilibrium continuation
        (``periodic=False``).
    kwargs_1D_lc_cont
        Extra kwargs merged into the 1D LC continuation triggered by a GH
        point.
    precision
        Decimal places used to dedupe codim-2 points along the curve.
    codim2_types
        Tuple of codim-2 bifurcation types to recurse on. Default
        ``('ZH', 'GH', 'BT')``.
    kwargs
        Remaining keyword arguments forwarded to ``pyauto_instance.run``.

    Returns
    -------
    dict
        Mapping of continuation name to ``ODESystem`` continuation key for
        every run performed (the initial continuations plus every
        successful recursive sub-run).
    """
    found_codim2: dict = {bf: [] for bf in codim2_types}
    continuations: dict = dict()
    name = kwargs.pop('name', f"{params[0]}/{params[1]}")
    param_cols = [_resolve_param_for_extract(pyauto_instance, p) for p in params]

    for p in starting_points:

        # initial 2-parameter codim-1 continuation
        kwargs_tmp = dict(kwargs)
        if periodic:
            kwargs_tmp.update({'ILP': 0, 'IPS': 2, 'ISW': 2, 'ISP': 2,
                               'ICP': list(params) + [11]})
            if kwargs_2D_lc_cont:
                kwargs_tmp.update(kwargs_2D_lc_cont)
        else:
            kwargs_tmp.update({'ILP': 0, 'IPS': 1, 'ISW': 2, 'ISP': 2,
                               'ICP': list(params)})
            if kwargs_2D_cont:
                kwargs_tmp.update(kwargs_2D_cont)

        name_tmp = f"{name}:{p}"
        try:
            sols, cont = pyauto_instance.run(
                starting_point=p, origin=origin, name=name_tmp,
                bidirectional=True, **kwargs_tmp,
            )
        except Exception as exc:
            warnings.warn(
                f"codim2_search: initial continuation from {p!r} failed "
                f"({type(exc).__name__}: {exc}); skipping this starting point.",
                UserWarning, stacklevel=2,
            )
            continue
        continuations[name_tmp] = cont

        if recursion >= max_recursion_depth:
            continue

        # walk the curve for codim-2 points
        try:
            curve, _ = pyauto_instance.extract(
                ['bifurcation'] + param_cols, cont=cont,
            )
        except KeyError as exc:
            warnings.warn(
                f"codim2_search: could not extract bifurcation summary from "
                f"{name_tmp!r} ({exc}); skipping codim-2 recursion on this branch.",
                UserWarning, stacklevel=2,
            )
            continue

        for bf, p1, p2 in curve.itertuples(index=False, name=None):
            param_pos = (round(float(p1), precision),
                         round(float(p2), precision))
            bf_s = str(bf)
            for codim2_type in codim2_types:
                if codim2_type not in bf_s:
                    continue
                if param_pos in found_codim2[codim2_type]:
                    continue
                found_codim2[codim2_type].append(param_pos)
                idx = len(found_codim2[codim2_type])

                sub = _recurse_codim2(
                    codim2_type=codim2_type, pyauto=pyauto_instance,
                    origin=cont, idx=idx, params=params, name=name_tmp,
                    recursion=recursion,
                    max_recursion_depth=max_recursion_depth,
                    kwargs_1D_lc_cont=kwargs_1D_lc_cont,
                    base_kwargs=kwargs,
                )
                continuations.update(sub)
                break  # don't apply multiple handlers to the same point

    return continuations


def _recurse_codim2(codim2_type: str, pyauto: ODESystem, origin: Any, idx: int,
                    params: list, name: str, recursion: int,
                    max_recursion_depth: int, kwargs_1D_lc_cont: dict,
                    base_kwargs: dict) -> dict:
    """Dispatch codim-2 recursion to the type-specific handler."""
    if codim2_type == 'ZH':
        return _recurse_zh(pyauto, origin, idx, params, name,
                            recursion, max_recursion_depth, base_kwargs)
    if codim2_type == 'GH':
        return _recurse_gh(pyauto, origin, idx, params, name,
                            recursion, max_recursion_depth,
                            kwargs_1D_lc_cont, base_kwargs)
    if codim2_type == 'BT':
        return _recurse_bt(pyauto, origin, idx, params, name,
                            recursion, max_recursion_depth, base_kwargs)
    return {}


def _recurse_zh(pyauto: ODESystem, origin: Any, idx: int, params: list, name: str,
                recursion: int, max_recursion_depth: int, base_kwargs: dict) -> dict:
    """ZH: 1D equilibrium continuation in params[0] stopping at the nearest
    fold or Hopf, then recursive 2D continuation of that codim-1 point.
    """
    kwargs_tmp = dict(base_kwargs)
    kwargs_tmp.update({
        'ILP': 1, 'IPS': 1, 'ISW': 1, 'ISP': 2,
        'ICP': params[0], 'STOP': ['LP1', 'HB1'],
    })
    try:
        s_tmp, c_tmp = pyauto.run(
            starting_point=f'ZH{idx}', origin=origin,
            bidirectional=True, **kwargs_tmp,
        )
    except Exception as exc:
        warnings.warn(
            f"codim2_search: ZH{idx} recursion failed "
            f"({type(exc).__name__}: {exc}); skipping.",
            UserWarning, stacklevel=3,
        )
        return {}

    codim1_bifs, _ = pyauto.extract(['bifurcation'], cont=c_tmp)
    bif_values = codim1_bifs['bifurcation']
    if _bif_series_contains(bif_values, 'LP'):
        next_p = 'LP1'
    elif _bif_series_contains(bif_values, 'HB'):
        next_p = 'HB1'
    else:
        return {}
    sub_name = f"{name}/ZH{idx}"
    return codim2_search(
        params=params, starting_points=[next_p], origin=c_tmp,
        pyauto_instance=pyauto, recursion=recursion + 1,
        max_recursion_depth=max_recursion_depth, periodic=False,
        name=sub_name, **base_kwargs,
    )


def _recurse_gh(pyauto: ODESystem, origin: Any, idx: int, params: list, name: str,
                recursion: int, max_recursion_depth: int,
                kwargs_1D_lc_cont: dict, base_kwargs: dict) -> dict:
    """GH (Bautin): switch to LC continuation in params[0] + PAR(11).

    GH lies on a Hopf curve and marks a change of criticality
    (super- ↔ subcritical Hopf). The local LC family changes accordingly;
    following it in 1 parameter is the standard exploratory move. If a
    fold-of-cycle (LP) appears on the LC branch, recurse on it as a 2-param
    LC continuation.
    """
    kwargs_tmp = dict(base_kwargs)
    kwargs_tmp.update({
        'ILP': 1, 'IPS': 2, 'ISW': -1, 'ISP': 2,
        'ICP': [params[0], 11], 'STOP': ['LP1'],
    })
    if kwargs_1D_lc_cont:
        kwargs_tmp.update(kwargs_1D_lc_cont)
    sub_name = f"{name}/GH{idx}"
    try:
        _, c_tmp = pyauto.run(
            starting_point=f'GH{idx}', origin=origin, name=sub_name,
            bidirectional=True, **kwargs_tmp,
        )
    except Exception as exc:
        warnings.warn(
            f"codim2_search: GH{idx} LC continuation failed "
            f"({type(exc).__name__}: {exc}); skipping. Pass "
            f"kwargs_1D_lc_cont={{...}} to tune the LC continuation, or handle "
            f"this point manually.",
            UserWarning, stacklevel=3,
        )
        return {}
    continuations = {sub_name: c_tmp}

    if recursion + 1 < max_recursion_depth:
        try:
            codim1_bifs, _ = pyauto.extract(['bifurcation'], cont=c_tmp)
        except KeyError:
            return continuations
        if _bif_series_contains(codim1_bifs['bifurcation'], 'LP'):
            continuations.update(codim2_search(
                params=params, starting_points=['LP1'], origin=c_tmp,
                pyauto_instance=pyauto, recursion=recursion + 1,
                max_recursion_depth=max_recursion_depth, periodic=True,
                name=sub_name, **base_kwargs,
            ))
    return continuations


def _recurse_bt(pyauto: ODESystem, origin: Any, idx: int, params: list, name: str,
                recursion: int, max_recursion_depth: int, base_kwargs: dict) -> dict:
    """BT (Bogdanov-Takens): 1D equilibrium continuation from BT stopping at
    the nearest Hopf, then 2D continuation of the resulting Hopf curve.

    A homoclinic curve also emerges from BT; this is **not** followed here
    because it requires HomCont parameters (IPS=9). Users who want that
    can take the returned 2D-Hopf continuation and start a homoclinic
    continuation manually following the auto-07p documentation.
    """
    kwargs_tmp = dict(base_kwargs)
    kwargs_tmp.update({
        'ILP': 1, 'IPS': 1, 'ISW': 1, 'ISP': 2,
        'ICP': params[0], 'STOP': ['HB1'],
    })
    try:
        _, c_tmp = pyauto.run(
            starting_point=f'BT{idx}', origin=origin,
            bidirectional=True, **kwargs_tmp,
        )
    except Exception as exc:
        warnings.warn(
            f"codim2_search: BT{idx} equilibrium continuation failed "
            f"({type(exc).__name__}: {exc}); skipping. Note that homoclinic "
            f"continuation from BT is not auto-handled — see auto-07p docs "
            f"for IPS=9 / HomCont if needed.",
            UserWarning, stacklevel=3,
        )
        return {}

    codim1_bifs, _ = pyauto.extract(['bifurcation'], cont=c_tmp)
    if not _bif_series_contains(codim1_bifs['bifurcation'], 'HB'):
        return {}
    if recursion + 1 >= max_recursion_depth:
        return {}
    sub_name = f"{name}/BT{idx}"
    return codim2_search(
        params=params, starting_points=['HB1'], origin=c_tmp,
        pyauto_instance=pyauto, recursion=recursion + 1,
        max_recursion_depth=max_recursion_depth, periodic=False,
        name=sub_name, **base_kwargs,
    )
