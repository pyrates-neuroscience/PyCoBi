import re
from pathlib import Path
from typing import Any, Iterable, List, Optional, Union

import numpy as np
import pandas as pd


# auto-object related utility functions
#######################################


# ---------------------------------------------------------------------------
# Regexes for parsing auto-07p's per-point diagnostic text blocks.
#
# A typical block (one per point) looks like:
#
#     BR    PT  IT         PAR           L2-NORM
#      1     5   0       -4.56412E+00   1.82904E+00
#      ...
#      1     5         Hopf Function:   0.00000E+00
#      1     5         Eigenvalues  :   Stable:   2
#      1     5         Eigenvalue  1:  -2.12891E+00   0.00000E+00
#      1     5         Eigenvalue  2:  -5.17236E+00   0.00000E+00
#      1     5         Fold Function:   9.49088E-01
#      ...
#
# Equilibrium continuations emit "Eigenvalues" + "Eigenvalue N:"; limit-cycle
# continuations emit "Multipliers" + "Multiplier N:" with Floquet multipliers.
# ---------------------------------------------------------------------------

_STABLE_COUNT_RE = re.compile(
    r'(?:Eigenvalues|Multipliers)\s*:\s*Stable:\s*(\d+)',
)
# Capture (index, real, imag) per eigenvalue / multiplier line. The leading
# `\d+\s+-?\d+` consumes the BR/PT prefix; PT is signed in auto's output.
_EIG_LINE_RE = re.compile(
    r'^\s*\d+\s+-?\d+\s+(?:Eigenvalue|Multiplier)\s+(\d+):\s*(\S+)\s+(\S+)\s*$',
    re.MULTILINE,
)
_NO_CONVERGENCE = "NOTE:No converge"


def _normalize_icp(branch_id: Any, icp: Any) -> tuple:
    """Coerce auto-07p's ICP entry (int for a single-parameter continuation,
    iterable for multi-parameter) to a tuple."""
    icp = (icp,) if isinstance(icp, int) else tuple(icp)
    return branch_id, icp


def get_branch_info(branch: Any) -> tuple:
    """Extract the branch id and continuation parameter(s) from an auto-07p
    solution branch object.

    Three call shapes are supported, tried in order:

      1. `bifDiag` (auto's top-level multi-branch container) — read `.BR` /
         `.c['ICP']` from the first sub-branch.
      2. A single branch object indexed by `'BR'` with its own `.c` config.
      3. An IVP-style nested structure where the BR field is buried under an
         `'RG'` (regular point) label on the first sub-branch. The previous
         implementation walked the first 10 label indices and re-raised
         whichever ``KeyError`` happened last; now we iterate the keys
         explicitly and, on exhaustion, raise a single ``KeyError`` naming
         every key that was tried.

    Returns
    -------
    tuple
        ``(branch_id, icp)`` where ``icp`` is a tuple of one or more PAR
        indices.
    """
    # --- Path 1: bifDiag with .BR / .c on the first sub-branch.
    try:
        return _normalize_icp(branch[0].BR, branch[0].c['ICP'])
    except (AttributeError, ValueError, KeyError, TypeError, IndexError):
        pass

    # --- Path 2: single branch object with 'BR' as a key.
    try:
        return _normalize_icp(branch['BR'], branch.c['ICP'])
    except (AttributeError, KeyError, TypeError):
        pass

    # --- Path 3: IVP-style — find a labeled point on branch[0] whose nested
    # ``RG`` solution carries the BR field, and read it from there.
    try:
        icp = branch[0].c['ICP']
        labels = branch[0].labels.by_index
    except (AttributeError, KeyError, TypeError, IndexError) as e:
        raise AttributeError(
            f"get_branch_info: cannot navigate {type(branch).__name__!s}; "
            f"none of the supported branch shapes matched. Last error: {e}"
        ) from e

    tried = []
    for sol_key, sol_entry in labels.items():
        tried.append(sol_key)
        try:
            br_id = sol_entry['RG']['solution']['data']['BR']
        except (KeyError, TypeError):
            continue
        return _normalize_icp(br_id, icp)

    raise KeyError(
        f"get_branch_info: no labeled point on branch[0] carries an "
        f"RG solution with a 'BR' field. Tried {len(tried)} label "
        f"indices: {tried}."
    )


def get_solution_keys(branch: Any) -> List[str]:
    """Extract the names of all solutions on a branch.

    Parameters
    ----------
    branch
        Solution branch object as returned by `auto`.

    Returns
    -------
    List[str]
        List with solution names.
    """
    keys = []
    for idx in range(len(branch.data)):
        keys += [key for key, val in branch[idx].labels.by_index.items()
                 if val and 'solution' in tuple(val.values())[0]]
    return keys


def get_point_idx(diag: list, point: int) -> int:
    """Extract list idx of correct diagnostics string for continuation point with index `point`.

    Parameters
    ----------
    diag
        Diagnostics as stored by `auto` on solution objects.
    point
        Index of the solution of interest on the branch.

    Returns
    -------
    int
        Point index for `diag`.
    """

    idx = point
    while idx < len(diag)-1:

        diag_tmp = diag[idx]['Text']
        if "Location of special point" in diag_tmp and "Convergence" not in diag_tmp:
            idx += 1
        elif "NOTE:Retrying step" in diag_tmp or "Total Time " in diag_tmp:
            idx += 1
        else:

            # find index of line after first appearance of BR
            diag_tmp = diag_tmp.split('\n')
            idx_line = 1
            while idx_line < len(diag_tmp)-1:
                if 'BR' in diag_tmp[idx_line].split(' '):
                    break
                idx_line += 1
            diag_tmp = diag_tmp[idx_line+1]

            # find point number in text
            line_split = [d for d in diag_tmp.split(' ') if d != ""]
            if abs(int(line_split[1])) < point+1:
                idx += 1
            elif abs(int(line_split[1])) == point+1:
                return idx
            else:
                raise ValueError(f"Point with index {point+1} was not found on solution. Last auto output line that "
                                 f"was checked: \n {diag_tmp}")
    return idx


# extract system state information from auto objects
####################################################

def get_point_diagnostics(s: Any) -> str:
    """Return the Auto-07p string that contains diagnostic data for a particular solution.

    Parameters
    ----------
    s
        Solution object on the branch.

    Returns
    -------
    str
        String with solution diagnostics.
    """

    try:

        diag = s.b.branch.diagnostics.data
        point_idx = get_point_idx(diag, point=s.b.idx)
        diag = diag[point_idx]["Text"]

    except AttributeError as e:

        raise e
        # diag = branch[0].diagnostics.data
        # point_idx = get_point_idx(diag, point=point)
        # diag = diag[point_idx]['Text']

    return diag

def parse_point_diagnostics(s: Any, diag: Optional[str] = None) -> dict:
    """Parse the auto-07p per-point diagnostic text block into structured fields.

    Used as a one-shot parser by `get_solution_stability` and `get_solution_eigenvalues`
    so multi-metric summaries (`_create_summary` with stability + eigenvalues +
    lyapunov) only pay the regex cost once per point.

    Parameters
    ----------
    s
        Solution object as returned by auto-07p.
    diag
        Optional pre-fetched diagnostic text for `s` (skips the
        `get_point_diagnostics(s)` lookup). Pass this when the caller already
        has the text, e.g. when iterating multiple points and batching the
        fetch.

    Returns
    -------
    dict with keys:
        ``stable``      bool or None — True if the solution is stable. None if
                        no Eigenvalues/Multipliers line is present AND auto's
                        PT-sign convention can't be read from `s.b`.
        ``eigenvalues`` list[complex] — eigenvalues (equilibrium continuation)
                        or Floquet multipliers (limit-cycle continuation), in
                        the order auto-07p emits them. Empty if the point did
                        not converge or no spectrum is recorded.
        ``text``        str — the raw diagnostic block (kept for callers that
                        still need to grep for things this parser doesn't
                        surface, e.g. Hopf/Fold/BP function values).
    """
    if diag is None:
        diag = get_point_diagnostics(s)
    out = {'text': diag, 'stable': None, 'eigenvalues': []}

    if _NO_CONVERGENCE in diag:
        return out

    # --- stability ---
    m = _STABLE_COUNT_RE.search(diag)
    if m:
        try:
            ndim = s.data['NDIM']
        except (AttributeError, KeyError, TypeError):
            ndim = None
        if ndim is not None:
            out['stable'] = int(m.group(1)) >= ndim
    else:
        # No spectrum recorded — fall back to auto's native PT-sign convention
        # (auto encodes stability as the sign of the point index in fort.7).
        try:
            out['stable'] = s.b['PT'] > 0
        except (AttributeError, KeyError, TypeError):
            pass  # leave as None

    # --- eigenvalues / Floquet multipliers ---
    for m in _EIG_LINE_RE.finditer(diag):
        real, imag = m.group(2), m.group(3)
        try:
            out['eigenvalues'].append(complex(float(real), float(imag)))
        except ValueError:
            # malformed numeric in the diagnostic line — skip this entry but
            # keep parsing the rest
            pass

    return out


def get_solution_stability(s: Any) -> bool:
    """Return stability of a solution.

    Thin wrapper around `parse_point_diagnostics` — when extracting multiple
    diagnostic fields for the same point, call `parse_point_diagnostics` once
    and reuse its dict to avoid re-parsing the text.
    """
    return bool(parse_point_diagnostics(s)['stable'])


def get_solution_variables(s: Any, variables: list, extract_time: bool = False) -> list:
    """Extract state variable values (time series) from a steady-state (periodic) solution object.

    Parameters
    ----------
    s
        Solution object as returned by `auto`.
    variables
        Keys of the state variables to be extracted.
    extract_time
        If true, attempt to extract a time vector as well.

    Returns
    -------
    list
        List of variable values/time series of the provided solution.
    """
    if hasattr(s, 'b') and extract_time:
        s = s.b['solution']
        solutions = [s.indepvararray]
    else:
        solutions = []
    solutions = [s[v] for v in variables] + solutions
    return solutions


def get_solution_params(s: Any, params: list) -> list:
    """Extract parameter values from a solution object.

    Parameters
    ----------
    s
        Solution object as returned by `auto`.
    params
        Keys of the parameters to be extracted.

    Returns
    -------
    list
        List of parameter values.
    """
    if hasattr(s, 'b'):
        s = s.b['solution']
    return [s[p] for p in params]


def get_solution_eigenvalues(s: Any, branch: int = None, point: int = None) -> list:
    """Return eigenvalues (or Floquet multipliers, for periodic solutions) of a point.

    Thin wrapper around `parse_point_diagnostics`. The `branch` and `point`
    arguments are kept for backward compatibility but are no longer used —
    `parse_point_diagnostics` reads the diagnostic block via `s.b.branch` /
    `s.b.idx` directly.
    """
    return parse_point_diagnostics(s)['eigenvalues']


def get_lyapunov_exponents(eigenvals, period) -> list:
    """Calculate Lyapunov exponents from eigenvalues/floquet multipliers.

    Parameters
    ----------
    eigenvals
        Eigenvalue or floquet multiplier spectrum.
    period
        Period of the periodic solution, if `eigenvals` is a Floquet multiplier spectrum of a periodic solution.

    Returns
    -------
    list
        Local Lyapunov exponent spectrum.
    """
    return [np.real(np.log(ev)/period) if period else np.real(ev) for ev in eigenvals]


def fractal_dimension(lyapunov_exponents: list) -> float:
    """Calculates the fractal or information dimension of an attractor of a dynamical system from its lyapunov
    epxonents, according to the Kaplan-Yorke formula (Kaplan and Yorke, 1979).

    Parameters
    ----------
    lyapunov_exponents
        List containing the lyapunov spectrum of a solution of a dynamical system.

    Returns
    -------
    float
        Fractal dimension of the attractor of the system.
    """

    LEs = np.sort(lyapunov_exponents)[::-1]
    if np.sum(LEs) > 0:
        return len(LEs)
    k = 0
    for j in range(len(LEs)-1):
        k = j+1
        if np.sum(LEs[:k]) < 0:
            k -= 1
            break
    return k + np.sum(LEs[:k]) / np.abs(LEs[k])


def write_auto_dat(
    df: pd.DataFrame,
    path: Union[str, Path],
    *,
    columns: Optional[Iterable] = None,
    n_points: Optional[int] = None,
    period: Optional[float] = None,
    normalize_time: bool = True,
) -> Path:
    r"""Write a time-series ``DataFrame`` as an auto-07p ``.dat`` file.

    Auto-07p reads ``.dat`` files via the ``dat=`` argument of its python
    ``run()`` command to initialise periodic-orbit continuations from a
    user-supplied time series (one column of time + one column per state
    variable, whitespace-separated). This helper writes such a file from a
    ``pandas.DataFrame`` — most commonly the output of
    :meth:`pyrates.CircuitTemplate.run`, but any time-indexed numeric
    DataFrame works.

    Parameters
    ----------
    df
        Time-indexed DataFrame. The index is interpreted as time; columns hold
        the state variables. A ``MultiIndex`` on columns is fine — every
        column entry survives into the output file in the order encountered
        (or the order given by ``columns``).
    path
        Output file path. A ``.dat`` suffix is appended if missing.
    columns
        Optional subset of columns to write, in the order they should appear.
        Defaults to every column of ``df``.
    n_points
        If given, linearly interpolate the time series onto this many
        evenly-spaced points before writing. Useful when the simulation
        produced more samples than auto-07p needs, or to ensure a regular
        mesh from an adaptive solver. Defaults to writing every input row.
    period
        Length (in the same units as the index) of one period of the orbit
        the file represents. Only consulted when ``normalize_time=True``.
        Defaults to ``df.index[-1] - df.index[0]`` — i.e. assumes the
        DataFrame already covers exactly one period.
    normalize_time
        When True (default), shift the time column so it starts at 0 and
        divide by ``period`` so it ends at 1. Auto-07p's default convention
        is a normalised ``[0, 1]`` mesh and the orbit's actual period lives in
        ``PAR(11)``; turn this off only if the consuming workflow expects raw
        simulation time.

    Returns
    -------
    pathlib.Path
        The (possibly suffix-extended) path the file was written to.

    Notes
    -----
    Recipe for hooking the file into a periodic-orbit continuation:

    .. code-block:: python

        write_auto_dat(sim_df, 'orbit')                 # → orbit.dat
        ode = ODESystem.from_template(template, ..., auto_constants=('lc',))
        sols, cont = ode.run(c='lc', ICP=['eta', 11], dat='orbit', name='lc0')

    The ``dat='orbit'`` argument tells auto-07p to read ``orbit.dat`` as the
    starting time series; ``ICP=['eta', 11]`` puts the orbit's period in
    PAR(11) so it varies along the continuation.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"expected pandas.DataFrame, got {type(df).__name__}")
    if len(df) < 2:
        raise ValueError("DataFrame must have at least 2 rows to define a time series")

    cols = list(df.columns) if columns is None else list(columns)
    if not cols:
        raise ValueError("at least one column must be selected")
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"columns not found in DataFrame: {missing}")

    data = df[cols].to_numpy()
    if not np.issubdtype(data.dtype, np.number):
        raise TypeError(
            f"selected columns must be numeric; got dtype {data.dtype}"
        )
    data = data.astype(float, copy=False)
    time = np.asarray(df.index, dtype=float)

    # Optional interpolation onto a regular mesh of n_points samples — done
    # before time normalisation so the interpolation uses the original spacing.
    if n_points is not None:
        if n_points < 2:
            raise ValueError("n_points must be at least 2")
        t_uniform = np.linspace(time[0], time[-1], n_points)
        data = np.column_stack([
            np.interp(t_uniform, time, data[:, j]) for j in range(data.shape[1])
        ])
        time = t_uniform

    if normalize_time:
        T = period if period is not None else float(time[-1] - time[0])
        if T <= 0:
            raise ValueError(
                f"cannot normalize time: period={T} is not positive"
            )
        time = (time - time[0]) / T

    out_path = Path(path)
    if out_path.suffix != ".dat":
        out_path = out_path.with_suffix(out_path.suffix + ".dat") \
            if out_path.suffix else out_path.with_suffix(".dat")

    # Format: free-form whitespace-separated, matching auto-07p's own demos
    # (e.g. pen.dat). 14-digit scientific notation gives ~12 significant
    # figures — well within double precision and trivial for auto's
    # free-format Fortran reader to consume.
    fmt = "%23.14E"
    table = np.column_stack([time, data])
    np.savetxt(out_path, table, fmt=fmt)
    return out_path
