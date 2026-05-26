import re
from typing import List, Any, Optional

import numpy as np


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


def get_branch_info(branch: Any) -> tuple:
    """Extract the relevant data from a solution branch object.

    Parameters
    ----------
    branch
        Solution branch as returned by `auto`.

    Returns
    -------
    tuple
        Branch data, continuation parameter.
    """
    try:
        branch, icp = branch[0].BR, branch[0].c['ICP']
    except (AttributeError, ValueError):
        try:
            branch, icp = branch['BR'], branch.c['ICP']
        except AttributeError:
            icp = branch[0].c['ICP']
            i = 0
            while i < 10:
                try:
                    sol_key = list(branch[0].labels.by_index.keys())[i]
                    branch = branch[0].labels.by_index[sol_key]['RG']['solution']['data']['BR']
                    break
                except KeyError as e:
                    i += 1
                    if i == 10:
                        raise e
    icp = (icp,) if type(icp) is int else tuple(icp)
    return branch, icp


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
