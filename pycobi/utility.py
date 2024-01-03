from typing import List, Any
import numpy as np


# auto-object related utility functions
#######################################


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


def get_solution_stability(branch: Any, solution: Any, point: int) -> bool:
    """Return stability of a solution.

    Parameters
    ----------
    branch
        Solution branch object as returned by `auto`.
    solution
        Solution object on the branch.
    point
        Index of the solution on the branch.

    Returns
    -------
    bool
        True, if the solution is stable.
    """
    point_idx = get_point_idx(branch[0].diagnostics.data, point=point)
    diag = branch[0].diagnostics.data[point_idx]['Text']

    if "Eigenvalues" in diag:
        diag_line = "Eigenvalues  :   Stable:"
    elif "Multipliers" in diag:
        diag_line = "Multipliers:     Stable:"
    else:
        return solution.b['solution'].b['PT'] < 0
    idx = diag.find(diag_line) + len(diag_line)
    value = int(diag[idx:].split("\n")[0])
    target = solution.data['NDIM']
    return value >= target


def get_solution_variables(solution: Any, variables: list, extract_time: bool = False) -> list:
    """Extract state variable values (time series) from a steady-state (periodic) solution object.

    Parameters
    ----------
    solution
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
    if hasattr(solution, 'b') and extract_time:
        solution = solution.b['solution']
        solutions = [solution.indepvararray]
    else:
        solutions = []
    solutions = [solution[v] for v in variables] + solutions
    return solutions


def get_solution_params(solution: Any, params: list) -> list:
    """Extract parameter values from a solution object.

    Parameters
    ----------
    solution
        Solution object as returned by `auto`.
    params
        Keys of the parameters to be extracted.

    Returns
    -------
    list
        List of parameter values.
    """
    if hasattr(solution, 'b'):
        solution = solution.b['solution']
    return [solution[p] for p in params]


def get_solution_eigenvalues(solution: Any, branch: int, point: int) -> list:
    """Extracts eigenvalue spectrum from a solution point on a branch of solutions.

    Parameters
    ----------
    solution
        Solution object as returned by `auto`.
    branch
        Branch ID as assigned by `auto`.
    point
        Index of the solution on the branch.

    Returns
    -------
    list
        List of eigenvalues. If solution is periodic, the list contains floquet multipliers instead of eigenvalues.
    """

    eigenvals = []

    # extract point index from diagnostics
    point_idx = get_point_idx(solution[0].diagnostics, point=point)
    diag = solution[0].diagnostics.data[point_idx]['Text']

    if "NOTE:No converge" in diag:
        return eigenvals

    # find index of line in diagnostic output where eigenvalue information are starting
    idx = diag.find('Stable:')
    if not idx:
        return eigenvals
    diag = diag[idx:]
    diag_split = diag.split("\n")

    # check whether branch and point identifiers match the targets
    branch_str = f' {branch} '
    point_str = f' {point+1} '
    if branch_str in diag_split[1] and point_str in diag_split[1] and \
            ('Eigenvalue' in diag_split[1] or 'Multiplier' in diag_split[1]):

        # go through lines of system diagnostics and extract eigenvalues/floquet multipliers
        i = 1
        while i < len(diag_split):

            diag_tmp = diag_split[i]
            diag_tmp_split = [d for d in diag_tmp.split(' ') if d != ''][2:]

            # check whether line contains eigenvals or floquet mults. If not, stop while loop.
            if not diag_tmp_split:
                break
            if 'Eigenvalue' not in diag_tmp_split[0] and 'Multiplier' not in diag_tmp_split[0]:
                break

            # extract eigenvalues/floquet multipliers
            try:
                idx2 = diag_tmp_split.index(f"{i}")
            except ValueError:
                idx2 = diag_tmp_split.index(f"{i}:")
            real = float(diag_tmp_split[idx2+1])
            imag = float(diag_tmp_split[idx2+2])
            eigenvals.append(complex(real, imag))

            i += 1

    return eigenvals


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
