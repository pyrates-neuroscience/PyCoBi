from .pycobi import ODESystem
from typing import Union, Any
import numpy as np


def continue_period_doubling_bf(solution: dict, continuation: Union[str, int, Any], pyauto_instance: Any,
                                max_iter: int = 1000, iteration: int = 0, precision: int = 3, pds: list = None,
                                **kwargs) -> tuple:
    """Automatically continue a cascade of period doubling bifurcations. Returns the labels of the continuation and the
    pycobi instance they were run on.

    Parameters
    ----------
    solution
    continuation
    pyauto_instance
    max_iter
    iteration
    precision
    pds
    kwargs

    Returns
    -------
    tuple
    """
    if pds is None:
        pds = []
    solutions = []
    params = kwargs['ICP']
    i = 1
    name = f'pd_{iteration}'
    solutions.append(name)

    if iteration >= max_iter:
        return solutions, pyauto_instance

    for point, point_info in solution.items():

        if 'PD' in point_info['bifurcation']:

            s_tmp, cont = pyauto_instance.run(starting_point=f'PD{i}', name=name, origin=continuation,
                                              **kwargs)
            bfs = s_tmp.loc[:, ['bifurcation', f'PAR({params[0]})', f'PAR({params[1]})']]

            for bf, p1, p2 in bfs.itertuples(index=False, name=None):

                param_pos = np.round([p1, p2], decimals=precision)

                if "PD" in bf and not any([p[0] == param_pos[0] and p[1] == param_pos[1] for p in pds]):

                    pds.append(param_pos)
                    s_tmp2, pyauto_instance = continue_period_doubling_bf(solution=s_tmp, continuation=cont,
                                                                          pyauto_instance=pyauto_instance,
                                                                          iteration=iteration + 1,
                                                                          precision=precision, pds=pds,
                                                                          **kwargs)
                    solutions += s_tmp2
                    iteration += len(s_tmp2)
            i += 1

    return solutions, pyauto_instance


def codim2_search(params: list, starting_points: list, origin: Union[str, int, Any],
                  pyauto_instance: ODESystem, max_recursion_depth: int = 3, recursion: int = 0, periodic: bool = False,
                  kwargs_2D_lc_cont: dict = None, kwargs_2D_cont: dict = None,
                  precision=2, **kwargs) -> dict:
    """Performs automatic continuation of codim 1 bifurcation points in 2 parameters and searches for codimension 2
    bifurcations along the solution curves.

    Currently only Zero-Hopf (ZH) points trigger recursive sub-searches — the GH
    (generalised Hopf) and BT (Bogdanov-Takens) branches are not yet implemented;
    those bifurcation types are detected and ignored rather than recursed on.

    Parameters
    ----------
    params
    starting_points
    origin
    pyauto_instance
    max_recursion_depth
    recursion
    periodic
    kwargs_2D_lc_cont
    kwargs_2D_cont
    precision
    kwargs

    Returns
    -------

    """

    zhs = dict()
    continuations = dict()
    name = kwargs.pop('name', f"{params[0]}/{params[1]}")

    for p in starting_points:

        # continue curve of special solutions in 2 parameters
        kwargs_tmp = kwargs.copy()
        if periodic:
            kwargs_tmp.update({'ILP': 0, 'IPS': 2, 'ISW': 2, 'ISP': 2, 'ICP': list(params) + [11]})
            if kwargs_2D_lc_cont:
                kwargs_tmp.update(kwargs_2D_lc_cont)
        else:
            kwargs_tmp.update({'ILP': 0, 'IPS': 1, 'ISW': 2, 'ISP': 2, 'ICP': params})
            if kwargs_2D_cont:
                kwargs_tmp.update(kwargs_2D_cont)

        name_tmp = f"{name}:{p}"
        sols, cont = pyauto_instance.run(starting_point=p, origin=origin, name=name_tmp, bidirectional=True,
                                         **kwargs_tmp)
        continuations[name_tmp] = cont

        if recursion < max_recursion_depth:

            # get types of all solutions along curve
            codim2_bifs = sols.loc[:, ['bifurcation', f'PAR({params[0]})', f'PAR({params[1]})']]

            for bf, p1, p2 in codim2_bifs.itertuples(index=False, name=None):

                param_pos = np.round([p1, p2], decimals=precision)

                if "ZH" in bf and (p not in zhs or not any([p_tmp[0] == param_pos[0] and p_tmp[1] == param_pos[1]
                                                            for p_tmp in zhs[p]['pos']])):

                    if p not in zhs:
                        zhs[p] = {'count': 1, 'pos': [param_pos]}
                    else:
                        zhs[p]['count'] += 1
                        zhs[p]['pos'].append(param_pos)

                    # perform 1D continuation to find nearby fold bifurcation
                    kwargs_tmp = kwargs.copy()
                    kwargs_tmp.update({'ILP': 1, 'IPS': 1, 'ISW': 1, 'ISP': 2, 'ICP': params[0], 'STOP': ['LP1', 'HB1']
                                       })
                    s_tmp, c_tmp = pyauto_instance.run(starting_point=f"ZH{zhs[p]['count']}", origin=cont,
                                                       bidirectional=True, **kwargs_tmp)

                    codim1_bifs = s_tmp.loc[:, 'bifurcation']
                    if "LP" in codim1_bifs:
                        p_tmp = 'LP1'
                        name_tmp2 = f"{name_tmp}/ZH{zhs[p]['count']}"
                    elif "HB" in codim1_bifs:
                        p_tmp = 'HB1'
                        name_tmp2 = f"{name_tmp}/ZH{zhs[p]['count']}"
                    else:
                        continue

                    # perform 2D continuation of the fold or hopf bifurcation
                    continuations.update(codim2_search(params=params, starting_points=[p_tmp], origin=c_tmp,
                                                       pyauto_instance=pyauto_instance, recursion=recursion + 1,
                                                       max_recursion_depth=max_recursion_depth, periodic=False,
                                                       name=name_tmp2, **kwargs))

                # GH (generalised Hopf) and BT (Bogdanov-Takens) recursion
                # not yet implemented — those bifurcation types are detected
                # by auto-07p but currently ignored here.

    return continuations