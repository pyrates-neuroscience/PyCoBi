import os
import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame, MultiIndex, Series
from mpl_toolkits.mplot3d import Axes3D
from pyrates import CircuitTemplate, clear
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from typing import Union, Any, Optional
import pickle

from .utility import get_solution_stability, get_solution_keys, get_branch_info, get_solution_variables, \
    get_solution_params, get_solution_eigenvalues, get_lyapunov_exponents


class ODESystem:

    __slots__ = ["auto_solutions", "results", "_orig_dir", "dir", "_auto", "_last_cont", "_cont_num", "_results_map",
                 "_branches", "_bifurcation_styles", "_temp", "additional_attributes", "_eq", "_var_map",
                 "_var_map_inv"]

    def __init__(self, eq_file: str, working_dir: str = None, auto_dir: str = None, init_cont: bool = True,
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
            If true, an integration with respect to time will be performed, using the equation file provided via the
            keyword argument `e=<fname>` (a file named `<fname>.f90` should exist in `working_dir`)
            and the auto constants provided via the keyword argument `c=<fname>`
            (a file named `c.<fname>` should exist in `working_dir`).
        params
            Optional ordered list with names of all parameters in the model equations. Can be used to refer to model
            parameters.
        state_vars
            Optional ordered list that provides a name for each entry in the state vector of the model equations.
        kwargs
            Additional keyword arguments that will be provided to the `ODESystem.run` method for performing the time
            integration.
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

        # create a map that links variable/parameter indices to string-based keys
        self._var_map = {"t": {"cont": 14, "plot": "PAR(14)"}}
        self._var_map_inv = {}
        if params:
            for i, key in enumerate(params):
                self._var_map[key] = {"cont": i+1, "plot": f"PAR({i+1})"}
        if state_vars:
            for i, key in enumerate(state_vars):
                self._var_map[key] = {"cont": i+1, "plot": f"U({i+1})"}
        for key, val in self._var_map.items():
            self._var_map_inv[val["cont"]] = key
            self._var_map_inv[val["plot"]] = key

        # perform initial continuation in time to ensure convergence to steady-state solution
        if init_cont:
            _ = self.run(ICP=[14], **kwargs)

    def __getitem__(self, item):
        try:
            return self.results[item]
        except KeyError:
            return self.results[self._results_map[item]]

    @property
    def pyrates_template(self):
        return self._temp

    def close_session(self, clear_files: bool = False, **kwargs):
        if clear_files:
            clear(self._temp, **kwargs)
        os.chdir(self._orig_dir)

    @classmethod
    def from_yaml(cls, path: str, working_dir: str = None, auto_dir: str = None, init_cont: bool = True,
                  init_kwargs: dict = None, **kwargs):
        """

        Parameters
        ----------
        path
            Path to a YAML model definition file.
        working_dir
            Directory in which all the fortran equation and auto-07p constant files are saved.
        auto_dir
            Installation directory of auto-07p.
        init_cont
            If true, an integration with respect to time will be performed, using the equation file provided via the
            keyword argument `e=<fname>` (a file named `<fname>.f90` should exist in `working_dir`)
            and the auto constants provided via the keyword argument `c=<fname>`
            (a file named `c.<fname>` should exist in `working_dir`).
        init_kwargs
            Additional keyword arguments that will be provided to the `ODESystem.run` method for performing the time
            integration.
        kwargs
            Additional keyword arguments provided to the `pyrates.CircuitTemplate.get_run_func` method that is used to
            generate the fortran equation file and the auto constants file that will be used to initialize `ODESystem`.

        Returns
        -------
        ODESystem
            `ODESystem` instance.
        """

        # preparations
        func_name = kwargs.pop("func_name", "vector_field")
        file_name = kwargs.pop("file_name", "system_equations")
        file_name_full = f"{working_dir}/{file_name}" if working_dir else file_name
        dt = kwargs.pop("step_size", 1e-3)
        solver = kwargs.pop("solver", "scipy")
        if init_kwargs is None:
            init_kwargs = {}

        # initialize circuit template
        template = CircuitTemplate.from_yaml(path)

        # update circuit template variables
        if "node_vars" in kwargs:
            template.update_var(node_vars=kwargs.pop("node_vars"))
        if "edge_vars" in kwargs:
            template.update_var(edge_vars=kwargs.pop("edge_vars"))

        # generate fortran files
        _, _, params, state_vars = template.get_run_func(func_name, dt, file_name=file_name_full, backend="fortran",
                                                         float_precision="float64", auto=True, vectorize=False,
                                                         solver=solver, **kwargs)

        # initialize ODESystem
        return cls(working_dir=working_dir, auto_dir=auto_dir, init_cont=init_cont, e=file_name, c="ivp",
                   template=template, params=params[3:], state_vars=list(state_vars), eq_file=file_name_full,
                   **init_kwargs)

    @classmethod
    def from_file(cls, filename: str, auto_dir: str = None):
        """Load `ODESystem` from file using pickle.

        Parameters
        ----------
        filename
            Name of the file that contains summary and solution data of an `ODESystem` instance.
        auto_dir
            Installation directory of `auto-07p`.

        Returns
        -------
        ODESystem
            ODESystem instance containing all the attributes that are available as dictionary entries in `filename`.
        """
        pyauto_instance = cls('', auto_dir=auto_dir, init_cont=False)
        data = pickle.load(open(filename, 'rb'))
        for key, val in data.items():
            attr = getattr(pyauto_instance, key)
            if type(attr) is dict:
                attr.update(val)
            else:
                raise AttributeError(f'Attribute {key} already exists on this `ODESystem` instance.')
        return pyauto_instance

    def to_file(self, filename: str, results_only: bool = True, **kwargs) -> None:
        """Save continuation results on disc via pickle.

        Parameters
        ----------
        filename
            Name of the file to which the results should be saved.
        results_only
            If true, only the `PyCoBi` recordings will be saved.
        kwargs
            Additional data/information to be saved to file. Will be available under the key 'additional_attributes'.

        Returns
        -------
        None
        """

        if results_only:
            data = {'results': self.results, '_branches': self._branches, '_results_map': self._results_map}
        else:
            skip = ["dir", "_orig_dir"]
            data = {key: getattr(self, key) for key in self.__slots__ if key not in skip}
        data.update({'additional_attributes': kwargs})

        try:
            pickle.dump(data, open(filename, 'x'))
        except (FileExistsError, TypeError):
            pickle.dump(data, open(filename, 'wb'))

    def run(self, origin: Union[int, str, object] = None, starting_point: Union[str, int] = None, variables: list = None,
            params: list = None, get_stability: bool = True, get_period: bool = False, get_timeseries: bool = False,
            get_eigenvals: bool = False, get_lyapunov_exp: bool = False, bidirectional: bool = False, name: str = None,
            **auto_kwargs) -> tuple:
        """
        Wraps auto-07p command `run` and stores requested solution details on instance.

        Parameters
        ----------
        origin
            Key of the solution branch that contains the solution `starting_point`, from which the new continuation will
            be started.
        starting_point
            Key of the solution on the solution branch `origin`, from which the continuation procedure will be initiated.
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
        bidirectional
            If true, parameter continuation will be performed into both directions for a given continuation parameter.
        name
            Name, under which the resulting solution branch will be accessible for future continuations.
        auto_kwargs
            Additional keyword arguments to be passed to the auto command `run`.

        Returns
        -------
        tuple
            DataFrame with the results, auto solution branch object.
        """

        # auto call
        ###########

        # extract starting point of continuation
        if self._last_cont == 0:
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
        constant_file = auto_kwargs.pop('c', None)
        auto_kwargs = self._map_auto_kwargs(auto_kwargs)
        if constant_file:
            solution = self._call_auto(starting_point, origin, c=constant_file, **auto_kwargs)
            auto_kwargs['c'] = constant_file
        else:
            solution = self._call_auto(starting_point, origin, **auto_kwargs)

        # extract information from auto solution
        ########################################

        # extract branch and solution info
        new_branch, new_icp = get_branch_info(solution)
        new_points = get_solution_keys(solution)

        # get all passed variables and params
        _, solution_tmp = self.get_solution(point=new_points[0], cont=solution)
        if variables is None:
            variables = self._get_all_var_keys(solution_tmp)
        if params is None:
            try:
                params = self._get_all_param_keys(solution_tmp)
            except KeyError:
                n_params = auto_kwargs['NPAR']
                params = [f"PAR({i})" for i in range(1, n_params+1)]

        # extract continuation results
        if new_icp[0] == 14:
            get_stability = False
        summary = self._create_summary(solution=solution, points=new_points, variables=variables,
                                       params=params, timeseries=get_timeseries, stability=get_stability,
                                       period=get_period, eigenvals=get_eigenvals, lyapunov_exp=get_lyapunov_exp)

        # store solution and extracted information in pycobi
        ####################################################

        # merge auto solutions if necessary and create key for auto solution
        if new_branch in self._branches and origin in self._branches[new_branch] \
                and new_icp in self._branches[new_branch][origin]:

            # get key from old solution and merge with new solution
            solution_old = self.get_solution(origin)
            pyauto_key = solution_old.pycobi_key
            solution, summary = self.merge(pyauto_key, solution, summary, new_icp)

        elif name == 'bidirect:cont2' and not bidirectional and 'DS' in auto_kwargs and auto_kwargs['DS'] == '-':

            # get key from old solution and merge with new solution
            solution_old = self.auto_solutions[self._last_cont]
            pyauto_key = solution_old.pycobi_key
            solution, summary = self.merge(pyauto_key, solution, summary, new_icp)

        else:

            # create pycobi key for solution
            pyauto_key = self._cont_num + 1 if self._cont_num in self.auto_solutions else self._cont_num
            solution.pycobi_key = pyauto_key

        # set up dictionary fields in _branches for new solution
        if new_branch not in self._branches:
            self._branches[new_branch] = {pyauto_key: []}
        elif pyauto_key not in self._branches[new_branch]:
            self._branches[new_branch][pyauto_key] = []

        # store auto solution under unique pycobi cont
        self.auto_solutions[pyauto_key] = solution
        self._last_cont = pyauto_key
        self._branches[new_branch][pyauto_key].append(new_icp)

        self.results[pyauto_key] = summary
        self._cont_num = len(self.auto_solutions)
        if name and name != 'bidirect:cont2':
            self._results_map[name] = pyauto_key

        # if continuation should be bidirectional, call this method again with reversed continuation direction
        ######################################################################################################

        if bidirectional:
            auto_kwargs.pop('DS', None)
            _, solution = self.run(origin, starting_point, variables=variables, params=params,
                                   get_stability=get_stability, get_period=get_period, get_timeseries=get_timeseries,
                                   get_eigenvals=get_eigenvals, get_lyapunov_exp=get_lyapunov_exp, bidirectional=False,
                                   name='bidirect:cont2', DS='-', **auto_kwargs)

        return self.results[pyauto_key], solution

    def merge(self, key: int, cont, summary: DataFrame, icp: tuple):
        """Merges two solutions from two separate auto continuations.

        Parameters
        ----------
        key
            PyCoBi identifier under which the merged solution should be stored. Must be equal to identifier of first
            continuation.
        cont
            auto continuation object that should be merged with the continuation object under `key`.
        summary
            PyCoBi continuation summary that should be merged with continuation summary under `key`.
        icp
            Continuation parameter that was used in both continuations that are to be merged.
        """

        # merge solutions
        #################

        # call merge in auto
        solution = self._auto.merge(self.auto_solutions[key] + cont)
        solution.pycobi_key = key

        # store solution in pycobi
        self.auto_solutions[key] = solution
        self._last_cont = solution

        # merge results summaries
        #########################

        summary_old = self.results[key]

        # extract end points
        points_old, points_new = summary_old.index, summary.index
        start_old, start_new, end_old, end_new = points_old[0], points_new[0], points_old[-1], points_new[-1]

        # extract ICP values
        key = f"PAR({icp[0]})"
        if key in self._var_map_inv:
            key = self._var_map_inv[key]
        icp_old = summary_old.loc[end_old, key].values[0]
        icp_new = summary.loc[end_new, key].values[0]

        # connect starting points of both continuations and re-label points accordingly
        if icp_new < icp_old:
            conts_sorted = [summary, summary_old]
            keys_sorted = [points_new, points_old]
        else:
            conts_sorted = [summary_old, summary]
            keys_sorted = [points_old, points_new]
        end_point = keys_sorted[0][-1]

        # move points into combined summary
        summary_final = {}
        for p1, p2 in zip(keys_sorted[0], keys_sorted[0][::-1]):
            summary_final[p1] = _extract_merge_point(p2, conts_sorted[0])
        for p in keys_sorted[1]:
            summary_final[p+end_point] = _extract_merge_point(p, conts_sorted[1])

        # store updated summary
        self.results[key] = DataFrame(summary_final).T

        return solution, self.results[key]

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
        branch, icp = get_branch_info(cont)

        if point is None:
            return cont

        # extract solution point from continuation object and its solution type
        try:
            s = cont(point)
            solution_name = point[:2]
        except (AttributeError, KeyError, TypeError):
            for idx in range(len(cont.data)):
                if type(point) is int:
                    s = cont[idx].labels.by_index[point]
                else:
                    count = 1
                    for p in cont[idx].labels.by_index.values():
                        key = list(p)[0]
                        if key in point and int(point.replace(key, "")) == count:
                            s = p
                            break
                        elif key in point:
                            count += 1
                    else:
                        raise ValueError(f'Invalid point {point} for continuation with ICP={icp} on branch {branch}.')
                if s:
                    solution_name = list(s.keys())[0]
                    break
            else:
                raise ValueError(f'Invalid point {point} for continuation with ICP={icp} on branch {branch}.')
            if solution_name != 'No Label':
                s = s[solution_name]['solution']

        return solution_name, s

    def extract(self, keys: list, cont: Union[Any, str, int], point: Union[str, int] = None) -> DataFrame:
        """Extract properties from a solution.

        Parameters
        ----------
        keys
            Keys of the properties (e.g. state variable names, parameter names, ...).
        cont
            Key of the solution branch.
        point
            Key of the solution on the branch.

        Returns
        -------
        DataFrame
            Contains the requested properties of the solution.
        """
        summary = self.get_summary(cont, point=point)
        columns = [k for k, _ in list(summary.keys())]
        keys = [key if key in columns else self._var_map_inv[key] for key in keys]
        if point:
            return summary.loc[point, keys]
        return summary.loc[:, keys]

    def plot_continuation(self, x: str, y: str, cont: Union[Any, str, int], ax: plt.Axes = None,
                          force_axis_lim_update: bool = False, **kwargs) -> LineCollection:
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
            results = self.extract([x, y], cont=cont)
            results['stability'] = np.asarray([True] * len(results[x]))
            results['bifurcation'] = np.asarray(['RG'] * len(results[x]))
        else:
            results = self.extract([x, y, 'stability', 'bifurcation'], cont=cont)

        # plot bifurcation points
        bifurcation_point_kwargs = ['default_color', 'default_marker', 'default_size', 'custom_bf_styles',
                                    'ignore']
        kwargs_tmp = {key: kwargs.pop(key) for key in bifurcation_point_kwargs if key in kwargs}
        ax = self.plot_bifurcation_points(solution_types=results['bifurcation'], x_vals=results[x],
                                          y_vals=results[y], ax=ax, **kwargs_tmp)

        # set title variable if passed
        tvar = kwargs.pop('title_var', None)
        if tvar:
            tvar_results = self.extract([tvar], cont=cont)
            tval = tvar_results[tvar][0]
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

        return line_col

    def plot_trajectory(self, variables: Union[list, tuple], cont: Union[Any, str, int], point: Union[str, int] = None,
                        ax: plt.Axes = None, force_axis_lim_update: bool = False, cutoff: float = None, **kwargs
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
            If true, the axis limits of x and y axis will be updated after creating the line plots.
        cutoff
            Initial time to be disregarded for plotting.
        kwargs
            Additional keyword arguments that allow to control the appearance of the line plot.

        Returns
        -------
        LineCollection
            Line object that was created.
        """

        # extract information from branch solutions
        try:
            results = self.extract(list(variables) + ['stability'], cont=cont, point=point)
        except KeyError:
            results = self.extract(list(variables), cont=cont, point=point)
            results['stability'] = None

        # apply cutoff, if passed
        if cutoff:
            try:
                time = self.extract(['PAR(14)'], cont=cont, point=point)['PAR(14)']
            except KeyError:
                try:
                    time = self.extract(['time'], cont=cont, point=point)['time']
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
            line_col = self._get_3d_line_collection(x=x, y=y, z=z, stability=results['stability'], **kwargs)
            ax.add_collection3d(line_col)
            ax.autoscale()

            # cosmetics
            ax.tick_params(axis='both', which='major', pad=tick_pad)
            ax.set_xlabel(variables[0], labelpad=label_pad)
            ax.set_ylabel(variables[1], labelpad=label_pad)
            ax.set_zlabel(variables[2], labelpad=label_pad)
            self._update_axis_lims(ax, [x, y, z], padding=axislim_pad, force_update=force_axis_lim_update)

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
        ax
            Axis in which to plot the data. If not provided, a new figure will be created.
        linespecs
            Keyword arguments that control the appearance of the line created for each entry in `points`.
        kwargs
            Additional keyword arguments that control the appearance of the plot.

        Returns
        -------
        plt.Axes
            Axis object that contains the plotted timeseries.
        """

        # extract information from branch solutions
        if not points:
            points = ['RG']
            points_tmp = self.results[self._results_map[cont] if type(cont) is str else cont].keys()
            results_tmp = [self.extract([var] + ['time'], cont=cont, point=p) for p in points_tmp]
            results = [{key: [] for key in results_tmp[0].keys()}]
            for r in results_tmp:
                for key in r:
                    results[0][key].append(np.squeeze(r[key]))
            for key in results[0].keys():
                results[0][key] = np.asarray(results[0][key]).squeeze()
        else:
            results = [self.extract([var] + ['time'], cont=cont, point=p) for p in points]

        # create plot
        if ax is None:
            fig, ax = plt.subplots()

        # plot phase trajectory
        if not linespecs:
            linespecs = [dict() for _ in range(len(points))]
        for i in range(len(points)):
            time = results[i]['time']
            kwargs_tmp = kwargs.copy()
            kwargs_tmp.update(linespecs[i])
            line_col = self._get_line_collection(x=time, y=results[i][var], **kwargs_tmp)
            ax.add_collection(line_col)
        ax.autoscale()
        ax.legend(points)

        return ax

    def plot_bifurcation_points(self, solution_types: DataFrame, x_vals: DataFrame, y_vals: DataFrame, ax: plt.Axes,
                                default_color: str = 'k', default_marker: str = '*', default_size: float = 50,
                                ignore: list = None, custom_bf_styles: dict = None) -> plt.Axes:
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
        plt.Axes
            Axis object that contains the plotted special solutions.
        """

        if not ignore:
            ignore = []

        # set bifurcation styles
        if custom_bf_styles:
            for key, args in custom_bf_styles.items():
                self.update_bifurcation_style(key, **args)
        bf_styles = self._bifurcation_styles.copy()
        plt.sca(ax)

        # draw bifurcation points
        for bf, x, y in zip(solution_types.values, x_vals.values, y_vals.values):
            if bf not in "EPMXRG" and bf not in ignore:
                if bf in bf_styles:
                    m = bf_styles[bf]['marker']
                    c = bf_styles[bf]['color']
                else:
                    m = default_marker
                    c = default_color
                if y.shape and np.sum(y.shape) > 1:
                    plt.scatter(x, y.max(), s=default_size, marker=m, c=c)
                    plt.scatter(x, y.min(), s=default_size, marker=m, c=c)
                else:
                    plt.scatter(x, y, s=default_size, marker=m, c=c)
        return ax

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
                        timeseries: bool, stability: bool, period: bool, eigenvals: bool, lyapunov_exp: bool
                        ) -> DataFrame:
        """Creates summary of auto continuation and stores it in dictionary.

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

        Returns
        -------
        DataFrame
        """

        columns_2d, columns_1d, indices, data_2d, data_1d = [], [], [], [], []
        add_columns = True
        for point in points:

            data_2d_tmp = []
            data_1d_tmp = []

            # get solution
            solution_type, s = self.get_solution(cont=solution, point=point)

            if solution_type != 'No Label' and solution_type != 'MX':

                indices.append(point)

                # extract variables and params from solution
                var_vals = get_solution_variables(s, variables, timeseries)
                param_vals = get_solution_params(s, params)

                # store solution type
                data_1d_tmp.append(solution_type)
                if add_columns:
                    columns_1d.append('bifurcation')

                # store parameter and variable information
                branch, _ = get_branch_info(s)
                for var, val in zip(variables, var_vals):
                    for i, v in enumerate(val):
                        data_2d_tmp.append(v)
                        if add_columns:
                            columns_2d.append((var, i))
                for param, val in zip(params, param_vals):
                    data_1d_tmp.append(val)
                    if add_columns:
                        columns_1d.append(param)

                # store time information, if requested
                if len(var_vals) > len(variables) and timeseries:
                    data_1d_tmp.append(var_vals[-1])
                    if add_columns:
                        columns_1d.append('time')

                # store stability information
                if stability:
                    data_1d_tmp.append(get_solution_stability(solution, s, point))
                    if add_columns:
                        columns_1d.append('stability')

                if period or lyapunov_exp or eigenvals:
                    p = get_solution_params(s, ['PAR(11)'])[0]

                # store information about oscillation periods
                if period:
                    data_1d_tmp.append(p)
                    if add_columns:
                        columns_1d.append('period')

                # store information about local eigenvalues/lyapunov exponents
                if eigenvals or lyapunov_exp:
                    evs = get_solution_eigenvalues(solution, branch, point)
                    if eigenvals:
                        for i, v in enumerate(evs):
                            data_2d_tmp.append(evs)
                            if add_columns:
                                columns_2d.append(('eigenvalues', i))
                    if lyapunov_exp:
                        for i, p in enumerate(get_lyapunov_exponents(evs, p)):
                            data_2d_tmp.append(p)
                            if add_columns:
                                columns_2d.append(('lyapunov_exponents', i))

            data_2d.append(data_2d_tmp)
            data_1d.append(data_1d_tmp)
            add_columns = False

        # arrange data into DataFrame
        df = self._to_dataframe(data_2d, columns=columns_2d, index=indices)
        df2 = self._to_dataframe(data_1d, columns=columns_1d, index=indices)
        for i, key in enumerate(df2.columns.values):
            df[key] = df2.loc[:, key]
        return df

    def _call_auto(self, starting_point: Union[str, int], origin: Union[Any, dict], **auto_kwargs) -> Any:
        if starting_point:
            _, s = self.get_solution(point=starting_point, cont=origin)
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
                kwargs["ICP"] = self._var_map[val]["cont"]
            elif type(val) in [list, tuple]:
                kwargs["ICP"] = [self._var_map[v]["cont"] if type(v) is str else v for v in val]
            else:
                kwargs["ICP"] = val

        # handle the user-point parameter
        if "UZR" in kwargs:
            uzr_dict = kwargs.pop("UZR")
            kwargs["UZR"] = {self._var_map[key]["cont"] if type(key) is str else key: val for key, val in uzr_dict.items()}

        return kwargs

    def _to_dataframe(self, data: list, columns: Union[list, tuple], index: list) -> DataFrame:

        # map variable/parameter keys to string-based keys
        columns_new = []
        multi_idx = False
        for c in columns:
            if type(c) is tuple:
                col = (self._var_map_inv[c[0]] if c[0] in self._var_map_inv else c[0], c[1])
                multi_idx = True
            else:
                col = self._var_map_inv[c] if c in self._var_map_inv else c
            columns_new.append(col)
        if multi_idx:
            columns = MultiIndex.from_tuples(tuple(columns_new))
        else:
            columns = columns_new

        # create dataframe
        try:
            return DataFrame(data=data, columns=columns, index=index)
        except ValueError as e:
            if len(data) > len(index):
                return DataFrame(data=data[:-1], columns=columns, index=index)
            else:
                raise e

    @staticmethod
    def _get_all_var_keys(solution):
        return [f'U({i+1})' for i in range(solution['NDIM'])]

    @staticmethod
    def _get_all_param_keys(solution):
        return solution.PAR.coordnames

    def _start_from_solution(self, solution: Any) -> Any:
        diag = str(solution[0].diagnostics)
        sol_keys = get_solution_keys(solution)
        if 'Starting direction of the free parameter(s)' in diag and len(sol_keys) == 1 and \
                "EP" in list(solution[0].labels.by_index[sol_keys[0]])[0]:
            _, s = solution[0].labels.by_index.popitem()
            solution = self._auto.run(s['EP']['solution'])
        return solution

    @staticmethod
    def _get_line_collection(x, y, stability=None, line_style_stable='solid', line_style_unstable='dotted',
                             line_color_stable='k', line_color_unstable='k', **kwargs) -> LineCollection:
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

        # combine y and param vals
        try:
            x = np.reshape(x, (x.squeeze().shape[0], 1))
        except IndexError:
            pass
        if hasattr(y[0], "shape") and sum(y[0].shape) > 1:
            y = np.asarray([y[i] for i in range(y.shape[0])])
            y_max = np.reshape(y.max(axis=1), (y.shape[0], 1))
            y_min = np.reshape(y.min(axis=1), (y.shape[0], 1))
            y_min = np.append(x, y_min, axis=1)
            y = y_max
            add_min = True
        else:
            y = np.reshape(y, (y.squeeze().shape[0], 1))
            add_min = False
        y = np.append(x, y, axis=1)

        # if stability was passed, collect indices for stable line segments
        ###################################################################

        if stability is not None and np.sum(stability.shape) > 1:

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

        colors = kwargs.pop('colors', colors)
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

        if stability is not None and np.sum(stability.shape) > 1:

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
        line_col = Line3DCollection(segments=lines, linestyles=styles, **kwargs)

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
