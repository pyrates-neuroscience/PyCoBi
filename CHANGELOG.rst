Changelog
=========

0.8
---

0.8.4
~~~~~

- updates to `ODESystem.extract` and the plotting functions that use it: Users can now flexibly switch between using the parameter/variable indices when specifying the variables to plot, or whether they want to use the pyrates-style naming scheme
- added option to `ODESystem.run` to store only the minimum and maximum of each selected state variable for periodic solutions (set keyword argument `reduce_limit_cycle=True`)

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
