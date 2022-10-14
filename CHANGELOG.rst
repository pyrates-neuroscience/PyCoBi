Changelog
=========

0.7
---

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
