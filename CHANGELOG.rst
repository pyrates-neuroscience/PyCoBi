Changelog
=========

0.6
---

0.6.4
~~~~~

- added possibility to call pytest scripts with auto-07p directory as additional (last) argument
- updated gitignore

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
