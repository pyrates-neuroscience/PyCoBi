PyAuto
======

*PyAuto* is a Python interface to *Auto-07p* [1]. It still requires user-supplied Fortran files for parameter continuations,
but allows for a more intuitive usage of *Auto-07p* commands within Python scripts. It provides direct access to 
solutions, branches, and their properties (i.e. special solutions, eigenvalues, etc.) as well as a range of plotting 
functions to visualize bifurcation diagrams and solutions.

**REQUIREMENT:** You will have to install [auto-07p](https://github.com/auto-07p/auto-07p) on your machine and follow
these [installation instructions](https://github.com/auto-07p/auto-07p/tree/master/doc) for any of the examples below
to work.

**Use Example:** Use examples will be provided here soon. For now, have a look at the [this example](https://pyrates.readthedocs.io/en/latest/auto_analysis/continuation.html#sphx-glr-auto-analysis-continuation-py)
which demonstrates how to create the required Fortran files for *Auto-07p* via [PyRates](https://github.com/pyrates-neuroscience/PyRates)
and use them to run a 1D parameter continuation and bifurcation detection via *PyAuto*.

Installation
============

Before installation of PyAuto, it is mandatory that *Auto-07p* is installed first.

Development version (github)
----------------------------

Currently, the only supported installation option is to clone this repository and run the following line
from the directory in which the repository was cloned:
```
python setup.py install
```

References
==========
 
[1] E.J. Doedel, T.F. Fairgrieve, B. Sandstede, A.R. Champneys, Y.A. Kuznetsov and W. Xianjun (2007) *Auto-07p:
       Continuation and bifurcation software for ordinary differential equations.* Technical report,
       Department of Computer Science, Concordia University, Montreal, Quebec.
