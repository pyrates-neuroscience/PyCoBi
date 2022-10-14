*******************
General Information
*******************

Basic features of PyCoBi
-------------------------

- perform parameter continuations with automated bifurcation detection in Python
- supports first-order ordinary differential equation (ODE) systems
- wrapper to `Auto-07p <https://github.com/auto-07p/auto-07p>`_, i.e. comes with all the functionalities and meta parameters of :code:`Auto-07p`
- create your own fortran files (i.e. use the standard :code:`Auto-07p` user files) or generate Fortran files automatically via `PyRates <https://github.com/pyrates-neuroscience/PyRates>`_
- automatic Fortran file generation allows to define ODE systems via YAML files or in pure Python
- automatic merging of solution branches and other convenience features
- various plotting functionalities
- functionalities to save and load results

Reference
---------

If you use `PyCoBi`, please cite the most recent release:

Contact
-------

If you have questions, problems or suggestions regarding PyCoBi, please contact `Richard Gast <https://www.richardgast.me>`_.

Contribute
----------

PyCoBi is an open-source project that everyone is welcome to contribute to. Check out our `GitHub repository <https://github.com/pyrates-neuroscience/PyCoBi>`_
for all the source code, open issues etc. and send us a pull request, if you would like to contribute something to our software.

Useful links
------------

`PyCoBi` makes use of two essential Python tools:

- Frontend: `PyRates <https://github.com/pyrates-neuroscience/PyRates>`_
- Backend: `Auto-07p <https://github.com/auto-07p/auto-07p>`_

Each of these two tools comes with an extensive documentation that is complementary to the content covered on this documentation website.
