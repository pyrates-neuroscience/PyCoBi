*****************************
Installation and Requirements
*****************************

Prerequisites
-------------

`PyCoBi` has been build and tested for `Python >= 3.6`.
We recommend to use `Anaconda` to create a new python environment with `Python >= 3.6`.
In addition, `Auto-07p` has to be installed manually, since it is not uploaded at `PyPI` or `Anaconda` at the moment.
Below, you find standard installation instructions that should work on most systems.
If you run into problems, check out the specific installation instructions of `Auto-07p <https://github.com/auto-07p/auto-07p>`_
for your operating system. Additionally, you need to have `ninja-build` installed, which can be done via
:code:`apt install ninja-build` on Ubuntu.

Step 1:
~~~~~~~

Clone the *Auto-07p* github repository via the following command line call

:code:`git clone https://github.com/auto-07p/auto-07p`

Step 2:
~~~~~~~

Go to the directory that you cloned the *Auto-07p* repository into and run the configuration script from command line

:code:`./configure`

If this step fails (check the displayed output whether configuration was successful or not), one of the most likely
issues is that no fortran compiler is installed on your system. In that case, see the `Auto-07p documentation <https://github.com/auto-07p/auto-07p/doc>`_
for required prerequisites for your operating system (you will not need to install prerequisites for the plotting/user interface functionalities).

Step 3:
~~~~~~~

In the same directory, execute the installation of `Auto-07p` via the following command line call

:code:`make`

Step 4:
~~~~~~~

Finally, still in the same directory, install the Python package of `Auto-07p` via

:code:`python setup.py install`

It is important that you use the `Python` interpreter that you also plan to use for installation of `PyCoBi`.
After that, the installation instructions provided below should work independent of the operating system.

Dependencies
------------

`PyCoBi` has the following hard dependencies:

- `pyrates`
- `numpy`
- `matplotlib`
- `pandas`

Following the installation instructions below, these packages will be installed automatically, if not already installed within the `Python` environment you are using.

Installation
------------

`PyCoBi` can be installed via the `pip` command.  Simply run the following line from a terminal with the target Python
environment being activated:

.. code-block:: bash

   pip install pycobi


You can install optional (non-default) packages by specifying one or more options in brackets, e.g.:

.. code-block:: bash

   pip install pycobi[dev]


Currently, the only available option is `dev` (includes `pytest` and `bump2version`).

Alternatively, it is possible to clone this repository and run one of the following lines
from the directory in which the repository was cloned:

.. code-block:: bash

   python setup.py install

or

.. code-block:: bash

   pip install .[<options>]