************************************************
Parameter continuation and bifurcation detection
************************************************

In this tutorial, you will learn how to perform a 1D
`numerical parameter continuation <http://www.scholarpedia.org/article/Numerical_analysis#Numerical_solution_of_
differential_and_integral_equations>`_ in `PyCoBi` with automatic fold `bifurcation
<http://www.scholarpedia.org/article/Bifurcation>`_ detection.
Furthermore, you will learn how to plot a simple bifurcation diagram. Throughout this example, we will use
the quadratic integrate-and-fire population model [1]_, a detailed introduction of which is given in the model
introductions example gallery. The dynamic equations of this model read the following:

.. math::
    \tau \dot r &= \frac{\Delta}{\pi\tau} + 2 r v, \\
    \tau \dot v &= v^2 +\bar \eta + I(t) + J r \tau - (\pi r \tau)^2,

where :math:`r` is the average firing rate and :math:`v` is the average membrane potential of the QIF population.
It is governed by 4 parameters:

    - :math:`\tau` --> the population time constant
    - :math:`\bar \eta` --> the mean of a Lorentzian distribution over the neural excitability in the population
    - :math:`\Delta` --> the half-width at half maximum of the Lorentzian distribution over the neural excitability
    - :math:`J` --> the strength of the recurrent coupling inside the population

In this tutorial, we will demonstrate how to (1) , (2) perform a simple 1D parameter continuation in
:math:`\bar \eta`, and (3) plot the corresponding bifurcation diagram.
The latter has also been done in [1]_, so you can compare the resulting plot with the results reported by Montbrió et
al.

References
^^^^^^^^^^
.. [1] E. Montbrió, D. Pazó, A. Roxin (2015) *Macroscopic description for networks of spiking neurons.* Physical
       Review X, 5:021028, https://doi.org/10.1103/PhysRevX.5.021028.
.. [2] E.J. Doedel, T.F. Fairgrieve, B. Sandstede, A.R. Champneys, Y.A. Kuznetsov and W. Xianjun (2007) *Auto-07p:
       Continuation and bifurcation software for ordinary differential equations.* Technical report,
       Department of Computer Science, Concordia University, Montreal, Quebec.

Part 1: Creating an ODESystem Instance
--------------------------------------

In this first part, we will be concerned with how to create a model representation that is compatible with auto-07p,
which is the software that is used for parameter continuations and bifurcation analysis in PyCoBi [2]_.
To this end, we make use of the code generation software `PyRates <>`_, which allows us to generate the Fortran files
required to run auto-07p from `YAML` model definition files.
The QIF mean-field model comes pre-implemented with `PyRates`, so we can simply point to the ready-to-use `YAML` definition file.

.. code-block::

    from pycobi import ODESystem

    # path to YAML model definition
    model = "model_templates.neural_mass_models.qif.qif"

    # installation directory of auto-07p
    auto_dir = "~/projects/auto-07p"

    # ODESystem initialization
    ode, _ = ODESystem.from_yaml(model, auto_dir=auto_dir)

If you haven't manually configured your system environment variables such that the Python installation of `auto-07p` is
recognized by your Python interpreter, you can simply add the installation directory of `auto-07p` to the :code:`ODESystem` instantiation
(as above) and it will take care of these environment variables for you.
During the initialization of :code:`ODESystem`, an integration of the QIF model over time is automatically performed in
order to ensure that the system converged to a steady-state.
This is done, because it is required for parameter continuations that the model converged to an
`equilibrium solution <http://www.scholarpedia.org/article/Equilibrium>`_, i.e. that
:math:`\dot r = 0` and :math:`\dot v = 0` in our example.
You can check whether the system indeed converged to a steady-state solution by plotting the results of the time integration:

.. code-block::

    ode.plot_continuation("PAR(14)", "U(1)", cont=0)
    plt.show()

The above code plots the first state variable :math:`r` (corresponding to :code:`"U(1)"`) against time (corresponding to :code:`"PAR(14)"`).
Alternatively, you can look at the auto-07p output in the terminal, where you will see that the values for :code:`U(1)` and :code:`U(2)`
converged to certain values. These two values represent the values of our state variables :math:`r` and
:math:`v` for the steady-state solution that the QIF system converged to.

Part 2: Performing Parameter Continuations
------------------------------------------

In this part, we will demonstrate how to perform simple 1D parameter continuations via the :code:`ODESystem.run()` method.
Since our model converged to an equilibrium and we are now save to perform the continuation in our parameter of
interest: :math:`\bar \eta`. We can do this via the :code:`ODESystem.run` method:

.. code-block::

    eta_sols, eta_cont = ode.run(
        origin=0, starting_point='EP1', name='eta', bidirectional=True,
        ICP=4, RL0=-20.0, RL1=20.0, IPS=1, ILP=1, ISP=2, ISW=1, NTST=400,
        NCOL=4, IAD=3, IPLT=0, NBC=0, NINT=0, NMX=2000, NPR=10, MXBF=5, IID=2,
        ITMX=40, ITNW=40, NWTN=12, JAC=0, EPSL=1e-06, EPSU=1e-06, EPSS=1e-04,
        DS=1e-4, DSMIN=1e-8, DSMAX=5e-2, IADS=1, THL={}, THU={}, UZR={}, STOP={}
    )

In this call, we specified the full set of auto-07p constants. Don't worry, usually, you do not have to bother with
most of them. It is common practice to specify most of them in constants files that you would pass to :code:`ODESystem.run`
via the keyword argument :code:`c=name`. In such a case, you would specify all auto-07p
constants that do not change between calls to the :code:`.run()` method in a file with the name *c.name* and only
provide the constants that need to be altered between :code:`.run()` calls directly to the :code:`.run()` method.

Checking the terminal output of auto-07p, you will realize that the output in column *TY* shows *LP* for two of the
solutions we computed along our branch in :math:`\bar\eta`. These indicate the detection of *limit point* or
`fold bifurcations <http://www.scholarpedia.org/article/Saddle-node_bifurcation>`_. We can visualize the full
bifurcation diagram via the following call:

.. code-block::

    ode.plot_continuation('PAR(4)', 'U(1)', cont='eta')
    plt.show()

The curve in this plot represents the value of :math:`r` (y-axis) at the equilibrium solutions that exist for each
value of :math:`\bar \eta` (x-axis). A solid line indicates that the equilibrium is stable, whereas a dotted line
indicated that the equilibrium is unstable. The triangles mark the points at which auto-07p detected fold
bifurcations. At a fold bifurcation, the critical eigenvalue of the vector field defined by the right-hand sides of
our model's ODEs crosses the imaginary axis (i.e. its real part changes the sign). This indicates a change of
`stability <http://www.scholarpedia.org/article/Bifurcation_diagram>`_ of an equilibrium solution, which happens
because an unstable and a stable equilibrium approach and annihilate each other. This behavior can be read from the
plot as well. The solid and the dotted line approach each other towards the fold bifurcation marks. After they meet
at the fold bifurcation, both cease to exist.
