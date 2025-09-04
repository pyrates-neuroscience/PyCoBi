"""
1D and 2D Parameter Continuations
=================================

In this tutorial, you will learn how to perform 1D and 2D
`numerical parameter continuations <http://www.scholarpedia.org/article/Numerical_analysis#Numerical_solution_of_
differential_and_integral_equations>`_ in `PyCoBi` with automatic fold `bifurcation
<http://www.scholarpedia.org/article/Bifurcation>`_ detection.
Furthermore, you will learn how to plot a simple bifurcation diagram. As an example model, we use the
mean-field model of a population of globally coupled Izhikevich (IK) neurons.
The model equations were derived in [1]_ and are given by:

.. math::

    C \\dot r &= \\frac{\\Delta k^2 (v-v_r)}{\\pi C} + r (k (2v - v_r - v_t) - g r), \n
    C \\dot v &= k (v - v_r) (v - v_t) + \\eta + I(t) + g r (E - v) - u - \\pi C r (\\Delta + \\frac{\\pi C r}{k}), \n
    \\tau_u \\dot u &= b(v-v_r) - u + d \\tau_u r

where :math:`r` is the average firing rate, :math:`v` is the average membrane potential, and :math:`u` is a global
recovery variable [1]_. This mean-field model approximates the macroscopic dynamics of a spiking neural network
consisting of IK neurons with `Cauchy <https://en.wikipedia.org/wiki/Cauchy_distribution>`_
distributed spike thresholds, with center :math:`v_r` and width :math:`\\Delta`.

It has been demonstrated that a QIF population with excitatory coupling expresses a bi-stable regime for sufficiently
 strong coupling :math:`g` and weak heterogeneity :math:`\\Delta` [1]_.
In the example below, we will replicate these results via bifurcation analysis in `PyCoBi`.

References
^^^^^^^^^^

.. [1] R. Gast, Sara A. Solla, A. Kennedy (2023) *Macroscopic dynamics of neural networks with heterogeneous spiking
       thresholds.* Physical Review E, 107 (2): 024306, https://doi.org/10.1103/PhysRevE.107.024306.
.. [2] R. Gast, Sara A. Solla, A. Kennedy (2024) *Neural heterogeneity controls computations in spiking neural networks.*
       PNAS, 121 (3): e2311885121, https://doi.org/10.1073/pnas.2311885121.
"""

from pycobi import ODESystem
import numpy as np
import matplotlib.pyplot as plt

# %%
# Step 1: Model Initialization
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# In this first part, we will be concerned with how to create a model representation that is compatible with auto-07p,
# which is the software that is used for parameter continuations and bifurcation analysis in PyCoBi.
# To this end, we make use of the code generation software `PyRates <https://github.com/pyrates-neuroscience/PyRates>`_,
# which allows us to generate the Fortran files required to run auto-07p from `YAML` model definition files.
# The IK mean-field model comes pre-implemented with `PyRates`, so we can simply point to the ready-to-use `YAML` definition file.

model = "model_templates.neural_mass_models.ik.ik_theta"
ode = ODESystem.from_yaml(model, auto_dir="~/PycharmProjects/auto-07p", NMX=10000, DSMAX=0.1,
                          node_vars={"p/ik_theta_op/Delta": 0.1, "p/ik_theta_op/g": 20.0, "p/ik_theta_op/d": 0.0})

# %%
# If you haven't manually configured your system environment variables such that the Python installation of `auto-07p` is
# recognized by your Python interpreter, you can simply add the installation directory of `auto-07p` to the :code:`ODESystem` instantiation
# (as above) and it will take care of these environment variables for you.
# During the initialization of :code:`ODESystem`, an integration of the IK model over time is automatically performed in
# order to ensure that the system converged to a steady-state.
# This is done, because it is required for parameter continuations that the model converged to an
# `equilibrium solution <http://www.scholarpedia.org/article/Equilibrium>`_, i.e. that
# :math:`\dot r = 0`, :math:`\dot v = 0`, and :math:`\dot u = 0` in our example.
# The keyword arguments :code:`NMX` and :code:`DSMAX` control the maximum number of integration steps and maximum
# integration step-size, respectively. Finally, the argument :code:`node_vars` allows you to change some of the model parameters
# before the start of the continuation, to start in a desirable dynamic regime.
# You can check whether the system indeed converged to a steady-state solution by plotting the results of the time integration:

ode.plot_continuation("PAR(14)", "U(1)", cont=0)
plt.show()

# %%
# The above code plots the first state variable :math:`r` (corresponding to :code:`"U(1)"`) against time (corresponding to :code:`"PAR(14)"`).
# Alternatively, you can look at the auto-07p output in the terminal, where you will see that the values for :code:`U(1)`, :code:`U(2)`, etc.
# converged to certain values. These values represent the values of our state variables for the
# steady-state solution that the IK system converged to.
#
# Part 2: Performing 1D Parameter Continuations
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# In this part, we will demonstrate how to perform simple 1D parameter continuations via the :code:`ODESystem.run()` method.
# Since our model converged to an equilibrium and we are now save to perform the continuation in our parameter of
# interest: :math:`\bar \eta`. We can do this via the :code:`ODESystem.run` method:

sols_1d, cont_1d = ode.run(
        origin=0, starting_point='EP2', name='eta', bidirectional=True,
        ICP=9, RL0=-10.0, RL1=400.0, IPS=1, ILP=1, ISP=2, ISW=1, NTST=400,
        NCOL=4, IAD=3, IPLT=0, NBC=0, NINT=0, NMX=4000, NPR=10, MXBF=5, IID=2,
        ITMX=40, ITNW=40, NWTN=12, JAC=0, EPSL=1e-06, EPSU=1e-06, EPSS=1e-04,
        DS=1e-4, DSMIN=1e-8, DSMAX=0.1, IADS=1, THL={}, THU={}, UZR={}, STOP={}
    )

# %%
# In this call, we specified the full set of auto-07p constants. Don't worry, usually, you do not have to bother with
# most of them. It is common practice to specify most of them in constants files that you would pass to :code:`ODESystem.run`
# via the keyword argument :code:`c=file_name`. In such a case, you would specify all `auto-07p`
# constants that do not change between calls to the :code:`.run()` method in a file with the name *c.name* and only
# provide the constants that need to be altered between :code:`.run()` calls directly to the :code:`.run()` method.
#
# Checking the terminal output of auto-07p, you will realize that the output in column *TY* shows *LP* for two of the
# solutions we computed along our branch in :math:`\bar\eta`. These indicate the detection of *limit point* or
# `fold bifurcations <http://www.scholarpedia.org/article/Saddle-node_bifurcation>`_. We can visualize the full
# bifurcation diagram via the following call:

ode.plot_continuation('PAR(9)', 'U(1)', cont='eta')
plt.show()

# %%
# The curve in this plot represents the value of :math:`r` (y-axis) at the equilibrium solutions that exist for each
# value of :math:`\eta` (x-axis). A solid line indicates that the equilibrium is stable, whereas a dotted line
# indicated that the equilibrium is unstable. The triangles mark the points at which auto-07p detected
# `fold bifurcations <http://www.scholarpedia.org/article/Saddle-node_bifurcation>`_.
# At a fold bifurcation, the critical eigenvalue of the vector field defined by the right-hand sides of
# our model's ODEs crosses the imaginary axis (i.e. its real part changes the sign). This indicates a change of
# `stability <http://www.scholarpedia.org/article/Stability>`_ of an equilibrium solution, which happens
# because an unstable and a stable equilibrium approach and annihilate each other. This behavior can be read from the
# plot as well. The solid and the dotted line approach each other towards the fold bifurcation marks. After they meet
# at the fold bifurcation, both cease to exist.
#
# Part 3: Performing 2D Parameter Continuations
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Any special solution, such as bifurcation points, can be continued in two model parameters that are simultaneously varied.
# This procedure results in a curve of codimension one bifurcations in 2D parameter space. Solutions along such curves
# can undergo bifurcations themselves, which are then of codimension two. An overview over these concepts and their
# relationships and cross-links can be found `here <http://www.scholarpedia.org/article/Bifurcation>`_.
# Below, we compute the 2D curves of the fold bifurcations detected in Part 2 in :math:`\eta` and :math:`\Delta`:

ode.run(
    origin=cont_1d, starting_point='LP1', name='Delta/eta:lp1', bidirectional=True,
    ICP=[5, 9], RL0=0.0, RL1=5.0, IPS=1, ILP=0, ISP=2, ISW=2, NTST=400,
    NCOL=4, IAD=3, IPLT=0, NBC=0, NINT=0, NMX=4000, NPR=10, MXBF=5, IID=2,
    ITMX=40, ITNW=40, NWTN=12, JAC=0, EPSL=1e-06, EPSU=1e-06, EPSS=1e-04,
    DS=1e-4, DSMIN=1e-8, DSMAX=0.01, IADS=1, THL={}, THU={}, UZR={}, STOP={}
)
ode.run(
    origin=cont_1d, starting_point='LP2', name='Delta/eta:lp2', bidirectional=True,
    ICP=[5, 9], RL0=0.0, RL1=5.0, IPS=1, ILP=0, ISP=2, ISW=2, NTST=400,
    NCOL=4, IAD=3, IPLT=0, NBC=0, NINT=0, NMX=4000, NPR=10, MXBF=5, IID=2,
    ITMX=40, ITNW=40, NWTN=12, JAC=0, EPSL=1e-06, EPSU=1e-06, EPSS=1e-04,
    DS=1e-4, DSMIN=1e-8, DSMAX=0.01, IADS=1, THL={}, THU={}, UZR={}, STOP={}
)

# %%
# As you can see, we ran the 2D parameter continuation for each of the two fold bifurcations we detected in Part 2.
# Importantly, we set :code:`ILP=0`, because fold bifurcation detection is only relevant for 1D parameter continuations.
# Also, we set :code:`ISW=2` to indicate that this is a two-parameter continuation problem, and set :code:`ICP=[5, 9]` to
# indicate the two parameters to treat as free variables. The continuation results can be visualized as follows:

fig, ax = plt.subplots(figsize=(6, 6))
ode.plot_continuation('PAR(9)', 'PAR(5)', cont='Delta/eta:lp1', ax=ax)
ode.plot_continuation('PAR(9)', 'PAR(5)', cont='Delta/eta:lp2', ax=ax)
plt.show()

# %%
# As can be seen, the two fold bifurcations exist for a range of parameters. As :math:`\Delta` increases, the fold
# bifurcations approach each other, until they eventually annihilate each other in a
# `cusp bifurcation <http://www.scholarpedia.org/article/Cusp_bifurcation>`_. This, increasing the heterogeneity in the
# IK network via :math:`\Delta` decreases the non-linearity of the steady-state solution curve that we plotted as a
# function of :math:`\eta` in Part 2. For more results on the effects of heterogeneity on IK networks, see [2]_.

ode.close_session(clear_files=True)
