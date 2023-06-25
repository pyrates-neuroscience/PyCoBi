"""
Hopf Bifurcation and Limit Cycle Continuation
=============================================

Here, we will demonstrate how to perform branch switching at an `Andronov-Hopf bifurcation <http://www.scholarpedia.org/article/Andronov-Hopf_bifurcation>`_
 and continue a branch of periodic solutions.
As an example model, we use the mean-field model of a population of globally coupled quadratic integrate-and-fire (QIF)
neurons with spike-frequency adaptation (SFA). The model equations were derived in [1]_ and are given by:

.. math::

    \\tau \\dot r &= \\frac{\\Delta}{\\pi\\tau} + 2 r v, \n
    \\tau \\dot v &= v^2 +\\bar\\eta + I(t) + J r \\tau - a - (\\pi r \\tau)^2, \n
    \\tau_a \\dot a &= -a + \\alpha \\tau_a r,

where :math:`r` is the average firing rate, :math:`v` is the average membrane potential, and :math:`a` is the average
SFA of the QIF population [1]_.
The model dynamics are governed by 6 parameters:
    - :math:`\\tau` --> the population time constant,
    - :math:`\\bar \\eta` --> the mean of a Cauchy distribution over the neural excitability in the population,
    - :math:`\\Delta` --> the half-width at half maximum of the Cauchy distribution over the neural excitability,
    - :math:`J` --> the strength of the recurrent coupling inside the population,
    - :math:`\\tau_a` --> the SFA time scale,
    - :math:`\\alpha` --> the SFA rate.
This mean-field model is an exact representation of the macroscopic firing rate and membrane potential dynamics of a
spiking neural network consisting of QIF neurons with `Cauchy <https://en.wikipedia.org/wiki/Cauchy_distribution>`_
distributed background excitabilities.

It has been demonstrated that SFA induces synchronized oscillations in a population of excitatorily coupled QIF
neurons [1]_. In the example below, we will replicate these results via bifurcation analysis in `PyCoBi`.

References
^^^^^^^^^^

.. [1] R. Gast, H. Schmidt, T.R. Knösche (2020) *A Mean-Field Description of Bursting Dynamics in Spiking Neural
       Networks with Short-Term Adaptation.* Neural Computation 32 (9): 1615-1634, https://doi.org/10.1162/neco_a_01300.
"""

# %%
# Step 1: Model Initialization
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Let's begin by loading the QIF SFA mean-field model into `PyCoBi`:

from pycobi import ODESystem
import numpy as np

ode = ODESystem.from_yaml(
    "model_templates.neural_mass_models.qif.qif_sfa", auto_dir="~/PycharmProjects/auto-07p",
    node_vars={'p/qif_sfa_op/Delta': 2.0, 'p/qif_sfa_op/alpha': 1.0, 'p/qif_sfa_op/eta': 3.0},
    edge_vars=[('p/qif_sfa_op/r', 'p/qif_sfa_op/r_in', {'weight': 15.0*np.sqrt(2.0)})],
    NPR=100, NMX=30000
)

# %%
# In the code above, we adjusted some default parameters of the model in accordance with the parameters reported
# in [1]_. Also, an initial integration of the ODE system with respect to time was performed, such that the system
# was able to converge to a steady-state solution. Let's examine whether the system did indeed converge to a
# steady-state solution.

import matplotlib.pyplot as plt
ode.plot_continuation("t", "p/qif_sfa_op/r", cont=0)
plt.show()

# %%
# If the system didn't converge yet, try increasing the parameter :code:`NMX` provided to the :code:`ODESystem.from_yaml`
# method, which controls the number of time integration steps, until the system is assumed to have converged.

# %%
# Step 2: 1D parameter continuation of a steady-state solution
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# A 1D parameter continuation in :math:`\bar \eta` can be carried out as follows:

eta_sols, eta_cont = ode.run(
    origin=0, starting_point='EP', name='eta', bidirectional=True,
    ICP="p/qif_sfa_op/eta", RL0=-20.0, RL1=20.0, IPS=1, ILP=1, ISP=2, ISW=1, NTST=400,
    NCOL=4, IAD=3, IPLT=0, NBC=0, NINT=0, NMX=2000, NPR=10, MXBF=5, IID=2,
    ITMX=40, ITNW=40, NWTN=12, JAC=0, EPSL=1e-06, EPSU=1e-06, EPSS=1e-04,
    DS=1e-4, DSMIN=1e-8, DSMAX=5e-2, IADS=1, THL={}, THU={}, UZR={"p/qif_sfa_op/eta": 3.0}, STOP={}
)

# %%
# In the above code, we provided all the `auto-07p` meta parameters as keyword arguments. As an alternative, you can
# create constants file with name :code:`c.name` and provide it to :code:`ODESystem.run` via keyword argument
# :code:`c=name`. See the `auto-07p documentation <https://github.com/auto-07p/auto-07p/doc>`_ for a detailed explanation
# of each of these parameters.
# The following keyword arguments are specific to `PyCoBi` and deviate from standard `auto-07p` syntax:
# :code:`bidirectional = True` leads to automatic continuation of the solution branch into both directions of the
# continuation parameter :math:`\bar \eta`. :code:`origin` indicates the key of the solution branch, where a solution with
# the name :code:`starting_point='EP'` exists, which is used as the starting point of the parameter continuation.
# We can plot the results of the 1D parameter continuation, to examine the results.

ode.plot_continuation("p/qif_sfa_op/eta", "p/qif_sfa_op/r", cont="eta")
plt.show()

# %%
# The plot shows a 1D bifurcation diagram in the state variable :math:`r` as a function of the continuation parameter
# :math:`\bar \eta`. Solid lines indicate stable solutions of the system, whereas dotted lines indicate unstable
# solutions. We can see that the steady-state solution branch undergoes a number of bifurcations,
# marked by the green circles (`Hopf bifurcations <http://www.scholarpedia.org/article/Andronov-Hopf_bifurcation>`_)
# and grey triangles (`fold bifurcations <http://www.scholarpedia.org/article/Saddle-node_bifurcation>`_),
# at which the stability of the solution branch changes.

# %%
# Step 3: Branch switching at a Hopf bifurcation to the limit cycle solutions
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# We know from bifurcation theory that the Hopf bifurcation at :math:`\bar \eta \approx 4` is a subcritical Hopf
# bifurcation, and that an unstable limit cycle branch must emerge from that bifurcation. We can switch to this branch
# and continue the limit cycle solution curve as follows.

hopf_sols, hopf_cont = ode.run(
    origin=eta_cont, starting_point='HB2', name='eta_hopf',
    IPS=2, ISP=2, ISW=-1, UZR={"p/qif_sfa_op/eta": -2.0}
)

# %%
# Let's visualize the 1D bifurcation diagram again, this time with the limit cycle branch included:

fig, ax = plt.subplots()
ode.plot_continuation("p/qif_sfa_op/eta", "p/qif_sfa_op/r", cont="eta", ax=ax)
ode.plot_continuation("p/qif_sfa_op/eta", "p/qif_sfa_op/r", cont="eta_hopf", ax=ax)
plt.show()

# %%
# As we can see, the limit cycle solution curve that branches off the Hopf bifurcation is indeed unstable.
# However, it undergoes a fold-of-limit-cycle bifurcation, giving rise to a stable limit cycle solution branch.
# During continuation of the limit cycle branch, we set a user point at :math:`\bar \eta = 2.0` via
# :code:`UZR={4: -2.0}`. In the code below, we integrate the system equations with respect to time to show that the
# system dynamics are indeed governed by the stable limit cycle solution in that regime.

ode.run(origin=hopf_cont, starting_point="UZ1", c="ivp", name="lc")
ode.plot_continuation("t", "p/qif_sfa_op/r", cont="lc")
plt.show()

# %%
# This concludes our use example. As a final step, we clear all the temporary files we created and make sure all
# working directory changes that have happened in the background during file creation and parameter continuation are
# undone. This can be done via a single call:

ode.close_session(clear_files=True)

# %%
# If you want to keep the Fortran files for future parameter continuations etc., just set :code:`clear_files=False`.
