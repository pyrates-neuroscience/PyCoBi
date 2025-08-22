"""
Bifurcation analysis of the Jansen-Rit Model
============================================

Jansen-Rit model, a neural mass model of the dynamic interactions between 3 populations:

    - pyramidal cells (PCs)
    - excitatory interneurons (EINs)
    - inhibitory interneurons (IINs)

Originally, the model has been developed to describe the waxing-and-waning of EEG activity in the alpha frequency range
(8-12 Hz) in the visual cortex [1]_. In the past years, however, it has been used as a generic model to describe the
macroscopic electrophysiological activity within a cortical column [2]_. Membrane potential deflections that are caused
at the somata of a neural population by synaptic input are modeled by a convolution of the input with an alpha
kernel. This choice has been shown to reflect the dynamic process of polarization propagation from the synapse via the
dendritic tree to the soma [3]_. The convolution operation can be expressed via a second-order differential equation:

.. math::
        \\dot V &= I, \n
        \\dot I &= \\frac{H}{\\tau} m_{in} - \\frac{2 I}{\\tau} - \\frac{V}{\\tau^2},

where :math:`V` represents the average post-synaptic potential and :math:`H` and :math:`\\tau` are the efficacy and
the timescale of the synapse, respectively. As a second operator, the translation of the average membrane potential
deflection at the soma to the average firing of the population is given by a sigmoidal function:

.. math::
        m_{out} = S(V) = \\frac{m_{max}}{1 + e^{(r (V_{thr} - V))}}.

In this equation, :math:`m_{out}` and :math:`V` represent the average firing rate and membrane potential, respectively,
while :math:`m_{max}`, :math:`r` and :math:`V_{thr}` are constants defining the maximum firing rate, firing threshold
variance and average firing threshold within the modeled population, respectively.

By using the linearity of the convolution operation, the dynamic interactions between PCs, EINs and IINs can be
expressed via 6 coupled ordinary differential equations that are composed of the two operators defined above:

.. math::

        \\dot V_{pce} &= I_{pce}, \n
        \\dot I_{pce} &= \\frac{H_e}{\\tau_e} c_4 S(c_3 V_{in}) - \\frac{2 I_{pce}}{\\tau_e} - \\frac{V_{pce}}{\\tau_e^2}, \n
        \\dot V_{pci} &= I_{pci}, \n
        \\dot I_{pci} &= \\frac{H_i}{\\tau_i} c_2 S(c_1 V_{in}) - \\frac{2 I_{pci}}{\\tau_i} - \\frac{V_{pci}}{\\tau_i^2}, \n
        \\dot V_{in} &= I_{in}, \n
        \\dot I_{in} &= \\frac{H_e}{\\tau_e} S(V_{pce} - V_{pci}) - \\frac{2 I_{in}}{\\tau_e} - \\frac{V_{in}}{\\tau_e^2},

where :math:`V_{pce}`, :math:`V_{pci}`, :math:`V_{in}` are used to represent the average membrane potential deflection
caused by the excitatory synapses at the PC population, the inhibitory synapses at the PC population, and the excitatory
synapses at both interneuron populations, respectively.

Below, we will reproduce the 1D bifurcation diagram for this model that is depicted in Fig.2 of [2]_.

**References**

.. [1] B.H. Jansen & V.G. Rit (1995) *Electroencephalogram and visual evoked potential generation in a mathematical
       model of coupled cortical columns.* Biological Cybernetics, 73(4): 357-366.

.. [2] A. Spiegler, S.J. Kiebel, F.M. Atay, T.R. Knösche (2010) *Bifurcation analysis of neural mass models: Impact of
       extrinsic inputs and dendritic time constants.* NeuroImage, 52(3): 1041-1058,
       https://doi.org/10.1016/j.neuroimage.2009.12.081.

.. [3] P.A. Robinson, C.J. Rennie, J.J. Wright (1997) *Propagation and stability of waves of electrical activity in the
       cerebral cortex.* Physical Review E, 56(826), https://doi.org/10.1103/PhysRevE.56.826.
"""

# %% Libraries import:
from pycobi import ODESystem
import matplotlib.pyplot as plt
from pyrates.frontend import CircuitTemplate

# %% Initialize the Jansen-Rit model in PyRates:
jrc = CircuitTemplate.from_yaml("model_templates.neural_mass_models.jansenrit.JRC")

# update input variable to start in a high-activity state
jrc.update_var(node_vars={'pc/rpo_e_in/u': 400.0})

# set state variables close to steady-state value
jrc.update_var(node_vars={'ein/rpo_e/v': 1.762243e-02 }) # y0
jrc.update_var(node_vars={'pc/rpo_e_in/v': 3.052549e-02}) # y1
jrc.update_var(node_vars={'pc/rpo_i/v': -2.195986e-02})
jrc.update_var(node_vars={'iin/rpo_e/v': 4.405607e-03 })

# %% PyCobi model initialization:
jrc_auto = ODESystem.from_template(jrc, auto_dir="~/PycharmProjects/auto-07p", init_cont=False)

# %% Time continuation in PyCobi:
#   - `DS` defines the initial step-size of the time continuation (in ms)
#   - `DSMIN` defines the minimal step-size of the time continuation (in ms)
#   - `DSMAX` defines the maximal step-size of the time continuation (in ms)
#   - `NMX` defines the maximum number of continuation steps to perform
#   - `UZR={14: 1000.0}` tells auto-07p to create a user-specified marker when the parameter 14, which is the
#     default parameter field in auto-07p in which time is stored, reaches a value of 1000.0 (ms)
#   - `STOP={'UZ1'}` tells auto-07p to stop the continuation ones it hits the first user-specified marker

t_sols, t_cont = jrc_auto.run(
    c='ivp', name='time', DS=1e-4, DSMIN=1e-10, EPSL=1e-08, EPSU=1e-08, EPSS=1e-06,
    DSMAX=1e-3, NMX=50000, UZR={14: 2.0}, STOP={'UZ1'})

# %% Outputs of the time continuation:
#   - t_sols: pandas Dataframe with the summary of the results of the simulation
#   - t_cont: type 'branch' or 'bifDiag', an auto-07p object that can be used for subsequent parameter continuations

# %% plotting the solutions:
v_pce = t_sols["pc/rpo_e_in/v"]
v_pci = t_sols["pc/rpo_i/v"]
pc = v_pce + v_pci
t = t_sols["t"]
plt.plot(t,v_pce,label='V_PCE')
plt.plot(t,v_pci,label='V_PCI')
plt.plot(t,pc,label='PC')
plt.title("Initial state of Jansen-Rit model")
plt.legend()
plt.show()

# %% BIFURCATION ANALYSIS
# 1D parameter continuation in the input parameter 'u':

u_sols, u_cont = jrc_auto.run(
    origin=t_cont, 
    starting_point='UZ1', 
    name='u', # Name of continuation parameter (variable to vary) 
    bidirectional=True, #  Continue both forward and backward from starting point
    ICP='pc/rpo_e_in/u',
    RL0=-50.0, RL1=600.0, # Lower and upper bounds for continuation parameter(s)
    IPS=1, # Problem type: 1=Equilibrium, 2=Periodic orbit continuation
    ILP=1, # Detect limit points (folds): 1=Yes, 0=No 
    ISP=2, # Bifurcation detection level: 0=none, 1=some, 2=more 
    ISW=1, 
    NTST=400,
    NCOL=4, # recommended in the auto docs
    IAD=3, # recommended in the auto docs
    IPLT=0, 
    NBC=0, 
    NINT=0, 
    NMX=8000, 
    NPR=20, # Print/report frequency (every NPR steps)
    MXBF={}, # Maximum number of bifurcations to locate
    IID=2,
    ITMX=1000, # max.num. of iterations allowed in the accurate location of special solutions
    ITNW=40, 
    NWTN=12, 
    JAC=0, 
    EPSL=1e-07, # recommended in the auto docs
    EPSU=1e-07, # recommended in the auto docs
    EPSS=1e-05, # recommended in the auto docs (approx.100/1000 times EPSL,EPSU)
    DS=1e-03, # initial continuation step size (0.5)
    DSMIN=1e-8, # minimum allowed step size
    DSMAX=0.5, # maximum allowed step size (4)
    IADS=1, 
    THL={}, 
    THU={}, 
    UZR={}, 
    STOP={}
)

# plot the continuation results
jrc_auto.plot_continuation('pc/rpo_e_in/u', 'pc/rpo_e_in/v', cont='u')
plt.show()

# %% From the 1D continuation in 'u', we see that the steady-state solution undergoes 3 Hopf bifurcations and 2 Fold
# bifurcations as we decreased the input parameter 'u' from 'u=400'. Below, we switch onto the limit cycle branches
# emerging from each Hopf bifurcation and continue the limit cycles in 'u':

for i in range(1,4):
    jrc_auto.run(origin=u_cont, starting_point=f'HB{i}', name=f'u_lc{i}', IPS=2, ISP=2, ISW=-1, STOP=['BP2'])

# %% Now, we plot the full 1D bifurcation diagram, including the limit cycle continuations:

fig, ax = plt.subplots()
jrc_auto.plot_continuation('pc/rpo_e_in/u', 'pc/rpo_e_in/v', cont='u', ax=ax)
for i in range(1,4):
    jrc_auto.plot_continuation('pc/rpo_e_in/u', 'pc/rpo_e_in/v', cont=f'u_lc{i}', ax=ax, ignore=["UZ", "BP"],
                               line_color_stable="green")
plt.show()

# %% We find the same bistable oscillatory regime that is depicted in Fig.2 of [2]_, where a large-amplitude limit
# cycle co-exists with a small-amplitude limit cycle.
#
# This concludes our use example. As a final step, we clear all the temporary files we created and make sure all
# working directory changes that have happened in the background during file creation and parameter continuation are
# undone. This can be done via a single call:

jrc.clear()
jrc_auto.close_session(clear_files=True)
