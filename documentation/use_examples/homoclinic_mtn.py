r"""
Saddle-Node Homoclinic Continuation: the ``mtn`` Demo
======================================================

This gallery reproduces auto-07p's ``demos/mtn`` homoclinic-bifurcation
tutorial — M. Scheffer's predator-prey model with a saddle-node homoclinic
orbit — using PyCoBi's :meth:`ODESystem.continue_homoclinic` entry point.
It's the small, well-conditioned test case for PyCoBi's HomCont support
(auto-07p Ch. 20): a 2-D system with a single saddle-node, a pre-computed
near-homoclinic orbit shipped as ``mtn.dat`` next to this script, and
sign-flip detection of the PSI(15) / PSI(16) test functions that
identify *non-central homoclinic to saddle-node* bifurcations [1]_.

The model is:

.. math::

    \dot u &= 0.5\, u\, (1 - u/K) - w \cdot \psi_Z(u) + d K, \\
    \dot w &= 0.6\, \psi_Z(u)\, w - 0.15\, w - G_f\, \psi_F(w),

with the Holling-II functional responses

.. math::

    \psi_Z(u) = \frac{0.4\, u}{0.6 + u}, \qquad
    \psi_F(w) = \frac{w^2}{w^2 + H_z^2}.

The seed orbit lives at :math:`(K, G_f, d, H_z) = (6.0, 0.0673, 0.01, 0.5)`
with the saddle equilibrium at :math:`(u^\star, w^\star) = (5.7386, 0.5108)`
and a truncation interval :math:`\text{PAR}(11) = 1046.18`.  Continuing in
:math:`(K, G_f)` while monitoring the SNIC test functions, the first run
walks the saddle-node-homoclinic curve and flags a UZ point at
:math:`K \approx 6.61,\ G_f \approx 0.069` where :math:`\psi_{15} = 0`
(matching auto-07p's mtn.1 step 26 to all significant figures).

References
^^^^^^^^^^

.. [1] M. Scheffer (1991) *Should we expect strange attractors behind
   plankton dynamics — and if so, should we bother?*, J. Plankton Res. 13,
   1291-1305; reproduced as auto-07p's ``demos/mtn`` (Doedel et al.).
.. [2] AUTO manual, Chapter 20 (Homoclinic Bifurcation Analysis), §20.5
   (the PSI test-function table) and §20.6 (the ``mtn`` walkthrough).
"""

# %%
# Step 1: build the model + ship the seed alongside
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# The model is small enough to define inline with a PyRates
# :class:`OperatorTemplate` — no YAML required.  We rename a couple of
# parameters from the auto-07p Fortran (``GF`` → ``Gf``, ``HZ`` → ``Hz``)
# because ``GF`` collides with sympy's Galois-field constructor and the
# equation parser then fails to sympify the RHS.
#
# The starting orbit lives in ``mtn.dat`` next to this script — auto-07p
# pre-computed it for the demo; we copy the file verbatim into the working
# directory before the HomCont run so :meth:`continue_homoclinic` can read
# it (``ISTART=1`` requires the seed on disk).

import shutil
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from pyrates import CircuitTemplate, NodeTemplate, OperatorTemplate

from pycobi import ODESystem

here = Path(__file__).resolve().parent
work = here  # write the c.* / .f90 / output files alongside the script
shutil.copy(here / 'mtn.dat', work / 'mtn_seed.dat')

op = OperatorTemplate(
    name='mtn_op',
    equations=[
        "u' = 0.5*u*(1.0 - u/K) - w*0.4*u/(0.6 + u) + d*K",
        "w' = 0.6*0.4*u/(0.6 + u)*w - 0.15*w - Gf*w*w/(w*w + Hz*Hz)",
    ],
    variables={
        'u': 'output(5.738626)', 'w': 'variable(0.5108401)',
        'K': 6.0, 'Gf': 0.06729762, 'd': 0.01, 'Hz': 0.5,
    },
)
node = NodeTemplate(name='mtn_node', operators=[op])
circuit = CircuitTemplate(name='mtn', nodes={'p': node})

# %%
# Step 2: instantiate the ODESystem and request the 'hom' scenario
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# ``auto_constants=('hom',)`` generates a ``c.hom`` file pre-configured for
# auto-07p's ``IPS=9`` mode with sensible HomCont defaults.  PyRates emits
# the analytical Jacobian (``analytical_jacobian=True``); HomCont consumes
# it through the same DFDU path as the equilibrium and limit-cycle
# scenarios.

ode = ODESystem.from_template(
    template=circuit,
    working_dir=str(work),
    auto_dir="~/PycharmProjects/auto-07p",
    init_cont=False,
    analytical_jacobian=True,
    auto_constants=('hom',),
)

# %%
# Step 3: continue the saddle-node homoclinic curve (auto-07p's mtn.1)
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# This is the heart of the demo.  Mode-B :meth:`continue_homoclinic`
# (omit ``origin``, pass ``dat_basename``) skips the LC seed-extraction
# step and reads ``mtn_seed.dat`` straight off disk.  The ``PAR={...}``
# override supplies the seed orbit's parameter values, including
# ``PAR(11) = 1046.18`` (the truncation interval) and ``PAR(12)``,
# ``PAR(13)`` (the saddle's equilibrium values that ``IEQUIB=2`` projects
# the unstable manifold from).
#
# ``NUNSTAB = NSTAB = -1`` tells HomCont to derive the eigenvalue split
# from ``NDIM``; for this 2-D model it picks the saddle-node-compatible
# split internally.  ``IPSI=(15, 16)`` selects the *non-central
# saddle-node homoclinic* test functions; auto-07p stores them at
# ``PAR(35)``, ``PAR(36)``.  ``UZR={35: 0.0, 36: 0.0}`` flags zero
# crossings of either as a ``UZ`` label — that's where the saddle-node
# homoclinic codim-2 phenomenon happens.

sols_fwd, _ = ode.continue_homoclinic(
    dat_basename='mtn_seed',
    ICP=['K', 'Gf'],
    NUNSTAB=-1, NSTAB=-1,
    IEQUIB=2, ITWIST=0, IPSI=(15, 16),
    name='snhom_fwd',
    PAR={11: 1046.178, 12: 5.738626, 13: 0.5108401},
    NTST=100, NCOL=4, IAD=1, ISP=0, ILP=0,
    NMX=30, NPR=10, IID=2, ITMX=10, ITNW=7, NWTN=3,
    DS=0.01, DSMIN=0.001, DSMAX=0.05,
    UZR={35: 0.0, 36: 0.0},
)

print(f"\nForward HomCont branch ({len(sols_fwd)} points)")
print(f"  bifurcations: {dict(sols_fwd['bifurcation'].value_counts())}")
sel = sols_fwd[[c for c in sols_fwd.columns
                if (isinstance(c, tuple) and c[0] in
                    ('K', 'Gf', 'PAR(35)', 'PAR(36)', 'bifurcation'))]]
print(sel.to_string())
print(
    "Reference (auto-07p mtn.1, step 26):  K = 6.61046, Gf = 0.069325, "
    "PSI(15) = 5.129e-9 → UZ"
)

# %%
# Step 4: continue in the opposite direction (auto-07p's mtn.2)
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# The reverse-direction half of the saddle-node-homoclinic curve uses a
# negative ``DS`` so auto-07p walks toward smaller ``K`` instead.  The
# resulting branch picks up the other SNIC point — where the second PSI
# test function (PSI(16) = non-central homoclinic to saddle-node in the
# unstable manifold) flips sign.  We restart from the original seed
# rather than chaining off ``snhom_fwd`` because the ``.dat`` file is the
# cleanest entry point in PyCoBi's mode-B workflow.

sols_rev, _ = ode.continue_homoclinic(
    dat_basename='mtn_seed',
    ICP=['K', 'Gf'],
    NUNSTAB=-1, NSTAB=-1,
    IEQUIB=2, ITWIST=0, IPSI=(15, 16),
    name='snhom_rev',
    PAR={11: 1046.178, 12: 5.738626, 13: 0.5108401},
    NTST=50, NCOL=4, IAD=1, ISP=0, ILP=0,
    NMX=40, NPR=10, IID=2, ITMX=10, ITNW=7, NWTN=3,
    DS=-0.01, DSMIN=0.001, DSMAX=0.05,
    UZR={35: 0.0, 36: 0.0},
)

print(f"\nReverse HomCont branch ({len(sols_rev)} points)")
print(f"  bifurcations: {dict(sols_rev['bifurcation'].value_counts())}")

# %%
# Step 5: visualise the saddle-node-homoclinic curve in :math:`(K, G_f)`
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

def _col(df, name):
    """Pull a column out of a possibly-MultiIndexed PyCoBi summary."""
    for c in df.columns:
        if (isinstance(c, tuple) and c[0] == name) or c == name:
            return df[c].to_numpy(dtype=float)
    raise KeyError(name)

fig, ax = plt.subplots(figsize=(7, 5))
for sols, lbl, col in [(sols_fwd, 'forward', '#1F77B4'),
                       (sols_rev, 'reverse', '#FF7F0E')]:
    K = _col(sols, 'K')
    Gf = _col(sols, 'Gf')
    ax.plot(K, Gf, '-', color=col, lw=1.5, label=f'SN-homoclinic ({lbl})')

# Overlay each UZ — auto-07p flags these where the PSI test functions
# pass through zero, i.e. where the homoclinic touches the saddle-node
# line (the codim-2 SNIC point).
for sols, col in [(sols_fwd, '#1F77B4'), (sols_rev, '#FF7F0E')]:
    bif = sols[('bifurcation', '')] if isinstance(sols.columns, type(sols.columns)) \
        and (('bifurcation', '') in sols.columns) else sols['bifurcation']
    K = _col(sols, 'K')
    Gf = _col(sols, 'Gf')
    for i, b in enumerate(bif):
        if str(b).strip() == 'UZ':
            ax.scatter(K[i], Gf[i], marker='*', s=200, c=col,
                       edgecolor='k', linewidth=0.8, zorder=10,
                       label='SNIC (PSI=0)' if i == 0 else None)

ax.set_xlabel('K (carrying capacity)')
ax.set_ylabel(r'$G_f$ (predator gain)')
ax.set_title('Scheffer predator-prey: saddle-node homoclinic curve')
ax.legend(loc='best')
plt.tight_layout()
plt.show()
