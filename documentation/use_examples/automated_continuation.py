"""
Automated Codim-2 Search and Period-Doubling Cascades
=====================================================

`PyCoBi` ships two convenience helpers for the bookkeeping-heavy parts of a
bifurcation analysis:

* :func:`codim2_search` — given a list of codim-1 starting points (folds,
  Hopfs, period-doublings) on an existing continuation, runs a 2-parameter
  continuation of each, walks the resulting curves for codim-2 points, and
  recursively continues the codim-1 bifurcations that emerge from them.
  Supports recursive handling of zero-Hopf (ZH), generalised-Hopf (GH /
  Bautin), and Bogdanov-Takens (BT) points.
* :func:`continue_period_doubling_bf` — chases a cascade of period-doubling
  bifurcations in 2 parameters, recursing on every new PD encountered.
  Useful for tracing the boundaries of period-doubling routes to chaos.

The example below uses the bi-exponential QIF-SFA mean-field model (a QIF
population with spike-frequency adaptation implemented as a convolution
of the firing rate :math:`r` with a difference of exponentials, see
``qif_biexp_sfa.yaml`` next to this script). The interest in this model:

* At :math:`\\tau_r \\approx \\tau_d` (alpha-kernel-equivalent regime), the
  Hopf curve in the :math:`(\\bar\\eta,\\, \\tau_r)` plane carries a clear
  *generalised-Hopf* (GH) point and the interaction with the fold curve
  yields a Bogdanov-Takens (BT) scenario — ideal for :func:`codim2_search`.
* At :math:`\\tau_r \\to 0` (mono-exponential adaptation limit), the system
  enters a *period-doubling cascade* leading to chaos — ideal for
  :func:`continue_period_doubling_bf`.

References
^^^^^^^^^^

.. [1] R. Gast, H. Schmidt, T.R. Knösche (2020) *A Mean-Field Description
       of Bursting Dynamics in Spiking Neural Networks with Short-Term
       Adaptation.* Neural Computation 32 (9): 1615-1634.

.. [2] C. Scholl (2020) *Onset of Chaos in a QIF Mean-Field Model with
       Bi-Exponential Spike-Frequency Adaptation.* MPI CBS lab rotation
       report, supervised by R. Gast and T.R. Knösche.
"""

# %%
# Step 1: Load the model
# ^^^^^^^^^^^^^^^^^^^^^^
#
# The YAML next to this script defines the four-dimensional bi-exponential
# QIF-SFA system. Adaptation rate ``alpha = 0.05`` and rise/decay times
# ``tau_r = tau_d = 10`` recover the alpha-kernel limit (no period doubling,
# clean BT/GH structure). We start the IVP at :math:`\bar\eta = -3` —
# safely inside the steady-state region — and continue down into the
# oscillatory regime.

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from pycobi import ODESystem
from pycobi.automated_continuation import (
    codim2_search,
    continue_period_doubling_bf,
)

here = Path(__file__).resolve().parent
yaml_path = str(here / 'qif_biexp_sfa' / 'qif_biexp_sfa')

ode = ODESystem.from_yaml(
    yaml_path,
    auto_dir="~/PycharmProjects/auto-07p",
    node_vars={
        'p/qif_biexp_sfa_op/eta': -3.0,
        'p/qif_biexp_sfa_op/alpha': 0.05,
        'p/qif_biexp_sfa_op/Delta': 2.0,
        'p/qif_biexp_sfa_op/tau_r': 10.0,
        'p/qif_biexp_sfa_op/tau_d': 10.0,
    },
    edge_vars=[('p/qif_biexp_sfa_op/r', 'p/qif_biexp_sfa_op/r_in',
                {'weight': 15.0 * np.sqrt(2.0)})],
    init_cont=True, NPR=500, NMX=20000,
)

# %%
# Step 2: Locate codim-1 starting points
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# A bidirectional equilibrium continuation in :math:`\bar\eta` sweeps
# through the fold and Hopf bifurcations that bound the bursting regime.

eta_sols, eta_cont = ode.run(
    starting_point='EP2', name='eta_branch',
    ICP='p/qif_biexp_sfa_op/eta', bidirectional=True,
    RL0=-10.0, RL1=0.0,
    IPS=1, ILP=1, ISP=2, ISW=1, NTST=200, NCOL=4,
    NMX=500, NPR=20, DS=1e-3, DSMIN=1e-7, DSMAX=5e-2,
    ITMX=20, ITNW=20, NWTN=12,
)
print("bifurcations on the eta branch:")
print(eta_sols['bifurcation'].value_counts())

# %%
# The branch should report two Hopf points (HB1, HB2) flanking the
# oscillatory regime and two fold bifurcations (LP1, LP2). Quick look at
# the resulting 1D bifurcation diagram:

fig, ax = plt.subplots(figsize=(6, 4))
ode.plot_continuation('p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/r',
                       cont='eta_branch', ax=ax)
ax.set_title(r'1D continuation in $\bar\eta$ (steady state)')
plt.tight_layout()
plt.show()

# %%
# Step 3: Codim-2 BT/GH scenario in :math:`(\bar\eta,\, \tau_r)`
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# :func:`codim2_search` takes a list of codim-1 labels and runs a
# 2-parameter continuation of each. We start from the Hopf ``HB1`` we just
# found and continue it as a function of :math:`\bar\eta` and the
# adaptation rise time :math:`\tau_r`. ``codim2_types=('GH', 'BT')`` tells
# the helper to recurse on the codim-2 points expected for this model
# (the GH point on the Hopf curve and the BT point at the Hopf/fold
# intersection); ZH recursion is suppressed for clarity.

hopf_curves = codim2_search(
    params=['p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/tau_r'],
    starting_points=['HB1'],
    origin=eta_cont,
    pyauto_instance=ode,
    max_recursion_depth=1,
    codim2_types=('GH', 'BT'),
    RL0=-10.0, RL1=0.0,
    NMX=300, NPR=20,
    DS=1e-3, DSMIN=1e-7, DSMAX=5e-2,
    name='hopf_codim2',
)
print(f"\nhopf_codim2 continuations: {list(hopf_curves.keys())}")

# %%
# The returned dict maps continuation name → continuation key. ``hopf_codim2:HB1``
# is the 2-parameter Hopf curve itself; entries with ``/GH{i}`` or ``/BT{i}``
# suffixes are the sub-continuations the helper ran at each codim-2 point.
#
# At each GH point the helper switches to limit-cycle continuation in
# ``[params[0], PAR(11)]`` (PAR(11) being the period) to expose the change
# of Hopf criticality (super- vs sub-critical). At each BT point it runs a
# 1D equilibrium continuation stopping at the nearest Hopf, then continues
# that Hopf curve in 2 parameters. The homoclinic curve emerging from BT
# is *not* automatically followed — that requires auto-07p's HomCont
# package with ``IPS=9``, which is left as a manual follow-up if needed.
#
# We can also extend the search to the fold curve so the BT point at the
# Hopf/fold intersection is captured from both sides:

fold_curves = codim2_search(
    params=['p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/tau_r'],
    starting_points=['LP1'],
    origin=eta_cont,
    pyauto_instance=ode,
    max_recursion_depth=1,
    codim2_types=('BT',),  # GH is a Hopf-only phenomenon
    RL0=-10.0, RL1=0.0,
    NMX=300, NPR=20,
    DS=1e-3, DSMIN=1e-7, DSMAX=5e-2,
    name='fold_codim2',
)
print(f"fold_codim2 continuations: {list(fold_curves.keys())}")

# %%
# Step 4: Plot the codim-2 bifurcation diagram
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Overlay every continuation that was recorded. Stable vs unstable segments
# render via the standard convention (solid black for stable, gray dotted
# for unstable — new defaults in PyCoBi 1.0.0). Bifurcation markers along
# each curve label the codim-2 points the search may have triggered
# recursion on.

fig, ax = plt.subplots(figsize=(7, 5))
for name in hopf_curves:
    ode.plot_continuation(
        'p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/tau_r',
        cont=name, ax=ax,
    )
for name in fold_curves:
    ode.plot_continuation(
        'p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/tau_r',
        cont=name, ax=ax,
    )
ax.set_title(r'codim-2 bifurcation diagram in $(\bar\eta,\, \tau_r)$')
ax.set_xlabel(r'$\bar\eta$')
ax.set_ylabel(r'$\tau_r$')
plt.tight_layout()
plt.show()

# %%
# At small :math:`\tau_r` (bottom of the diagram), the Hopf curve loses
# stability and the limit cycle born at it should eventually undergo
# period-doubling. Let's track that explicitly in the next step.

# %%
# Step 5: Period-doubling cascades — `continue_period_doubling_bf`
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# The companion helper :func:`continue_period_doubling_bf` chases a cascade
# of period-doubling (PD) bifurcations in two parameters. It expects the
# solution dict from a limit-cycle continuation that already contains at
# least one PD point. Whether a model exhibits PD bifurcations is regime-
# specific — for the bi-exponential QIF-SFA used here the lab report
# documents a clean cascade at small :math:`\tau_r` (the mono-exponential
# adaptation limit), but reproducing it cleanly requires careful tuning of
# the limit-cycle continuation. We outline the recipe here without
# committing to numerical reproduction:
#
# .. code-block:: python
#
#     # 1. From a Hopf on the equilibrium branch, switch to the limit-
#     #    cycle family and continue it bidirectionally in eta:
#     lc_eta_sols, lc_eta_cont = ode.run(
#         origin=eta_cont, starting_point='HB1', name='lc_eta',
#         IPS=2, ISW=-1, ISP=2,
#         ICP='p/qif_biexp_sfa_op/eta',
#         UZR={4: -5.2},                # PAR index for eta (cf. c.* parnames)
#         NMX=400, NPR=10,
#         bidirectional=True, get_period=True,
#     )
#
#     # 2. From the labelled user point on the LC, continue the limit cycle
#     #    in tau_r down towards the mono-exponential limit:
#     lc_tau_sols, lc_tau_cont = ode.run(
#         origin=lc_eta_cont, starting_point='UZ1', name='lc_tau_r',
#         IPS=2, ISW=1, ISP=2,
#         ICP=['p/qif_biexp_sfa_op/tau_r', 11],  # PAR(11) = period
#         RL0=0.01, RL1=10.0,
#         NMX=400, NPR=5, DS=-1e-2, DSMAX=1e-1,
#     )
#
#     # 3. Chase the PD cascade in (eta, tau_r):
#     pd_names, _ = continue_period_doubling_bf(
#         solution=ode.results[ode.get_continuation('lc_tau_r').key],
#         continuation=lc_tau_cont,
#         pyauto_instance=ode,
#         max_iter=3,        # max recursion depth (new semantics in 1.0.0)
#         precision=3,       # decimals used to dedupe revisited PD points
#         ICP=['p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/tau_r'],
#         IPS=2, ISW=2, ISP=2,
#         NMX=200, NPR=20,
#     )
#
# ``max_iter`` was previously an iteration counter that only fired once on
# entry; it now bounds the *recursion depth* of the cascade — the variable
# that actually prevents runaway recursion on chaotic models. Names of
# successive continuations follow the ``pd_d{depth}_n{i}`` convention so
# parallel sub-cascades don't collide on shared labels.
#
# Sub-runs that fail inside :func:`ODESystem.run` are surfaced as
# :class:`UserWarning` and the cascade continues with the remaining PD
# points rather than aborting — useful when chasing a cascade that
# eventually enters a chaotic regime where individual continuations may
# fail to converge.

# %%
# Step 6: Failure modes and what to expect on unfamiliar models
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Both helpers wrap a series of nested ``ODESystem.run`` calls. If any
# individual sub-run raises inside auto-07p (a common occurrence for
# low-quality parameter regimes or aggressive step sizes), the failure
# surfaces as a :class:`UserWarning` and the search continues with the
# remaining starting points rather than aborting. Read the warnings
# carefully: they cite the auto-07p exception type and message, the
# sub-run label that failed, and the kwargs hook (``kwargs_1D_lc_cont``,
# ``kwargs_2D_cont``, ...) you can use to override the default constants
# for that path.
#
# A reasonable rule of thumb: if more than a few sub-runs fail, the
# parameter regime is probably unsuitable for the default solver settings.
# Either tighten ``DS`` / ``DSMAX`` via the relevant ``kwargs_*`` hook, or
# break the search up into smaller scopes (narrower ``codim2_types``,
# explicit ``starting_points``) so the model-specific tuning lives in the
# user code rather than buried inside the helper.

# %%
# Step 7: Clean up
# ^^^^^^^^^^^^^^^^

ode.close_session(clear_files=True)
