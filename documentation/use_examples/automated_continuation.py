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

This example uses the QIF mean-field model with bi-exponential
spike-frequency adaptation (``qif_biexp_sfa.yaml`` next to this script).
The bi-exponential kernel is a strict generalisation of the alpha-kernel
QIF-SFA from :ref:`Hopf Bifurcation and Limit Cycle Continuation`:
both kernels satisfy
:math:`\\tau_a^2 A'' + 2 \\tau_a A' + A = \\alpha r \\tau_a`
when ``tau_r == tau_d == tau_a``, so the codim-2 structure in the
:math:`(\\bar\\eta,\\, \\Delta)` plane (generalised-Hopf, Bogdanov-Takens,
cusps) is identical to the alpha-kernel case at that parameter point —
verified by side-by-side bifurcation analysis. Picking the bi-exponential
model now means we can also explore the period-doubling regime that opens
up when ``tau_r`` is taken much smaller than ``tau_d`` (Section 2.1 below).

References
^^^^^^^^^^

.. [1] R. Gast, H. Schmidt, T.R. Knösche (2020) *A Mean-Field Description
       of Bursting Dynamics in Spiking Neural Networks with Short-Term
       Adaptation.* Neural Computation 32 (9): 1615-1634.
"""

# %%
# Step 1: Load the model
# ^^^^^^^^^^^^^^^^^^^^^^
#
# Bi-exponential QIF-SFA from the co-located YAML, with rise and decay
# adaptation time constants both set to 10 — the alpha-kernel-equivalent
# parameter point. Coupling ``J = 15 sqrt(2)``, ``alpha = 1.0``,
# ``Delta = 2.0`` to match the QIF-SFA example.

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
        'p/qif_biexp_sfa_op/Delta': 2.0,
        'p/qif_biexp_sfa_op/alpha': 0.8,
        # Start at tau_r = 11 (slightly above tau_d=10) and continue down in
        # tau_r before any eta work — see Step 2 below.
        'p/qif_biexp_sfa_op/tau_r': 11.0,
        'p/qif_biexp_sfa_op/tau_d': 10.0,
        'p/qif_biexp_sfa_op/eta': -8.0,
    },
    edge_vars=[('p/qif_biexp_sfa_op/r', 'p/qif_biexp_sfa_op/r_in',
                {'weight': 15.0 * np.sqrt(2.0)})],
    init_cont=True, NPR=100, NMX=30000,
)

# Auto-07p doesn't ship a built-in style for the cusp (CP) marker; add one
# so it shows up in the codim-2 diagram below alongside BT and GH.
ode.update_bifurcation_style('CP', marker='d', color='#7F4FBF')

# %%
# Step 2: Pre-scan in :math:`\\tau_r` to anchor the two regimes
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# We'll cover two parameter regimes of this model in one example:
#
# 1. :math:`\\tau_r = \\tau_d = 10` (alpha-kernel-equivalent): clean codim-2
#    structure in :math:`(\\bar\\eta,\\, \\Delta)` — BT / GH / CP.
# 2. :math:`\\tau_r = 0.3 \\ll \\tau_d`: period-doubling cascade on the LC
#    born at the upper Hopf.
#
# Rather than re-instantiate the model twice, we continue the steady state
# in :math:`\\tau_r` from the initial :math:`\\tau_r = 11` down to small
# values, planting user points at :math:`\\tau_r = \\tau_d = 10` (UZ1) and
# at :math:`\\tau_r = 0.3` (UZ2). Each UZ then serves as the starting point
# for an :math:`\\bar\\eta` continuation in its own regime.

tau_sols, tau_cont = ode.run(
    starting_point='EP2', name='tau_r_branch',
    ICP='p/qif_biexp_sfa_op/tau_r',
    IPS=1, ILP=1, ISP=2, ISW=1, NTST=400, NCOL=4,
    NMX=5000, NPR=100, DS=-1e-3, DSMIN=1e-9, DSMAX=5e-2,
    UZR={'p/qif_biexp_sfa_op/tau_r': [10.0, 0.3]},
    UZSTOP={'p/qif_biexp_sfa_op/tau_r': [0.1, 11.0]},
)
print(f"tau_r pre-scan bifurcations: {dict(tau_sols['bifurcation'].value_counts())}")

# %%
# Step 3: 1D :math:`\\bar\\eta` scan at :math:`\\tau_r = \\tau_d = 10` (UZ1)
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# A bidirectional equilibrium continuation in :math:`\bar\eta` sweeps
# through the fold and Hopf bifurcations we'll feed into
# :func:`codim2_search` next.

eta_sols, eta_cont = ode.run(
    origin='tau_r_branch', starting_point='UZ1', name='eta_branch',
    ICP='p/qif_biexp_sfa_op/eta', bidirectional=True,
    RL0=-8.0, RL1=2.0,
    IPS=1, ILP=1, ISP=2, ISW=1, NTST=400, NCOL=4,
    NMX=2000, NPR=10,
    DS=1e-4, DSMIN=1e-8, DSMAX=5e-2,
    ITMX=40, ITNW=40, NWTN=12,
)
print("bifurcations on the eta branch:")
print(eta_sols['bifurcation'].value_counts())

ode.plot_continuation('p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/r', cont='eta_branch')
plt.title(r'1D continuation in $\bar\eta$ (steady state)')
plt.tight_layout()
plt.show()

# %%
# The 1D scan should report two Hopf bifurcations (``HB1``, ``HB2``) and
# two fold bifurcations (``LP1``, ``LP2``) flanking the bistable / Hopf
# regime. Those are the codim-1 starting points for the codim-2 search.

# %%
# Step 4: Codim-2 fold and Hopf curves in :math:`(\bar\eta,\, \Delta)`
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# :func:`codim2_search` takes a list of codim-1 labels on ``origin`` and
# runs a 2-parameter continuation of each. We run three separate searches
# so each codim-1 starting point lands in its own continuation key and we
# can colour the resulting curves distinctly on the bifurcation diagram
# below.
#
# ``max_recursion_depth=0`` disables the recursive sub-continuations the
# helper would otherwise run at each detected codim-2 point — keeps the
# output dict flat (one continuation per starting label) and the figure
# tractable. Pass ``max_recursion_depth=1`` (or more) if you want the
# helper to also explore the codim-1 branches emerging from each codim-2
# point.
#
# ``codim2_search`` defaults to ``get_stability=False`` for 2-parameter
# codim-1 continuations: per-point stability flags on a *curve of
# bifurcations* aren't meaningful (the equilibrium is degenerate by
# construction along the whole curve) and would otherwise toggle on
# numerical noise.

# We run two unidirectional continuations per codim-1 starting point —
# `DS=+1e-3` (forward) and `DS=-1e-3` (reverse) — rather than a single
# bidirectional one. Bidirectional continuation can get stuck bouncing at
# the cusp adjacent to LP1 / LP2; explicit per-direction calls give us
# clean traces in both directions for every starting point and let us
# stop each direction at a sensible Delta boundary via ``UZSTOP``.
#
# ``NPR=10`` records points every 10 continuation steps — dense enough
# that the resulting curves render as real lines (rather than a handful
# of sparsely sampled points) in the bifurcation diagram below.

shared = dict(
    pyauto_instance=ode,
    params=['p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/Delta'],
    origin=eta_cont,
    max_recursion_depth=0,
    NMX=1500, NPR=10,
    DSMIN=1e-9, DSMAX=5e-2,
    RL0=-15.0, RL1=5.0,
    bidirectional=False,
    UZSTOP={'p/qif_biexp_sfa_op/Delta': [0.0, 4.0]},
)

# (starting_point, bifurcation type) — bifurcation type is used below to
# pick the curve colour and the codim-1 markers to ignore.
codim1_points = [
    ('LP1', 'fold'),
    ('LP2', 'fold'),
    ('HB1', 'Hopf'),
    ('HB2', 'Hopf'),
]

# Run forward + reverse from each codim-1 point; collect (key, type)
# pairs for the plot loop further down.
codim2_curves = []  # list of (continuation key, 'fold' | 'Hopf')
for sp, bif_type in codim1_points:
    for ds in (1e-3, -1e-3):
        try:
            result = codim2_search(
                starting_points=[sp], DS=ds,
                name=f'{sp}_ds{"pos" if ds > 0 else "neg"}',
                **shared,
            )
            codim2_curves.append((list(result.keys())[0], bif_type))
        except Exception as exc:
            print(f"{sp} DS={ds:+g}: skipped ({type(exc).__name__}: {exc})")

# Quick summary of what was detected per curve.
print("\ncodim-2 curves recorded:")
for key, bif_type in codim2_curves:
    bif_counts = ode.get_summary(key)['bifurcation'].value_counts()
    print(f"  {key} ({bif_type}): {dict(bif_counts)}")

# %%
# Across the eight codim-2 continuations you should see ``CP`` (cusp)
# points on the fold curves emerging from ``LP1`` / ``LP2``, ``BT``
# (Bogdanov-Takens) where the fold curve and the Hopf curve from ``HB1``
# meet, and ``GH`` (generalised Hopf / Bautin) on the Hopf curve from
# ``HB2``. These are the three codim-2 phenomena most relevant for this
# model in the :math:`(\bar\eta,\, \Delta)` plane.

# %%
# Step 5: Plot the codim-2 bifurcation diagram
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# All fold curves share one colour (blue), all Hopf curves another
# (orange) — so the legend collapses to two entries even though we drew
# eight continuations. Individual codim-1 starting points are identified
# by labelled stars rather than separate legend entries.
# ``ignore=['LP', 'HB', 'UZ']`` strips both the codim-1-self markers
# (every point on a fold curve is an ``LP`` by construction, every point
# on a Hopf curve an ``HB``; drawing them would clutter the diagram) and
# user-defined ``UZ`` stop points which are not bifurcations.

CURVE_COLORS = {'fold': '#1F77B4', 'Hopf': '#FF7F0E'}

fig, ax = plt.subplots(figsize=(7, 5))
labels_used: set = set()
for key, bif_type in codim2_curves:
    color = CURVE_COLORS[bif_type]
    label = f'{bif_type} curve' if bif_type not in labels_used else None
    labels_used.add(bif_type)
    ode.plot_continuation(
        'p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/Delta', cont=key, ax=ax,
        line_color_stable=color, line_color_unstable=color,
        line_style_stable='solid', line_style_unstable='solid',
        bifurcation_legend=False, get_stability=False,
        ignore=['LP', 'HB', 'UZ'], label=label,
    )

# Mark the original 1D HB / LP locations with labelled stars so the
# reader can trace each codim-2 curve back to its starting codim-1
# bifurcation visually (without polluting the legend).
for sp, bif_type in codim1_points:
    sol, _, _ = ode.get_solution(point=sp, cont='eta_branch')
    eta_p = float(sol['eta'])
    delta_p = float(sol['Delta'])
    ax.scatter(eta_p, delta_p, marker='*', s=120,
                c=CURVE_COLORS[bif_type],
                edgecolor='k', linewidth=0.5, zorder=10)
    ax.annotate(sp, (eta_p, delta_p), xytext=(5, 5),
                 textcoords='offset points', fontsize=9)

# Restrict the y-axis to the physically meaningful region Delta >= 0
# (Delta is a half-width-at-half-maximum and can't be negative). The
# codim-2 continuations may briefly cross into Delta < 0 to detect the
# BT point cleanly, but those segments are unphysical and hidden here.
ax.set_xlim(-6.0, 2.0)
ax.set_ylim(0.0, 2.5)
ax.set_xlabel(r'$\bar\eta$')
ax.set_ylabel(r'$\Delta$')
ax.set_title(r'codim-2 bifurcation diagram in $(\bar\eta,\, \Delta)$')
ax.legend(loc='best')
plt.tight_layout()
plt.show()

# %%
# Reading the diagram:
#
# * Blue curves are fold (saddle-node) manifolds traced from ``LP1`` and
#   ``LP2``, both directions each. The cusp (``CP``) where these
#   collide marks the closing of the bistable wedge.
# * Orange curves are Hopf manifolds traced from ``HB1`` and ``HB2``,
#   both directions each. The ``BT`` (Bogdanov-Takens) point is where
#   the Hopf curve from ``HB1`` meets the fold manifold — two fold
#   branches pass through BT and one Hopf branch terminates at it
#   (a homoclinic curve also emerges, but tracking it requires
#   auto-07p's HomCont package, ``IPS=9``, and is left as a manual
#   follow-up).
# * The ``GH`` (generalised-Hopf / Bautin) points on the Hopf curve from
#   ``HB2`` are where the first Lyapunov coefficient changes sign —
#   supercritical Hopfs on one side, subcritical on the other.

# %%
# Step 6: Period-doubling cascade at :math:`\\tau_r \\ll \\tau_d` (UZ2)
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Now we shift to the second regime planted in Step 2: the user point
# ``UZ2`` at :math:`\\tau_r = 0.3 \\ll \\tau_d = 10`. Here the
# adaptation kernel is far from the alpha-kernel limit and the model's
# bifurcation structure changes — in particular, the limit cycle born
# at the upper Hopf undergoes a period-doubling.
#
# Step 6a: a fresh :math:`\\bar\\eta` scan from ``UZ2`` finds the codim-1
# Hopfs / folds at this new :math:`\\tau_r`.

eta_lo_sols, eta_lo_cont = ode.run(
    origin='tau_r_branch', starting_point='UZ2', name='eta_branch_tau03',
    ICP='p/qif_biexp_sfa_op/eta', bidirectional=True,
    RL0=-8.0, RL1=2.0,
    IPS=1, ILP=1, ISP=2, ISW=1, NTST=400, NCOL=4,
    NMX=2000, NPR=10,
    DS=1e-4, DSMIN=1e-8, DSMAX=5e-2,
    ITMX=40, ITNW=40, NWTN=12,
)
print(f"\neta scan at tau_r=0.3: {dict(eta_lo_sols['bifurcation'].value_counts())}")

# Step 6b: branch-switch to the LC at HB2 (the upper Hopf) and continue
# in :math:`\\bar\\eta`. PD points appear along this LC family when
# :math:`\\tau_r \\ll \\tau_d`.
lc_sols, lc_cont = ode.run(
    origin='eta_branch_tau03', starting_point='HB2', name='lc_pd',
    IPS=2, ISP=2, ISW=-1,
    ICP=['p/qif_biexp_sfa_op/eta', 11],
    NMX=2000, NPR=20, DS=1e-3, DSMIN=1e-9, DSMAX=5e-2,
    bidirectional=True, get_period=True,
)
print(f"LC bifurcations: {dict(lc_sols['bifurcation'].value_counts())}")

# %%
# The ``PD`` count in the LC summary tells you the cascade is in play.
# :func:`continue_period_doubling_bf` traces each PD point as a codim-1
# manifold in :math:`(\\bar\\eta,\\, \\tau_r)`, recursing through the
# cascade up to ``max_iter`` levels deep.
#
# (``max_iter`` was previously a per-call iteration counter that only
# fired once on entry; since PyCoBi 1.0.0 it bounds the *recursion depth*
# of the cascade — the variable that actually prevents runaway recursion.
# Names of successive continuations follow the ``pd_d{depth}_n{i}``
# convention so parallel sub-cascades don't collide on shared labels.)
#
# Sub-runs that fail inside :func:`ODESystem.run` surface as
# :class:`UserWarning` and the cascade continues with the remaining PD
# points rather than aborting — useful when chasing a cascade that
# eventually enters a chaotic regime where individual continuations may
# fail to converge.

pd_continuation_data = ode.results[ode.get_continuation('lc_pd').key]
pd_names, _ = continue_period_doubling_bf(
    solution=pd_continuation_data,
    continuation=lc_cont,
    pyauto_instance=ode,
    max_iter=2, precision=3,
    ICP=['p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/tau_r'],
    IPS=2, ISW=2, ISP=2,
    NMX=400, NPR=20, DS=1e-3, DSMIN=1e-9, DSMAX=5e-2,
    RL0=-8.0, RL1=2.0,
    UZSTOP={'p/qif_biexp_sfa_op/tau_r': [0.05, 12.0]},
)
print(f"PD continuations: {pd_names}")

# %%
# Plot the PD continuations in :math:`(\\bar\\eta,\\, \\tau_r)`, overlaid
# with the user points ``UZ1 = (any eta, tau_r=10)`` and
# ``UZ2 = (any eta, tau_r=0.3)`` as horizontal reference lines so the
# reader can see where the codim-2 work (Steps 3-5) lives versus where
# the period-doubling cascade emerges (this step).

fig, ax = plt.subplots(figsize=(7, 5))
labels_seen: set = set()
for pd_name in pd_names:
    label = 'PD curve' if 'PD' not in labels_seen else None
    labels_seen.add('PD')
    ode.plot_continuation(
        'p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/tau_r', cont=pd_name,
        ax=ax,
        line_color_stable='#D62728', line_color_unstable='#D62728',
        line_style_stable='solid', line_style_unstable='solid',
        bifurcation_legend=False, get_stability=False,
        ignore=['PD', 'UZ'], label=label,
    )
ax.axhline(10.0, color='0.7', linewidth=0.5, linestyle='--',
            label=r'$\tau_r = \tau_d$ (Step 3-5 regime)')
ax.axhline(0.3,  color='0.7', linewidth=0.5, linestyle=':',
            label=r'$\tau_r = 0.3$ (PD-cascade slice)')
ax.set_xlabel(r'$\bar\eta$')
ax.set_ylabel(r'$\tau_r$')
ax.set_title(r'period-doubling curve in $(\bar\eta,\, \tau_r)$')
ax.legend(loc='best')
plt.tight_layout()
plt.show()

# %%
# Step 7: Failure modes and what to expect on unfamiliar models
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Both helpers wrap a series of nested ``ODESystem.run`` calls. If any
# individual sub-run raises inside auto-07p (common on low-quality
# parameter regimes or aggressive step sizes), the failure surfaces as a
# :class:`UserWarning` and the search continues with the remaining starting
# points rather than aborting. Read the warnings carefully: they cite the
# auto-07p exception type and message, the sub-run label that failed, and
# the kwargs hook (``kwargs_1D_lc_cont``, ``kwargs_2D_cont``, ...) you can
# use to override the default constants for that path.

# %%
# Step 8: Clean up
# ^^^^^^^^^^^^^^^^

ode.close_session(clear_files=True)
