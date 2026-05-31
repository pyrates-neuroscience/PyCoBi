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
#    structure in :math:`(\\bar\\eta,\\, \\Delta)` — BT / GH / CP, no PD.
# 2. :math:`\\tau_r = 0.1 \\ll \\tau_d`: same codim-2 backbone plus a
#    period-doubling cascade on the LC born at the upper Hopf.
#
# Rather than re-instantiate the model twice, we continue the steady state
# in :math:`\\tau_r` from the initial :math:`\\tau_r = 11` down to small
# values, planting user points at :math:`\\tau_r = \\tau_d = 10` (UZ1) and
# at :math:`\\tau_r = 0.1` (UZ2). Each UZ then serves as the starting point
# for an :math:`\\bar\\eta` continuation in its own regime.

tau_sols, tau_cont = ode.run(
    starting_point='EP2', name='tau_r_branch',
    ICP='p/qif_biexp_sfa_op/tau_r',
    IPS=1, ILP=1, ISP=2, ISW=1, NTST=400, NCOL=4,
    NMX=5000, NPR=100, DS=-1e-3, DSMIN=1e-9, DSMAX=5e-2,
    UZR={'p/qif_biexp_sfa_op/tau_r': [10.0, 0.1]},
    UZSTOP={'p/qif_biexp_sfa_op/tau_r': [0.05, 11.0]},
)
print(f"tau_r pre-scan bifurcations: {dict(tau_sols['bifurcation'].value_counts())}")

# %%
# Shared bookkeeping
# ^^^^^^^^^^^^^^^^^^
#
# We run the same analysis pipeline at two values of :math:`\\tau_r`
# (Steps 3-4 and 5-6 below), so it pays to factor the codim-1 / codim-2
# bookkeeping into a small helper. ``codim1_points`` lists the four
# starting points we expect on every :math:`\\bar\\eta` scan (the model
# has two folds and two Hopfs in the parameter regime we're exploring);
# ``CURVE_COLORS`` keys colours by bifurcation type so the codim-2
# diagrams collapse their legends to two entries even though we draw
# eight curves.

codim1_points = [
    ('LP1', 'fold'),
    ('LP2', 'fold'),
    ('HB1', 'Hopf'),
    ('HB2', 'Hopf'),
]
CURVE_COLORS = {'fold': '#1F77B4', 'Hopf': '#FF7F0E'}
PD_COLOR = '#D62728'


def run_codim2_in_delta(origin_name, name_prefix):
    """Run the 8-call codim-2 search (LP1/LP2/HB1/HB2 × DS=±1e-3) in
    :math:`(\\bar\\eta,\\, \\Delta)` and return ``[(key, bif_type)]``.
    Forward + reverse continuations are kept separate (rather than one
    bidirectional call) because bidirectional gets stuck bouncing at
    the cusp adjacent to LP1 / LP2.
    """
    shared_kwargs = dict(
        pyauto_instance=ode,
        params=['p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/Delta'],
        origin=origin_name,
        max_recursion_depth=0,
        NMX=1500, NPR=10,
        DSMIN=1e-9, DSMAX=5e-2,
        RL0=-15.0, RL1=5.0,
        bidirectional=False,
        UZSTOP={'p/qif_biexp_sfa_op/Delta': [0.0, 4.0]},
    )
    curves = []
    for sp, bif_type in codim1_points:
        for ds in (1e-3, -1e-3):
            try:
                result = codim2_search(
                    starting_points=[sp], DS=ds,
                    name=f'{name_prefix}_{sp}_{"pos" if ds > 0 else "neg"}',
                    **shared_kwargs,
                )
                curves.append((list(result.keys())[0], bif_type))
            except Exception as exc:
                print(f"{name_prefix} {sp} DS={ds:+g}: skipped "
                      f"({type(exc).__name__}: {exc})")
    return curves


def plot_2d_diagram(ax, codim2_curves, eta_origin, *, pd_names=(),
                     title, xlim=(-6.0, 2.0), ylim=(0.0, 2.5)):
    """Plot a 2D codim-2 diagram in :math:`(\\bar\\eta,\\, \\Delta)` from
    the curves produced by :func:`run_codim2_in_delta`. Optionally overlay
    one or more PD continuations (rendered in red). Fold curves all share
    the blue colour, Hopf curves all share orange — so the legend has at
    most three entries (fold curve, Hopf curve, PD curve)."""
    labels_used = set()
    for key, bif_type in codim2_curves:
        color = CURVE_COLORS[bif_type]
        label = f'{bif_type} curve' if bif_type not in labels_used else None
        labels_used.add(bif_type)
        ode.plot_continuation(
            'p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/Delta',
            cont=key, ax=ax,
            line_color_stable=color, line_color_unstable=color,
            line_style_stable='solid', line_style_unstable='solid',
            bifurcation_legend=False, get_stability=False,
            ignore=['LP', 'HB', 'UZ'], label=label,
        )
    for pd_name in pd_names:
        label = 'PD curve' if 'PD' not in labels_used else None
        labels_used.add('PD')
        ode.plot_continuation(
            'p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/Delta',
            cont=pd_name, ax=ax,
            line_color_stable=PD_COLOR, line_color_unstable=PD_COLOR,
            line_style_stable='solid', line_style_unstable='solid',
            bifurcation_legend=False, get_stability=False,
            ignore=['PD', 'UZ'], label=label,
        )
    for sp, bif_type in codim1_points:
        try:
            sol, _, _ = ode.get_solution(point=sp, cont=eta_origin)
            eta_p = float(sol['eta'])
            delta_p = float(sol['Delta'])
            ax.scatter(eta_p, delta_p, marker='*', s=120,
                        c=CURVE_COLORS[bif_type],
                        edgecolor='k', linewidth=0.5, zorder=10)
            ax.annotate(sp, (eta_p, delta_p), xytext=(5, 5),
                         textcoords='offset points', fontsize=9)
        except Exception:
            pass
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r'$\bar\eta$')
    ax.set_ylabel(r'$\Delta$')
    ax.set_title(title)
    ax.legend(loc='best')


def plot_1d_diagram(ax, eta_cont_name, lc_cont_name, *, title,
                     xlim=(-8.0, 2.0)):
    """Plot a 1D bifurcation diagram in :math:`(\\bar\\eta,\\, r)` with the
    equilibrium branch (blue) and the limit cycle branch (orange) overlaid.
    Codim-1 markers (HB / LP / PD) appear automatically; UZ / BP / EP are
    filtered out as cosmetic noise."""
    ode.plot_continuation(
        'p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/r',
        cont=eta_cont_name, ax=ax,
        line_color_stable='#1F77B4', line_color_unstable='#1F77B4',
        bifurcation_legend=False, label='equilibrium branch',
        ignore=['UZ', 'BP', 'EP'],
    )
    ode.plot_continuation(
        'p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/r',
        cont=lc_cont_name, ax=ax,
        line_color_stable='#FF7F0E', line_color_unstable='#FF7F0E',
        bifurcation_legend=False, label='limit cycle branch',
        ignore=['UZ', 'BP', 'EP'],
    )
    ax.set_xlim(*xlim)
    ax.set_xlabel(r'$\bar\eta$')
    ax.set_ylabel(r'$r$')
    ax.set_title(title)
    ax.legend(loc='best')

# %%
# Step 3: Continuation + codim-2 search at :math:`\\tau_r = \\tau_d = 10`
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# This is the alpha-kernel-equivalent regime — clean codim-2 structure
# (BT / GH / CP), no period doublings. A bidirectional :math:`\\bar\\eta`
# scan from ``UZ1`` exposes the codim-1 starting points; we then
# branch-switch to the limit cycle born at ``HB2`` and run
# :func:`codim2_search` from each codim-1 point.

eta_sols, eta_cont = ode.run(
    origin='tau_r_branch', starting_point='UZ1', name='eta_branch',
    ICP='p/qif_biexp_sfa_op/eta', bidirectional=True,
    RL0=-8.0, RL1=2.0,
    IPS=1, ILP=1, ISP=2, ISW=1, NTST=400, NCOL=4,
    NMX=2000, NPR=10, DS=1e-4, DSMIN=1e-8, DSMAX=5e-2,
    ITMX=40, ITNW=40, NWTN=12,
)
print("\neta scan (tau_r=10):", dict(eta_sols['bifurcation'].value_counts()))

lc_sols, lc_cont = ode.run(
    origin='eta_branch', starting_point='HB2', name='lc_branch',
    IPS=2, ISP=2, ISW=-1,
    ICP=['p/qif_biexp_sfa_op/eta', 11],
    NMX=400, NPR=10, DS=1e-3, DSMIN=1e-9, DSMAX=5e-2,
    bidirectional=True, get_period=True,
)
print("LC (tau_r=10):", dict(lc_sols['bifurcation'].value_counts()))

codim2_curves = run_codim2_in_delta('eta_branch', name_prefix='delta10')
print("codim-2 curves at tau_r=10:")
for key, bif_type in codim2_curves:
    print(f"  {key} ({bif_type}): "
          f"{dict(ode.get_summary(key)['bifurcation'].value_counts())}")

# %%
# Step 4: Figures at :math:`\\tau_r = 10`
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Two figures: a 1D bifurcation diagram in :math:`(\\bar\\eta,\\, r)` with
# the equilibrium and limit-cycle branches overlaid, and a 2D codim-2
# diagram in :math:`(\\bar\\eta,\\, \\Delta)` showing the BT / CP / GH
# structure.

fig, ax = plt.subplots(figsize=(7, 4.5))
plot_1d_diagram(ax, 'eta_branch', 'lc_branch',
                 title=r'1D bifurcation diagram, $\tau_r = 10$')
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=(7, 5))
plot_2d_diagram(ax, codim2_curves, 'eta_branch',
                 title=r'codim-2 bifurcation diagram, $\tau_r = 10$')
plt.tight_layout()
plt.show()

# %%
# Reading the 2D diagram:
#
# * Blue fold curves traced from ``LP1`` / ``LP2``, both directions
#   each. The cusp ``CP`` where these collide marks the closing of the
#   bistable wedge.
# * Orange Hopf curves traced from ``HB1`` / ``HB2``. The ``BT``
#   (Bogdanov-Takens) point is where the Hopf curve from ``HB1`` meets
#   the fold manifold — two fold branches pass through BT and one Hopf
#   branch terminates at it (a homoclinic curve also emerges from BT,
#   but tracking it requires auto-07p's HomCont package ``IPS=9`` and
#   is left as a manual follow-up).
# * The ``GH`` (generalised-Hopf / Bautin) point on the Hopf curve from
#   ``HB2`` is where the first Lyapunov coefficient changes sign —
#   supercritical Hopfs on one side, subcritical on the other.

# %%
# Step 5: Continuation + codim-2 + PD cascade at :math:`\\tau_r = 0.1`
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Same pipeline at the second user point ``UZ2``, with the addition of
# a :func:`continue_period_doubling_bf` call after the LC continuation —
# at :math:`\\tau_r = 0.1 \\ll \\tau_d` the LC born at ``HB2`` carries
# at least one period-doubling point.

eta_lo_sols, eta_lo_cont = ode.run(
    origin='tau_r_branch', starting_point='UZ2', name='eta_branch_lo',
    ICP='p/qif_biexp_sfa_op/eta', bidirectional=True,
    RL0=-8.0, RL1=2.0,
    IPS=1, ILP=1, ISP=2, ISW=1, NTST=400, NCOL=4,
    NMX=2000, NPR=10, DS=1e-4, DSMIN=1e-8, DSMAX=5e-2,
    ITMX=40, ITNW=40, NWTN=12,
)
print("\neta scan (tau_r=0.1):", dict(eta_lo_sols['bifurcation'].value_counts()))

lc_lo_sols, lc_lo_cont = ode.run(
    origin='eta_branch_lo', starting_point='HB2', name='lc_branch_lo',
    IPS=2, ISP=2, ISW=-1,
    ICP=['p/qif_biexp_sfa_op/eta', 11],
    NMX=2000, NPR=20, DS=1e-3, DSMIN=1e-9, DSMAX=5e-2,
    bidirectional=True, get_period=True,
)
print("LC (tau_r=0.1):", dict(lc_lo_sols['bifurcation'].value_counts()))

codim2_curves_lo = run_codim2_in_delta('eta_branch_lo', name_prefix='delta01')
print("codim-2 curves at tau_r=0.1:")
for key, bif_type in codim2_curves_lo:
    print(f"  {key} ({bif_type}): "
          f"{dict(ode.get_summary(key)['bifurcation'].value_counts())}")

# Trace the PD locus in (eta, Delta) so it can be overlaid on the same
# 2D diagram as the codim-2 curves above.
pd_names, _ = continue_period_doubling_bf(
    solution=ode.results[ode.get_continuation('lc_branch_lo').key],
    continuation=lc_lo_cont, pyauto_instance=ode,
    max_iter=2, precision=3,
    ICP=['p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/Delta'],
    IPS=2, ISW=2, ISP=2,
    NMX=400, NPR=20, DS=1e-3, DSMIN=1e-9, DSMAX=5e-2,
    RL0=-8.0, RL1=2.0,
    UZSTOP={'p/qif_biexp_sfa_op/Delta': [0.0, 5.0]},
)
print("PD continuations:", pd_names)

# %%
# Step 6: Figures at :math:`\\tau_r = 0.1`
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Same recipe as Step 4 (1D + 2D in :math:`(\\bar\\eta,\\, \\Delta)`),
# now with the PD point visible on the limit-cycle branch and the PD
# locus overlaid in red on the 2D diagram.

fig, ax = plt.subplots(figsize=(7, 4.5))
plot_1d_diagram(ax, 'eta_branch_lo', 'lc_branch_lo',
                 title=r'1D bifurcation diagram, $\tau_r = 0.1$')
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=(7, 5))
plot_2d_diagram(ax, codim2_curves_lo, 'eta_branch_lo',
                 pd_names=pd_names,
                 title=r'codim-2 + PD bifurcation diagram, $\tau_r = 0.1$')
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
