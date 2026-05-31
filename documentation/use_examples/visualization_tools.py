"""
Visualization Tools: A Tour of PyCoBi's Plot Methods
====================================================

PyCoBi exposes five visualization methods on top of matplotlib, each
suited to a different aspect of a bifurcation analysis:

* :meth:`ODESystem.plot_continuation` — the 1D / 2D bifurcation
  diagram.  Stable / unstable segmentation, codim-1 marker overlays,
  configurable colours, optional bifurcation-marker legend.
* :meth:`ODESystem.plot_continuation_grid` — lay out multiple
  continuations as a grid of subplots in one call. Per-panel titles,
  panel labels (the (a)/(b)/(c) annotations of published figures),
  shared kwargs.
* :meth:`ODESystem.plot_bifurcation_points` — the low-level marker
  helper that ``plot_continuation`` builds on. Useful when you want to
  overlay custom markers onto an arbitrary axis (e.g. to annotate a
  result coming from a *non*-continuation analysis).
* :meth:`ODESystem.plot_timeseries` — :math:`r(t)` over time for one or
  more labelled limit-cycle points on a continuation.
* :meth:`ODESystem.plot_trajectory` — 2D or 3D phase-space trajectory.
  3D plots get an optional colorbar that surfaces which projected
  coordinate the colormap encodes.

Plus the per-marker style hook:

* :meth:`ODESystem.update_bifurcation_style` — change the marker /
  colour PyCoBi uses for a given bifurcation type (BT, GH, CP, …).

This example uses the bi-exponential QIF-SFA mean-field model
(``qif_biexp_sfa.yaml`` next to this script) at the alpha-kernel-
equivalent parameter point :math:`\\tau_r = \\tau_d = 10` — same model
and parameter regime as :ref:`Automated Codim-2 Search and
Period-Doubling Cascades`, so we can reuse its continuations and focus
the prose here on the visualization side.

References
^^^^^^^^^^

.. [1] R. Gast, H. Schmidt, T.R. Knösche (2020) *A Mean-Field Description
       of Bursting Dynamics in Spiking Neural Networks with Short-Term
       Adaptation.* Neural Computation 32 (9): 1615-1634.
"""

# %%
# Step 1: Load the model and run the continuations we'll visualize
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Same model / parameters as the codim-2 example, plus a limit-cycle
# continuation in :math:`\\bar\\eta` from the upper Hopf so we have
# something interesting to overlay in the 1D diagram and to render in
# phase space. ``get_timeseries=True`` on the LC run is what
# :meth:`plot_timeseries` needs further down: it records the per-period
# :math:`r(t)` sampling on every labelled LC point.

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from pycobi import ODESystem

here = Path(__file__).resolve().parent
yaml_path = str(here / 'qif_biexp_sfa' / 'qif_biexp_sfa')

ode = ODESystem.from_yaml(
    yaml_path,
    auto_dir="~/PycharmProjects/auto-07p",
    node_vars={
        'p/qif_biexp_sfa_op/Delta': 2.0,
        'p/qif_biexp_sfa_op/alpha': 0.8,
        'p/qif_biexp_sfa_op/tau_r': 10.0,
        'p/qif_biexp_sfa_op/tau_d': 10.0,
        'p/qif_biexp_sfa_op/eta': -8.0,
    },
    edge_vars=[('p/qif_biexp_sfa_op/r', 'p/qif_biexp_sfa_op/r_in',
                {'weight': 15.0 * np.sqrt(2.0)})],
    init_cont=True, NPR=100, NMX=30000,
)

# Equilibrium continuation in eta — exposes 2 HBs and 2 LPs.
eta_sols, eta_cont = ode.run(
    starting_point='EP2', name='eta_branch',
    ICP='p/qif_biexp_sfa_op/eta', bidirectional=True,
    RL0=-8.0, RL1=2.0,
    IPS=1, ILP=1, ISP=2, ISW=1, NTST=400, NCOL=4,
    NMX=2000, NPR=10, DS=1e-4, DSMIN=1e-8, DSMAX=5e-2,
    ITMX=40, ITNW=40, NWTN=12,
)

# Limit-cycle continuation from HB2. Two flags worth singling out:
#
# * ``get_timeseries=True`` records the per-period time sampling on every
#   labelled LC point — required by :meth:`plot_timeseries`.
# * ``reduce_limit_cycle=False`` stores the *full* per-period state-vector
#   samples (rather than the default ``(min, max)`` envelope). This costs
#   more memory but is what :meth:`plot_timeseries` and the
#   point-specific :meth:`plot_trajectory` need to render the LC orbit —
#   the envelope-only form only contains two points per LC.
#   :meth:`plot_continuation` still works on the full-sample form by
#   internally extracting ``(min, max)`` per row.
#
# (We deliberately *don't* set ``bidirectional=True`` here: bidirectional
# continuation re-builds the summary at the merge step with the default
# ``reduce_limit_cycle=True``, which would clobber the full per-period
# samples this example needs. The forward LC family is enough to drive
# all of the figures below.)
lc_sols, lc_cont = ode.run(
    origin='eta_branch', starting_point='HB2', name='lc_branch',
    IPS=2, ISP=2, ISW=-1,
    ICP=['p/qif_biexp_sfa_op/eta', 11],
    NMX=400, NPR=10, DS=1e-3, DSMIN=1e-9, DSMAX=5e-2,
    get_period=True, get_timeseries=True,
    reduce_limit_cycle=False,
    STOP=["BP3", "LP5"],
)
print("eta-branch:", dict(eta_sols['bifurcation'].value_counts()))
print("lc-branch :", dict(lc_sols['bifurcation'].value_counts()))

# Pin a consistent palette across the example so the same line / marker
# colours mean the same thing from figure to figure.
EQ_COLOR = '#1F77B4'   # blue  — equilibrium branches
LC_COLOR = '#FF7F0E'   # orange — limit-cycle branches
HL_COLOR = '#D62728'   # red   — highlighted / annotated trajectories


# %%
# Step 2: :meth:`plot_continuation` — a clean 1D bifurcation diagram
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# The bread-and-butter method. Stable segments render solid, unstable
# segments dotted-gray (1.0.0 default — see Step 8 for how to revert).
# Codim-1 bifurcation markers come for free via
# :meth:`plot_bifurcation_points`; ``bifurcation_legend=True`` (the
# default) adds them to the matplotlib legend if any are present.

fig, ax = plt.subplots(figsize=(7, 4.5))
ode.plot_continuation(
    'p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/r',
    cont='eta_branch', ax=ax,
    line_color_stable=EQ_COLOR, line_color_unstable=EQ_COLOR,
    label='equilibrium branch',
)
ode.plot_continuation(
    'p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/r',
    cont='lc_branch', ax=ax,
    line_color_stable=LC_COLOR, line_color_unstable=LC_COLOR,
    ignore=['UZ', 'BP', 'EP'], label='limit cycle branch',
)
ax.set_xlim(-8.0, 2.0)
ax.set_xlabel(r'$\bar\eta$')
ax.set_ylabel(r'$r$')
ax.set_title('plot_continuation — equilibrium + limit cycle')
ax.legend(loc='best')
plt.tight_layout()
plt.show()

# %%
# Two patterns to notice in the call:
#
# * The same axis hosts both calls — passing ``ax=ax`` is how you
#   overlay multiple continuations on one figure.
# * ``line_color_stable`` and ``line_color_unstable`` are independent
#   knobs. Setting both to the same colour forces a single hue per
#   branch; the linestyle (solid / dotted) still distinguishes stable
#   vs unstable.
# * ``ignore=['UZ', 'BP', 'EP']`` suppresses the markers that aren't
#   genuine bifurcations — user-defined stop points, start / end of
#   branch, etc. — so only the HBs / LPs / PDs annotate the curve.

# %%
# Step 3: :meth:`plot_continuation_grid` — multi-panel comparison
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# ``plot_continuation_grid`` takes a list of subplot specs (one dict
# per panel). Each spec supplies the ``x`` / ``y`` / ``cont`` triple
# plus any keyword arguments the user wants to forward to
# :meth:`plot_continuation` for that specific panel (``ignore``,
# ``line_color_*``, ``get_stability``, …).  Per-panel ``title`` and
# ``panel_label`` keys are extracted before forwarding — the latter
# draws the (a)/(b)/(c) annotation in the upper-left corner of each
# subplot, matching the convention of most published bifurcation
# figures.

plots = [
    {'x': 'p/qif_biexp_sfa_op/eta', 'y': 'r',  'cont': 'eta_branch',
     'title': r'equilibria — $r$ vs. $\bar\eta$', 'panel_label': '(a)',
     'line_color_stable': EQ_COLOR, 'line_color_unstable': EQ_COLOR},
    {'x': 'p/qif_biexp_sfa_op/eta', 'y': 'A',  'cont': 'eta_branch',
     'title': r'equilibria — adaptation $A$', 'panel_label': '(b)',
     'line_color_stable': EQ_COLOR, 'line_color_unstable': EQ_COLOR},
    {'x': 'p/qif_biexp_sfa_op/eta', 'y': 'r',  'cont': 'lc_branch',
     'title': r'limit cycle — $r$ vs. $\bar\eta$', 'panel_label': '(c)',
     'line_color_stable': LC_COLOR, 'line_color_unstable': LC_COLOR,
     'ignore': ['UZ', 'BP', 'EP']},
    {'x': 't', 'y': 'r', 'cont': 0,
     'title': r'IVP — $r$ vs. $t$', 'panel_label': '(d)',
     'get_stability': False, 'line_color_stable': HL_COLOR},
]
fig, axes, line_cols = ode.plot_continuation_grid(
    plots, ncols=2, figsize=(10, 7),
    bifurcation_legend=False,
)
fig.suptitle('plot_continuation_grid — four views of the same model',
              y=1.02, fontsize='large')
plt.show()

# %%
# Notes on the figure:
#
# * Panels (a) and (b) plot the same equilibrium scan against different
#   y-axis variables; codim-1 markers (HB / LP) annotate each curve
#   independently.
# * Panel (c) renders the limit cycle as a (min, max) envelope — for an
#   LC ``plot_continuation`` packs each period into two traces
#   bracketing the oscillation amplitude.
# * Panel (d) draws the time-domain IVP that ``init_cont=True`` ran at
#   model construction. With ``get_stability=False`` we tell it not to
#   try to colour the trace by stability (the IVP doesn't record any).
# * The empty position in the 2x2 grid would be deleted automatically
#   if we'd passed three plots instead of four.

# %%
# Step 4: :meth:`plot_bifurcation_points` — overlay markers manually
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# ``plot_continuation`` calls ``plot_bifurcation_points`` internally to
# draw the codim-1 markers. You can also call it yourself when you want
# markers on an axis whose curves came from somewhere else — e.g. to
# annotate the original 1D HB / LP locations on a *codim-2* diagram so
# the reader can trace each curve back to its starting point. Below we
# render the equilibrium curve manually (the codim-1 markers are
# suppressed via ``bifurcation_legend=False, ignore=['LP', 'HB']``) and
# then re-draw just the HB / LP markers with custom colours on top.

eta_extracted, vmap = ode.extract(
    ['p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/r', 'bifurcation'],
    cont='eta_branch',
)
eta_vals  = eta_extracted[vmap['p/qif_biexp_sfa_op/eta']]
r_vals    = eta_extracted[vmap['p/qif_biexp_sfa_op/r']]
bif_types = eta_extracted['bifurcation']

fig, ax = plt.subplots(figsize=(7, 4.5))
ode.plot_continuation(
    'p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/r',
    cont='eta_branch', ax=ax,
    line_color_stable=EQ_COLOR, line_color_unstable=EQ_COLOR,
    bifurcation_legend=False, ignore=['LP', 'HB', 'EP'],
)
# Re-draw HB / LP markers in custom colours that match the rest of the
# example. plot_bifurcation_points takes the type / x / y arrays
# directly and writes to whichever axis we pass.
ode.plot_bifurcation_points(
    solution_types=bif_types, x_vals=eta_vals, y_vals=r_vals, ax=ax,
    custom_bf_styles={
        'HB': {'marker': 'o', 'color': '#2CA02C'},
        'LP': {'marker': 'D', 'color': '#9467BD'},
    },
)
ax.set_xlim(-8.0, 2.0)
ax.set_xlabel(r'$\bar\eta$')
ax.set_ylabel(r'$r$')
ax.set_title('plot_bifurcation_points — custom markers on top of an existing curve')
ax.legend(loc='best')
plt.tight_layout()
plt.show()

# %%
# ``custom_bf_styles={ <type>: { 'marker': ..., 'color': ... }, … }`` is
# the per-call override; for a persistent change (the styles stick on
# the ``ODESystem`` instance) use :meth:`update_bifurcation_style` —
# see Step 8.

# %%
# Step 5: :meth:`plot_timeseries` — :math:`r(t)` along the LC family
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# ``plot_timeseries`` shows how one state variable evolves *within* a
# limit-cycle period, at one or more labelled LC points. Pass
# ``points=[<point_id>, …]`` to pick specific LC labels — useful for
# comparing how the period and waveform change as you move along the
# branch. Below we grab three labels from the LC family (start, middle,
# end) and overlay their :math:`r(t)` traces on one axis.

lc_labels = list(ode.get_summary('lc_branch').index)
pick = [lc_labels[1], lc_labels[len(lc_labels) // 2], lc_labels[-2]]
print(f"LC labels picked: {pick}")

fig, ax = plt.subplots(figsize=(7, 4.5))
ode.plot_timeseries(
    var='r', cont='lc_branch', points=pick, ax=ax,
    linespecs=[
        {'colors': EQ_COLOR},
        {'colors': LC_COLOR},
        {'colors': HL_COLOR},
    ],
)
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$r$')
ax.set_title(r'plot_timeseries — $r(t)$ at three LC labels')
plt.tight_layout()
plt.show()

# %%
# ``linespecs`` is a per-trace dict of LineCollection kwargs that gets
# merged into the shared ``**kwargs`` for each point — handy when you
# want to set per-trace colour / style. Pass ``points=None`` (the
# default) to plot every labelled LC point on the branch.

# %%
# Step 6: :meth:`plot_trajectory` — 2D and 3D phase space
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# ``plot_trajectory`` draws a parametric curve through the model's
# state space. Pass two state variables for a 2D projection, three for
# a 3D plot. The 3D version accepts ``colorbar=True`` to attach a
# colorbar to the figure surfacing which projected coordinate the
# Line3DCollection's colormap encodes (default: the x-axis variable).

mid_lc = lc_labels[len(lc_labels) // 2]

# 2D phase plot — (r, v) projection of one LC period.
fig, ax = plt.subplots(figsize=(6, 4.5))
ode.plot_trajectory(
    variables=['r', 'v'], cont='lc_branch', point=mid_lc, ax=ax,
    colors=HL_COLOR,
)
ax.set_xlabel(r'$r$')
ax.set_ylabel(r'$v$')
ax.set_title(r'plot_trajectory — 2D LC orbit in $(r,\, v)$')
plt.tight_layout()
plt.show()

# 3D phase plot — same LC in (r, v, A) with a colorbar over the period.
fig = plt.figure(figsize=(7, 5.5))
ax = fig.add_subplot(111, projection='3d')
ode.plot_trajectory(
    variables=['r', 'v', 'A'], cont='lc_branch', point=mid_lc, ax=ax,
    colorbar=True, colorbar_label=r'$r$ (firing rate)',
)
ax.set_title(r'plot_trajectory — 3D LC orbit in $(r,\, v,\, A)$ with colorbar')
plt.tight_layout()
plt.show()

# %%
# The colormap mapping is controlled by ``array='x'`` / ``'y'`` / ``'z'``
# in the 3D variant — the default ``'x'`` colours the curve by the
# first variable (here ``r``), so a glance at the colorbar tells you
# where on the orbit each colour band sits. ``colorbar_label`` overrides
# the default colorbar label (otherwise the variable name is used).

# %%
# Step 7: :meth:`update_bifurcation_style` — persistent marker overrides
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Auto-07p detects more codim-1 / codim-2 types than PyCoBi ships
# default styles for (e.g. ``CP`` cusps don't have a built-in marker).
# :meth:`update_bifurcation_style` registers a marker / colour for a
# type on the active ``ODESystem``; every subsequent
# :meth:`plot_continuation` / :meth:`plot_bifurcation_points` call uses
# the registered style:
#
# .. code-block:: python
#
#     ode.update_bifurcation_style('CP', marker='d', color='#7F4FBF')
#     ode.update_bifurcation_style('PD', marker='h', color='#D62728')
#
# Unlike the per-call ``custom_bf_styles=`` override in Step 4, this
# sticks for the lifetime of the ``ode`` instance. It's the right
# extension point when you're generating a series of figures (e.g. the
# codim-2 diagrams in :ref:`Automated Codim-2 Search and
# Period-Doubling Cascades`) and want them all to share a custom legend
# style for CP / GH / BT.

# %%
# Step 8: Default-style note — gray for unstable
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Pre-1.0 PyCoBi rendered both stable and unstable line segments as
# solid black, distinguishable only by the linestyle. 1.0.0 changes the
# default unstable colour to ``'gray'`` so the two regimes are
# visually distinguishable at a glance. The default applies wherever
# the line-collection builder is called internally —
# :meth:`plot_continuation`, :meth:`plot_trajectory`,
# :meth:`plot_timeseries`, and :meth:`plot_continuation_grid`.
#
# To revert to the legacy all-black appearance, pass
# ``line_color_unstable='k'`` through ``**kwargs`` on any of the plot
# methods:
#
# .. code-block:: python
#
#     ode.plot_continuation(
#         'p/qif_biexp_sfa_op/eta', 'r', cont='eta_branch',
#         line_color_unstable='k',  # legacy 0.x default
#     )

# %%
# Step 9: Clean up
# ^^^^^^^^^^^^^^^^

ode.close_session(clear_files=True)
