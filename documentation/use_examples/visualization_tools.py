"""
Visualization Tools: Grid Layouts and 3D Phase Trajectories
===========================================================

PyCoBi 1.0.0 introduced two visualization additions and one default-style
change:

* :meth:`ODESystem.plot_continuation_grid` — lay out multiple continuations
  as a grid of subplots in a single figure call. Each subplot is described
  by an ``{x, y, cont, ...}`` dict; per-subplot ``title`` / ``panel_label``
  keys and shared keyword arguments cover most figure-quality layouts
  without dropping back to manual ``plt.subplots`` boilerplate.
* :meth:`ODESystem.plot_trajectory(colorbar=True)` — for 3D phase
  trajectories, attach a colorbar that surfaces the color-mapped projected
  coordinate (default: the x-axis variable) so the reader can quickly read
  off where on the trajectory each colour band sits.
* Default ``line_color_unstable`` changed from ``'k'`` to ``'gray'`` so
  stable and unstable line segments are visually distinguishable without
  requiring users to pass an explicit override. User-supplied overrides
  still take precedence; the line-style distinction (``'solid'`` vs
  ``'dotted'``) is unchanged.

This example uses the bi-exponential QIF-SFA model (see
:ref:`Automated Codim-2 Search and Period-Doubling Cascades` for its
equations) so the figures actually have something interesting to show:
multiple Hopf and fold bifurcations to label, a limit-cycle branch to
plot, and a 4-state phase space we can project into 3D.
"""

# %%
# Step 1: Load the model and run a few continuations
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# We'll need several continuations in `ode` before the grid plot has
# anything interesting to lay out. Pull the same starting parameters as
# the codim-2 example.

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

# Equilibrium continuation in eta.
ode.run(
    starting_point='EP2', name='eta_eq',
    ICP='p/qif_biexp_sfa_op/eta', bidirectional=True,
    RL0=-10.0, RL1=0.0,
    IPS=1, ILP=1, ISP=2, ISW=1, NTST=200,
    NMX=500, NPR=20, DS=1e-3, DSMIN=1e-7, DSMAX=5e-2,
)

# Limit-cycle continuation from the upper Hopf.
ode.run(
    origin='eta_eq', starting_point='HB2', name='lc',
    IPS=2, ISW=-1, ISP=2,
    ICP='p/qif_biexp_sfa_op/eta',
    NMX=300, NPR=10, DS=1e-3, DSMIN=1e-7, DSMAX=5e-2,
    bidirectional=True, get_period=True,
)

# %%
# Step 2: `plot_continuation_grid` — multi-panel layout
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# `plot_continuation_grid` takes a list of subplot specs. Each spec is a
# dict containing ``x``, ``y``, ``cont`` (forwarded to
# :meth:`plot_continuation`) plus the optional ``title`` and
# ``panel_label`` keys. Per-subplot keys override the shared kwargs passed
# at the helper level — useful when most panels share an option (e.g.
# ``bifurcation_legend=False``) but one needs to be different.

plots = [
    {'x': 'p/qif_biexp_sfa_op/eta', 'y': 'r',
     'cont': 'eta_eq', 'title': r'equilibria in $\bar\eta$',
     'panel_label': '(a)'},
    {'x': 'p/qif_biexp_sfa_op/eta', 'y': 'r',
     'cont': 'lc', 'title': r'limit cycle in $\bar\eta$',
     'panel_label': '(b)', 'ignore': ['UZ', 'BP']},
    {'x': 'p/qif_biexp_sfa_op/eta', 'y': 'A',
     'cont': 'eta_eq', 'title': r'adaptation $A$ vs. $\bar\eta$',
     'panel_label': '(c)'},
    {'x': 't', 'y': 'r',
     'cont': 0, 'title': 'IVP convergence (t vs r)',
     'panel_label': '(d)', 'get_stability': False},
]
fig, axes, line_cols = ode.plot_continuation_grid(
    plots, ncols=2, figsize=(10, 7),
    bifurcation_legend=False,
)
fig.suptitle('bi-exp QIF-SFA — four views', y=1.02, fontsize='large')
plt.show()

# %%
# A few things worth noting from the figure:
#
# * Panels (a) and (c) carry codim-1 bifurcation markers (folds, Hopfs)
#   placed by :meth:`plot_bifurcation_points`; the gray ``'dotted'``
#   segments mark unstable branches under the new 1.0.0 default colour
#   convention.
# * Panel (b) overlays the limit-cycle envelope — for limit cycles
#   :meth:`plot_continuation` automatically packs each period into a
#   (min, max) envelope, drawn as a pair of lines.
# * The trailing position in the 2×2 grid would be empty if we had passed
#   only three plot specs; :meth:`plot_continuation_grid` deletes those
#   unused axes so ``tight_layout`` lays out the panels cleanly.
#
# The same recipe scales — pass ``ncols=3`` and a longer list to compare
# six or nine continuations side by side.

# %%
# Step 3: `plot_trajectory(colorbar=True)` — 3D phase with colorbar
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# For a 3-element ``variables`` argument, :meth:`plot_trajectory` draws
# the phase-space curve in 3D with a colourmap encoding one of the
# projected coordinates (``array='x'`` by default — i.e. the first
# variable). ``colorbar=True`` attaches a colorbar to the figure that
# surfaces this mapping; ``colorbar_label`` sets a custom label
# (otherwise the variable name is used).

# Pass ``point=None`` so :meth:`extract` returns the full DataFrame of
# the IVP — each row is one labelled time point and the (r, v, A) columns
# trace the orbit through 3D phase space. Here we use the initial IVP
# that the `init_cont=True` constructor ran (cont=0) — the trajectory
# starts at the initial condition and converges to the steady state at
# the model's starting :math:`\bar\eta = -3`.
fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111, projection='3d')
ode.plot_trajectory(
    variables=['r', 'v', 'A'],
    cont=0,
    ax=ax, colorbar=True,
)
ax.set_title(r'IVP phase trajectory in $(r,\, v,\, A)$')
plt.tight_layout()
plt.show()

# %%
# The colorbar surfaces how the curve traverses the firing-rate
# coordinate — useful for reading the direction of flow at a glance,
# since matplotlib's default 3D rendering can otherwise make it hard to
# tell which loop comes first in time. To colour by a different
# coordinate, pass ``array='y'`` or ``array='z'`` through ``**kwargs``;
# the colorbar label updates accordingly.
#
# When ``cutoff`` is set, the helper discards trajectory data below that
# time threshold so initial transients don't squash the visible range.

# %%
# Step 4: Default-style change — gray for unstable
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Pre-1.0 PyCoBi rendered both stable and unstable line segments as solid
# black, distinguishable only by the linestyle (``'solid'`` vs
# ``'dotted'``). 1.0.0 changes the default unstable colour to ``'gray'``
# so the two regimes are visually distinguishable at a glance. The
# default applies wherever :meth:`_get_line_collection` is called
# internally — :meth:`plot_continuation`, :meth:`plot_trajectory`,
# :meth:`plot_timeseries`, and :meth:`plot_continuation_grid`.
#
# To revert to the legacy all-black appearance, pass
# ``line_color_unstable='k'`` through ``**kwargs`` on any of the plot
# methods:
#
# .. code-block:: python
#
#     ode.plot_continuation(
#         'p/qif_biexp_sfa_op/eta', 'r', cont='eta_eq',
#         line_color_unstable='k',  # legacy 0.x default
#     )

# %%
# Step 5: Clean up
# ^^^^^^^^^^^^^^^^

ode.close_session(clear_files=True)
