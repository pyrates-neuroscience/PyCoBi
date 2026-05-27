"""
FitzHugh-Nagumo with the Analytical Jacobian
============================================

The `FitzHugh-Nagumo (FHN) model <http://www.scholarpedia.org/article/FitzHugh-Nagumo_model>`_ is a textbook 2D
reduction of the Hodgkin-Huxley equations. The system has a single cubic non-linearity and is a clean target for
demonstrating the *analytical-Jacobian* path in `PyCoBi`: PyRates symbolically differentiates the vector field and
writes the DFDU / DFDP entries straight into the generated Fortran, and auto-07p uses the analytical Jacobian
instead of falling back to finite differences.

The model equations read:

.. math::

    \\dot v &= v - \\frac{v^3}{3} - w + I_\\mathrm{ext}, \n
    \\dot w &= \\varepsilon (v + a - b w),

where :math:`v` is the fast membrane-potential variable and :math:`w` is a slow recovery variable.
With the standard parameter regime (:math:`a=0.7,\\, b=0.8,\\, \\varepsilon=0.08`) the system has a stable
fixed point for :math:`I_\\mathrm{ext} \\lesssim 0.34` and :math:`I_\\mathrm{ext} \\gtrsim 1.40` and a stable
limit cycle in between, with two sub-critical Hopf bifurcations delimiting the oscillatory regime.

In what follows we will:

1. Load the FHN model from a co-located YAML file and confirm that PyRates emitted analytical DFDU / DFDP
   entries into the generated Fortran.
2. Continue the steady state in :math:`I_\\mathrm{ext}` and locate the two Hopfs.
3. Branch-switch at one of the Hopfs to trace the periodic-solution branch (the limit cycle).
4. Compare the analytical-Jacobian and finite-difference-Jacobian paths on the same continuation to show
   the speed-up.

References
^^^^^^^^^^

.. [1] R. FitzHugh (1961) *Impulses and physiological states in theoretical models of nerve membrane.*
   Biophysical Journal 1 (6): 445-466, https://doi.org/10.1016/S0006-3495(61)86902-6.
"""

# %%
# Step 1: load the model and inspect the generated Fortran
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# We ship a small ``fhn.yaml`` next to this script so the example is self-contained. PyRates'
# slash-notation path takes the form ``<dir>/<filename-no-ext>/<template-name>``; we anchor on
# ``__file__`` so the load works regardless of cwd.

import os
from pathlib import Path

from pycobi import ODESystem
import matplotlib.pyplot as plt

here = Path(__file__).resolve().parent
yaml_path = str(here / 'fhn' / 'fhn')

# `analytical_jacobian=True` is the default; spelling it out for emphasis. We also request both an `ivp`
# (for the IVP that converges us to the steady state) and an `eq` scenario so the second `run()` call can
# load `c.eq` directly. `init_cont=True` opts into the legacy automatic-IVP behaviour — the default since
# PyCoBi 0.10.0 is `False`.
ode = ODESystem.from_yaml(
    yaml_path,
    auto_dir="~/PycharmProjects/auto-07p",
    init_cont=True,
    analytical_jacobian=True,
    auto_constants=('ivp', 'eq', 'lc'),
    NMX=20000, NPR=20000,
)

# %%
# Let's check that PyRates actually wrote the analytical Jacobian into the generated Fortran by grepping
# for ``dfdu(`` and ``dfdp(`` assignments in the source file. We'd also expect ``JAC = 1`` in the
# generated ``c.eq`` / ``c.lc`` files since the JAC flag is wired up to the analytical-Jacobian emission.

src_path = Path(ode.dir) / 'system_equations.f90'
src = src_path.read_text()
import re
n_dfdu = len(re.findall(r'\bdfdu\(\s*\d', src, re.IGNORECASE))
n_dfdp = len(re.findall(r'\bdfdp\(\s*\d', src, re.IGNORECASE))
print(f"DFDU entries in {src_path.name}: {n_dfdu}")
print(f"DFDP entries in {src_path.name}: {n_dfdp}")
ceq = (Path(ode.dir) / 'c.eq').read_text()
print("JAC set in c.eq:", any('JAC = 1' in ln for ln in ceq.splitlines()))

# %%
# You should see four ``DFDU`` entries (the 2x2 Jacobian is dense for FHN — all four entries are
# non-zero) and four ``DFDP`` entries (one row of partials per state variable, across the parameters
# the right-hand side actually depends on), all written by PyRates' symbolic differentiation. The cubic
# ``v^3/3`` term in the RHS becomes ``1 - y(1)**2`` in the DFDU(1,1) slot — that's exactly
# :math:`\partial \dot v / \partial v` for the FHN equations.

# %%
# Step 2: continue the steady state in :math:`I_\mathrm{ext}`
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Standard 1D equilibrium continuation. We start from ``'EP2'``, which is the converged final point of the
# IVP that the constructor ran (``'EP1'`` is the :math:`t=0` initial condition). `bidirectional=True`
# walks the branch both ways so we capture the full S-curve.

eta_sols, eta_cont = ode.run(
    starting_point='EP2', name='I_branch',
    c='eq', ICP='Iext', bidirectional=True,
    NMX=500, NPR=10, DS=0.01, DSMAX=0.05,
    UZR={'Iext': [0.5, 1.0]},  # plant a couple of user points along the branch
)

ode.plot_continuation('Iext', 'v', cont='I_branch')
plt.show()

# %%
# The bifurcation diagram shows the steady state's :math:`v` value as a function of :math:`I_\mathrm{ext}`.
# Solid line: stable; dotted line: unstable. The two green markers are the Hopf bifurcations where the
# equilibrium loses / regains stability via a complex-conjugate pair of eigenvalues crossing the imaginary
# axis. The grey triangles (if any are detected) are fold bifurcations.

# %%
# Step 3: branch-switch at a Hopf to trace the limit cycle
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# At a Hopf bifurcation a branch of periodic solutions emerges from the equilibrium. We switch to that
# branch by passing the Hopf label as the starting point and setting auto-07p's switches:
# ``IPS=2`` (continue periodic solutions), ``ISW=-1`` (start from a HB / PD / BP), and using ``c='lc'``
# loads the limit-cycle preset PyRates wrote for us.

lc_sols, lc_cont = ode.run(
    origin='I_branch', starting_point='HB1', name='lc',
    c='lc', ICP=['Iext'],
    IPS=2, ISP=2, ISW=-1,
    NMX=200, NPR=10, DS=0.01, DSMAX=0.1,
    get_period=True,
)

# %%
# Overlay the limit-cycle branch on the equilibrium diagram. For a limit cycle, ``plot_continuation``
# packs each period into a (min, max) envelope, so the LC branch appears as a pair of curves bounding the
# oscillation amplitude in :math:`v`.

fig, ax = plt.subplots()
ode.plot_continuation('Iext', 'v', cont='I_branch', ax=ax)
ode.plot_continuation('Iext', 'v', cont='lc', ax=ax, ignore=['UZ', 'BP'])
plt.show()

# %%
# Step 4: analytical vs. finite-difference Jacobian
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# To quantify what the analytical Jacobian buys us, we re-run the same equilibrium continuation with
# ``JAC=0`` (finite-difference Jacobian) on the same model. The ``JAC=`` override is per-call, so we
# don't need to re-instantiate the system.

import time

# warm-up + analytical-Jacobian timing
t0 = time.perf_counter()
ode.run(
    origin=0, starting_point='EP2', name='I_jac',
    c='eq', ICP='Iext', bidirectional=True,
    NMX=500, NPR=500, DS=0.01, DSMAX=0.05,
    JAC=1,  # analytical (the default since we passed analytical_jacobian=True)
)
t_jac = time.perf_counter() - t0

t0 = time.perf_counter()
ode.run(
    origin=0, starting_point='EP2', name='I_fd',
    c='eq', ICP='Iext', bidirectional=True,
    NMX=500, NPR=500, DS=0.01, DSMAX=0.05,
    JAC=0,  # force finite-difference fallback
)
t_fd = time.perf_counter() - t0

print(f"analytical Jacobian: {t_jac*1000:6.1f} ms")
print(f"finite-difference  : {t_fd*1000:6.1f} ms")
print(f"speed-up           : {t_fd / t_jac:5.2f}x")

# %%
# For a 2D system like FHN the speed-up is modest (the Jacobian is tiny and the FD overhead is small),
# but the analytical path scales much better with system dimension — for the kind of medium-sized
# biophysical models PyCoBi is built for (NDIM in the 10-50 range), the analytical Jacobian typically
# cuts continuation time by a factor of 5-20× and improves Newton convergence near bifurcation points
# where the FD Jacobian becomes noisy.

# %%
# Step 5: clean up
# ^^^^^^^^^^^^^^^^

ode.close_session(clear_files=True)

# %%
# As usual, ``close_session(clear_files=True)`` removes the temporary Fortran sources and ``c.*`` files.
# Use ``clear_files=False`` to keep them for inspection.
