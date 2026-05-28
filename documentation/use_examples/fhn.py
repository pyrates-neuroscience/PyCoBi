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
fixed point for :math:`I_\\mathrm{ext} \\lesssim 0.34` and :math:`I_\\mathrm{ext} \\gtrsim 1.42` and a stable
limit cycle in between, bounded by two Hopf bifurcations.

In what follows we will:

1. Load the FHN model from a co-located YAML file and confirm that PyRates emitted analytical DFDU / DFDP
   entries into the generated Fortran.
2. Continue the steady state in :math:`I_\\mathrm{ext}` and locate the two Hopfs.
3. Branch-switch at HB1 to trace the periodic-solution branch (the limit cycle) and overlay it on the
   equilibrium diagram.
4. Compare the analytical and finite-difference Jacobian paths on the limit-cycle continuation — where the
   Newton step solves a larger BVP system than on the equilibrium branch, so the analytical path has more
   room to pay off than it does on a plain equilibrium scan.

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

import re
import time
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
n_dfdu = len(re.findall(r'\bdfdu\(\s*\d', src, re.IGNORECASE))
n_dfdp = len(re.findall(r'\bdfdp\(\s*\d', src, re.IGNORECASE))
print(f"DFDU entries in {src_path.name}: {n_dfdu}")
print(f"DFDP entries in {src_path.name}: {n_dfdp}")
ceq = (Path(ode.dir) / 'c.eq').read_text()
clc = (Path(ode.dir) / 'c.lc').read_text()
print("JAC = 1 in c.eq:", any('JAC = 1' in ln for ln in ceq.splitlines()))
print("JAC = 1 in c.lc:", any('JAC = 1' in ln for ln in clc.splitlines()))

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

eq_sols, eq_cont = ode.run(
    starting_point='EP2', name='I_branch',
    c='eq', ICP='Iext', bidirectional=True,
    NMX=500, NPR=10, DS=0.01, DSMAX=0.05,
)
print("\nEquilibrium-branch bifurcation counts:")
print(eq_sols['bifurcation'].value_counts())
print(eq_sols.loc[eq_sols['bifurcation'] == 'HB', ['bifurcation', 'Iext']])

ode.plot_continuation('Iext', 'v', cont='I_branch')
plt.title('FHN: equilibria as a function of $I_\\mathrm{ext}$')
plt.show()

# %%
# The bifurcation diagram shows the steady state's :math:`v` value as a function of :math:`I_\mathrm{ext}`.
# Solid line: stable; dotted gray: unstable. The green circles mark the two Hopf bifurcations
# (:math:`\mathrm{HB}_1 \approx 0.331` and :math:`\mathrm{HB}_2 \approx 1.419`) where the equilibrium
# loses / regains stability via a complex-conjugate pair of eigenvalues crossing the imaginary axis.

# %%
# Step 3: branch-switch at HB1 and continue the limit cycle
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# At a Hopf bifurcation a branch of periodic solutions emerges from the equilibrium. We switch onto that
# branch by passing the Hopf label as ``starting_point`` and setting auto-07p's switches:
# ``IPS=2`` (continue periodic solutions), ``ISW=-1`` (start from HB / PD / BP), and ``c='lc'`` which
# loads the limit-cycle preset PyRates wrote for us. ``ICP=['Iext', 11]`` declares Iext as the
# continuation parameter plus ``PAR(11)`` for the period (auto's convention for LC continuations).
#
# For FHN with these parameters the LC born at HB1 is *sub-critical*: a small-amplitude unstable cycle
# emerges, undergoes a fold-of-cycle (LP) at :math:`I_\mathrm{ext} \approx 0.324`, becomes stable, and
# then the stable branch sweeps across the oscillatory regime up to HB2. We use a larger ``NMX`` than
# the equilibrium scan to cover that whole sweep.

lc_sols, lc_cont = ode.run(
    origin='I_branch', starting_point='HB1', name='lc',
    c='lc', ICP=['Iext', 11], IPS=2, ISP=2, ISW=-1,
    NMX=5000, NPR=500, DS=0.01, DSMAX=0.5,
    get_period=True,
)
print("\nLimit-cycle branch bifurcation counts:")
print(lc_sols['bifurcation'].value_counts())

# %%
# Overlay the LC branch on the equilibrium diagram. For a limit cycle, ``plot_continuation`` packs each
# period into a (min, max) envelope, so the LC branch appears as a pair of curves bounding the
# oscillation amplitude in :math:`v`.

fig, ax = plt.subplots()
ode.plot_continuation('Iext', 'v', cont='I_branch', ax=ax)
ode.plot_continuation('Iext', 'v', cont='lc', ax=ax, ignore=['UZ', 'BP'])
ax.set_title('FHN: equilibrium + limit cycle in $I_\\mathrm{ext}$')
plt.show()

# %%
# Step 4: analytical vs. finite-difference Jacobian
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# The ``JAC=`` override is per-call, so we can re-run the same LC continuation with the FD fallback
# without re-instantiating the system. We time both on the limit-cycle branch because that's where the
# Newton step solves a much bigger BVP system than on the equilibrium scan (auto's BVP system has
# dimension ``NTST × NCOL × NDIM`` per Newton iteration — for the c.lc defaults of ``NTST=50``,
# ``NCOL=4``, and ``NDIM=2`` that's a 400-dim Jacobian per iteration, versus a 2x2 on the equilibrium
# branch).

t0 = time.perf_counter()
ode.run(
    origin='I_branch', starting_point='HB1', name='lc_jac',
    c='lc', ICP=['Iext', 11], IPS=2, ISP=2, ISW=-1,
    NMX=5000, NPR=5000, DS=0.01, DSMAX=0.5,
    JAC=1,  # analytical (the default since we passed analytical_jacobian=True)
)
t_jac = time.perf_counter() - t0

t0 = time.perf_counter()
ode.run(
    origin='I_branch', starting_point='HB1', name='lc_fd',
    c='lc', ICP=['Iext', 11], IPS=2, ISP=2, ISW=-1,
    NMX=5000, NPR=5000, DS=0.01, DSMAX=0.5,
    JAC=0,  # force finite-difference fallback
)
t_fd = time.perf_counter() - t0

print(f"\nLC continuation, analytical Jacobian: {t_jac:5.2f}s")
print(f"LC continuation, finite-difference   : {t_fd:5.2f}s")
print(f"speed-up                              : {t_fd / t_jac:5.2f}x")

# %%
# For a 2D system the speed-up is modest (the symbolic 2x2 Jacobian only saves us a handful of RHS
# evaluations per Newton iteration), but the analytical path scales much better with state-vector
# dimension. The finite-difference Jacobian needs ``NDIM`` RHS evaluations to fill each Newton column
# — :math:`O(\mathrm{NDIM}^2)` total — while the analytical path fills the matrix in a single
# symbolic sweep regardless of NDIM. For the medium-sized biophysical models PyCoBi was built for
# (NDIM in the 10–50 range), the same comparison typically cuts continuation time by an order of
# magnitude and improves Newton convergence near fold-of-cycle / period-doubling points where the
# finite-difference Jacobian becomes noisy.

# %%
# Step 5: clean up
# ^^^^^^^^^^^^^^^^

ode.close_session(clear_files=True)

# %%
# As usual, ``close_session(clear_files=True)`` removes the temporary Fortran sources and ``c.*`` files.
# Use ``clear_files=False`` to keep them for inspection.
