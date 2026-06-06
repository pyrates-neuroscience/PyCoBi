r"""
Boundary-Value Problems: the Gelfand-Bratu Fold
================================================

This example demonstrates PyCoBi's support for two-point boundary-value problems
(BVPs) — auto-07p's ``IPS=4`` mode. Two new ``ODESystem.from_template`` kwargs
let you express boundary and integral residuals directly in PyRates-name
space, and PyCoBi auto-derives ``NBC`` / ``NINT`` from the list length:

.. code-block:: python

    boundary_conditions=['u0_u', 'u1_u'],        # u(0) = 0, u(1) = 0
    integral_constraints=['u_u - par_intval'],   # ∫u dt = intval

For escape-hatch cases (PDE spatial BCs, custom Jacobians, anything the DSL
can't reach) the same call accepts ``bcnd_fortran=...`` / ``icnd_fortran=...``
plus an explicit ``nbc`` / ``nint``, and the raw block is dropped verbatim into
the generated ``BCND`` / ``ICND`` subroutines.

The model is the Gelfand-Bratu BVP [1]_:

.. math::

    u''(t) + \lambda\, e^{u(t)} = 0,  \qquad u(0) = u(1) = 0,
    \qquad \int_0^1 u\, \mathrm{d} t = c.

In first-order form,

.. math::

    u' = v, \qquad v' = -\lambda\, e^{u}.

The integral constraint pins a free direction so auto-07p can continue jointly
in :math:`(\lambda, c)`: starting at :math:`u \equiv 0,\ \lambda = c = 0`, the
Newton step grows :math:`u`, walks :math:`\lambda` up to the *Bratu fold*
:math:`\lambda^\star \approx 3.513830719` [2]_, then turns the branch back
toward :math:`\lambda \to 0` along the upper sheet where :math:`u` keeps
growing. This S-curve is the classical signature of the problem; beyond
:math:`\lambda^\star` no solution exists.

References
^^^^^^^^^^

.. [1] I. M. Gelfand (1963). *Some problems in the theory of quasilinear
   equations.* American Mathematical Society Translations 29, 295-381.
.. [2] D. A. Frank-Kamenetskii (1969). *Diffusion and Heat Exchange in
   Chemical Kinetics.* Plenum Press.
"""

# %%
# Step 1: build the model + the BVP residuals
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# The PyRates ``OperatorTemplate`` carries the first-order ODE.  We declare
# ``intval`` as a model parameter so the integral constraint can fix it.

import re
from pathlib import Path

import matplotlib.pyplot as plt
from pyrates import CircuitTemplate, NodeTemplate, OperatorTemplate

from pycobi import ODESystem

op = OperatorTemplate(
    name="bratu_op",
    equations=[
        "u' = v",
        "v' = -lam*exp(u)",
    ],
    variables={
        "u": "output(0.0)", "v": "variable(0.0)",
        "lam": 0.0, "intval": 0.0,
    },
)
node = NodeTemplate(name="bratu_node", operators=[op])
circuit = CircuitTemplate(name="bratu", nodes={"p": node})

# %%
# The boundary-condition DSL: each entry is a residual that should equal 0 at
# the solution. ``u0_<var>`` denotes the state at :math:`t=0`; ``u1_<var>``
# the state at :math:`t=1`. The integral-constraint DSL adds the prefixes
# ``u_<var>`` (mesh-point state), ``uold_<var>``, ``udot_<var>``,
# ``upold_<var>``; ``par_<param>`` reaches model parameters from either DSL.
# ``NBC = 2`` and ``NINT = 1`` are auto-derived from the list lengths.

ode = ODESystem.from_template(
    template=circuit,
    auto_dir="~/PycharmProjects/auto-07p",
    working_dir=str(Path(__file__).resolve().parent),
    init_cont=False,
    analytical_jacobian=True,
    auto_constants=("bvp",),
    boundary_conditions=["u0_u", "u1_u"],          # u(0) = 0, u(1) = 0
    integral_constraints=["u_u - par_intval"],     # ∫u dt = intval
    NTST=20, NCOL=4,
    NMX=400, NPR=200,
)

# %%
# Confirm the generated Fortran: ``BCND`` carries the two Dirichlet residuals,
# ``ICND`` carries the integral constraint, and ``c.bvp`` reflects the
# auto-derived ``NBC`` / ``NINT``.

src = (Path(ode.dir) / "system_equations.f90").read_text()
cbvp = (Path(ode.dir) / "c.bvp").read_text()
bcnd = re.search(r"subroutine bcnd[\s\S]*?end subroutine bcnd", src).group(0)
icnd = re.search(r"subroutine icnd[\s\S]*?end subroutine icnd", src).group(0)
print(bcnd)
print()
print(icnd)
print()
for ln in cbvp.splitlines():
    if any(key in ln for key in ("IPS", "NBC", "NINT", "JAC", "parnames")):
        print(ln)

# %%
# Step 2: continue the Bratu branch through the fold
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Start auto-07p directly from ``STPNT`` (the trivial solution at
# :math:`\lambda = 0,\ c = 0`) and continue jointly in :math:`(\lambda,\, c)`.
# The non-trivial profile :math:`u(t) = -\frac{\lambda}{2} t(t-1) + O(\lambda^2)`
# grows immediately as :math:`\lambda` increases — :math:`c = \int_0^1 u\,
# \mathrm{d}t` is its summary observable.  Auto-07p turns the branch at the
# Bratu fold (labelled ``LP1``) at :math:`\lambda^\star \approx 3.513830719`,
# then continues on the upper sheet where :math:`c` keeps growing as
# :math:`\lambda \to 0`.

sols, cont = ode.run(
    c="bvp",
    ICP=["lam", "intval"], name="bratu_branch",
    RL0=0.0, RL1=10.0,
    DS=0.2, DSMAX=0.2,
)
print("\nBifurcations encountered:")
print(sols["bifurcation"].value_counts())
lp_rows = sols.loc[sols["bifurcation"] == "LP", ["bifurcation", "lam", "intval"]]
if len(lp_rows):
    lam_star = float(lp_rows["lam"].iloc[0])
    print(f"\nBratu fold detected at lambda* = {lam_star:.6f}")
    print(f"(reference value: 3.513830719…   error: {abs(lam_star - 3.513830719):.2e})")

# %%
# Step 3: visualise the S-curve
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Plot :math:`c = \int_0^1 u\,\mathrm{d}t` versus :math:`\lambda`.  The branch
# leaves the origin, walks up to :math:`\lambda^\star`, and turns back —
# everything to the right of the fold is unreachable.

fig, ax = plt.subplots(figsize=(7, 5))
ode.plot_continuation("lam", "intval", cont="bratu_branch", ax=ax,
                      bifurcation_legend=True, default_size=8)
ax.axvline(3.513830719, color="grey", linestyle=":", linewidth=1,
           label=r"$\lambda^\star = 3.5138$")
ax.set_xlabel(r"$\lambda$")
ax.set_ylabel(r"$c = \int_0^1 u\,\mathrm{d}t$")
ax.set_title("Gelfand-Bratu BVP: branch and fold")
ax.legend()
plt.tight_layout()
plt.show()

# %%
# Step 5: when the DSL isn't enough — the raw-Fortran escape hatch
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# For BCs / ICs that the DSL can't express — PDE spatial discretisations,
# parameter-dependent residuals that need custom helper functions, full
# analytical Jacobian rows on the BC residuals — pass ``bcnd_fortran=`` (or
# ``icnd_fortran=``) plus an explicit ``nbc`` / ``nint``. Auto-07p's standard
# indexing applies: ``fb(i) = ...`` fills the i-th boundary residual,
# ``u0(j)`` / ``u1(j)`` reads the j-th state at the left / right boundary,
# and ``args(k)`` reads the k-th parameter (matching the order PyRates emits
# into ``parnames``). The same applies to ``icnd_fortran``, with the
# additional ``uold(j)`` / ``udot(j)`` / ``upold(j)`` arrays available inside
# the subroutine.
#
# Sketch — the same Dirichlet BC expressed as raw Fortran:
#
# .. code-block:: python
#
#     ode_raw = ODESystem.from_template(
#         template=circuit,
#         auto_constants=("bvp",),
#         bcnd_fortran=(
#             "fb(1) = u0(1)\n"   # u(0) = 0
#             "fb(2) = u1(1)\n"   # u(1) = 0
#         ),
#         nbc=2,
#         init_cont=False,
#         ...
#     )
#
# The DSL and the raw escape hatch are mutually exclusive on each routine;
# raw blocks must be paired with the corresponding ``nbc`` / ``nint`` so the
# emitted ``c.bvp`` matches the number of residuals the body fills.
