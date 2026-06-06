r"""
Homoclinic Continuation in a QIF Mean-Field Model
===================================================

This example demonstrates PyCoBi's support for auto-07p's HomCont extension
(:math:`IPS=9`, AUTO manual Ch. 20) via three new entry points:

* The ``'hom'`` scenario in ``auto_constants=(...)`` generates a ``c.hom``
  file pre-configured with HomCont's eight specific constants (``NUNSTAB``,
  ``NSTAB``, ``IEQUIB``, ``ITWIST``, ``ISTART``, ``IREV``, ``IFIXED``,
  ``IPSI``).  These flow through ``from_template(...)`` /
  ``run(...)`` as ordinary kwargs.
* :meth:`ODESystem.extract_orbit_to_dat` writes a labelled limit-cycle
  profile from auto-07p's solution object straight into the ``.dat`` format
  HomCont reads when ``ISTART=1``.
* :meth:`ODESystem.continue_homoclinic` orchestrates the full pipeline:
  extract the near-homoclinic orbit, seed a HomCont continuation in two
  parameters, append the requested PSI test-function PARs
  (``PAR(20 + IPSI[j])``) to ``ICP`` so they land in the summary, then
  post-process the result to flag every PSI zero-crossing as a custom
  ``'SNIC'`` bifurcation in the bifurcation column.

The model is the bi-exponential QIF mean field with spike-frequency
adaptation from the :ref:`Automated Codim-2 Search` example.  At
:math:`\tau_r = \tau_d = 10` (the alpha-kernel limit, Gast 2020 [1]_)
the upper Hopf (``HB2``) gives birth to a stable limit cycle whose period
diverges as :math:`\bar\eta` is pushed toward a saddle-loop homoclinic.
We pick the near-homoclinic LC up at its endpoint and continue the
homoclinic curve in :math:`(\bar\eta,\, \Delta)`.

References
^^^^^^^^^^

.. [1] R. Gast, H. Schmidt, T.R. Knösche (2020) *A Mean-Field Description of
   Bursting Dynamics in Spiking Neural Networks with Short-Term Adaptation.*
   Neural Computation 32 (9): 1615-1634.
"""

# %%
# Step 1: load the QIF-SFA model + run the equilibrium continuation
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Same setup as ``automated_continuation.py``: alpha-kernel-equivalent
# regime, four-dimensional state ``(r, v, A, B)``.  ``auto_constants``
# now also requests ``'hom'`` so the ``c.hom`` file is ready for use.

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
    auto_constants=('ivp', 'eq', 'lc', 'hom'),
)

# Equilibrium continuation in eta — locates HB1 / HB2 / LP1 / LP2.
eta_sols, _ = ode.run(
    starting_point='EP2', name='eta_branch',
    c='eq', ICP='p/qif_biexp_sfa_op/eta', bidirectional=True,
    RL0=-8.0, RL1=2.0,
    NMX=2000, NPR=10, DS=1e-4, DSMIN=1e-8, DSMAX=5e-2,
    ITMX=40, ITNW=40, NWTN=12, NTST=400, NCOL=4,
)
print("eta scan:", dict(eta_sols['bifurcation'].value_counts()))

# %%
# Step 2: continue the limit cycle from HB2 toward the homoclinic
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# The Hopf at the upper edge of the spiking regime gives birth to a stable
# limit cycle.  As :math:`\bar\eta` is continued in the positive direction,
# the orbit's period grows: the LC is approaching a homoclinic to a saddle
# equilibrium.  We push the continuation hard (``NMX=5000``, very small
# ``DSMIN``) to reach a high-period near-homoclinic profile that HomCont
# can take over from.

lc_sols, _ = ode.run(
    starting_point='HB2', name='lc_branch',
    c='lc', origin='eta_branch',
    ICP=['p/qif_biexp_sfa_op/eta', 11],
    NMX=5000, NPR=10, DS=1e-3, DSMIN=1e-12, DSMAX=5e-2,
    bidirectional=False, get_period=True,
)
periods = np.asarray(lc_sols[('PAR(11)', '')], dtype=float)
print(f"LC bifurcations: {dict(lc_sols['bifurcation'].value_counts())}")
print(f"LC period grows from {periods.min():.1f} to {periods.max():.1f}")

# Quick 1D bifurcation diagram of the eta scan + LC envelope.
fig, ax = plt.subplots(figsize=(8, 4))
ode.plot_continuation('p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/r',
                      cont='eta_branch', ax=ax,
                      line_color_stable='#1F77B4', line_color_unstable='#1F77B4',
                      bifurcation_legend=False, ignore=['UZ', 'BP', 'EP'],
                      label='equilibrium')
ode.plot_continuation('p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/r',
                      cont='lc_branch', ax=ax,
                      line_color_stable='#FF7F0E', line_color_unstable='#FF7F0E',
                      bifurcation_legend=False, ignore=['UZ', 'BP', 'EP'],
                      label='limit cycle')
ax.set_xlabel(r'$\bar\eta$')
ax.set_ylabel(r'$r$')
ax.set_title('1D bifurcation diagram — LC born at HB2 approaches a homoclinic')
ax.legend(loc='best')
plt.tight_layout()
plt.show()

# %%
# Step 3: extract the near-homoclinic orbit and write it as a .dat file
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Auto-07p's HomCont uses an existing orbit as the initial guess; the
# orbit lives on disk as a whitespace-separated ``.dat`` file of (time,
# state, state, ...) columns.  :meth:`extract_orbit_to_dat` reads the
# per-mesh state values straight off auto-07p's solution object and writes
# the ``.dat`` for us — no need to have run the LC continuation with
# ``get_timeseries=True``.
#
# Inspecting the seed profile confirms it has the characteristic
# "slow near saddle, fast excursion" shape of a near-homoclinic orbit.

seed_path = ode.extract_orbit_to_dat(
    cont='lc_branch', point='EP1', path='lc_seed', n_points=201,
)
data = np.loadtxt(seed_path)
print(f"seed written to {seed_path.name}: shape {data.shape}")
print(f"r profile: min={data[:,1].min():.3f}, max={data[:,1].max():.3f}, "
      f"median={np.median(data[:,1]):.3f}")

fig, ax = plt.subplots(figsize=(8, 3))
ax.plot(data[:, 0], data[:, 1], color='#FF7F0E', lw=1.5)
ax.set_xlabel(r'$t / T$')
ax.set_ylabel(r'$r(t)$')
ax.set_title(f'Near-homoclinic orbit profile (period $\\approx${periods.max():.0f})')
plt.tight_layout()
plt.show()

# %%
# Step 4: continue the homoclinic curve in :math:`(\bar\eta,\, \Delta)`
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# :meth:`continue_homoclinic` handles four things in one call:
#
# 1. Pulls the seed orbit out of ``lc_branch`` via :meth:`extract_orbit_to_dat`
#    (reusing the file we just inspected — pass ``dat_basename='lc_seed'``).
# 2. Reads the orbit's parameter values off the auto solution object and
#    forwards them via ``PAR={...}`` so HomCont's Newton step starts at the
#    correct point in parameter space (without this, ``STPNT``'s YAML
#    defaults clobber :math:`\bar\eta = -4.6` back to :math:`-8`).
# 3. Appends ``PAR(20 + IPSI[j])`` for each chosen PSI test function to
#    ``ICP`` so the test-function values are recorded along the branch.
# 4. After the run, scans those PAR columns for sign changes and writes
#    ``'SNIC'`` into the bifurcation column at each crossing — the
#    non-central-homoclinic-to-saddle-node detection that motivates this
#    whole machinery (AUTO §20.5; ``HOMCONT_PSI_NAMES[15]`` /
#    ``[16]`` for the meaning of those two specifically).
#
# .. note::
#    HomCont convergence is sensitive to (a) the seed orbit being close
#    enough to a true homoclinic, (b) the eigenvalue-split kwargs
#    ``NUNSTAB`` / ``NSTAB`` matching the saddle's manifold structure,
#    and (c) the ``IEQUIB`` flag.  For the QIF-SFA model at this
#    parameter regime the saddle is 4-dimensional with one unstable +
#    three stable directions (``NUNSTAB=1, NSTAB=3``).  If you adapt
#    this recipe to a different model, expect to spend some time tuning
#    these — auto-07p's HomCont demos (``demos/cir``, ``demos/she``)
#    are the canonical references.

try:
    snic_sols, _ = ode.continue_homoclinic(
        origin='lc_branch', starting_point='EP1',
        ICP=['p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/Delta'],
        NUNSTAB=1, NSTAB=3,
        IEQUIB=1, ITWIST=0,
        IPSI=(1, 4, 15, 16),     # saddle-loop test funcs + SNIC indicators
        n_points=201,
        NMX=500, NPR=10, DSMAX=0.1, DS=0.01,
        RL0=-15.0, RL1=5.0,
        UZSTOP={'p/qif_biexp_sfa_op/Delta': [0.01, 4.0]},
        name='homoclinic',
    )
    bifs = dict(snic_sols['bifurcation'].value_counts())
    print(f"HomCont curve: {bifs}")
    n_snic = bifs.get('SNIC', 0)
    if n_snic:
        print(f"\n{n_snic} non-central-homoclinic-to-saddle-node point(s) detected:")
        print(snic_sols[snic_sols['bifurcation'] == 'SNIC'][
            [('p/qif_biexp_sfa_op/eta', ''), ('p/qif_biexp_sfa_op/Delta', '')]
        ])
    hom_curve_available = len(snic_sols) > 2
except Exception as exc:
    print(f"HomCont continuation did not converge: "
          f"{type(exc).__name__}: {exc}")
    print("This is expected when the seed orbit isn't yet close enough to a "
          "true homoclinic, or when (NUNSTAB, NSTAB, IEQUIB) don't match the "
          "saddle structure of the model.  Push NMX higher on the LC step "
          "above to get a higher-period seed, or try different eigenvalue "
          "splits — see the AUTO demos/cir tutorial for a worked example "
          "of HomCont parameter tuning.")
    hom_curve_available = False

# %%
# Step 5: overlay the homoclinic curve on the 2D bifurcation diagram
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# When the HomCont continuation produced a meaningful curve, plot it in the
# :math:`(\bar\eta,\, \Delta)` plane next to the codim-1 fold and Hopf
# curves a codim-2 search would produce.  Any ``'SNIC'`` labels show up
# as red stars — they pin the SNIC bifurcation points along the
# homoclinic locus.

if hom_curve_available:
    fig, ax = plt.subplots(figsize=(8, 5))
    ode.plot_continuation(
        'p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/Delta',
        cont='homoclinic', ax=ax,
        line_color_stable='#D62728', line_color_unstable='#D62728',
        line_style_stable='solid', line_style_unstable='dashed',
        bifurcation_legend=False, get_stability=False,
        ignore=['UZ', 'EP', 'RG'],
        label='homoclinic curve',
    )
    ax.set_xlabel(r'$\bar\eta$')
    ax.set_ylabel(r'$\Delta$')
    ax.set_title('Homoclinic curve in the $(\\bar\\eta,\\, \\Delta)$ plane')
    ax.legend(loc='best')
    plt.tight_layout()
    plt.show()

# %%
# Reference table — HomCont PSI test functions
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# The IPSI list above selects which "PSIHO" test functions auto-07p computes
# along the homoclinic continuation.  Their values land in
# ``PAR(20 + IPSI[j])`` and zero crossings flag codim-2 phenomena on the
# homoclinic curve.  PyCoBi keeps a reference table so users picking
# ``IPSI=[...]`` have the meanings close at hand:

for j, desc in ODESystem.HOMCONT_PSI_NAMES.items():
    marker = " ←" if j in (15, 16) else ""
    print(f"  PSI({j:>2}) → PAR({20+j}): {desc}{marker}")
