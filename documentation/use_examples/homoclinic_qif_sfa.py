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
        'p/qif_biexp_sfa_op/alpha': 0.5,
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
# UZR pins :math:`\bar\eta = -5.39406` (the homoclinic point we will
# discover from the LC's BP cluster in Step 2) so the saddle equilibrium's
# state values are available as labelled UZ solutions later in Step 4.
eta_sols, _ = ode.run(
    starting_point='EP2', name='eta_branch',
    c='eq', ICP='p/qif_biexp_sfa_op/eta', bidirectional=True,
    RL0=-8.0, RL1=2.0,
    NMX=2000, NPR=10, DS=1e-4, DSMIN=1e-8, DSMAX=5e-2,
    ITMX=40, ITNW=40, NWTN=12, NTST=400, NCOL=4,
    UZR={'p/qif_biexp_sfa_op/eta': [-5.39406]},
)
print("eta scan:", dict(eta_sols['bifurcation'].value_counts()))

# %%
# Step 2: continue the limit cycle from HB2 toward the homoclinic
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# The Hopf at the upper edge of the spiking regime gives birth to a limit
# cycle whose period grows as :math:`\bar\eta` walks the branch toward a
# saddle-loop homoclinic.  At the canonical Gast-2020 Fig. 5 parameter
# point (:math:`\alpha=1.0,\ \Delta=2.0,\ J=15\sqrt{\Delta}`) the LC has
# exactly two folds (``LP1`` near the upper turning point in :math:`\bar
# \eta`, ``LP2`` near the homoclinic side) and *no* period doublings.
#
# .. note::
#    For limit-cycle continuations that span large periods, the auto-07p
#    Newton tolerances (``EPSL``, ``EPSU``) and the bifurcation-detection
#    tolerance (``EPSS``) need to be tightened relative to the global
#    defaults — otherwise the test functions used to flag ``LP`` /
#    ``PD`` / ``BP`` along the LC become dominated by numerical noise and
#    auto-07p reports double-digit *spurious* bifurcations (the tell-tale
#    symptom: an ``LP`` is detected but the ``stability`` column doesn't
#    flip around it; ``PD``s appear despite the analytical solution
#    having no period doubling).
#
#    PyCoBi's ``'lc'`` scenario ships with ``EPSL = EPSU = 1e-7`` and
#    ``EPSS = 1e-5``.  At those defaults this exact LC reports ~24 LPs
#    and 4 PDs.  Tightening to ``1e-9`` / ``1e-7`` with ``NTST = 400``
#    (the values below) gives **exactly 2 LPs and 0 PDs** — matching
#    Gast-2020 Fig. 5 — at the cost of a slower continuation.  Add
#    ``STOP=['LP3']`` so auto-07p doesn't waste steps in the noisy
#    near-homoclinic tail past the second fold.

lc_sols, _ = ode.run(
    starting_point='HB2', name='lc_branch',
    c='lc', origin='eta_branch',
    ICP=['p/qif_biexp_sfa_op/eta', 11],
    NMX=8000, NPR=20, DS=1e-3, DSMIN=1e-12, DSMAX=5e-2,
    EPSL=1e-9, EPSU=1e-9, EPSS=1e-7,
    NTST=400, NCOL=4,
    bidirectional=False, get_period=True,
    STOP=['LP3'],
)
periods = np.asarray(lc_sols[('PAR(11)', '')], dtype=float)
print(f"LC bifurcations: {dict(lc_sols['bifurcation'].value_counts())}")
print(f"LC period grows from {periods.min():.1f} to {periods.max():.1f}")
# Show the LPs and confirm they're real (stability flips) — reproduces
# the two-fold structure of Gast-2020 Fig. 5.
_bif = lc_sols[('bifurcation', '')].values
_stab = lc_sols[('stability', '')].values
_eta_col = ('eta', '') if ('eta', '') in lc_sols.columns \
    else ('p/qif_biexp_sfa_op/eta', '')
_eta = lc_sols[_eta_col].values
for i, b in enumerate(_bif):
    if str(b).strip() == 'LP':
        flip = (_stab[i-1] != _stab[i]) if i > 0 else None
        print(f"  LP @ eta={float(_eta[i]):+.5f}, period={periods[i]:7.2f}, "
              f"stability flip: {bool(flip)}")

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
# Step 4: find the saddle equilibrium (the homoclinic's target)
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# A homoclinic orbit by definition leaves and returns to a saddle equilibrium.
# At :math:`\bar\eta \approx -5.394` the QIF-SFA has three equilibria: a
# stable down-state, an unstable middle saddle, and a stable up-state.
# We pull the equilibrium values straight off the UZ labels we planted on
# the ``eta_branch`` continuation in Step 1 — the saddle is identified as
# the unstable middle equilibrium whose ``r`` coordinate matches the LC's
# minimum ``r`` (the orbit's slow-point).

saddle_state = None
for uz_label in ('UZ1', 'UZ2', 'UZ3'):
    try:
        s, _, _ = ode.get_solution(cont='eta_branch', point=uz_label)
        if hasattr(s, 'b') and isinstance(getattr(s, 'b', None), dict):
            s = s.b['solution']
        coords = {c: float(np.asarray(s[c]).ravel()[0]) for c in s.coordnames}
        # the saddle is the unstable middle equilibrium; r matches LC's min r
        if abs(coords['r'] - data[:, 1].min()) < 0.05:
            saddle_state = coords
            print(f"saddle equilibrium at {uz_label}: r={coords['r']:.4f}, "
                  f"v={coords['v']:.4f}, A={coords['A']:.4f}, B={coords['B']:.4f}")
            break
    except (KeyError, AttributeError):
        pass

# %%
# Step 5: continue the homoclinic curve in :math:`(\bar\eta,\, \Delta)`
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# :meth:`continue_homoclinic` now wraps the full LC-to-HomCont seed
# preparation pipeline:
#
# 1. ``extract_orbit_to_dat`` reads the labelled LC profile straight off
#    auto-07p's solution object.
# 2. ``phase_shift_to=saddle_state`` re-rolls the orbit so the point
#    closest to the saddle sits at :math:`t = 0` — without this, the LC's
#    arbitrary phase choice leaves HomCont's stable / unstable manifold
#    projections misaligned and Newton MXs at step 2.
# 3. ``saddle_state`` also populates ``PAR(11 + i)`` with the saddle's
#    coordinates (one per state variable, in uname order).  Combined with
#    ``IEQUIB = 0``, this pins HomCont's equilibrium at the precomputed
#    saddle rather than asking Newton to find it.
# 4. The seed's parameter values (from the auto solution's PAR vector)
#    are forwarded via ``PAR={...}`` so the homoclinic starts at the
#    correct :math:`(\bar\eta, \Delta, \alpha, \dots)`.
# 5. After the run, ``_flag_psi_zero_crossings`` scans ``PAR(35) = PSI(15)``
#    and ``PAR(36) = PSI(16)`` for sign changes and marks each as
#    ``'SNIC'`` (non-central homoclinic to saddle-node).
#
# .. note::
#    Eigenvalue-split kwargs (``NUNSTAB``, ``NSTAB``) still need to match
#    the model's saddle structure; for this 4-D QIF-SFA the saddle has one
#    unstable + three stable directions, hence ``NUNSTAB=1, NSTAB=3``.
#    Adapting the recipe to a different model means picking these from
#    the local eigenvalue spectrum of the saddle.

snic_sols, _ = ode.continue_homoclinic(
    origin='lc_branch', starting_point='EP1',
    ICP=['p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/Delta'],
    NUNSTAB=1, NSTAB=3,
    IEQUIB=0, ITWIST=0,
    IPSI=(15, 16),
    saddle_state=saddle_state,
    n_points=201,
    NTST=100, NCOL=4, IAD=1, ISP=0, ILP=0,
    NMX=400, NPR=10, IID=2, ITMX=10, ITNW=7, NWTN=3,
    DS=0.001, DSMIN=1e-5, DSMAX=0.05,
    UZSTOP={'p/qif_biexp_sfa_op/Delta': [0.01, 4.0]},
    name='homoclinic',
)
bifs = dict(snic_sols['bifurcation'].value_counts())
print(f"\nHomCont curve in (eta, Delta): {bifs}, {len(snic_sols)} points")
n_snic = bifs.get('SNIC', 0)
if n_snic:
    snic_col = ('bifurcation', '')
    print(f"\n{n_snic} non-central-homoclinic-to-saddle-node point(s) detected:")
    snic_rows = snic_sols[snic_sols[snic_col] == 'SNIC']
    eta_col = [c for c in snic_sols.columns
               if isinstance(c, tuple) and 'eta' in c[0]][0]
    delta_col = [c for c in snic_sols.columns
                 if isinstance(c, tuple) and 'Delta' in c[0]][0]
    for _, row in snic_rows.iterrows():
        print(f"  eta = {float(row[eta_col]):+.5f}, Delta = {float(row[delta_col]):.4f}")
hom_curve_available = len(snic_sols) > 2

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
