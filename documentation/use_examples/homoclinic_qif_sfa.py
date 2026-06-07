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
  post-process the result to flag every PSI zero-crossing as a
  ``'SNIC'`` bifurcation in the bifurcation column.

  .. note::
     The ``'SNIC'`` label is the historical PyCoBi name and *not* the
     standard codim-1 SNIC bifurcation (saddle-node on invariant cycle).
     With the default ``IPSI=(15, 16)`` we are detecting the **codim-2
     non-central SNIC** of Nechyporenko, Ashwin & Tsaneva-Atanasova 2026
     (arXiv:2412.12298): a homoclinic orbit to a saddle-node equilibrium
     whose return trajectory enters the saddle-node along its stable
     manifold rather than along the central (zero-eigenvalue) direction.
     This is a *vertex* terminating a curve of codim-1 SNIC bifurcations,
     not the codim-1 SNIC itself.

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
.. [2] K. Nechyporenko, P. Ashwin, K. Tsaneva-Atanasova (2026) *A novel
   route to oscillations via non-central SNICeroclinic bifurcation:
   unfolding the separatrix loop between a saddle-node and a saddle.*
   SIAM J. Appl. Dyn. Syst., to appear.  arXiv:2412.12298.
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
# UZR pins :math:`\bar\eta = -5.394` (a value slightly different from the
# LC's homoclinic terminus at :math:`\bar\eta \approx -5.39358`, see the
# note in Step 4) so the saddle equilibrium's state values are available
# as labelled UZ solutions later when seeding HomCont.
eta_sols, _ = ode.run(
    starting_point='EP2', name='eta_branch',
    c='eq', ICP='p/qif_biexp_sfa_op/eta', bidirectional=True,
    RL0=-8.0, RL1=2.0,
    NMX=2000, NPR=10, DS=1e-4, DSMIN=1e-8, DSMAX=5e-2,
    ITMX=40, ITNW=40, NWTN=12, NTST=400, NCOL=4,
    UZR={'p/qif_biexp_sfa_op/eta': [-5.394]},
)
print("eta scan:", dict(eta_sols['bifurcation'].value_counts()))

# %%
# Step 2: continue the limit cycle from HB2 toward the homoclinic
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# The Hopf at the upper edge of the spiking regime gives birth to an
# unstable limit cycle whose period grows as :math:`\bar\eta` walks the
# branch toward a saddle-loop homoclinic.  At :math:`\alpha=0.5,\ \Delta=2.0,
# \ J=15\sqrt{\Delta},\ \tau_r=\tau_d=10` the LC reaches the homoclinic
# at :math:`\bar\eta \approx -5.394` with a period of ~388.
#
# .. note::
#    For high-period limit-cycle continuations approaching a homoclinic,
#    auto-07p's Newton tolerances and bifurcation-detection tolerance
#    need to be tighter than the global defaults — otherwise the test
#    functions used to flag ``LP`` / ``PD`` / ``BP`` along the LC become
#    dominated by numerical noise and auto-07p reports spurious
#    bifurcations.  PyCoBi's ``'lc'`` scenario ships with
#    ``EPSL = EPSU = 1e-7`` and ``EPSS = 1e-5``.  For this LC the
#    progression is:
#
#    * Defaults (``NTST=50, EPSL=1e-7``): ~24 LPs, ~8 BPs, ~4 PDs.
#    * ``NTST=400, EPSL=1e-9``: ~24 LPs drop to <10, but 3 PDs remain.
#    * ``NTST=800, EPSL=1e-9`` (the values below): all PDs disappear.
#
#    The remaining ~800 BP labels concentrate at the homoclinic
#    :math:`\bar\eta` and are **not** numerical noise — they're auto-07p
#    detecting the LC branch nearly intersecting the saddle equilibrium's
#    branch, which is exactly what a saddle-loop homoclinic *is*.  We
#    keep them and ignore them in the 1D plot.

lc_sols, _ = ode.run(
    starting_point='HB2', name='lc_branch',
    c='lc', origin='eta_branch',
    ICP=['p/qif_biexp_sfa_op/eta', 11],
    NMX=8000, NPR=20, DS=1e-3, DSMIN=1e-12, DSMAX=5e-2,
    EPSL=1e-9, EPSU=1e-9, EPSS=1e-7,
    NTST=800, NCOL=4,
    bidirectional=False, get_period=True,
)
periods = np.asarray(lc_sols[('PAR(11)', '')], dtype=float)
print(f"LC bifurcations (raw): {dict(lc_sols['bifurcation'].value_counts())}")
print(f"LC period grows from {periods.min():.1f} to {periods.max():.1f}")

# Heuristically tag the homoclinic terminus as 'HC' on the LC summary:
# the max-to-min period ratio is ~22 here, well above the default 5x
# threshold for the heuristic to fire.
hc_idx = ode.label_homoclinic_terminus('lc_branch')
if hc_idx is not None:
    _eta_col = ('eta', '') if ('eta', '') in lc_sols.columns \
        else ('p/qif_biexp_sfa_op/eta', '')
    print(f"HC marker placed at row {hc_idx}: "
          f"eta = {float(lc_sols[_eta_col].iloc[hc_idx]):+.5f}, "
          f"period = {periods[hc_idx]:.1f}")

# 1D bifurcation diagram of the eta scan + LC envelope + HC marker.
fig, ax = plt.subplots(figsize=(9, 4))
ode.plot_continuation('p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/r',
                      cont='eta_branch', ax=ax,
                      line_color_stable='#1F77B4', line_color_unstable='#1F77B4',
                      bifurcation_legend=False,
                      ignore=['UZ', 'BP', 'EP'],
                      label='equilibrium')
ode.plot_continuation('p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/r',
                      cont='lc_branch', ax=ax,
                      line_color_stable='#FF7F0E', line_color_unstable='#FF7F0E',
                      bifurcation_legend=False,
                      ignore=['UZ', 'BP', 'EP', 'RG'],
                      label='limit cycle')
# Plot the HC marker directly — plot_continuation's bifurcation table
# doesn't know about our custom 'HC' label yet, so we drop a marker
# manually at the labelled row.
if hc_idx is not None:
    eta_col = ('eta', '') if ('eta', '') in lc_sols.columns \
        else ('p/qif_biexp_sfa_op/eta', '')
    r_col = [c for c in lc_sols.columns
             if isinstance(c, tuple) and c[0] == 'r'][0]
    eta_hc = float(lc_sols[eta_col].iloc[hc_idx])
    r_hc = lc_sols[r_col].iloc[hc_idx]
    if hasattr(r_hc, '__len__'):
        r_hc = float(np.max(r_hc))
    else:
        r_hc = float(r_hc)
    ax.scatter([eta_hc], [r_hc], marker='X', s=120,
               c='#D62728', edgecolor='k', linewidth=0.8, zorder=10,
               label=f'HC (period ≈ {periods[hc_idx]:.0f})')
ax.set_xlabel(r'$\bar\eta$')
ax.set_ylabel(r'$r$')
ax.set_title('1D bifurcation diagram — LC born at HB2 terminates at HC')
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
# the unstable middle equilibrium whose ``r`` coordinate is between the LC
# orbit's min and median ``r`` (i.e. the "slow-point" the orbit hovers near).
#
# .. note::
#    The UZR :math:`\bar\eta = -5.394` doesn't *exactly* match the LC's
#    EP1 :math:`\bar\eta = -5.39358` — and that small mismatch is
#    intentional.  Reading the saddle off the LC's orbit slow-point puts
#    us AT the codim-2 non-central SNIC vertex (PSI(15) ≈ 0), where
#    HomCont has no direction to advance from.  A saddle taken from the
#    eta-branch at a nominally-nearby :math:`\bar\eta` gives a slightly
#    off-vertex starting orbit (PSI(15) of order 0.3), and HomCont can
#    walk the codim-1 homoclinic curve forward.

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
#    ``'SNIC'`` in the bifurcation column.  Per the docstring note above:
#    these are codim-2 **non-central SNIC** bifurcations (Nechyporenko et
#    al. 2026, arXiv:2412.12298) — endpoints of a curve of standard SNIC
#    bifurcations — not the codim-1 SNIC itself.
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
    print(f"\n{n_snic} non-central SNIC bifurcation(s) detected on the "
          f"homoclinic curve (Nechyporenko et al. 2026, arXiv:2412.12298):")
    snic_rows = snic_sols[snic_sols[snic_col] == 'SNIC']
    eta_col = [c for c in snic_sols.columns
               if isinstance(c, tuple) and 'eta' in c[0]][0]
    delta_col = [c for c in snic_sols.columns
                 if isinstance(c, tuple) and 'Delta' in c[0]][0]
    for _, row in snic_rows.iterrows():
        print(f"  eta = {float(row[eta_col]):+.5f}, Delta = {float(row[delta_col]):.4f}")
hom_curve_available = len(snic_sols) > 2

# %%
# Step 6: codim-2 continuation of HB / LP curves in :math:`(\bar\eta,\, \Delta)`
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# PyCoBi's :func:`codim2_search` wraps the standard auto-07p 2-parameter
# continuation of codim-1 starting points.  Tracking the four equilibrium
# bifurcations we located on ``eta_branch`` (``HB1``, ``HB2``, ``LP1``,
# ``LP2``) in :math:`(\bar\eta,\, \Delta)` gives the codim-1 backbone we
# need to plot alongside the homoclinic locus from Step 5.
#
# .. note::
#    ``NPR=5`` is deliberately tighter than the typical ``NPR=20`` you'd
#    use for a 1-parameter scan.  ``NPR`` is auto-07p's stored-point
#    spacing: with the default 20, the first stored row on a codim-2
#    branch lies ~20 continuation steps (up to ~1 unit of arclength at
#    ``DSMAX=5e-2``) past the LP/HB starting point on ``eta_branch``,
#    leaving a visible gap in the :math:`(\bar\eta,\, \Delta)` plot
#    between each codim-1 source on ``eta_branch`` and the start of its
#    own codim-2 curve.  Lowering to ``NPR=5`` shrinks the first-step
#    offset to ~0.01 units — visually contiguous on the plot.

from pycobi.automated_continuation import codim2_search

codim2_curves = []
for sp in ('LP1', 'LP2', 'HB1', 'HB2'):
    bif_type = 'fold' if sp.startswith('LP') else 'Hopf'
    for ds in (1e-3, -1e-3):
        try:
            res = codim2_search(
                pyauto_instance=ode,
                params=['p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/Delta'],
                starting_points=[sp],
                origin='eta_branch',
                name=f'codim2_{sp}_{"pos" if ds > 0 else "neg"}',
                max_recursion_depth=0,
                NMX=3000, NPR=5, DSMIN=1e-9, DSMAX=5e-2, DS=ds,
                RL0=-15.0, RL1=5.0,
                bidirectional=False,
                UZSTOP={'p/qif_biexp_sfa_op/Delta': [0.0, 4.0]},
            )
            curve_name = next(iter(res))
            codim2_curves.append((curve_name, bif_type))
        except Exception as exc:
            print(f"  codim-2 {sp} DS={ds:+g} skipped: "
                  f"{type(exc).__name__}: {str(exc)[:80]}")

# %%
# Step 7: 2D bifurcation diagram with all the codim-1 curves + homoclinic
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Stack everything in the :math:`(\bar\eta,\, \Delta)` plane: fold curves
# (blue), Hopf curves (orange), homoclinic curve (red), and any detected
# non-central SNIC points (red stars).  The result is the complete
# codim-1 bifurcation portrait of the QIF-SFA model in this slice.

FOLD_COLOR = '#1F77B4'
HOPF_COLOR = '#FF7F0E'
HOM_COLOR  = '#D62728'

fig, ax = plt.subplots(figsize=(9, 6))
fold_labelled, hopf_labelled = False, False
for curve_name, bif_type in codim2_curves:
    color = FOLD_COLOR if bif_type == 'fold' else HOPF_COLOR
    if bif_type == 'fold':
        label = None if fold_labelled else 'fold of equilibria'
        fold_labelled = True
    else:
        label = None if hopf_labelled else 'Hopf'
        hopf_labelled = True
    try:
        ode.plot_continuation(
            'p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/Delta',
            cont=curve_name, ax=ax,
            line_color_stable=color, line_color_unstable=color,
            line_style_stable='solid', line_style_unstable='solid',
            bifurcation_legend=False, get_stability=False,
            ignore=['LP', 'HB', 'UZ', 'EP'],
            label=label,
        )
    except (KeyError, ValueError):
        pass

if hom_curve_available:
    ode.plot_continuation(
        'p/qif_biexp_sfa_op/eta', 'p/qif_biexp_sfa_op/Delta',
        cont='homoclinic', ax=ax,
        line_color_stable=HOM_COLOR, line_color_unstable=HOM_COLOR,
        line_style_stable='solid', line_style_unstable='dashed',
        bifurcation_legend=False, get_stability=False,
        ignore=['UZ', 'EP', 'RG'],
        label='homoclinic',
    )
    # Plant explicit markers at the codim-2 non-central SNIC points
    # (the PSI(15)/(16) zero-crossings flagged as 'SNIC' in Step 5).
    snic_col = ('bifurcation', '')
    eta_col = [c for c in snic_sols.columns
               if isinstance(c, tuple) and 'eta' in c[0]][0]
    delta_col = [c for c in snic_sols.columns
                 if isinstance(c, tuple) and 'Delta' in c[0]][0]
    snic_rows = snic_sols[snic_sols[snic_col] == 'SNIC']
    if len(snic_rows):
        ax.scatter(snic_rows[eta_col].to_numpy(dtype=float),
                   snic_rows[delta_col].to_numpy(dtype=float),
                   marker='*', s=200, c=HOM_COLOR,
                   edgecolor='k', linewidth=0.8, zorder=10,
                   label='non-central SNIC')

ax.set_xlabel(r'$\bar\eta$')
ax.set_ylabel(r'$\Delta$')
ax.set_xlim(-12.0, 5.0)
ax.set_ylim(0.0, 4.0)
ax.set_title('Complete codim-1 bifurcation portrait in '
             r'$(\bar\eta,\, \Delta)$')
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
