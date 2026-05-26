"""Test suite for the ODESystem class.
"""

# imports
import os
import re

from pycobi import ODESystem
from pyrates import clear
from pytest import fixture, approx
import numpy as np


# Matches an analytical-Jacobian assignment like ``dfdu(1,2) = ...`` (or DFDU);
# does NOT match the array declaration ``dfdu(ndim,ndim)`` in the func signature.
_DFDU_ASSIGN = re.compile(r'\bdfdu\(\s*\d', re.IGNORECASE)
_DFDP_ASSIGN = re.compile(r'\bdfdp\(\s*\d', re.IGNORECASE)

# meta infos
__author__ = "Richard Gast"
__status__ = "Development"

# Utility
#########


def setup_module():
    print("\n")
    print("==========================")
    print("| Test Suite : ODESystem |")
    print("==========================")


@fixture(scope="session")  # type: str
def auto_dir(pytestconfig) -> str:
    return pytestconfig.getoption("auto_dir")


# test accuracy
accuracy = 1e-4

# tests
#######


def test_1_1_init(auto_dir):
    """Testing the different instantiation options of the `ODESystem` class.
    """

    # initialize ODESystem the default way
    ode1 = ODESystem(eq_file='qif_eq', working_dir="resources", auto_dir=auto_dir, init_cont=False)
    ode1.close_session()

    # initialize ODESystem from YAML file
    model = "model_templates.neural_mass_models.qif.qif"
    ode2 = ODESystem.from_yaml(model, init_cont=False, file_name='qif_eq2', func_name="qif2", auto_dir=auto_dir)
    ode2.close_session(clear_files=True)

    # initialize ODESystem from YAML file with different parameters
    ode3 = ODESystem.from_yaml(model, init_cont=False, file_name='qif_eq3', func_name="qif3", auto_dir=auto_dir,
                               node_vars={'p/qif_op/eta': 2.0})
    ode3.close_session(clear_files=True)

    # these tests should pass
    assert isinstance(ode1, ODESystem)
    assert isinstance(ode2, ODESystem)
    assert isinstance(ode3, ODESystem)
    assert ode1.dir != ode2.dir
    assert ode2.dir == ode3.dir


def test_1_2_run(auto_dir):
    """Testing the run method for running auto commands via `ODESystem`.
    """

    # initialize ODESystem the default way
    ode1 = ODESystem(eq_file="qif_eq", working_dir="resources", auto_dir=auto_dir, init_cont=True, c="ivp", NPR=100)
    ode1.close_session()

    # initialize ODESystem from YAML file
    model = "model_templates.neural_mass_models.qif.qif"
    ode2 = ODESystem.from_yaml(model, init_cont=True, file_name='qif_eq2', func_name="qif2", auto_dir=auto_dir,
                               NPR=100, NMX=5000)
    ode2.close_session(clear_files=True)

    # initialize ODESystem from YAML file with different parameters
    ode3 = ODESystem.from_yaml(model, init_cont=True, file_name='qif_eq3', func_name="qif3", auto_dir=auto_dir,
                               node_vars={'p/qif_op/eta': 2.0}, NPR=100, NMX=5000)
    ode3.close_session(clear_files=True)

    print(ode2[0].columns.values)

    # PyRates >= 1.1 emits `unames` into c.* so auto exposes state variables under
    # their bare local name ('r') rather than 'U(1)'. ode1 uses the hand-written
    # qif_eq.f90 (no unames) → column 'U(1)'; ode2/ode3 are from_yaml → column 'r'.
    assert (ode1[0]["U(1)"] - ode2[0]["r"]).sum()[0] == approx(0.0, rel=accuracy, abs=accuracy)
    assert abs((ode2[0]["r"] - ode3[0]["r"]).sum()[0]) > 0


def test_2_1_jacobian_parity(auto_dir, tmp_path):
    """Equilibrium continuation should agree between PyRates' analytical Jacobian (JAC=1)
    and auto-07p's finite-difference Jacobian (JAC=0).

    This pins the PyCoBi <-> PyRates contract: if PyRates' DFDU/DFDP emission goes wrong
    (sign flip, missing parameter dependency, indexing bug), the analytical-Jac curve will
    diverge from the FD curve, and this test catches it.
    """
    model = "model_templates.neural_mass_models.qif.qif"

    # ---- analytical Jacobian (default) ----
    ode_jac = ODESystem.from_yaml(
        model,
        working_dir=str(tmp_path), auto_dir=auto_dir,
        file_name='qif_jac_mod', func_name='qif_jac',
        init_cont=True, analytical_jacobian=True,
        auto_constants=('ivp', 'eq'),
        NPR=100, NMX=5000,
    )

    # structural assertions: DFDU/DFDP entries emitted and JAC=1 in constants file
    f90_jac = (tmp_path / 'qif_jac_mod.f90').read_text()
    assert _DFDU_ASSIGN.search(f90_jac), "analytical_jacobian=True should emit dfdu(i,j) = ... entries"
    assert _DFDP_ASSIGN.search(f90_jac), "analytical_jacobian=True should emit dfdp(i,j) = ... entries"
    ceq_jac = (tmp_path / 'c.eq').read_text()
    assert 'JAC = 1' in ceq_jac, f"expected JAC = 1 in c.eq, got:\n{ceq_jac}"

    # equilibrium continuation in eta, starting from the IVP's *final* EP
    # (EP1 is the t=0 initial condition; EP2 is the converged steady state).
    sols_jac, _ = ode_jac.run(
        c='eq', ICP='eta', name='eta_jac',
        starting_point='EP2',
        DS=0.05, DSMAX=0.1, NMX=30,
    )
    ode_jac.close_session(clear_files=True)

    # ---- finite-difference Jacobian ----
    ode_fd = ODESystem.from_yaml(
        model,
        working_dir=str(tmp_path), auto_dir=auto_dir,
        file_name='qif_fd_mod', func_name='qif_fd',
        init_cont=True, analytical_jacobian=False,
        auto_constants=('ivp', 'eq'),
        NPR=100, NMX=5000,
    )

    f90_fd = (tmp_path / 'qif_fd_mod.f90').read_text()
    assert not _DFDU_ASSIGN.search(f90_fd), "analytical_jacobian=False should not emit dfdu(i,j) = ... entries"
    assert not _DFDP_ASSIGN.search(f90_fd), "analytical_jacobian=False should not emit dfdp(i,j) = ... entries"
    ceq_fd = (tmp_path / 'c.eq').read_text()
    assert 'JAC = 0' in ceq_fd, f"expected JAC = 0 in c.eq, got:\n{ceq_fd}"

    sols_fd, _ = ode_fd.run(
        c='eq', ICP='eta', name='eta_fd',
        starting_point='EP2',
        DS=0.05, DSMAX=0.1, NMX=30,
    )
    ode_fd.close_session(clear_files=True)

    # ---- numerical parity: r(eta) curves should match ----
    # PyRates emits bare local names via auto's `parnames`/`unames`, so the
    # DataFrame columns are 'eta' / 'r' (not 'p/qif_op/eta' / 'p/qif_op/r').
    eta_jac = np.asarray(sols_jac['eta'].values).squeeze()
    eta_fd = np.asarray(sols_fd['eta'].values).squeeze()
    r_jac = np.asarray(sols_jac['r'].values).squeeze()
    r_fd = np.asarray(sols_fd['r'].values).squeeze()

    # sort by eta (defensive — auto's stepper may visit points out of order at boundaries)
    order_jac = np.argsort(eta_jac)
    order_fd = np.argsort(eta_fd)
    eta_jac, r_jac = eta_jac[order_jac], r_jac[order_jac]
    eta_fd, r_fd = eta_fd[order_fd], r_fd[order_fd]

    # interpolate onto common eta grid and compare
    lo = max(eta_jac.min(), eta_fd.min())
    hi = min(eta_jac.max(), eta_fd.max())
    assert hi > lo, "no overlapping eta range between analytical and FD continuations"
    eta_grid = np.linspace(lo, hi, 20)
    r_jac_i = np.interp(eta_grid, eta_jac, r_jac)
    r_fd_i = np.interp(eta_grid, eta_fd, r_fd)
    assert np.allclose(r_jac_i, r_fd_i, atol=1e-3, rtol=1e-3), (
        f"analytical vs FD Jacobian continuations diverge:\n"
        f"  eta_grid = {eta_grid}\n  r_jac    = {r_jac_i}\n  r_fd     = {r_fd_i}"
    )


def test_3_1_pycobi_vs_auto_ops(auto_dir, tmp_path):
    """End-to-end parity vs auto-07p's own analytical-Jacobian demo.

    auto-07p ships the `ops` demo (a FitzHugh-Nagumo-like 3D ODE used in
    "optimization of periodic solutions") with a hand-derived DFDU/DFDP
    and JAC=1. We run PyCoBi over the demo's native Fortran source as
    the reference, then build the same model from a PyRates YAML and run
    the same equilibrium continuation in PAR(3) ('p3').

    The two paths agreeing means: (1) PyRates' symbolic DFDU/DFDP entries
    match the hand-derived ones; (2) PyCoBi's PyRates-driven setup yields
    the same auto-07p output as a raw demo run.

    We chose `ops` over the slightly simpler `abcb` because abcb's RHS
    uses `exp(x3)` and PyRates' equation parser currently leaves `exp`
    as an undefined function — `sympy.diff` then returns an unevaluated
    `Derivative(exp(x3), x3)` that gfortran refuses to compile.
    """
    import shutil

    demo_src = os.path.join(auto_dir, 'demos', 'ops')

    # ---- reference: auto-07p's native ops demo, driven through PyCoBi ----
    demo_dir = tmp_path / 'demo'
    demo_dir.mkdir()
    for fn in ('ops.f90', 'c.ops'):
        shutil.copy(os.path.join(demo_src, fn), demo_dir)
    ode_demo = ODESystem(
        eq_file='ops', working_dir=str(demo_dir),
        auto_dir=auto_dir, init_cont=False,
    )
    # `parnames={}, unames={}` clears auto-07p's runner cache: by design (see
    # auto/runAUTO.py config()) auto retains parnames/unames across run() calls
    # within a process and merges them with new c.* files. Without this, any
    # prior test that used a PyRates-generated model with unames/parnames would
    # leak its column names into ode_demo's DataFrame.
    sols_demo, _ = ode_demo.run(c='ops', name='ops_demo', parnames={}, unames={})
    ode_demo.close_session()

    # ---- mirror: same model through PyCoBi/PyRates with the analytical Jac ----
    yaml_dir = tmp_path / 'yaml'
    yaml_dir.mkdir()
    # PyRates' slash-notation path: <dir>/<filename-no-ext>/<template-name>.
    # Anchor on this test file's directory rather than cwd — preceding tests
    # leave cwd inside their own tmp_path via ODESystem._orig_dir handling.
    test_dir = os.path.dirname(os.path.abspath(__file__))
    yaml_template = os.path.join(test_dir, 'resources', 'ops', 'ops')
    ode_yaml = ODESystem.from_yaml(
        yaml_template,
        working_dir=str(yaml_dir), auto_dir=auto_dir,
        file_name='ops_mod', func_name='ops_fn',
        init_cont=False, analytical_jacobian=True,
        auto_constants=('eq',),
        # mirror the demo's c.ops auto-07p constants:
        NTST=15, NCOL=4, ISP=1, ILP=1, MXBF=10,
        NMX=25, NPR=500, IID=2, ITMX=8, ITNW=5, NWTN=3,
        EPSL=1e-7, EPSU=1e-7, EPSS=1e-4,
        DS=0.01, DSMIN=1e-3, DSMAX=0.5, IADS=1,
        # Stop when p3 reaches 0.95; the demo uses UZSTOP={3: 0.95} (literal
        # PAR(3)) because it doesn't declare parnames, but PyRates writes
        # `parnames` into c.eq and may not preserve the YAML's parameter order
        # (e.g. p3 ends up at PAR(4) for this model). Use the named form so
        # auto-07p resolves p3 via parnames regardless of its PAR index.
        UZSTOP={'p3': 0.95},
    )
    sols_yaml, _ = ode_yaml.run(c='eq', ICP='p3', name='ops_yaml')

    # ---- the analytical Jacobian path actually fired ----
    # (read these BEFORE close_session(clear_files=True) wipes the build dir)
    f90 = (yaml_dir / 'ops_mod.f90').read_text()
    ceq = (yaml_dir / 'c.eq').read_text()
    assert _DFDU_ASSIGN.search(f90)
    assert _DFDP_ASSIGN.search(f90)
    assert 'JAC = 1' in ceq

    ode_yaml.close_session(clear_files=True)

    # ---- numerical parity ----
    # demo uses default auto-07p names (U(1), PAR(3)); PyRates emits the user
    # names via unames/parnames (x1, p3).
    p3_demo = np.asarray(sols_demo['PAR(3)'].values).squeeze()
    x1_demo = np.asarray(sols_demo['U(1)'].values).squeeze()
    p3_yaml = np.asarray(sols_yaml['p3'].values).squeeze()
    x1_yaml = np.asarray(sols_yaml['x1'].values).squeeze()

    order = np.argsort(p3_demo)
    p3_demo, x1_demo = p3_demo[order], x1_demo[order]
    order = np.argsort(p3_yaml)
    p3_yaml, x1_yaml = p3_yaml[order], x1_yaml[order]

    lo = max(p3_demo.min(), p3_yaml.min())
    hi = min(p3_demo.max(), p3_yaml.max())
    assert hi - lo > 0.02, (
        f"insufficient overlap between PyCoBi/auto curves: p3 in [{lo}, {hi}]"
    )
    grid = np.linspace(lo, hi, 20)
    x1_demo_i = np.interp(grid, p3_demo, x1_demo)
    x1_yaml_i = np.interp(grid, p3_yaml, x1_yaml)
    assert np.allclose(x1_demo_i, x1_yaml_i, atol=1e-4, rtol=1e-4), (
        f"PyCoBi/PyRates ops curve diverges from auto-07p demo:\n"
        f"  p3_grid    = {grid}\n  x1_demo    = {x1_demo_i}\n  x1_yaml    = {x1_yaml_i}"
    )
