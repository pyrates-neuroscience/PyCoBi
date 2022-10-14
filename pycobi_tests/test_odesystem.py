"""Test suite for the ODESystem class.
"""

# imports
from pycobi import ODESystem
from pyrates import clear
from pytest import fixture, approx
import numpy as np

# meta infos
__author__ = "Richard Gast"
__status__ = "Development"

# Utility
#########


def setup_module():
    print("\n")
    print("============================")
    print("| Test Suite : Input layer |")
    print("============================")


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
    ode1 = ODESystem(working_dir="resources", auto_dir=auto_dir, init_cont=False)
    ode1.close_session()

    # initialize ODESystem from YAML file
    model = "model_templates.neural_mass_models.qif.qif"
    ode2, t1 = ODESystem.from_yaml(model, init_cont=False, file_name='qif_eq2', func_name="qif2", auto_dir=auto_dir)
    ode2.close_session()
    clear(t1)

    # initialize ODESystem from YAML file with different parameters
    ode3, t2 = ODESystem.from_yaml(model, init_cont=False, file_name='qif_eq3', func_name="qif3", auto_dir=auto_dir,
                                   node_vars={'p/qif_op/eta': 2.0})
    ode3.close_session()
    clear(t2)

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
    ode1 = ODESystem(working_dir="resources", auto_dir=auto_dir, init_cont=True, e="qif_eq", c="ivp", NPR=100)
    ode1.close_session()

    # initialize ODESystem from YAML file
    model = "model_templates.neural_mass_models.qif.qif"
    ode2, t1 = ODESystem.from_yaml(model, init_cont=True, file_name='qif_eq2', func_name="qif2", auto_dir=auto_dir,
                                   NPR=100, NMX=5000)
    ode2.close_session()
    clear(t1)

    # initialize ODESystem from YAML file with different parameters
    ode3, t2 = ODESystem.from_yaml(model, init_cont=True, file_name='qif_eq3', func_name="qif3", auto_dir=auto_dir,
                                   node_vars={'p/qif_op/eta': 2.0}, NPR=100, NMX=5000)
    ode3.close_session()
    clear(t2)

    # these tests should pass
    assert (ode1[0].loc[:, "U(1)"] - ode2[0].loc[:, "U(1)"]).sum()[0] == approx(0.0, rel=accuracy, abs=accuracy)
    assert abs((ode2[0].loc[:, "U(1)"] - ode3[0].loc[:, "U(1)"]).sum()[0]) > 0
