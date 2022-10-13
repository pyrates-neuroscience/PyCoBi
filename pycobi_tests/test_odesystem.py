"""Test suite for the ODESystem class.
"""

# imports
from pycobi import ODESystem
import pytest
import numpy as np
import sys

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


# test accuracy
accuracy = 1e-4

# auto installation directory
auto_dir = sys.argv[-1] if len(sys.argv) > 2 else "~/PycharmProjects/auto-07p"

# tests
#######


def test_1_1_init():
    """Testing the different instantiation options of the `ODESystem` class.
    """

    # parameters
    model = "model_templates.neural_mass_models.qif.qif"

    ode = ODESystem.from_yaml(model, init_cont=False, file_name='qif_equations', func_name="qif", auto_dir=auto_dir)
    assert type(ode.results) is dict
