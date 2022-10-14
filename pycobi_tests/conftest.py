"""Config file for pytest for defining additional parser options.
"""


def pytest_addoption(parser):
    parser.addoption("--auto_dir", action="store", default="/home/circleci/project/auto-07p")
