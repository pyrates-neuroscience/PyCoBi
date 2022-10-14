
# -*- coding: utf-8 -*-
#
#
# PyCoBi software framework for parameter continuation and automated bifurcation detection. See also:
# https://github.com/pyrates-neuroscience/PyCoBi
#
# Copyright (C) 2021-2022, Richard Gast.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>
#
# CITATION:
#
# Richard Gast et al., in preparation.

"""Python package for parameter continuations and bifurcation analysis via Auto-07p in Python.
"""

__author__ = "Richard Gast"
__status__ = "Development"
__version__ = "0.7.0"

from .pycobi import ODESystem, get_from_solutions, continue_period_doubling_bf
