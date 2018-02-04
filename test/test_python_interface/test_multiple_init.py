#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ===============================================================================
# @file   test_multiple_init.py
#
# @author Fabian Barras <fabian.barras@epfl.ch>
#
# @date creation: Tue Jan 5 2016
#
# @brief  Testing multiple initialize calls through Python
#
# @section LICENSE
#
# Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# ===============================================================================
from __future__ import print_function
import sys
import os
import akantu as aka

os.system('gmsh -order 2 -2 -o mesh_dcb_2d.msh mesh_dcb_2d.geo')

aka.initialize('input_test.dat')

print('First initialisation')
mesh = aka.Mesh(2)
mesh.read('mesh_dcb_2d.msh')
model = aka.SolidMechanicsModel(mesh)
model.initFull(aka.SolidMechanicsModelOptions(aka._static))
del model
del mesh

print('Second initialisation')
mesh = aka.Mesh(2)
mesh.read('mesh_dcb_2d.msh')
model = aka.SolidMechanicsModel(mesh)
model.initFull(aka.SolidMechanicsModelOptions(aka._static))
del model
del mesh

print('All right')
