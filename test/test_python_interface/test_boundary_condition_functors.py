#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ===============================================================================
# @file   test_boundary_condition_functors.py
#
# @author Lucas Frérot <lucas.frerot@epfl.ch>
#
# @date creation: Thu Feb 11 2016
#
# @brief  Testing the apply*BC functions in python interface
#
# @section LICENSE
#
# Copyright (©) 2016 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

import numpy as np
import akantu as aka

######################################################################
# Boundary conditions founctors
######################################################################
class FixedValue:
  """Functor for Dirichlet boundary conditions"""
  def __init__(self, value, axis):
    self.value = value
    axis_dict = {'x':0, 'y':1, 'z':2}
    self.axis = axis_dict[axis]

  def operator(self, node, flags, primal, coord):
    primal[self.axis] = self.value
    flags[self.axis]  = True

#---------------------------------------------------------------------

class FromStress:
  """Functor for Neumann boundary conditions"""
  def __init__(self, stress):
    self.stress = stress

  def operator(self, quad_point, dual, coord, normals):
    dual[:] = np.dot(self.stress,normals)

######################################################################

def main():
    aka.initialize("input_test.dat")

    mesh = aka.Mesh(2)
    mesh.read('mesh_dcb_2d.msh')
    mesh.createGroupsFromStringMeshData("physical_names")

    model = aka.SolidMechanicsModel(mesh, 2)
    model.initFull()

    model.applyDirichletBC(FixedValue(0.0, 'x'), "edge")

    stress = np.array([[1, 0],
                       [0, 0]])


    blocked_nodes = mesh.getElementGroup("edge").getNodes().flatten()
    boundary = model.getBlockedDOFs()

    # Testing that nodes are correctly blocked
    for n in blocked_nodes:
        if not boundary[n, 0]:
            return -1

    boundary.fill(False)

    model.applyNeumannBC(FromStress(stress), "edge")
    force = model.getForce()

    # Checking that nodes have a force in the correct direction
    for n in blocked_nodes:
        if not force[n, 0] > 0:
            return -1

    aka.finalize()
    return 0

if __name__ == "__main__":
    sys.exit(main())
