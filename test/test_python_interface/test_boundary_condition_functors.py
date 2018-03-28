#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
__author__ = "Lucas Frérot"
__copyright__ = "Copyright (C) 2016-2018, EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__credits__ = ["Lucas Frérot"]
__license__ = "L-GPLv3"
__maintainer__ = "Lucas Frérot"
__email__ = "lucas.frerot@epfl.ch"
# ------------------------------------------------------------------------------

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
    aka.parseInput("input_test.dat")

    mesh = aka.Mesh(2)
    mesh.read('mesh_dcb_2d.msh')

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

    return 0

if __name__ == "__main__":
    sys.exit(main())
