#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
__author__ = "Fabian Barras"
__copyright__ = "Copyright (C) 2016-2018, EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__credits__ = ["Fabian Barras"]
__license__ = "L-GPLv3"
__maintainer__ = "Fabian Barras"
__email__ = "fabian.barras@epfl.ch"
# ------------------------------------------------------------------------------

from __future__ import print_function
import sys
import os
import akantu as aka

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
