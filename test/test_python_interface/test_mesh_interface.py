#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ===============================================================================
# @file   test_mesh_interface.py
#
# @author Lucas Frérot <lucas.frerot@epfl.ch>
#
# @date creation: Mon Feb 01 2016
#
# @brief  Testing the Mesh python interface
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
import akantu as aka

def main():
    aka.initialize()

    mesh = aka.Mesh(2)
    mesh.read('mesh_dcb_2d.msh')

    # Tests the getNbElement() function
    if mesh.getNbElement(aka._quadrangle_8) != mesh.getNbElementByDimension(2):
        print("Number of elements wrong, should be {}".format(mesh.getNbElementByDimension(2)))
        return -1

    # TODO test the other functions in Mesh
    aka.finalize()
    return 0

if __name__ == "__main__":
    sys.exit(main())
