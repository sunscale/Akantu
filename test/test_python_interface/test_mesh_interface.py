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

import sys
import os
import akantu as aka

def main():
    os.system('gmsh -order 2 -2 -o mesh_dcb_2d.msh mesh_dcb_2d.geo')
    print ""

    aka.initialize()

    mesh = aka.Mesh(2)
    mesh.read('mesh_dcb_2d.msh')

    print mesh.getNbElement(aka._segment_2)
    print mesh.getNbElement(aka._segment_3)
    print mesh.getNbElement(aka._triangle_3)
    print mesh.getNbElement(aka._triangle_6)
    print mesh.getNbElement(aka._quadrangle_4)
    print mesh.getNbElement(aka._quadrangle_8)
    print ""
    print mesh.getNbElement(1)
    print mesh.getNbElement(2)
    # Tests the getNbElement() function
    if mesh.getNbElement(aka._quadrangle_8) != mesh.getNbElement(2):
        print "Number of elements wrong, should be {}".format(mesh.getNbElement(2))
        return -1

    aka.finalize()
    return 0

if __name__ == "__main__":
    sys.exit(main())
