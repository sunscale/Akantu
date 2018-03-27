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
import akantu as aka

def main():
    mesh = aka.Mesh(2)
    mesh.read('mesh_dcb_2d.msh')

    # Tests the getNbElement() function
    if mesh.getNbElement(aka._quadrangle_8) \
       != mesh.getNbElementByDimension(2):
        raise Exception("Number of elements wrong, should be"\
                        " {}".format(mesh.getNbElementByDimension(2)))
        return -1

    # TODO test the other functions in Mesh
    return 0

if __name__ == "__main__":
    sys.exit(main())
