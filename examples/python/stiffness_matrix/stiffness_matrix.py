#!/usr/bin/env python3

import akantu
import numpy as np



def getStiffnessMatrix(material_file, mesh_file, traction):
    akantu.parseInput(material_file)
    spatial_dimension = 2

    ################################################################
    # Initialization
    ################################################################
    mesh = akantu.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = akantu.SolidMechanicsModel(mesh)
    model.initFull(akantu._static)
    model.assembleStiffnessMatrix()
    K = model.getDOFManager().getMatrix('K')
    stiff = akantu.AkantuSparseMatrix(K).toarray()
    return stiff

################################################################
# main
################################################################


def main():

    import os
    mesh_file = 'plate.msh'
    # if mesh was not created the calls gmsh to generate it
    if not os.path.isfile(mesh_file):
        import subprocess
        ret = subprocess.call(
            'gmsh -format msh2 -2 plate.geo {0}'.format(mesh_file), shell=True)
        if not ret == 0:
            raise Exception(
                'execution of GMSH failed: do you have it installed ?')

    material_file = 'material.dat'

    traction = 1.
    mat = getStiffnessMatrix(material_file, mesh_file, traction)
    print(mat)

################################################################
if __name__ == "__main__":
    main()
