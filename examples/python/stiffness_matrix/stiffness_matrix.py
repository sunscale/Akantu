#!/usr/bin/env python3
import akantu as aka


def getStiffnessMatrix(material_file, mesh_file, traction):
    aka.parseInput(material_file)
    spatial_dimension = 2

    # --------------------------------------------------------------------------
    # Initialization
    # --------------------------------------------------------------------------
    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = aka.SolidMechanicsModel(mesh)
    model.initFull(_analysis_method=aka._static)

    model.assembleStiffnessMatrix()
    K = model.getDOFManager().getMatrix('K')
    stiff = aka.AkantuSparseMatrix(K).toarray()

    return stiff


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    mesh_file = 'plate.msh'
    material_file = 'material.dat'

    traction = 1.
    mat = getStiffnessMatrix(material_file, mesh_file, traction)
    print(mat)


# --------------------------------------------------------------------------
if __name__ == "__main__":
    main()
