#!/bin/env python
import pytest
import numpy as np
import py11_akantu_test_common as aka


def test_array_size():
    ptr, array = aka.createArray(1000, 3)
    assert array.shape == (1000, 3)


def test_array_nocopy():
    ptr, array = aka.createArray(1000, 3)
    through_python = aka.getRawPointerArray(array)
    assert(ptr == through_python)


def test_modify_array():
    ptr, array = aka.createArray(3, 3)
    array[0, :] = (1., 2., 3.)
    array2 = aka.getArray(ptr)
    assert(np.linalg.norm(array-array2) < 1e-15)

    for i in [1, 2, 3]:
        ptr, array = aka.createArray(10000, i)
        array[:, :] = np.random.random((10000, i))
        array2 = aka.getArray(ptr)
        assert(np.linalg.norm(array-array2) < 1e-15)


def test_array_copy():
    ptr, array = aka.createArray(1000, 3)
    array2 = aka.copyArray(ptr)
    ptr2 = aka.getRawPointerArray(array2)
    assert(ptr != ptr2)


def test_vector_size():
    ptr, vector = aka.createVector(3)
    assert vector.shape == (3,)


def test_vector_nocopy():
    ptr, vector = aka.createVector(3)
    through_python = aka.getRawPointerVector(vector)
    assert(ptr == through_python)


def test_modify_vector():
    ptr, vector = aka.createVector(3)
    vector[:] = (1., 2., 3.)
    vector2 = aka.getVector(ptr)
    assert(np.linalg.norm(vector-vector2) < 1e-15)

    for i in np.arange(1, 10):
        ptr, vector = aka.createVector(i)
        vector[:] = np.random.random(i)
        vector2 = aka.getVector(ptr)
        assert(np.linalg.norm(vector-vector2) < 1e-15)


def test_vector_copy():
    ptr, vector = aka.createVector(1000)
    vector2 = aka.copyVector(ptr)
    ptr2 = aka.getRawPointerVector(vector2)
    assert(ptr != ptr2)


def test_matrix_size():
    ptr, matrix = aka.createMatrix(3, 2)
    assert matrix.shape == (3, 2)


def test_matrix_nocopy():
    ptr, matrix = aka.createMatrix(3, 2)
    through_python = aka.getRawPointerMatrix(matrix)
    assert(ptr == through_python)


def test_modify_matrix():
    ptr, matrix = aka.createMatrix(2, 3)
    matrix[0, :] = (1., 2., 3.)
    matrix2 = aka.getMatrix(ptr)
    assert(np.linalg.norm(matrix-matrix2) < 1e-15)

    for i in np.arange(1, 10):
        for j in np.arange(1, 10):
            ptr, matrix = aka.createMatrix(i, j)
            matrix[:, :] = np.random.random((i, j))
            matrix2 = aka.getMatrix(ptr)
            assert(np.linalg.norm(matrix-matrix2) < 1e-15)


def test_matrix_copy():
    ptr, matrix = aka.createMatrix(10, 3)
    matrix2 = aka.copyMatrix(ptr)
    ptr2 = aka.getRawPointerMatrix(matrix2)
    assert(ptr != ptr2)


def test_multiple_init():
    aka.parseInput("elastic.dat")
    dcb_mesh = 'mesh_dcb_2d.msh'

    print('First initialisation')
    mesh = aka.Mesh(2)
    mesh.read(dcb_mesh)
    model = aka.SolidMechanicsModel(mesh)
    model.initFull(aka.SolidMechanicsModelOptions(aka._static))
    del model
    del mesh

    print('Second initialisation')
    mesh = aka.Mesh(2)
    mesh.read(dcb_mesh)
    model = aka.SolidMechanicsModel(mesh)
    model.initFull(aka.SolidMechanicsModelOptions(aka._static))
    del model
    del mesh

    print('All right')


def test_boundary_condition_functors():

    class FixedValue(aka.DirichletFunctor):
        def __init__(self, value, axis):
            super().__init__(axis)
            self.value = value
            self.axis = int(axis)

        def __call__(self, node, flags, primal, coord):
            primal[self.axis] = self.value
            flags[self.axis] = True

    class FromStress(aka.NeumannFunctor):
        def __init__(self, stress):
            super().__init__()
            self.stress = stress

        def __call__(self, quad_point, dual, coord, normals):
            dual[:] = np.dot(self.stress, normals)

    aka.parseInput("elastic.dat")

    mesh = aka.Mesh(2)
    mesh.read("mesh_dcb_2d.msh")

    model = aka.SolidMechanicsModel(mesh, 2)
    model.initFull()

    model.applyBC(FixedValue(0.0, aka._x), "edge")

    stress = np.array([[1, 0],
                       [0, 0]])

    blocked_nodes = \
        mesh.getElementGroup("edge").getNodeGroup().getNodes().flatten()
    boundary = model.getBlockedDOFs()

    # Testing that nodes are correctly blocked
    for n in blocked_nodes:
        assert boundary[n, 0]

    boundary.fill(False)

    model.applyBC(FromStress(stress), "edge")
    force = model.getExternalForce()

    # Checking that nodes have a force in the correct direction
    for n in blocked_nodes:
        assert force[n, 0] > 0

    return 0


def test_mesh_interface():
    mesh = aka.Mesh(2)
    mesh.read("mesh_dcb_2d.msh")

    # Tests the getNbElement() function
    if mesh.getNbElement(aka._quadrangle_8) != mesh.getNbElement(2):
        raise Exception("Number of elements wrong: "
                        " {0} != {1}".format(
                            mesh.getNbElement(aka._quadrangle_8),
                            mesh.getNbElement(2)))


def test_heat_transfer():
    mesh = aka.Mesh(2)
    model = aka.HeatTransferModel(mesh)
    print(aka._explicit_lumped_mass)
    model.initFull(aka._explicit_lumped_mass)


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
