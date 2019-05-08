# import pybind11
import pytest
import os
import subprocess
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


@pytest.fixture
def dcb_mesh():

    geo_file = 'mesh_dcb_2d.geo'
    mesh_file = 'mesh_dcb_2d.msh'

    if os.path.exists(mesh_file):
        return mesh_file

    open(geo_file, 'w').write("""
//Mesh size
dx = 57.8e-5;

Point(1) = {0,0,0,dx};
Point(2) = {0,0.00055,0,dx};
Point(3) = {0,-0.00055,0,dx};
Point(4) = {57.8e-3,0,0,dx};
Point(5) = {57.8e-3,0.00055,0,dx};
Point(6) = {57.8e-3,-0.00055,0,dx};
Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 4};
Line(4) = {1, 4};
Line(5) = {1, 3};
Line(6) = {6, 4};
Line(7) = {3, 6};
Line Loop(8) = {2, 3, -4, 1};
Plane Surface(9) = {-8};
Line Loop(10) = {5, 7, 6, -4};
Plane Surface(11) = {10};
Physical Surface("bulk") = {9,11};
Physical Line("coh") = {4};
Physical Line("edge") = {1};
Transfinite Surface "*";
Recombine Surface "*";
Mesh.SecondOrderIncomplete = 1;
    """)

    ret = subprocess.call(
        'gmsh -format msh2 -2 {0} -o {1}'.format(
            geo_file, mesh_file), shell=True)
    if not ret == 0:
        raise Exception(
            'execution of GMSH failed: do you have it installed ?')
    return mesh_file


@pytest.fixture
def elastic_material():
    mat_file = 'elastic.dat'
    open(mat_file, 'w').write("""
    material elastic [
         name    = bulk
         rho     = 2500
         nu      = 0.22
         E       = 71e9
]
""")
    return mat_file


def test_multiple_init(dcb_mesh, elastic_material):

    aka.parseInput(elastic_material)

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


def test_boundary_condition_functors(dcb_mesh, elastic_material):

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

    aka.parseInput(elastic_material)

    mesh = aka.Mesh(2)
    mesh.read(dcb_mesh)

    model = aka.SolidMechanicsModel(mesh, 2)
    model.initFull()

    model.applyBC(FixedValue(0.0, aka._x), "edge")

    stress = np.array([[1, 0],
                       [0, 0]])

    blocked_nodes = mesh.getElementGroup("edge").getNodes().flatten()
    boundary = model.getBlockedDOFs()

    # Testing that nodes are correctly blocked
    for n in blocked_nodes:
        if not boundary[n, 0]:
            return -1

    boundary.fill(False)

    model.applyBC(FromStress(stress), "edge")
    force = model.getForce()

    # Checking that nodes have a force in the correct direction
    for n in blocked_nodes:
        if not force[n, 0] > 0:
            return -1

    return 0
