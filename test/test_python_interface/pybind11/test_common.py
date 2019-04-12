# import pybind11
import numpy as np
import py11_akantu_test_common as pyaka_test


def test_array_size():
    ptr, array = pyaka_test.createArray(1000, 3)
    assert array.shape == (1000, 3)


def test_array_nocopy():
    ptr, array = pyaka_test.createArray(1000, 3)
    through_python = pyaka_test.getRawPointerArray(array)
    assert(ptr == through_python)


def test_modify_array():
    ptr, array = pyaka_test.createArray(3, 3)
    array[0, :] = (1., 2., 3.)
    array2 = pyaka_test.getArray(ptr)
    assert(np.linalg.norm(array-array2) < 1e-15)

    for i in [1, 2, 3]:
        ptr, array = pyaka_test.createArray(10000, i)
        array[:, :] = np.random.random((10000, i))
        array2 = pyaka_test.getArray(ptr)
        assert(np.linalg.norm(array-array2) < 1e-15)


def test_array_copy():
    ptr, array = pyaka_test.createArray(1000, 3)
    array2 = pyaka_test.copyArray(ptr)
    ptr2 = pyaka_test.getRawPointerArray(array2)
    assert(ptr != ptr2)


def test_vector_size():
    ptr, vector = pyaka_test.createVector(3)
    assert vector.shape == (3,)


def test_vector_nocopy():
    ptr, vector = pyaka_test.createVector(3)
    through_python = pyaka_test.getRawPointerVector(vector)
    assert(ptr == through_python)


def test_modify_vector():
    ptr, vector = pyaka_test.createVector(3)
    vector[:] = (1., 2., 3.)
    vector2 = pyaka_test.getVector(ptr)
    assert(np.linalg.norm(vector-vector2) < 1e-15)

    for i in np.arange(1, 10):
        ptr, vector = pyaka_test.createVector(i)
        vector[:] = np.random.random(i)
        vector2 = pyaka_test.getVector(ptr)
        assert(np.linalg.norm(vector-vector2) < 1e-15)


def test_vector_copy():
    ptr, vector = pyaka_test.createVector(1000)
    vector2 = pyaka_test.copyVector(ptr)
    ptr2 = pyaka_test.getRawPointerVector(vector2)
    assert(ptr != ptr2)


def test_matrix_size():
    ptr, matrix = pyaka_test.createMatrix(3, 2)
    assert matrix.shape == (3, 2)


def test_matrix_nocopy():
    ptr, matrix = pyaka_test.createMatrix(3, 2)
    through_python = pyaka_test.getRawPointerMatrix(matrix)
    assert(ptr == through_python)


def test_modify_matrix():
    ptr, matrix = pyaka_test.createMatrix(2, 3)
    matrix[0, :] = (1., 2., 3.)
    matrix2 = pyaka_test.getMatrix(ptr)
    assert(np.linalg.norm(matrix-matrix2) < 1e-15)

    for i in np.arange(1, 10):
        for j in np.arange(1, 10):
            ptr, matrix = pyaka_test.createMatrix(i, j)
            matrix[:, :] = np.random.random((i, j))
            matrix2 = pyaka_test.getMatrix(ptr)
            assert(np.linalg.norm(matrix-matrix2) < 1e-15)


def test_matrix_copy():
    ptr, matrix = pyaka_test.createMatrix(10, 3)
    matrix2 = pyaka_test.copyMatrix(ptr)
    ptr2 = pyaka_test.getRawPointerMatrix(matrix2)
    assert(ptr != ptr2)
