# import pybind11
import numpy as np
import py11_akantu_test_common as pyaka_test


def test_array_size():
    ptr, array = pyaka_test.createData(1000, 3)
    assert array.shape == (1000, 3)


def test_array_nocopy():
    ptr, array = pyaka_test.createData(1000, 3)
    through_python = pyaka_test.getRawPointer(array)
    assert(ptr == through_python)


def test_modify_values():
    ptr, array = pyaka_test.createData(3, 3)
    array[0, :] = (1., 2., 3.)
    array2 = pyaka_test.getData(ptr)
    assert(np.linalg.norm(array-array2) < 1e-15)

    for i in [1, 2, 3]:
        ptr, array = pyaka_test.createData(10000, i)
        array[:, :] = np.random.random((10000, i))
        array2 = pyaka_test.getData(ptr)
        assert(np.linalg.norm(array-array2) < 1e-15)


def test_array_copy():
    ptr, array = pyaka_test.createData(1000, 3)
    array2 = pyaka_test.getCopyData(ptr)
    ptr2 = pyaka_test.getRawPointer(array2)
    assert(ptr != ptr2)
