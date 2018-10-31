# import pybind11
import py11_akantu as pyaka
import py11_akantu_test_common as pyaka_test


def test_array():
    array = pyaka_test.getData()
    assert array.shape == (1000, 3)
