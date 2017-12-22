/* -------------------------------------------------------------------------- */
#include "pybind11_akantu.hh"
/* -------------------------------------------------------------------------- */
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */

namespace py = pybind11;

namespace akantu {

PYBIND11_MODULE(aka_test, m) {
  m.doc() = "Module for the tests of akantu";

  py::class_<ArrayProxy<Real>>(m, "ArrayProxy", py::buffer_protocol())
      .def_buffer([](ArrayProxy<Real> & a) {
        return py::buffer_info(
            a.storage(),                           /* Pointer to buffer */
            sizeof(Real),                          /* Size of one scalar */
            py::format_descriptor<Real>::format(), /* Python struct-style
                                                       format descriptor */
            2,                                     /* Number of dimensions */
            {a.size(), a.getNbComponent()},        /* Buffer dimensions */
            {sizeof(Real) *
                 a.getNbComponent(), /* Strides (in bytes) for each index */
             sizeof(Real)});
      })
      .def(py::init([](py::array & n) {
        /* Request a buffer descriptor from Python */
        py::buffer_info info = n.request();

        /* Some sanity checks ... */
        if (info.format != py::format_descriptor<Real>::format())
          throw std::runtime_error(
              "Incompatible format: expected a double array!");

        if (info.ndim != 2)
          throw std::runtime_error("Incompatible buffer dimension!");

        return std::make_unique<ArrayProxy<Real>>(static_cast<Real *>(info.ptr),
                                                  info.shape[0], info.shape[1]);
      }));
}
} // namespace akantu
