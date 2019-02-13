/**
 * @file   pybind11_akantu.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Dec 22 2017
 * @date last modification: Thu Dec 28 2017
 *
 * @brief  tool to wrap akantu classes in python
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

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

  py::class_<MatrixProxy<Real>>(m, "Matrix", py::buffer_protocol())
      .def_buffer([](MatrixProxy<Real> & a) {
        return py::buffer_info(
            a.storage(),                           /* Pointer to buffer */
            sizeof(Real),                          /* Size of one scalar */
            py::format_descriptor<Real>::format(), /* Python struct-style
                                                       format descriptor */
            2,                                     /* Number of dimensions */
            {a.size(0), a.size(1)},                /* Buffer dimensions */
            {sizeof(Real), a.size(0) * sizeof(Real)}
            /* Strides (in bytes) for each index */
            );
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

        return std::make_unique<MatrixProxy<Real>>(
            static_cast<Real *>(info.ptr), info.shape[0], info.shape[1]);
      }));
} // Module aka test
} // namespace akantu
