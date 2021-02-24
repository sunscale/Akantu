/* -------------------------------------------------------------------------- */
#include "py_solver.hh"
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <model.hh>
#include <non_linear_solver.hh>
#include <sparse_matrix_aij.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void register_solvers(py::module & mod) {
  py::class_<SparseMatrix>(mod, "SparseMatrix")
      .def("getMatrixType", &SparseMatrix::getMatrixType)
      .def("size", &SparseMatrix::size)
      .def("zero", &SparseMatrix::zero)
      .def("saveProfile", &SparseMatrix::saveProfile)
      .def("saveMatrix", &SparseMatrix::saveMatrix)
      .def(
          "add", [](SparseMatrix & self, UInt i, UInt j) { self.add(i, j); },
          "Add entry in the profile")
      .def(
          "add",
          [](SparseMatrix & self, UInt i, UInt j, Real value) {
            self.add(i, j, value);
          },
          "Add the value to the matrix")
      .def("__call__", [](const SparseMatrix & self, UInt i, UInt j) {
        return self(i, j);
      });

  py::class_<SparseMatrixAIJ, SparseMatrix>(mod, "SparseMatrixAIJ")
      .def("getIRN", &SparseMatrixAIJ::getIRN)
      .def("getJCN", &SparseMatrixAIJ::getJCN)
      .def("getA", &SparseMatrixAIJ::getA);

  py::class_<SolverVector>(mod, "SolverVector");
}

} // namespace akantu
