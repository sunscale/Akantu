/**
 * @file   python_user_defined_bc.cc
 *
 * @brief  python-user-defined boundary condition example
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <iostream>
#include <pybind11/embed.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

using namespace akantu;

class SineBoundary : public BC::Dirichlet::DirichletFunctor {
public:
  SineBoundary(Real amplitude, Real phase) {
    py_module = py::module::import("boundary_condition");
    py_sin_boundary = py_module.attr("SinBoundary")(amplitude, phase);
  }

public:
  inline void operator()(__attribute__((unused)) UInt node,
                         Vector<bool> & flags, Vector<Real> & primal,
                         const Vector<Real> & coord) const {
    py_sin_boundary.attr("compute")(primal, coord, flags);
  }

protected:
  py::object py_sin_boundary;
  py::module py_module;
};

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  py::scoped_interpreter guard{};

  UInt spatial_dimension = 2;

  Mesh mesh(spatial_dimension);
  mesh.read("fine_mesh.msh");

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull();

  /// boundary conditions
  Vector<Real> traction(2, 0.2);
  SineBoundary sin_boundary(.2, 10.);
  
  model.applyBC(sin_boundary, "Fixed_x");
  model.applyBC(BC::Dirichlet::FixedValue(0., _y), "Fixed_y");
  model.applyBC(BC::Neumann::FromTraction(traction), "Traction");

  // output a paraview file with the boundary conditions
  model.setBaseName("plate");
  model.addDumpFieldVector("displacement");
  model.addDumpFieldVector("external_force");
  model.addDumpField("blocked_dofs");
  model.dump();

  finalize();
  return EXIT_SUCCESS;
}
