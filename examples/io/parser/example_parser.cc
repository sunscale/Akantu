/**
 * @file   example_parser.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Mon Dec 14 2015
 * @date last modification: Mon Jan 18 2016
 *
 * @brief  Example on how to parse input text file
 *
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
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char * argv[]) {

  // Precise in initialize the name of the text input file to parse.
  initialize("input_file.dat", argc, argv);

  // Get the user ParserSection.
  const ParserSection & usersect = getUserParser();

  // getParameterValue() allows to extract data associated to a given parameter
  // name
  // and cast it in the desired type set as template paramter.
  Mesh mesh(usersect.getParameterValue<UInt>("spatial_dimension"));
  mesh.read(usersect.getParameterValue<std::string>("mesh_file"));

  // getParameter() can be used with variable declaration (destination type is
  // explicitly known).
  Int max_iter = usersect.getParameter("max_nb_iterations");
  Real precision = usersect.getParameter("precision");

  // Following NumPy convention, data can be interpreted as Vector or Matrix
  // structures.
  Matrix<Real> eigen_stress = usersect.getParameter("stress");

  SolidMechanicsModel model(mesh);
  model.initFull(SolidMechanicsModelOptions(_static));

  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x),
                usersect.getParameterValue<std::string>("outter_crust"));
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _y),
                usersect.getParameterValue<std::string>("outter_crust"));
  model.applyBC(BC::Neumann::FromStress(eigen_stress),
                usersect.getParameterValue<std::string>("inner_holes"));

  model.setDirectory("./paraview");
  model.setBaseName("swiss_cheese");
  model.addDumpFieldVector("displacement");

  auto & solver = model.getNonLinearSolver();
  solver.set("max_iterations", max_iter);
  solver.set("threshold", precision);

  model.solveStep();

  model.dump();

  finalize();

  return EXIT_SUCCESS;
}
