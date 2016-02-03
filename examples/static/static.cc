/**
 * @file   static.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jan 18 2016
 *
 * @brief  This code refers to the implicit static example from the user manual
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
/* -------------------------------------------------------------------------- */

using namespace akantu;

#define bar_length 0.01
#define bar_height 0.01

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  const UInt spatial_dimension = 2;

  Mesh mesh(spatial_dimension);
  mesh.read("square.msh");
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull(SolidMechanicsModelOptions(_static));

  model.setBaseName("static");
  model.addDumpFieldVector("displacement");
  model.addDumpField("force");
  model.addDumpField("residual");
  model.addDumpField("grad_u");

  /// Dirichlet boundary conditions
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "Fixed_x");
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "Fixed_y");
  model.applyBC(BC::Dirichlet::FixedValue(0.0001, _y), "Traction");
  model.dump();

  model.assembleStiffnessMatrix();
  model.getStiffnessMatrix().saveMatrix("stiffness.mtx");

  bool converged =
      model.solveStep<_scm_newton_raphson_tangent_modified, _scc_increment>(
          1e-4, 2);

  if (!converged)
    AKANTU_DEBUG_ERROR("Did not converged in 1 step");

  model.dump();

  finalize();

  return EXIT_SUCCESS;
}
