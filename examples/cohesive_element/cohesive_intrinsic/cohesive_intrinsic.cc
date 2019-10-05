/**
 * @file   cohesive_intrinsic.cc
 *
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Mon Jan 18 2016
 *
 * @brief  Test for cohesive elements
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
#include "element_group.hh"
#include "mesh_iterators.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

static void updateDisplacement(SolidMechanicsModelCohesive &,
                               const ElementGroup &, Real);

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  const UInt spatial_dimension = 2;
  const UInt max_steps = 350;

  Mesh mesh(spatial_dimension);
  mesh.read("triangle.msh");

  SolidMechanicsModelCohesive model(mesh);
  model.getElementInserter().setLimit(_x, -0.26, -0.24);

  /// model initialization
  model.initFull(_analysis_method = _explicit_lumped_mass,
                 _is_extrinsic = false);

  Real time_step = model.getStableTimeStep() * 0.8;
  model.setTimeStep(time_step);
  std::cout << "Time step: " << time_step << std::endl;

  Array<bool> & boundary = model.getBlockedDOFs();

  UInt nb_nodes = mesh.getNbNodes();

  /// boundary conditions
  for (UInt dim = 0; dim < spatial_dimension; ++dim) {
    for (UInt n = 0; n < nb_nodes; ++n) {
      boundary(n, dim) = true;
    }
  }

  model.setBaseName("intrinsic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpField("external_force");
  model.addDumpField("internal_force");
  model.dump();

  /// update displacement
  auto && elements = mesh.createElementGroup("diplacement");
  Vector<Real> barycenter(spatial_dimension);

  for_each_element(mesh,
                   [&](auto && el) {
                     mesh.getBarycenter(el, barycenter);
                     if (barycenter(_x) > -0.25)
                       elements.add(el, true);
                   },
                   _element_kind = _ek_regular);

  Real increment = 0.01;

  updateDisplacement(model, elements, increment);

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {
    model.solveStep();

    updateDisplacement(model, elements, increment);
    if (s % 1 == 0) {
      model.dump();
      std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }
  }

  Real Ed = model.getEnergy("dissipated");
  Real Edt = 2 * sqrt(2);

  std::cout << Ed << " " << Edt << std::endl;

  if (Ed < Edt * 0.999 || Ed > Edt * 1.001 || std::isnan(Ed)) {
    std::cout << "The dissipated energy is incorrect" << std::endl;
    return EXIT_FAILURE;
  }

  finalize();

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
static void updateDisplacement(SolidMechanicsModelCohesive & model,
                               const ElementGroup & group, Real increment) {
  Array<Real> & displacement = model.getDisplacement();

  for (auto && node : group.getNodeGroup().getNodes()) {
    displacement(node, 0) += increment;
  }
}
