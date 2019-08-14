/**
 * @file   cohesive_extrinsic_implicit.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jan 12 2016
 * @date last modification: Mon Jan 18 2016
 *
 * @brief  Example for extrinsic cohesive elements in implicit
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
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  debug::setDebugLevel(dblError);

  const UInt spatial_dimension = 2;
  const UInt max_steps = 20;
  const Real final_opening = 1e-4;

  Mesh mesh(spatial_dimension);
  mesh.read("dcb_2d.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull(SolidMechanicsModelCohesiveOptions(_static, true));

  //  CohesiveElementInserter inserter(mesh);
  model.limitInsertion(_y, -0.000001, 0.000001);
  model.updateAutomaticInsertion();

  Real eps = 1e-11;
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & position = mesh.getNodes();
  Array<Real> & displacement = model.getDisplacement();

  /// boundary conditions
  const Vector<Real> & lower = mesh.getLowerBounds();
  const Vector<Real> & upper = mesh.getUpperBounds();
  const Real left = lower[0];
  const Real right = upper[0];

  for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
    if (std::abs(position(n, 0) - left) < eps) {
      boundary(n, 1) = true;
      boundary(n, 0) = true;
    }
    if (std::abs(position(n, 0) - right) < eps && position(n, 1) < 0.0)
      boundary(n, 1) = true;
    if (std::abs(position(n, 0) - right) < eps && position(n, 1) > 0.0)
      boundary(n, 1) = true;
  }

  model.setBaseName("extr_impl");
  model.addDumpFieldVector("displacement");
  model.addDumpField("external_force");
  model.addDumpField("internal_force");
  model.addDumpField("stress");
  model.addDumpField("partitions");
  model.dump();

  // Dumping cohesive elements
  model.setBaseNameToDumper("cohesive elements", "cohe_elem_extr_impl");
  model.addDumpFieldVectorToDumper("cohesive elements", "displacement");
  model.addDumpFieldToDumper("cohesive elements", "damage");
  model.dump("cohesive elements");

  //  model.updateResidual();

  Real increment = final_opening / max_steps;
  Real tolerance = 1e-13;
  Real error;
  bool load_reduction = false;
  Real tol_increase_factor = 1.0e8;

  /// Main loop
  for (UInt nstep = 0; nstep < max_steps; ++nstep) {
    std::cout << "step no.  " << nstep << std::endl;

    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      if (std::abs(position(n, 0) - right) < eps && position(n, 1) > 0.0)
        displacement(n, 1) += increment;

      if (std::abs(position(n, 0) - right) < eps && position(n, 1) < 0.0)
        displacement(n, 1) -= increment;
    }

    model.solveStepCohesive<_scm_newton_raphson_tangent,
                            SolveConvergenceCriteria::_increment>(
        tolerance, error, 25, load_reduction, tol_increase_factor);

    // If convergence has not been reached, the load is reduced and
    // the incremental step is solved again.
    while (!load_reduction && error > tolerance) {
      load_reduction = true;

      std::cout << "LOAD STEP REDUCTION" << std::endl;
      increment = increment / 2.0;

      for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
        if (std::abs(position(n, 0) - right) < eps && position(n, 1) > 0.0)
          displacement(n, 1) -= increment;

        if (std::abs(position(n, 0) - right) < eps && position(n, 1) < 0.0)
          displacement(n, 1) += increment;
      }

      UInt nb_cohesive_elements =
          mesh.getNbElement(spatial_dimension, _not_ghost, _ek_cohesive);

      model.solveStepCohesive<_scm_newton_raphson_tangent,
                              SolveConvergenceCriteria::_increment>(
          tolerance, error, 25, load_reduction, tol_increase_factor);

      UInt new_nb_cohesive_elements =
          mesh.getNbElement(spatial_dimension, _not_ghost, _ek_cohesive);

      UInt nb_cohe[2];
      nb_cohe[0] = nb_cohesive_elements;
      nb_cohe[1] = new_nb_cohesive_elements;

      // Every time a new cohesive element is introduced, the variable
      // load_reduction is set to false, so that it is possible to
      // further iterate in the loop of load reduction. If no new
      // cohesive elements are introduced, usually there is no gain in
      // further reducing the load, even if convergence is not reached
      if (nb_cohe[0] == nb_cohe[1])
        load_reduction = true;
      else
        load_reduction = false;
    }

    model.dump();
    model.dump("cohesive elements");

    UInt nb_cohe_elems[1];
    nb_cohe_elems[0] =
        mesh.getNbElement(spatial_dimension, _not_ghost, _ek_cohesive);
    std::cout << "No. of cohesive elements: " << nb_cohe_elems[0] << std::endl;
  }

  finalize();

  return EXIT_SUCCESS;
}
