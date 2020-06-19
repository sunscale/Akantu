/**
 * @file   test_cohesive_buildfragments.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Tue Aug 22 2017
 *
 * @brief  Test for cohesive elements
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "fragment_manager.hh"
#include "material_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>
#include <limits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  Math::setTolerance(1e-14);

  const UInt spatial_dimension = 2;
  const UInt max_steps = 200;
  Real strain_rate = 1.e5;
  ElementType type = _quadrangle_4;

  Real L = 0.03;
  Real theoretical_mass = L * L / 20. * 2500;

  ElementType type_facet = Mesh::getFacetType(type);
  ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);

  Mesh mesh(spatial_dimension);
  mesh.read("mesh.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull(_analysis_method = _explicit_lumped_mass,
                 _is_extrinsic = true);

  Real time_step = model.getStableTimeStep() * 0.05;
  model.setTimeStep(time_step);
  //  std::cout << "Time step: " << time_step << std::endl;

  Real disp_increment = strain_rate * L / 2. * time_step;
  model.assembleMassLumped();

  Array<Real> & velocity = model.getVelocity();
  const Array<Real> & position = mesh.getNodes();
  UInt nb_nodes = mesh.getNbNodes();

  /// initial conditions
  for (UInt n = 0; n < nb_nodes; ++n)
    velocity(n, 0) = strain_rate * position(n, 0);

  /// boundary conditions
  model.applyBC(BC::Dirichlet::FixedValue(0, _x), "Left_side");
  model.applyBC(BC::Dirichlet::FixedValue(0, _x), "Right_side");

  UInt cohesive_index = 1;

  UInt nb_quad_per_facet =
      model.getFEEngine("FacetsFEEngine").getNbIntegrationPoints(type_facet);
  MaterialCohesive & mat_cohesive =
      dynamic_cast<MaterialCohesive &>(model.getMaterial(cohesive_index));
  const Array<Real> & damage = mat_cohesive.getDamage(type_cohesive);

  FragmentManager fragment_manager(model, false);
  const Array<Real> & fragment_mass = fragment_manager.getMass();

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {
    model.checkCohesiveStress();
    model.solveStep();

    /// apply boundary conditions
    model.applyBC(BC::Dirichlet::IncrementValue(-disp_increment, _x),
                  "Left_side");
    model.applyBC(BC::Dirichlet::IncrementValue(disp_increment, _x),
                  "Right_side");

    if (s % 1 == 0) {
      //      model.dump();
      std::cout << "passing step " << s << "/" << max_steps << std::endl;

      fragment_manager.computeAllData();

      /// check number of fragments
      UInt nb_fragment_num = fragment_manager.getNbFragment();

      UInt nb_cohesive_elements = mesh.getNbElement(type_cohesive);

      UInt nb_fragment = 1;
      for (UInt el = 0; el < nb_cohesive_elements; ++el) {
        UInt q = 0;
        while (q < nb_quad_per_facet &&
               Math::are_float_equal(damage(el * nb_quad_per_facet + q), 1))
          ++q;

        if (q == nb_quad_per_facet) {
          ++nb_fragment;
        }
      }

      if (nb_fragment != nb_fragment_num) {
        std::cout << "The number of fragments is wrong!" << std::endl;
        return EXIT_FAILURE;
      }

      /// check mass computation
      Real total_mass = 0.;
      for (UInt frag = 0; frag < nb_fragment_num; ++frag) {
        total_mass += fragment_mass(frag);
      }

      if (!Math::are_float_equal(theoretical_mass, total_mass)) {
        std::cout << "The fragments' mass is wrong!" << std::endl;
        return EXIT_FAILURE;
      }
    }
  }

  model.dump();

  /// check velocities
  UInt nb_fragment = fragment_manager.getNbFragment();
  const Array<Real> & fragment_velocity = fragment_manager.getVelocity();
  const Array<Real> & fragment_center = fragment_manager.getCenterOfMass();

  Real fragment_length = L / nb_fragment;
  Real initial_position = -L / 2. + fragment_length / 2.;

  for (UInt frag = 0; frag < nb_fragment; ++frag) {
    Real theoretical_center = initial_position + fragment_length * frag;

    if (!Math::are_float_equal(fragment_center(frag, 0), theoretical_center)) {
      std::cout << "The fragments' center is wrong!" << std::endl;
      return EXIT_FAILURE;
    }

    Real initial_vel = fragment_center(frag, 0) * strain_rate;

    Math::setTolerance(100);

    if (!Math::are_float_equal(fragment_velocity(frag), initial_vel)) {
      std::cout << "The fragments' velocity is wrong!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  finalize();

  std::cout << "OK: test_cohesive_buildfragments was passed!" << std::endl;
  return EXIT_SUCCESS;
}
