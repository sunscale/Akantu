/**
 * @file   test_cohesive_intrinsic_quadrangle.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Mon Dec 18 2017
 *
 * @brief  Intrinsic cohesive elements' test for quadrangles
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

/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>
#include <limits>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

static void updateDisplacement(SolidMechanicsModelCohesive &, Array<UInt> &,
                               ElementType, Real);

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  const UInt spatial_dimension = 2;
  const UInt max_steps = 350;

  const ElementType type = _quadrangle_4;

  Mesh mesh(spatial_dimension);
  mesh.read("quadrangle.msh");

  // debug::setDebugLevel(dblDump);
  // std::cout << mesh << std::endl;
  // debug::setDebugLevel(dblWarning);

  SolidMechanicsModelCohesive model(mesh);
  model.getElementInserter().setLimit(_x, -0.01, 0.01);
  /// model initialization
  model.initFull();

  Real time_step = model.getStableTimeStep() * 0.8;
  model.setTimeStep(time_step);
  model.assembleMassLumped();

  Array<bool> & boundary = model.getBlockedDOFs();
  //  const Array<Real> & residual = model.getResidual();

  UInt nb_nodes = mesh.getNbNodes();
  UInt nb_element = mesh.getNbElement(type);

  /// boundary conditions
  for (UInt dim = 0; dim < spatial_dimension; ++dim) {
    for (UInt n = 0; n < nb_nodes; ++n) {
      boundary(n, dim) = true;
    }
  }

  model.assembleInternalForces();

  model.setBaseName("intrinsic_quadrangle");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("internal_force");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpField("external_force");

  model.setBaseNameToDumper("cohesive elements",
                            "cohesive_elements_quadrangle");
  model.addDumpFieldVectorToDumper("cohesive elements", "displacement");
  model.addDumpFieldToDumper("cohesive elements", "damage");

  model.dump();
  model.dump("cohesive elements");

  /// update displacement
  Array<UInt> elements;
  Vector<Real> bary(spatial_dimension);
  for (UInt el = 0; el < nb_element; ++el) {
    mesh.getBarycenter({type, el, _not_ghost}, bary);
    if (bary(_x) > 0.)
      elements.push_back(el);
  }

  Real increment = 0.01;

  updateDisplacement(model, elements, type, increment);

  // for (UInt n = 0; n < nb_nodes; ++n) {
  //   if (position(n, 1) + displacement(n, 1) > 0) {
  //     if (position(n, 0) == 0) {
  // 	displacement(n, 1) -= 0.25;
  //     }
  //     if (position(n, 0) == 1) {
  // 	displacement(n, 1) += 0.25;
  //     }
  //   }
  // }

  // std::ofstream edis("edis.txt");
  // std::ofstream erev("erev.txt");

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {
    model.solveStep();

    updateDisplacement(model, elements, type, increment);

    if (s % 1 == 0) {
      model.dump();
      model.dump("cohesive elements");
      std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }

    // // update displacement
    // for (UInt n = 0; n < nb_nodes; ++n) {
    //   if (position(n, 1) + displacement(n, 1) > 0) {
    // 	displacement(n, 0) -= 0.01;
    //   }
    // }

    //    Real Ed = dynamic_cast<MaterialCohesive&>
    //    (model.getMaterial(1)).getDissipatedEnergy();
    //    Real Er = dynamic_cast<MaterialCohesive&>
    //    (model.getMaterial(1)).getReversibleEnergy();

    // edis << s << " "
    // 	 << Ed << std::endl;

    // erev << s << " "
    // 	 << Er << std::endl;
  }

  // edis.close();
  // erev.close();

  Real Ed = model.getEnergy("dissipated");

  Real Edt = 1;

  std::cout << Ed << " " << Edt << std::endl;

  if (Ed < Edt * 0.999 || Ed > Edt * 1.001) {
    std::cout << "The dissipated energy is incorrect" << std::endl;
    return EXIT_FAILURE;
  }

  finalize();

  std::cout << "OK: test_cohesive_intrinsic_quadrangle was passed!"
            << std::endl;
  return EXIT_SUCCESS;
}

static void updateDisplacement(SolidMechanicsModelCohesive & model,
                               Array<UInt> & elements, ElementType type,
                               Real increment) {

  Mesh & mesh = model.getFEEngine().getMesh();
  UInt nb_element = elements.size();
  UInt nb_nodes = mesh.getNbNodes();
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);

  const Array<UInt> & connectivity = mesh.getConnectivity(type);
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> update(nb_nodes);
  update.clear();

  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt node = connectivity(elements(el), n);
      if (!update(node)) {
        displacement(node, 0) += increment;
        //	displacement(node, 1) += increment;
        update(node) = true;
      }
    }
  }
}
