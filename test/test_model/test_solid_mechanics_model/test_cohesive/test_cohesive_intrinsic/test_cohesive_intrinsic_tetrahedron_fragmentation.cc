/**
 * @file   test_cohesive_intrinsic_tetrahedron_fragmentation.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Oct 09 2013
 * @date last modification: Mon Dec 18 2017
 *
 * @brief  Test for cohesive elements
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <fstream>
#include <iostream>
#include <limits>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  //  debug::setDebugLevel(dblDump);
  ElementType type = _tetrahedron_10;

  const UInt spatial_dimension = 3;
  const UInt max_steps = 100;

  Mesh mesh(spatial_dimension);
  mesh.read("tetrahedron_full.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull();

  Real time_step = model.getStableTimeStep() * 0.8;
  model.setTimeStep(time_step);
  //  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();

  model.assembleInternalForces();

  model.setBaseName("intrinsic_tetrahedron_fragmentation");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("internal_force");
  model.addDumpField("stress");
  model.addDumpField("grad_u");

  model.setBaseNameToDumper("cohesive elements",
                            "cohesive_elements_tetrahedron_fragmentation");
  model.addDumpFieldVectorToDumper("cohesive elements", "displacement");
  model.addDumpFieldToDumper("cohesive elements", "damage");

  model.dump();
  model.dump("cohesive elements");

  /// update displacement
  UInt nb_element = mesh.getNbElement(type);
  UInt nb_nodes = mesh.getNbNodes();
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  Vector<Real> bary(spatial_dimension);

  const Array<UInt> & connectivity = mesh.getConnectivity(type);
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> update(nb_nodes);

  for (UInt s = 0; s < max_steps; ++s) {
    Real increment = s / 10.;
    update.clear();

    for (UInt el = 0; el < nb_element; ++el) {
      mesh.getBarycenter({type, el, _not_ghost}, bary);
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
        UInt node = connectivity(el, n);
        if (!update(node)) {
          for (UInt dim = 0; dim < spatial_dimension; ++dim) {
            displacement(node, dim) = increment * bary(dim);
            update(node) = true;
          }
        }
      }
    }

    if (s % 10 == 0) {
      model.dump();
      model.dump("cohesive elements");
    }
  }

  if (nb_nodes != nb_element * Mesh::getNbNodesPerElement(type)) {
    std::cout << "Wrong number of nodes" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  finalize();
  std::cout << "OK: test_cohesive_intrinsic_tetrahedron was passed!"
            << std::endl;
  return EXIT_SUCCESS;
}
