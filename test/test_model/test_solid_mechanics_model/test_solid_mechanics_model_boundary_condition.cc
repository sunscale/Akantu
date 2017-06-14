/**
 * @file   test_solid_mechanics_model_boundary_condition.cc
 *
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Thu Feb 11 2016
 *
 * @brief  Test of the boundary condition functors and PBC
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

#include <iostream>
#include "aka_common.hh"
#include "solid_mechanics_model.hh"

using namespace akantu;

/* -------------------------------------------------------------------------- */

int main(int argc, char* argv[])
{
    UInt spatial_dimension(3);

    initialize("material.dat", argc, argv);

    Mesh mesh(spatial_dimension, "mesh_names");
    mesh.read("cube1.msh");

    const Array<UInt> & nodes = mesh.getElementGroup("Bottom").getNodes();
    Array<UInt>::const_scalar_iterator
      n_it  = nodes.begin(),
      n_end = nodes.end();

    SolidMechanicsModel model(mesh);
    model.setPBC(1, 0, 0);
    model.initFull(SolidMechanicsModelOptions(_static));

    Array<Real> & force = model.getForce();

    /// Testing FromTraction functor
    Real traction_ptr[] = {0, 0, 1};
    Vector<Real> surface_traction(traction_ptr, spatial_dimension);
    model.applyBC(BC::Neumann::FromTraction(surface_traction), "Bottom");

    Real total_force = 0;
    for (; n_it != n_end ; ++n_it) {
      if (!model.isPBCSlaveNode(*n_it) && !(force(*n_it, 2) > 0)) {
        std::cout << "FromTraction" << std::endl;
        return EXIT_FAILURE;
      }
      total_force += force(*n_it, 2);
    }

    Math::setTolerance(1e-14);
    if (!Math::are_float_equal(total_force, 1)) {
      std::cout << "Force balance " << total_force << " != 1" << std::endl;
      return EXIT_FAILURE;
    }

    // Copy force vector
    Array<Real> force_traction = force;
    force.set(0.);

/* -------------------------------------------------------------------------- */

    /// Testing FromHigherDim functor
    Real stress_ptr[] = {0, 0, 0,
                         0, 0, 0,
                         0, 0, 1};
    Matrix<Real> surface_stress(stress_ptr, spatial_dimension, spatial_dimension);
    model.applyBC(BC::Neumann::FromHigherDim(surface_stress), "Bottom");

    n_it = nodes.begin();
    for (; n_it != n_end ; ++n_it) {
      if (!Math::are_float_equal(force(*n_it, 2), force_traction(*n_it, 2))) {
        std::cout << "FromHigherDim" << std::endl;
        return EXIT_FAILURE;
      }
    }

/* -------------------------------------------------------------------------- */

    // Testing the periodic boundary conditions
    const Array<UInt> & xmin_nodes = mesh.getElementGroup("XMin").getNodes();
    const Array<UInt> & xmax_nodes = mesh.getElementGroup("XMax").getNodes();
    const Array<bool> & boundary = model.getBlockedDOFs();

    // Checking boundary on master and slave nodes
    n_it = xmin_nodes.begin();
    for (; n_it != xmin_nodes.end() ; ++n_it) {
      if (!model.isPBCSlaveNode(*n_it)) {
        for (UInt i = 0 ; i < spatial_dimension ; i++) {
          if (boundary(*n_it, i)) {
            std::cout << "PBC XMin : boundary on master node" << std::endl;
            return EXIT_FAILURE;
          }
        }
      } else {
        for (UInt i = 0 ; i < spatial_dimension ; i++) {
          if (!boundary(*n_it, i)) {
            std::cout << "PBC XMin : no boundary on slave node" << std::endl;
            return EXIT_FAILURE;
          }
        }
      }
    }

    n_it = xmax_nodes.begin();
    for (; n_it != xmax_nodes.end() ; ++n_it) {
      if (!model.isPBCSlaveNode(*n_it)) {
        for (UInt i = 0 ; i < spatial_dimension ; i++) {
          if (boundary(*n_it, i)) {
            std::cout << "PBC XMax : boundary on master node" << std::endl;
            return EXIT_FAILURE;
          }
        }
      } else {
        for (UInt i = 0 ; i < spatial_dimension ; i++) {
          if (!boundary(*n_it, i)) {
            std::cout << "PBC XMax : no boundary on slave node" << std::endl;
            return EXIT_FAILURE;
          }
        }
      }
    }

/* -------------------------------------------------------------------------- */

    /// Testing dirichlet BC functor
    model.applyBC(BC::Dirichlet::FixedValue(13.0, _x), "Bottom");
    model.applyBC(BC::Dirichlet::FixedValue(13.0, _y), "Bottom");
    model.applyBC(BC::Dirichlet::FixedValue(13.0, _z), "Bottom");

    Array<Real> & displacement = model.getDisplacement();

    n_it = nodes.begin();
    for (; n_it != n_end ; ++n_it) {
      for (UInt i = 0 ; i < spatial_dimension ; i++) {
        if (!boundary(*n_it, i) ||
            std::abs(displacement(*n_it, i) - 13.0) > Math::getTolerance()) {
          std::cout << "FixedValue" << std::endl;
          return EXIT_FAILURE;
        }
      }
    }

/* -------------------------------------------------------------------------- */
    model.assembleStiffnessMatrix();

    finalize();
    return EXIT_SUCCESS;
}
