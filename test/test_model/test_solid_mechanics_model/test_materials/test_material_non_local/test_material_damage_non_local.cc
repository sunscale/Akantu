/**
 * @file   test_material_damage_non_local.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Clement Roux <clement.roux@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  test for non-local damage materials on a 2D plate with a section gap
 * the sample should break at the notch
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  debug::setDebugLevel(dblWarning);

  akantu::initialize("material_damage_non_local.dat", argc, argv);
  UInt max_steps = 1100;

  const UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("mesh_section_gap.msh");
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  SolidMechanicsModel model(mesh);

  /// model initialization

  model.initFull();

  Real time_step = model.getStableTimeStep();
  model.setTimeStep(time_step / 2.5);

  std::cout << model << std::endl;

  model.applyBC(BC::Dirichlet::FixedValue(0.0), "Fixed");

  // Boundary condition (Neumann)
  Matrix<Real> stress(2, 2);
  stress.eye(5e8);
  model.applyBC(BC::Neumann::FromHigherDim(stress), "Traction");

  /*model.setBaseName("damage_non_local");
  model.addDumpFieldVector("displacement");
  model.addDumpField("mass"        );
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpFieldVector("external_force");
  model.addDumpFieldVector("internal_force");
  model.addDumpField("damage"      );
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );
  model.dump();*/

  for (UInt s = 0; s < max_steps; ++s) {
    model.solveStep();

    // if(s % 100 == 0) std::cout << "Step " << s+1 << "/" << max_steps
    // <<std::endl;  if(s % 100 == 0) model.dump();
  }

  const Vector<Real> & lower_bounds = mesh.getLowerBounds();
  const Vector<Real> & upper_bounds = mesh.getUpperBounds();
  Real L = upper_bounds(0) - lower_bounds(0);
  Real H = upper_bounds(1) - lower_bounds(1);

  const ElementTypeMapArray<UInt> & filter =
      model.getMaterial(0).getElementFilter();
  ElementTypeMapArray<UInt>::type_iterator it =
      filter.firstType(spatial_dimension);
  ElementTypeMapArray<UInt>::type_iterator end =
      filter.lastType(spatial_dimension);
  Vector<Real> barycenter(spatial_dimension);
  for (; it != end; ++it) {
    UInt nb_elem = mesh.getNbElement(*it);
    const UInt nb_gp = model.getFEEngine().getNbIntegrationPoints(*it);
    Array<Real> & material_damage_array =
        model.getMaterial(0).getArray<Real>("damage", *it);
    UInt cpt = 0;
    for (UInt nel = 0; nel < nb_elem; ++nel) {
      mesh.getBarycenter(nel, *it, barycenter.storage());
      if ((std::abs(barycenter(0) - (L / 2) + 0.025) < 0.025) &&
          (std::abs(barycenter(1) - (H / 2) + 0.025) < 0.025)) {
        if (material_damage_array(cpt, 0) < 0.9) {
          //	  std::cout << "barycenter(0) = " << barycenter(0) << ",
          //barycenter(1) = "
          //	    << barycenter(1)  <<std::endl;
          return EXIT_FAILURE;
        }
      }
      cpt += nb_gp;
    }
  }

  akantu::finalize();
  return EXIT_SUCCESS;
}
