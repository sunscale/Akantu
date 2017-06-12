/**
 * @file   test_local_material.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Clement Roux <clement.roux@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  test of the class SolidMechanicsModel with custom local damage on a
 * notched plate
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
#include "local_material_damage.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  akantu::initialize("material.dat", argc, argv);
  UInt max_steps = 1100;

  const UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("mesh_section_gap.msh");
  
  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull(SolidMechanicsModelOptions(_explicit_lumped_mass, true));
  model.registerNewCustomMaterials<LocalMaterialDamage>("local_damage");
  model.initMaterials();

  Real time_step = model.getStableTimeStep();
  model.setTimeStep(time_step / 2.5);

  model.assembleMassLumped();

  std::cout << model << std::endl;

  /// Dirichlet boundary conditions
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "Fixed");
  // model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "Fixed");

  // Boundary condition (Neumann)
  Matrix<Real> stress(2, 2);
  stress.eye(7e5);
  model.applyBC(BC::Neumann::FromHigherDim(stress), "Traction");

  for (UInt s = 0; s < max_steps; ++s) {
    if (s < 100) {
      // Boundary condition (Neumann)
      stress.eye(7e5);
      model.applyBC(BC::Neumann::FromHigherDim(stress), "Traction");
    }

    model.solveStep();
  }

  const Vector<Real> & lower_bounds = mesh.getLowerBounds();
  const Vector<Real> & upper_bounds = mesh.getUpperBounds();
  Real L = upper_bounds(0) - lower_bounds(0);

  const ElementTypeMapArray<UInt> & filter =
      model.getMaterial(0).getElementFilter();
  ElementTypeMapArray<UInt>::type_iterator it =
      filter.firstType(spatial_dimension);
  ElementTypeMapArray<UInt>::type_iterator end =
      filter.lastType(spatial_dimension);
  Vector<Real> barycenter(spatial_dimension);
  bool is_fully_damaged = false;
  for (; it != end; ++it) {
    UInt nb_elem = mesh.getNbElement(*it);
    const UInt nb_gp = model.getFEEngine().getNbIntegrationPoints(*it);
    Array<Real> & material_damage_array =
        model.getMaterial(0).getArray<Real>("damage", *it);
    UInt cpt = 0;
    for (UInt nel = 0; nel < nb_elem; ++nel) {
      if (material_damage_array(cpt, 0) > 0.9) {
        is_fully_damaged = true;
        mesh.getBarycenter(nel, *it, barycenter.storage());
        if ((std::abs(barycenter(0) - (L / 2)) < (L / 10))) {
          return EXIT_FAILURE;
        }
      }
      cpt += nb_gp;
    }
  }
  if (!is_fully_damaged)
    return EXIT_FAILURE;

  akantu::finalize();
  return EXIT_SUCCESS;
}
