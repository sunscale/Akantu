/**
 * @file   test_weight_computation.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Sep 23 16:30:12 2015
 *
 * @brief  test for non-local averaging of strain
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "test_material.hh"
#include "non_local_manager.hh"
#include "non_local_neighborhood.hh"
#include "dumper_paraview.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  akantu::initialize("material_weight_computation.dat", argc, argv);

  // some configuration variables
  const UInt spatial_dimension = 2;
  ElementType element_type = _quadrangle_4;
  GhostType ghost_type = _not_ghost;

  // mesh creation and read
  Mesh mesh(spatial_dimension);
  mesh.read("plate.msh");

  /// model creation
  SolidMechanicsModel  model(mesh);
 
  /// model initialization changed to use our material
  model.initFull(SolidMechanicsModelOptions(_static, true));
  model.registerNewCustomMaterials< TestMaterial<spatial_dimension> >("test_material");
  model.initMaterials();
  /// dump material index in paraview
  model.addDumpField("material_index");
  model.dump();

  /// apply constant strain field everywhere in the plate
  Matrix<Real> applied_strain(spatial_dimension, spatial_dimension);
  applied_strain.clear();
  for (UInt i = 0; i < spatial_dimension; ++i)
    applied_strain(i,i) = 2.;

  Array<Real> & grad_u = const_cast<Array<Real> &>(model.getMaterial(0).getInternal<Real>("grad_u")(element_type, ghost_type));
  Array<Real>::iterator< Matrix<Real> > grad_u_it = grad_u.begin(spatial_dimension, spatial_dimension);
  Array<Real>::iterator< Matrix<Real> > grad_u_end = grad_u.end(spatial_dimension, spatial_dimension);
  for (; grad_u_it != grad_u_end; ++grad_u_it) 
    (*grad_u_it) += applied_strain;

  /// compute the non-local strains
  model.getNonLocalManager().averageInternals(ghost_type);
  
  return EXIT_SUCCESS;
}
