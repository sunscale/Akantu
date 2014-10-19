/**
 * @file   test_solid_mechanics_model_material_eigenstrain.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Feb 10 2014
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  test the internal field prestrain
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

using namespace akantu;

Real alpha [3][4] = { { 0.01, 0.02, 0.03, 0.04 },
		      { 0.05, 0.06, 0.07, 0.08 },
		      { 0.09, 0.10, 0.11, 0.12 } };

/* -------------------------------------------------------------------------- */
template<ElementType type, bool is_plane_strain>
static Matrix<Real> prescribed_strain() {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  Matrix<Real> strain(spatial_dimension, spatial_dimension);

  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      strain(i, j) = alpha[i][j + 1];
    }
  }
  return strain;
}

template<ElementType type, bool is_plane_strain>
static Matrix<Real> prescribed_stress(Matrix<Real> prescribed_eigenstrain) {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  Matrix<Real> stress(spatial_dimension, spatial_dimension);

  //plane strain in 2d
  Matrix<Real> strain(spatial_dimension, spatial_dimension);
  Matrix<Real> pstrain;
  pstrain = prescribed_strain<type, is_plane_strain>();
  Real nu = 0.3;
  Real E  = 2.1e11;
  Real trace = 0;

  /// symetric part of the strain tensor
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      strain(i,j) = 0.5 * (pstrain(i, j) + pstrain(j, i));

  // elastic strain is equal to elastic strain minus the eigenstrain
  strain -= prescribed_eigenstrain;
  for (UInt i = 0; i < spatial_dimension; ++i) trace += strain(i,i);

  Real lambda   = nu * E / ((1 + nu) * (1 - 2*nu));
  Real mu       = E / (2 * (1 + nu));

  if(!is_plane_strain) {
    std::cout << "toto" << std::endl;
    lambda = nu * E / (1 - nu*nu);
  }

  if(spatial_dimension == 1) {
    stress(0, 0) =  E * strain(0, 0);
  } else {
    for (UInt i = 0; i < spatial_dimension; ++i)
      for (UInt j = 0; j < spatial_dimension; ++j) {
	stress(i, j) =  (i == j)*lambda*trace + 2*mu*strain(i, j);
      }
  }

  return stress;
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  initialize("material_elastic_plane_strain.dat", argc, argv);

  UInt dim = 3;
  const ElementType element_type = _tetrahedron_4;
  const bool plane_strain = true;
  Matrix<Real> prescribed_eigenstrain(dim, dim);
  prescribed_eigenstrain.clear();
  for (UInt i = 0; i < dim; ++i) {
    for (UInt j = 0; j < dim; ++j)
      prescribed_eigenstrain(i,j) += 0.1;
  }


  /// load mesh
  Mesh my_mesh(dim);

  std::stringstream filename; filename << "cube_3d_tet_4.msh";
  my_mesh.read(filename.str());

  /// declaration of model
  SolidMechanicsModel  my_model(my_mesh);
  /// model initialization
  my_model.initFull(SolidMechanicsModelOptions(_static));

  const Array<Real> & coordinates = my_mesh.getNodes();
  Array<Real> & displacement = my_model.getDisplacement();
  Array<bool> & boundary = my_model.getBlockedDOFs();
  MeshUtils::buildFacets(my_mesh);

  my_mesh.createBoundaryGroupFromGeometry();

  // Loop over (Sub)Boundar(ies)
  for(GroupManager::const_element_group_iterator it(my_mesh.element_group_begin());
      it != my_mesh.element_group_end(); ++it) {
    for(ElementGroup::const_node_iterator nodes_it(it->second->node_begin());
	nodes_it!= it->second->node_end(); ++nodes_it) {
      UInt n(*nodes_it);
      std::cout << "Node " << *nodes_it << std::endl;
      for (UInt i = 0; i < dim; ++i) {
	displacement(n, i) = alpha[i][0];
	for (UInt j = 0; j < dim; ++j) {
	  displacement(n, i) += alpha[i][j + 1] * coordinates(n, j);
	}
	boundary(n, i) = true;
      }
    }
  }

  /* ------------------------------------------------------------------------ */
  /* Apply eigenstrain in each element                                          */
  /* ------------------------------------------------------------------------ */


  Array<Real> & eigenstrain_vect = const_cast<Array<Real> &>(my_model.getMaterial(0).getInternal("eigenstrain")(element_type));
  Array<Real>::iterator< Matrix<Real> > eigenstrain_it = eigenstrain_vect.begin(dim, dim);
  Array<Real>::iterator< Matrix<Real> > eigenstrain_end = eigenstrain_vect.end(dim, dim);

  for (; eigenstrain_it != eigenstrain_end; ++eigenstrain_it) {
    for (UInt i = 0; i < dim; ++i)
      for (UInt j = 0; j < dim; ++j)
	(*eigenstrain_it)(i,j) += prescribed_eigenstrain(i,j);
  }
  /* ------------------------------------------------------------------------ */
  /* Static solve                                                             */
  /* ------------------------------------------------------------------------ */
  my_model.solveStep<_scm_newton_raphson_tangent_modified, _scc_residual>(2e-4, 2);

  /* ------------------------------------------------------------------------ */
  /* Checks                                                                   */
  /* ------------------------------------------------------------------------ */


  Array<Real> & stress_vect = const_cast<Array<Real> &>(my_model.getMaterial(0).getStress(element_type));

  Array<Real>::iterator< Matrix<Real> > stress_it = stress_vect.begin(dim, dim);
  Array<Real>::iterator< Matrix<Real> > stress_end = stress_vect.end(dim, dim);

  Matrix<Real> presc_stress;
  presc_stress = prescribed_stress<element_type, plane_strain>(prescribed_eigenstrain);

  Real stress_tolerance = 1e-13;

  for (; stress_it != stress_end; ++stress_it) {
    Matrix<Real> & stress = *stress_it;
    Matrix<Real> diff(dim, dim);

    diff  = stress;
    diff -= presc_stress;
    Real stress_error = diff.norm<L_inf>() / stress.norm<L_inf>();

    if(stress_error > stress_tolerance) {
      std::cerr << "stress error: " << stress_error << " > " << stress_tolerance << std::endl;
      std::cerr << "stress: " << stress << std::endl
		<< "prescribed stress: " << presc_stress << std::endl;
      return EXIT_FAILURE;
    } else {
      std::cerr << "stress error: " << stress_error << " < " << stress_tolerance << std::endl;
    }


  }


  finalize();

  return EXIT_SUCCESS;
}
