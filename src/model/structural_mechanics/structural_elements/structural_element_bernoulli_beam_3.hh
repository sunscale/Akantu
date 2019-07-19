/**
 * @file   structural_element_bernoulli_beam_3.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Wed Oct 11 2017
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Specific functions for bernoulli beam 3d
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef __AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_BEAM_3_HH__
#define __AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_BEAM_3_HH__

#include "structural_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <>
inline void StructuralMechanicsModel::assembleMass<_bernoulli_beam_3>() {
  AKANTU_DEBUG_IN();
#if 0
  GhostType ghost_type = _not_ghost;
  ElementType type = _bernoulli_beam_3;
  MyFEEngineType & fem = getFEEngineClass<MyFEEngineType>();
  UInt nb_element = getFEEngine().getMesh().getNbElement(type);
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = getFEEngine().getNbIntegrationPoints(type);
  UInt nb_fields_to_interpolate =
      getTangentStiffnessVoigtSize<_bernoulli_beam_3>();

  UInt nt_n_field_size = nb_degree_of_freedom * nb_nodes_per_element;

  Array<Real> * n =
      new Array<Real>(nb_element * nb_quadrature_points,
                      nb_fields_to_interpolate * nt_n_field_size, "N");
  n->clear();
  Array<Real> * rho_field =
      new Array<Real>(nb_element * nb_quadrature_points, "Rho");
  rho_field->clear();
  computeRho(*rho_field, type, _not_ghost);

  /* --------------------------------------------------------------------------
   */
  fem.computeShapesMatrix(type, nb_degree_of_freedom, nb_nodes_per_element, n,
                          0, 0, 0, true, ghost_type); // Ni ui ->         u

  fem.computeShapesMatrix(type, nb_degree_of_freedom, nb_nodes_per_element, n,
                          1, 1, 1, true, ghost_type); // Mi vi ->         v
  fem.computeShapesMatrix(type, nb_degree_of_freedom, nb_nodes_per_element, n,
                          2, 5, 1, true, ghost_type); // Li Theta_z_i ->  v

  fem.computeShapesMatrix(type, nb_degree_of_freedom, nb_nodes_per_element, n,
                          1, 2, 2, true, ghost_type); // Mi wi ->        w
  fem.computeShapesMatrix(type, nb_degree_of_freedom, nb_nodes_per_element, n,
                          2, 4, 2, false, ghost_type); // -Li Theta_y_i ->  w

  fem.computeShapesMatrix(type, nb_degree_of_freedom, nb_nodes_per_element, n,
                          0, 3, 3, true, ghost_type); // Ni Theta_x_i->Theta_x
  /* --------------------------------------------------------------------------
   */

  fem.assembleFieldMatrix(*rho_field, nb_degree_of_freedom, *mass_matrix, n,
                          rotation_matrix, type, ghost_type);

  delete n;
  delete rho_field;
#endif
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <>
void StructuralMechanicsModel::computeRotationMatrix<_bernoulli_beam_3>(
    Array<Real> & rotations) {
  ElementType type = _bernoulli_beam_3;
  Mesh & mesh = getFEEngine().getMesh();
  UInt nb_element = mesh.getNbElement(type);

  auto n_it = mesh.getNormals(type).begin(spatial_dimension);
  Array<UInt>::iterator<Vector<UInt>> connec_it =
      mesh.getConnectivity(type).begin(2);
  auto nodes_it = mesh.getNodes().begin(spatial_dimension);

  Matrix<Real> Pe(spatial_dimension, spatial_dimension);
  Matrix<Real> Pg(spatial_dimension, spatial_dimension);
  Matrix<Real> inv_Pg(spatial_dimension, spatial_dimension);
  Vector<Real> x_n(spatial_dimension); // x vect n

  Array<Real>::matrix_iterator R_it =
      rotations.begin(nb_degree_of_freedom, nb_degree_of_freedom);

  for (UInt e = 0; e < nb_element; ++e, ++n_it, ++connec_it, ++R_it) {
    Vector<Real> & n = *n_it;
    Matrix<Real> & R = *R_it;
    Vector<UInt> & connec = *connec_it;

    Vector<Real> x;
    x = nodes_it[connec(1)]; // X2
    Vector<Real> y;
    y = nodes_it[connec(0)]; // X1

    Real l = x.distance(y);
    x -= y; // X2 - X1
    x_n.crossProduct(x, n);

    Pe.eye();
    Pe(0, 0) *= l;
    Pe(1, 1) *= -l;

    Pg(0, 0) = x(0);
    Pg(0, 1) = x_n(0);
    Pg(0, 2) = n(0);
    Pg(1, 0) = x(1);
    Pg(1, 1) = x_n(1);
    Pg(1, 2) = n(1);
    Pg(2, 0) = x(2);
    Pg(2, 1) = x_n(2);
    Pg(2, 2) = n(2);

    inv_Pg.inverse(Pg);

    Pe *= inv_Pg;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = 0; j < spatial_dimension; ++j) {
        R(i, j) = Pe(i, j);
        R(i + spatial_dimension, j + spatial_dimension) = Pe(i, j);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <>
void StructuralMechanicsModel::computeTangentModuli<_bernoulli_beam_3>(
    Array<Real> & tangent_moduli) {
  UInt nb_element = getFEEngine().getMesh().getNbElement(_bernoulli_beam_3);
  UInt nb_quadrature_points =
      getFEEngine().getNbIntegrationPoints(_bernoulli_beam_3);
  UInt tangent_size = 4;

  tangent_moduli.clear();
  Array<Real>::matrix_iterator D_it =
      tangent_moduli.begin(tangent_size, tangent_size);

  for (UInt e = 0; e < nb_element; ++e) {
    UInt mat = element_material(_bernoulli_beam_3, _not_ghost)(e);
    Real E = materials[mat].E;
    Real A = materials[mat].A;
    Real Iz = materials[mat].Iz;
    Real Iy = materials[mat].Iy;
    Real GJ = materials[mat].GJ;
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++D_it) {
      Matrix<Real> & D = *D_it;
      D(0, 0) = E * A;
      D(1, 1) = E * Iz;
      D(2, 2) = E * Iy;
      D(3, 3) = GJ;
    }
  }
}

} // namespace akantu

#endif /* __AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_BEAM_3_HH__ */
