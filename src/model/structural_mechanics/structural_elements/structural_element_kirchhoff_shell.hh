/**
 * @file   structural_element_bernoulli_kirchhoff_shell.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation  Tue Sep 19 2017
 *
 * @brief Specific functions for bernoulli kirchhoff shell
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
#ifndef __AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_KIRCHHOFF_SHELL_HH__
#define __AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_KIRCHHOFF_SHELL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <>
inline void StructuralMechanicsModel::assembleMass<_discrete_kirchhoff_triangle_18>() {

  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
void StructuralMechanicsModel::computeRotationMatrix<_discrete_kirchhoff_triangle_18>(
    Array<Real> & rotations) {
  ElementType type = _discrete_kirchhoff_triangle_18;
  Mesh & mesh = getFEEngine().getMesh();
  UInt nb_element = mesh.getNbElement(type);

  Array<UInt>::iterator<Vector<UInt> > connec_it =
      mesh.getConnectivity(type).begin(3);
  Array<Real>::vector_iterator nodes_it =
      mesh.getNodes().begin(spatial_dimension);

  Matrix<Real> Pe(spatial_dimension, spatial_dimension);
  Matrix<Real> Pg(spatial_dimension, spatial_dimension);
  Matrix<Real> inv_Pg(spatial_dimension, spatial_dimension);

  Array<Real>::matrix_iterator R_it =
      rotations.begin(nb_degree_of_freedom, nb_degree_of_freedom);

  for (UInt e = 0; e < nb_element; ++e, ++connec_it, ++R_it) {

    Pe.eye();

    Matrix<Real> & R = *R_it;
    Vector<UInt> & connec = *connec_it;

    Vector<Real> x2;
    x2 = nodes_it[connec(1)]; // X2
    Vector<Real> x1;
    x1 = nodes_it[connec(0)]; // X1
    Vector<Real> x3;
    x3 = nodes_it[connec(2)]; // X3

    Vector<Real> Pg_col_1 = x2 - x1;

    Vector<Real> Pg_col_2 = x3 - x1;

    Vector<Real> Pg_col_3(spatial_dimension);
    Pg_col_3.crossProduct(Pg_col_1, Pg_col_2);

    for (UInt i = 0; i < spatial_dimension; ++i) {
      Pg(i, 0) = Pg_col_1(i);
      Pg(i, 1) = Pg_col_2(i);
      Pg(i, 2) = Pg_col_3(i);
    }

    inv_Pg.inverse(Pg);
    // Pe *= inv_Pg;
    Pe.eye();

    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = 0; j < spatial_dimension; ++j) {
        R(i, j) = Pe(i, j);
        R(i + spatial_dimension, j + spatial_dimension) = Pe(i, j);
      }
    }
  }
}

}  // akantu

#endif /* __AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_DISCRETE_KIRCHHOFF_TRIANGLE_18_HH__ */
