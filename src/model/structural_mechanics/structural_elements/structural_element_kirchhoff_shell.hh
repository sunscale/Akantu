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

inline UInt
StructuralMechanicsModel::getTangentStiffnessVoigtSize<_bernoulli_kirchhoff_shell>() {
  return 6;
}

/* -------------------------------------------------------------------------- */
template <>
inline void StructuralMechanicsModel::assembleMass<_kirchhoff_shell>() {

  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
void StructuralMechanicsModel::computeRotationMatrix<_kirchhoff_shell>(
    Array<Real> & rotations) {
  ElementType type = _kirchhoff_shell;
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

/* -------------------------------------------------------------------------- */
template <>
void StructuralMechanicsModel::computeTangentModuli<_bernoulli_beam_3>(
    Array<Real> & tangent_moduli) {
  UInt nb_element = getFEEngine().getMesh().getNbElement(_bernoulli_beam_3);
  UInt nb_quadrature_points =
      getFEEngine().getNbIntegrationPoints(_bernoulli_beam_3);
  UInt tangent_size = 4;

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

/* -------------------------------------------------------------------------- */
template <>
void StructuralMechanicsModel::transferBMatrixToSymVoigtBMatrix<
    _kirchhoff_shell>(Array<Real> & b, __attribute__((unused)) bool local) {
  MyFEEngineType & fem = getFEEngineClass<MyFEEngineType>();

  UInt nb_element = getFEEngine().getMesh().getNbElement(_kirchhoff_shell);
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(_kirchhoff_shell);
  UInt nb_quadrature_points =
      getFEEngine().getNbIntegrationPoints(_kirchhoff_shell);

  Array<Real>::const_matrix_iterator shape_Np =
      fem.getShapesDerivatives(_kirchhoff_shell, _not_ghost, 0)
          .begin(2, nb_nodes_per_element);
  Array<Real>::const_matrix_iterator shape_Nx1p =
      fem.getShapesDerivatives(_kirchhoff_shell, _not_ghost, 1)
          .begin(2, nb_nodes_per_element);
  Array<Real>::const_matrix_iterator shape_Nx2p =
      fem.getShapesDerivatives(_kirchhoff_shell, _not_ghost, 2)
          .begin(2, nb_nodes_per_element);
  Array<Real>::const_matrix_iterator shape_Nx3p =
      fem.getShapesDerivatives(_kirchhoff_shell, _not_ghost, 3)
          .begin(2, nb_nodes_per_element);
  Array<Real>::const_matrix_iterator shape_Ny1p =
      fem.getShapesDerivatives(_kirchhoff_shell, _not_ghost, 4)
          .begin(2, nb_nodes_per_element);
  Array<Real>::const_matrix_iterator shape_Ny2p =
      fem.getShapesDerivatives(_kirchhoff_shell, _not_ghost, 5)
          .begin(2, nb_nodes_per_element);
  Array<Real>::const_matrix_iterator shape_Ny3p =
      fem.getShapesDerivatives(_kirchhoff_shell, _not_ghost, 6)
          .begin(2, nb_nodes_per_element);

  UInt tangent_size = getTangentStiffnessVoigtSize<_kirchhoff_shell>();
  UInt bt_d_b_size = nb_nodes_per_element * nb_degree_of_freedom;

  b.clear();

  Array<Real>::matrix_iterator B_it = b.begin(tangent_size, bt_d_b_size);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      Matrix<Real> & B = *B_it;

      const Matrix<Real> & Np = *shape_Np;
      const Matrix<Real> & Nx1p = *shape_Nx1p;
      const Matrix<Real> & Nx2p = *shape_Nx2p;
      const Matrix<Real> & Nx3p = *shape_Nx3p;
      const Matrix<Real> & Ny1p = *shape_Ny1p;
      const Matrix<Real> & Ny2p = *shape_Ny2p;
      const Matrix<Real> & Ny3p = *shape_Ny3p;

      B(0, 0) = Np(0, 0);
      B(0, 6) = Np(0, 1);
      B(0, 12) = Np(0, 2);

      B(1, 1) = Np(1, 0);
      B(1, 7) = Np(1, 1);
      B(1, 13) = Np(1, 2);

      B(2, 0) = Np(1, 0);
      B(2, 1) = Np(0, 0);
      B(2, 6) = Np(1, 1);
      B(2, 7) = Np(0, 1);
      B(2, 12) = Np(1, 2);
      B(2, 13) = Np(0, 2);

      B(3, 2) = Nx1p(0, 0);
      B(3, 3) = Nx2p(0, 0);
      B(3, 4) = Nx3p(0, 0);
      B(3, 8) = Nx1p(0, 1);
      B(3, 9) = Nx2p(0, 1);
      B(3, 10) = Nx3p(0, 1);
      B(3, 14) = Nx1p(0, 2);
      B(3, 15) = Nx2p(0, 2);
      B(3, 16) = Nx3p(0, 2);

      B(4, 2) = Ny1p(1, 0);
      B(4, 3) = Ny2p(1, 0);
      B(4, 4) = Ny3p(1, 0);
      B(4, 8) = Ny1p(1, 1);
      B(4, 9) = Ny2p(1, 1);
      B(4, 10) = Ny3p(1, 1);
      B(4, 14) = Ny1p(1, 2);
      B(4, 15) = Ny2p(1, 2);
      B(4, 16) = Ny3p(1, 2);

      B(5, 2) = Nx1p(1, 0) + Ny1p(0, 0);
      B(5, 3) = Nx2p(1, 0) + Ny2p(0, 0);
      B(5, 4) = Nx3p(1, 0) + Ny3p(0, 0);
      B(5, 8) = Nx1p(1, 1) + Ny1p(0, 1);
      B(5, 9) = Nx2p(1, 1) + Ny2p(0, 1);
      B(5, 10) = Nx3p(1, 1) + Ny3p(0, 1);
      B(5, 14) = Nx1p(1, 2) + Ny1p(0, 2);
      B(5, 15) = Nx2p(1, 2) + Ny2p(0, 2);
      B(5, 16) = Nx3p(1, 2) + Ny3p(0, 2);

      ++B_it;
      ++shape_Np;
      ++shape_Nx1p;
      ++shape_Nx2p;
      ++shape_Nx3p;
      ++shape_Ny1p;
      ++shape_Ny2p;
      ++shape_Ny3p;
    }
  }
}

}  // akantu

#endif /* __AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_KIRCHHOFF_SHELL_HH__ */
