/**
 * @file   structural_element_bernoulli_beam_2.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation  Tue Sep 19 2017
 *
 * @brief Specific functions for bernoulli beam 2d
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
#ifndef __AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_BEAM_2_HH__
#define __AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_BEAM_2_HH__

namespace akantu {

inline UInt
StructuralMechanicsModel::getTangentStiffnessVoigtSize<_bernoulli_beam_2>() {
  return 2;
}

/* -------------------------------------------------------------------------- */
template <>
inline void StructuralMechanicsModel::assembleMass<_bernoulli_beam_2>() {
  AKANTU_DEBUG_IN();

  GhostType ghost_type = _not_ghost;
  ElementType type = _bernoulli_beam_2;
  MyFEEngineType & fem = getFEEngineClass<MyFEEngineType>();
  UInt nb_element = getFEEngine().getMesh().getNbElement(type);
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = getFEEngine().getNbIntegrationPoints(type);
  UInt nb_fields_to_interpolate =
      getTangentStiffnessVoigtSize<_bernoulli_beam_2>();

  UInt nt_n_field_size = nb_degree_of_freedom * nb_nodes_per_element;

  Array<Real> * n =
      new Array<Real>(nb_element * nb_quadrature_points,
                      nb_fields_to_interpolate * nt_n_field_size, "N");
  n->clear();
  Array<Real> * rho_field =
      new Array<Real>(nb_element * nb_quadrature_points, "Rho");
  rho_field->clear();
  computeRho(*rho_field, type, _not_ghost);

  bool sign = true;

  /* ------------------------------------------------------------------------ */
  fem.computeShapesMatrix(type, nb_degree_of_freedom, nb_nodes_per_element, n,
                          0, 0, 0, sign, ghost_type); // Ni ui -> u
  fem.computeShapesMatrix(type, nb_degree_of_freedom, nb_nodes_per_element, n,
                          1, 1, 1, sign, ghost_type); // Mi vi -> v
  fem.computeShapesMatrix(type, nb_degree_of_freedom, nb_nodes_per_element, n,
                          2, 2, 1, sign, ghost_type); // Li Theta_i -> v
  /* ------------------------------------------------------------------------ */
  fem.assembleFieldMatrix(*rho_field, nb_degree_of_freedom, *mass_matrix, n,
                          rotation_matrix, type, ghost_type);

  delete n;
  delete rho_field;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <>
void StructuralMechanicsModel::computeRotationMatrix<_bernoulli_beam_2>(
    Array<Real> & rotations) {

  ElementType type = _bernoulli_beam_2;
  Mesh & mesh = getFEEngine().getMesh();
  UInt nb_element = mesh.getNbElement(type);

  Array<UInt>::iterator<Vector<UInt> > connec_it =
      mesh.getConnectivity(type).begin(2);
  Array<Real>::vector_iterator nodes_it =
      mesh.getNodes().begin(spatial_dimension);
  Array<Real>::matrix_iterator R_it =
      rotations.begin(nb_degree_of_freedom, nb_degree_of_freedom);

  for (UInt e = 0; e < nb_element; ++e, ++R_it, ++connec_it) {
    Matrix<Real> & R = *R_it;
    Vector<UInt> & connec = *connec_it;

    Vector<Real> x2;
    x2 = nodes_it[connec(1)]; // X2
    Vector<Real> x1;
    x1 = nodes_it[connec(0)]; // X1

    Real le = x1.distance(x2);
    Real c = (x2(0) - x1(0)) / le;
    Real s = (x2(1) - x1(1)) / le;

    /// Definition of the rotation matrix
    R(0, 0) = c;
    R(0, 1) = s;
    R(0, 2) = 0.;
    R(1, 0) = -s;
    R(1, 1) = c;
    R(1, 2) = 0.;
    R(2, 0) = 0.;
    R(2, 1) = 0.;
    R(2, 2) = 1.;
  }
}

/* -------------------------------------------------------------------------- */
template <>
void StructuralMechanicsModel::computeTangentModuli<_bernoulli_beam_2>(
    Array<Real> & tangent_moduli) {
  UInt nb_element = getFEEngine().getMesh().getNbElement(_bernoulli_beam_2);
  UInt nb_quadrature_points =
      getFEEngine().getNbIntegrationPoints(_bernoulli_beam_2);
  UInt tangent_size = 2;

  Array<Real>::matrix_iterator D_it =
      tangent_moduli.begin(tangent_size, tangent_size);
  Array<UInt> & el_mat = element_material(_bernoulli_beam_2, _not_ghost);

  for (UInt e = 0; e < nb_element; ++e) {
    UInt mat = el_mat(e);
    Real E = materials[mat].E;
    Real A = materials[mat].A;
    Real I = materials[mat].I;
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++D_it) {
      Matrix<Real> & D = *D_it;
      D(0, 0) = E * A;
      D(1, 1) = E * I;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <>
void StructuralMechanicsModel::transferBMatrixToSymVoigtBMatrix<
    _bernoulli_beam_2>(Array<Real> & b, bool) {
  UInt nb_element = getFEEngine().getMesh().getNbElement(_bernoulli_beam_2);
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(_bernoulli_beam_2);
  UInt nb_quadrature_points =
      getFEEngine().getNbIntegrationPoints(_bernoulli_beam_2);

  MyFEEngineType & fem = getFEEngineClass<MyFEEngineType>();
  Array<Real>::const_vector_iterator shape_Np =
      fem.getShapesDerivatives(_bernoulli_beam_2, _not_ghost, 0)
          .begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Mpp =
      fem.getShapesDerivatives(_bernoulli_beam_2, _not_ghost, 1)
          .begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Lpp =
      fem.getShapesDerivatives(_bernoulli_beam_2, _not_ghost, 2)
          .begin(nb_nodes_per_element);

  UInt tangent_size = getTangentStiffnessVoigtSize<_bernoulli_beam_2>();
  UInt bt_d_b_size = nb_nodes_per_element * nb_degree_of_freedom;
  b.clear();
  Array<Real>::matrix_iterator B_it = b.begin(tangent_size, bt_d_b_size);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      Matrix<Real> & B = *B_it;
      const Vector<Real> & Np = *shape_Np;
      const Vector<Real> & Lpp = *shape_Lpp;
      const Vector<Real> & Mpp = *shape_Mpp;

      B(0, 0) = Np(0);
      B(0, 3) = Np(1);

      B(1, 1) = Mpp(0);
      B(1, 2) = Lpp(0);
      B(1, 4) = Mpp(1);
      B(1, 5) = Lpp(1);

      ++B_it;
      ++shape_Np;
      ++shape_Mpp;
      ++shape_Lpp;
    }

    // ++R_it;
  }
}

}  // akantu

#endif /* __AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_BEAM_2_HH__ */
