/**
 * @file   structural_mechanics_model_inline_impl.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  Implementation of inline functions of StructuralMechanicsModel
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
#include "structural_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_STRUCTURAL_MECHANICS_MODEL_INLINE_IMPL_CC__
#define __AKANTU_STRUCTURAL_MECHANICS_MODEL_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline UInt StructuralMechanicsModel::getTangentStiffnessVoigtSize() {
  AKANTU_DEBUG_TO_IMPLEMENT();
  return 0;
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void StructuralMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  UInt nb_element = getFEEngine().getMesh().getNbElement(type);
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = getFEEngine().getNbIntegrationPoints(type);

  UInt tangent_size = getTangentStiffnessVoigtSize<type>();

  Array<Real> tangent_moduli(nb_element * nb_quadrature_points,
                             tangent_size * tangent_size, 0.,
                             "tangent_stiffness_matrix");

  computeTangentModuli<type>(tangent_moduli);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  UInt bt_d_b_size = nb_degree_of_freedom * nb_nodes_per_element;

  Array<Real> bt_d_b(nb_element * nb_quadrature_points,
                     bt_d_b_size * bt_d_b_size, "B^t*D*B");

  Array<Real> b(nb_element * nb_quadrature_points, tangent_size * bt_d_b_size,
                "B");

  transferBMatrixToSymVoigtBMatrix<type>(b);

  Matrix<Real> Bt_D(bt_d_b_size, tangent_size);
  Matrix<Real> BT(tangent_size, bt_d_b_size);

  auto B = b.begin(tangent_size, bt_d_b_size);
  auto D = tangent_moduli.begin(tangent_size, tangent_size);
  auto Bt_D_B = bt_d_b.begin(bt_d_b_size, bt_d_b_size);
  auto T = rotation_matrix(type).begin(bt_d_b_size, bt_d_b_size);

  for (UInt e = 0; e < nb_element; ++e, ++T) {
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++B, ++D, ++Bt_D_B) {
      BT.mul<false, false>(*B, *T);
      Bt_D.mul<true, false>(BT, *D);
      Bt_D_B->mul<false, false>(Bt_D, BT);
    }
  }

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  Array<Real> int_bt_d_b(nb_element, bt_d_b_size * bt_d_b_size, "int_B^t*D*B");

  getFEEngine().integrate(bt_d_b, int_bt_d_b, bt_d_b_size * bt_d_b_size, type);

  getDOFManager().assembleElementalMatricesToMatrix("K", "displacement",
                                                    int_bt_d_b, type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void StructuralMechanicsModel::computeTangentModuli(Array<Real> &) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void StructuralMechanicsModel::transferBMatrixToSymVoigtBMatrix(Array<Real> &,
                                                                bool) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void StructuralMechanicsModel::computeStressOnQuad() {
  AKANTU_DEBUG_IN();

  Array<Real> & sigma = stress(type, _not_ghost);

  sigma.clear();
  const Mesh & mesh = getFEEngine().getMesh();

  UInt nb_element = mesh.getNbElement(type);
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = getFEEngine().getNbIntegrationPoints(type);

  UInt tangent_size = getTangentStiffnessVoigtSize<type>();

  Array<Real> * tangent_moduli =
      new Array<Real>(nb_element * nb_quadrature_points,
                      tangent_size * tangent_size, "tangent_stiffness_matrix");
  tangent_moduli->clear();

  computeTangentModuli<type>(*tangent_moduli);

  /// compute DB
  UInt d_b_size = nb_degree_of_freedom * nb_nodes_per_element;

  Array<Real> * d_b = new Array<Real>(nb_element * nb_quadrature_points,
                                      d_b_size * tangent_size, "D*B");

  Array<Real> * b = new Array<Real>(nb_element * nb_quadrature_points,
                                    tangent_size * d_b_size, "B");

  transferBMatrixToSymVoigtBMatrix<type>(*b);

  Array<Real>::matrix_iterator B = b->begin(tangent_size, d_b_size);
  Array<Real>::matrix_iterator D =
      tangent_moduli->begin(tangent_size, tangent_size);
  Array<Real>::matrix_iterator D_B = d_b->begin(tangent_size, d_b_size);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++B, ++D, ++D_B) {
      D_B->mul<false, false>(*D, *B);
    }
  }

  delete b;
  delete tangent_moduli;

  /// compute DBu
  D_B = d_b->begin(tangent_size, d_b_size);
  Array<Real>::iterator<Vector<Real>> DBu = sigma.begin(tangent_size);
  Vector<Real> ul(d_b_size);

  Array<Real> u_el(0, d_b_size);
  FEEngine::extractNodalToElementField(mesh, *displacement_rotation, u_el,
                                       type);

  Array<Real>::vector_iterator ug = u_el.begin(d_b_size);
  Array<Real>::matrix_iterator T =
      rotation_matrix(type).begin(d_b_size, d_b_size);

  for (UInt e = 0; e < nb_element; ++e, ++T, ++ug) {
    ul.mul<false>(*T, *ug);
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++D_B, ++DBu) {
      DBu->mul<false>(*D_B, ul);
    }
  }

  delete d_b;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void StructuralMechanicsModel::computeForcesByLocalTractionArray(
    const Array<Real> & tractions) {
  AKANTU_DEBUG_IN();

  UInt nb_element = getFEEngine().getMesh().getNbElement(type);
  UInt nb_nodes_per_element =
      getFEEngine().getMesh().getNbNodesPerElement(type);
  UInt nb_quad = getFEEngine().getNbIntegrationPoints(type);

  // check dimension match
  AKANTU_DEBUG_ASSERT(
      Mesh::getSpatialDimension(type) == getFEEngine().getElementDimension(),
      "element type dimension does not match the dimension of boundaries : "
          << getFEEngine().getElementDimension()
          << " != " << Mesh::getSpatialDimension(type));

  // check size of the vector
  AKANTU_DEBUG_ASSERT(
      tractions.size() == nb_quad * nb_element,
      "the size of the vector should be the total number of quadrature points");

  // check number of components
  AKANTU_DEBUG_ASSERT(tractions.getNbComponent() == nb_degree_of_freedom,
                      "the number of components should be the spatial "
                      "dimension of the problem");

  Array<Real> Nvoigt(nb_element * nb_quad, nb_degree_of_freedom *
                                               nb_degree_of_freedom *
                                               nb_nodes_per_element);
  transferNMatrixToSymVoigtNMatrix<type>(Nvoigt);

  auto N_it = Nvoigt.begin(nb_degree_of_freedom,
                           nb_degree_of_freedom * nb_nodes_per_element);
  auto T_it =
      rotation_matrix(type).begin(nb_degree_of_freedom * nb_nodes_per_element,
                                  nb_degree_of_freedom * nb_nodes_per_element);
  auto te_it = tractions.begin(nb_degree_of_freedom);

  Array<Real> funct(nb_element * nb_quad,
                    nb_degree_of_freedom * nb_nodes_per_element, 0.);
  auto Fe_it = funct.begin(nb_degree_of_freedom * nb_nodes_per_element);

  Vector<Real> fe(nb_degree_of_freedom * nb_nodes_per_element);
  for (UInt e = 0; e < nb_element; ++e, ++T_it) {
    const Matrix<Real> & T = *T_it;
    for (UInt q = 0; q < nb_quad; ++q, ++N_it, ++te_it, ++Fe_it) {
      const Matrix<Real> & N = *N_it;
      const Vector<Real> & te = *te_it;
      Vector<Real> & Fe = *Fe_it;

      // compute N^t tl
      fe.mul<true>(N, te);
      // turn N^t tl back in the global referential
      Fe.mul<true>(T, fe);
    }
  }

  // allocate the vector that will contain the integrated values
  Array<Real> int_funct(nb_element, nb_degree_of_freedom * nb_nodes_per_element,
                        std::string(id + ":integral_boundary"));

  // do the integration
  getFEEngine().integrate(funct, int_funct,
                          nb_degree_of_freedom * nb_nodes_per_element, type);

  // assemble the result into force vector
  getDOFManager().assembleElementalArrayLocalArray(
      int_funct, *internal_force_momentum, type, _not_ghost);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void StructuralMechanicsModel::computeForcesByGlobalTractionArray(
    const Array<Real> & traction_global) {
  AKANTU_DEBUG_IN();
  UInt nb_element = getFEEngine().getMesh().getNbElement(type);
  UInt nb_quad = getFEEngine().getNbIntegrationPoints(type);
  UInt nb_nodes_per_element =
      getFEEngine().getMesh().getNbNodesPerElement(type);

  std::stringstream name;
  name << id << ":structuralmechanics:imposed_linear_load";
  Array<Real> traction_local(nb_element * nb_quad, nb_degree_of_freedom,
                             name.str());

  auto T_it =
      rotation_matrix(type).begin(nb_degree_of_freedom * nb_nodes_per_element,
                                  nb_degree_of_freedom * nb_nodes_per_element);

  auto Te_it = traction_global.begin(nb_degree_of_freedom);
  auto te_it = traction_local.begin(nb_degree_of_freedom);

  Matrix<Real> R(nb_degree_of_freedom, nb_degree_of_freedom);
  for (UInt e = 0; e < nb_element; ++e, ++T_it) {
    const Matrix<Real> & T = *T_it;
    for (UInt i = 0; i < nb_degree_of_freedom; ++i)
      for (UInt j = 0; j < nb_degree_of_freedom; ++j)
        R(i, j) = T(i, j);

    for (UInt q = 0; q < nb_quad; ++q, ++Te_it, ++te_it) {
      const Vector<Real> & Te = *Te_it;
      Vector<Real> & te = *te_it;
      // turn the traction in the local referential
      te.mul<false>(R, Te);
    }
  }

  computeForcesByLocalTractionArray<type>(traction_local);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * @param myf pointer  to a function that fills a  vector/tensor with respect to
 * passed coordinates
 */
// template <ElementType type>
// inline void StructuralMechanicsModel::computeForcesFromFunction(
//     BoundaryFunction myf, BoundaryFunctionType function_type) {
//   /** function type is
//    ** _bft_forces : linear load is given
//    ** _bft_stress : stress function is given -> Not already done for this kind
//    *of model
//    */

//   std::stringstream name;
//   name << id << ":structuralmechanics:imposed_linear_load";
//   Array<Real> lin_load(0, nb_degree_of_freedom, name.str());
//   name.clear();

//   UInt offset = nb_degree_of_freedom;

//   // prepare the loop over element types
//   UInt nb_quad = getFEEngine().getNbIntegrationPoints(type);
//   UInt nb_element = getFEEngine().getMesh().getNbElement(type);

//   name.clear();
//   name << id << ":structuralmechanics:quad_coords";
//   Array<Real> quad_coords(nb_element * nb_quad, spatial_dimension,
//                           "quad_coords");

//   getFEEngineClass<MyFEEngineType>()
//       .getShapeFunctions()
//       .interpolateOnIntegrationPoints<type>(mesh.getNodes(), quad_coords,
//                                             spatial_dimension);
//   getFEEngineClass<MyFEEngineType>()
//       .getShapeFunctions()
//       .interpolateOnIntegrationPoints<type>(mesh.getNodes(), quad_coords,
//                                             spatial_dimension, _not_ghost,
//                                             empty_filter, true, 0, 1, 1);
//   if (spatial_dimension == 3)
//     getFEEngineClass<MyFEEngineType>()
//         .getShapeFunctions()
//         .interpolateOnIntegrationPoints<type>(
//             getFEEngine().getMesh().getNodes(), quad_coords, spatial_dimension,
//             _not_ghost, empty_filter, true, 0, 2, 2);
//   lin_load.resize(nb_element * nb_quad);
//   Real * imposed_val = lin_load.storage();

//   /// sigma/load on each quadrature points
//   Real * qcoord = quad_coords.storage();
//   for (UInt el = 0; el < nb_element; ++el) {
//     for (UInt q = 0; q < nb_quad; ++q) {
//       myf(qcoord, imposed_val, NULL, 0);
//       imposed_val += offset;
//       qcoord += spatial_dimension;
//     }
//   }

//   switch (function_type) {
//   case _bft_traction_local:
//     computeForcesByLocalTractionArray<type>(lin_load);
//     break;
//   case _bft_traction:
//     computeForcesByGlobalTractionArray<type>(lin_load);
//     break;
//   default:
//     break;
//   }
// }

} // namespace akantu

#endif /* __AKANTU_STRUCTURAL_MECHANICS_MODEL_INLINE_IMPL_CC__ */
