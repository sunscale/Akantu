/**
 * @file   structural_mechanics_model_inline_impl.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Implementation of inline functions of StructuralMechanicsModel
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "structural_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_STRUCTURAL_MECHANICS_MODEL_INLINE_IMPL_HH_
#define AKANTU_STRUCTURAL_MECHANICS_MODEL_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <ElementType type>
void StructuralMechanicsModel::computeTangentModuli(
    Array<Real> & /*tangent_moduli*/) {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void StructuralMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  auto nb_element = getFEEngine().getMesh().getNbElement(type);
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_quadrature_points = getFEEngine().getNbIntegrationPoints(type);
  auto tangent_size = ElementClass<type>::getNbStressComponents();

  auto tangent_moduli = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, tangent_size * tangent_size,
      "tangent_stiffness_matrix");
  computeTangentModuli<type>(*tangent_moduli);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  UInt bt_d_b_size = nb_degree_of_freedom * nb_nodes_per_element;

  auto bt_d_b = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, bt_d_b_size * bt_d_b_size, "B^t*D*B");

  const auto & b = getFEEngine().getShapesDerivatives(type);

  Matrix<Real> BtD(bt_d_b_size, tangent_size);

  for (auto && tuple :
       zip(make_view(b, tangent_size, bt_d_b_size),
           make_view(*tangent_moduli, tangent_size, tangent_size),
           make_view(*bt_d_b, bt_d_b_size, bt_d_b_size))) {
    auto & B = std::get<0>(tuple);
    auto & D = std::get<1>(tuple);
    auto & BtDB = std::get<2>(tuple);
    BtD.mul<true, false>(B, D);
    BtDB.template mul<false, false>(BtD, B);
  }

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  auto int_bt_d_b = std::make_unique<Array<Real>>(
      nb_element, bt_d_b_size * bt_d_b_size, "int_B^t*D*B");

  getFEEngine().integrate(*bt_d_b, *int_bt_d_b, bt_d_b_size * bt_d_b_size,
                          type);

  getDOFManager().assembleElementalMatricesToMatrix(
      "K", "displacement", *int_bt_d_b, type, _not_ghost, _symmetric);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void StructuralMechanicsModel::computeStressOnQuad() {
  AKANTU_DEBUG_IN();

  auto & sigma = stress(type, _not_ghost);

  auto nb_element = mesh.getNbElement(type);
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_quadrature_points = getFEEngine().getNbIntegrationPoints(type);
  auto tangent_size = ElementClass<type>::getNbStressComponents();

  auto tangent_moduli = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, tangent_size * tangent_size,
      "tangent_stiffness_matrix");

  computeTangentModuli<type>(*tangent_moduli);

  /// compute DB
  auto d_b_size = nb_degree_of_freedom * nb_nodes_per_element;

  auto d_b = std::make_unique<Array<Real>>(nb_element * nb_quadrature_points,
                                           d_b_size * tangent_size, "D*B");

  const auto & b = getFEEngine().getShapesDerivatives(type);

  auto B = b.begin(tangent_size, d_b_size);
  auto D = tangent_moduli->begin(tangent_size, tangent_size);
  auto D_B = d_b->begin(tangent_size, d_b_size);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++B, ++D, ++D_B) {
      D_B->template mul<false, false>(*D, *B);
    }
  }

  /// compute DBu
  D_B = d_b->begin(tangent_size, d_b_size);
  auto DBu = sigma.begin(tangent_size);

  Array<Real> u_el(0, d_b_size);
  FEEngine::extractNodalToElementField(mesh, *displacement_rotation, u_el,
                                       type);

  auto ug = u_el.begin(d_b_size);

  // No need to rotate because B is post-multiplied
  for (UInt e = 0; e < nb_element; ++e, ++ug) {
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++D_B, ++DBu) {
      DBu->template mul<false>(*D_B, *ug);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void StructuralMechanicsModel::computeForcesByLocalTractionArray(
    const Array<Real> & tractions) {
  AKANTU_DEBUG_IN();

  auto nb_element = getFEEngine().getMesh().getNbElement(type);
  auto nb_nodes_per_element =
      getFEEngine().getMesh().getNbNodesPerElement(type);
  auto nb_quad = getFEEngine().getNbIntegrationPoints(type);

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

  Array<Real> Ntbs(nb_element * nb_quad,
                   nb_degree_of_freedom * nb_nodes_per_element);
  Array<Real> TtNtbs(nb_element * nb_quad,
                     nb_degree_of_freedom * nb_nodes_per_element);

  auto & fem = getFEEngine();
  fem.computeNtb(tractions, Ntbs, type);

  auto T_it =
      rotation_matrix(type).begin(nb_degree_of_freedom * nb_nodes_per_element,
                                  nb_degree_of_freedom * nb_nodes_per_element);
  auto Ntb_it = Ntbs.begin(nb_degree_of_freedom * nb_nodes_per_element);
  auto TtNtb_it = TtNtbs.begin(nb_degree_of_freedom * nb_nodes_per_element);

  for (UInt e = 0; e < nb_element; ++e, ++T_it) {
    const auto & T = *T_it;
    for (UInt q = 0; q < nb_quad; ++q, ++Ntb_it, ++TtNtb_it) {
      const auto & Ntb = *Ntb_it;
      auto & TtNtb = *TtNtb_it;

      // turn N^t tl back in the global referential
      TtNtb.template mul<true>(T, Ntb);
    }
  }

  // allocate the vector that will contain the integrated values
  auto name = id + std::to_string(type) + ":integral_boundary";
  Array<Real> int_funct(nb_element, nb_degree_of_freedom * nb_nodes_per_element,
                        name);

  // do the integration
  getFEEngine().integrate(TtNtbs, int_funct,
                          nb_degree_of_freedom * nb_nodes_per_element, type);

  // assemble the result into force vector
  getDOFManager().assembleElementalArrayLocalArray(int_funct, *external_force,
                                                   type, _not_ghost, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void StructuralMechanicsModel::computeForcesByGlobalTractionArray(
    const Array<Real> & traction_global) {
  AKANTU_DEBUG_IN();

  UInt nb_element = mesh.getNbElement(type);
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
    const auto & T = *T_it;
    for (UInt i = 0; i < nb_degree_of_freedom; ++i) {
      for (UInt j = 0; j < nb_degree_of_freedom; ++j) {
        R(i, j) = T(i, j);
      }
    }

    for (UInt q = 0; q < nb_quad; ++q, ++Te_it, ++te_it) {
      const auto & Te = *Te_it;
      auto & te = *te_it;
      // turn the traction in the local referential
      te.template mul<false>(R, Te);
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
#if 0
template <ElementType type>
inline void StructuralMechanicsModel::computeForcesFromFunction(
    BoundaryFunction myf, BoundaryFunctionType function_type) {
  /** function type is
   ** _bft_forces : linear load is given
   ** _bft_stress : stress function is given -> Not already done for this kind
   *of model
   */

  std::stringstream name;
  name << id << ":structuralmechanics:imposed_linear_load";
  Array<Real> lin_load(0, nb_degree_of_freedom, name.str());
  name.zero();

  UInt offset = nb_degree_of_freedom;

  // prepare the loop over element types
  UInt nb_quad = getFEEngine().getNbIntegrationPoints(type);
  UInt nb_element = getFEEngine().getMesh().getNbElement(type);

  name.zero();
  name << id << ":structuralmechanics:quad_coords";
  Array<Real> quad_coords(nb_element * nb_quad, spatial_dimension,
                          "quad_coords");

  getFEEngineClass<MyFEEngineType>()
      .getShapeFunctions()
      .interpolateOnIntegrationPoints<type>(getFEEngine().getMesh().getNodes(),
                                            quad_coords, spatial_dimension);
  getFEEngineClass<MyFEEngineType>()
      .getShapeFunctions()
      .interpolateOnIntegrationPoints<type>(
          getFEEngine().getMesh().getNodes(), quad_coords, spatial_dimension,
          _not_ghost, empty_filter, true, 0, 1, 1);
  if (spatial_dimension == 3)
    getFEEngineClass<MyFEEngineType>()
        .getShapeFunctions()
        .interpolateOnIntegrationPoints<type>(
            getFEEngine().getMesh().getNodes(), quad_coords, spatial_dimension,
            _not_ghost, empty_filter, true, 0, 2, 2);
  lin_load.resize(nb_element * nb_quad);
  Real * imposed_val = lin_load.storage();

  /// sigma/load on each quadrature points
  Real * qcoord = quad_coords.storage();
  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt q = 0; q < nb_quad; ++q) {
      myf(qcoord, imposed_val, NULL, 0);
      imposed_val += offset;
      qcoord += spatial_dimension;
    }
  }

  switch (function_type) {
  case _bft_traction_local:
    computeForcesByLocalTractionArray<type>(lin_load);
    break;
  case _bft_traction:
    computeForcesByGlobalTractionArray<type>(lin_load);
    break;
  default:
    break;
  }
}
#endif
} // namespace akantu

#endif /* AKANTU_STRUCTURAL_MECHANICS_MODEL_INLINE_IMPL_HH_ */
