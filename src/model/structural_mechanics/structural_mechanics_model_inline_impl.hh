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
inline UInt StructuralMechanicsModel::addMaterial(StructuralMaterial & material,
                                                  const ID & name) {

  const auto material_index = materials.size();

  auto material_name = name;
  if (name.empty()) {
    material_name = "material_" + std::to_string(material_index);
  }

  if (materials_names_to_id.find(material_name) !=
      materials_names_to_id.end()) {
    AKANTU_EXCEPTION("The material " << material_name
                                     << " already exists in the model " << id);
  }

  AKANTU_DEBUG_ASSERT(material_index <=
                          (::std::size_t)::std::numeric_limits<UInt>::max(),
                      "Can not represent the material ID");

  materials_names_to_id[material_name] = material_index;
  materials.push_back(material); // add the material, might cause
                                 // reallocation.

  return UInt(material_index);
}

/* -------------------------------------------------------------------------- */
inline const StructuralMaterial &
StructuralMechanicsModel::getMaterialByElement(const Element & element) const {
  return materials[element_material(element)];
}

/* -------------------------------------------------------------------------- */
inline const StructuralMaterial &
StructuralMechanicsModel::getMaterial(UInt material_index) const {
  return materials.at(material_index);
}

/* -------------------------------------------------------------------------- */
inline const StructuralMaterial &
StructuralMechanicsModel::getMaterial(const ID & name) const {
  auto it = materials_names_to_id.find(name);
  if (it == materials_names_to_id.end()) {
    AKANTU_EXCEPTION("The material " << name << " was not found in the model "
                                     << id);
  }

  return materials.at(it->second);
}

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
