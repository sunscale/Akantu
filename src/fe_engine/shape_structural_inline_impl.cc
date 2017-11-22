/**
 * @file   shape_structural_inline_impl.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  ShapeStructural inline implementation
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
#include "mesh_iterators.hh"
#include "shape_structural.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SHAPE_STRUCTURAL_INLINE_IMPL_CC__
#define __AKANTU_SHAPE_STRUCTURAL_INLINE_IMPL_CC__

namespace akantu {

template <ElementKind kind>
inline void ShapeStructural<kind>::initShapeFunctions(
    const Array<Real> & /* unused */, const Matrix<Real> & /* unused */,
    const ElementType & /* unused */, const GhostType & /* unused */) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
#define INIT_SHAPE_FUNCTIONS(type)                                             \
  setIntegrationPointsByType<type>(integration_points, ghost_type);            \
  precomputeShapesOnIntegrationPoints<type>(nodes, ghost_type);                \
  precomputeShapeDerivativesOnIntegrationPoints<type>(nodes, ghost_type);

template <>
inline void ShapeStructural<_ek_structural>::initShapeFunctions(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    const ElementType & type, const GhostType & ghost_type) {
  AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(INIT_SHAPE_FUNCTIONS);
}

#undef INIT_SHAPE_FUNCTIONS

/* -------------------------------------------------------------------------- */
template <>
template <ElementType type>
void ShapeStructural<_ek_structural>::precomputeShapesOnIntegrationPoints(
    const Array<Real> & /*nodes*/, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  const auto & natural_coords = integration_points(type, ghost_type);
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto spatial_dimension = mesh.getSpatialDimension();
  auto nb_points = integration_points(type, ghost_type).cols();
  auto nb_element = mesh.getNbElement(type, ghost_type);
  auto nb_dof = InterpolationElement<
      ElementClassProperty<type>::interpolation_type>::getNbDegreeOfFreedom();

  auto itp_type = FEEngine::getInterpolationType(type);
  if (not shapes.exists(itp_type, ghost_type)) {
    auto size_of_shapes = this->getShapeSize(type);
    this->shapes.alloc(0, size_of_shapes, itp_type, ghost_type);
  }

  auto & shapes_ = this->shapes(itp_type, ghost_type);
  shapes_.resize(nb_element * nb_points);

  auto shapes_it = shapes_.begin_reinterpret(
      nb_dof, nb_dof * nb_nodes_per_element, nb_points, nb_element);

  for (UInt elem = 0; elem < nb_element; ++elem, ++shapes_it) {
    auto & N = *shapes_it;
    ElementClass<type>::computeShapes(natural_coords, N);
  }

  AKANTU_DEBUG_OUT();
} // namespace akantu

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeStructural<kind>::precomputeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & /*nodes*/, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  const auto & natural_coords = integration_points(type, ghost_type);
  auto spatial_dimension = mesh.getSpatialDimension();
  auto natural_spatial_dimension =
      ElementClass<type>::getNaturalSpaceDimension();
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_points = natural_coords.cols();
  auto nb_dof = ElementClass<type>::getNbDegreeOfFreedom();
  auto nb_element = mesh.getNbElement(type, ghost_type);
  auto nb_stress_components = InterpolationProperty<
      ElementClassProperty<type>::interpolation_type>::nb_stress_components;

  auto itp_type = FEEngine::getInterpolationType(type);
  if (not this->shapes_derivatives.exists(itp_type, ghost_type)) {
    auto size_of_shapesd = this->getShapeDerivativesSize(type);
    this->shapes_derivatives.alloc(0, size_of_shapesd, itp_type, ghost_type);
  }

  auto & shapesd = this->shapes_derivatives(itp_type, ghost_type);
  shapesd.resize(nb_element * nb_points);

  for (auto && tuple :
       zip(make_view(x_el, spatial_dimension, nb_nodes_per_element),
           make_view(shapesd, natural_spatial_dimension,
                     nb_nodes_per_element * nb_dof, nb_points),
           make_view(rot_matrices, spatial_dimension, spatial_dimension))) {
    // compute shape derivatives
    auto & X = std::get<0>(tuple);
    auto & B = std::get<1>(tuple);
    auto & R = std::get<2>(tuple);

    // Rotate to local basis
    auto x = R * X;
    Matrix<Real> node_coords(natural_spatial_dimension, nb_nodes_per_element);

    // Extract relevant first lines
    for (UInt j = 0 ; j < node_coords.rows() ; ++j)
      for (UInt i = 0 ; i < node_coords.cols() ; ++i)
	node_coords(i, j) = x(i, j);

    ElementClass<type>::computeShapeDerivatives(natural_coords, B, node_coords);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeStructural<kind>::interpolateOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_uq, UInt nb_dof,
    const GhostType & ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(out_uq.getNbComponent() == nb_dof,
                      "The output array shape is not correct");

  auto itp_type = FEEngine::getInterpolationType(type);
  const auto & shapes_ = shapes(itp_type, ghost_type);

  auto nb_element = mesh.getNbElement(type, ghost_type);
  auto nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();
  auto nb_quad_points_per_element = integration_points(type, ghost_type).cols();

  Array<Real> u_el(0, nb_nodes_per_element * nb_dof);
  FEEngine::extractNodalToElementField<type>(mesh, in_u, u_el, ghost_type, filter_elements);

  auto nb_quad_points = nb_quad_points_per_element * u_el.size();
  out_uq.resize(nb_quad_points);

  auto out_it = out_uq.begin_reinterpret(nb_dof, 1, nb_quad_points_per_element,
                                         u_el.size());
  auto shapes_it =
      shapes_.begin_reinterpret(nb_dof, nb_dof * nb_nodes_per_element,
                                nb_quad_points_per_element, nb_element);
  auto u_it = u_el.begin_reinterpret(nb_dof * nb_nodes_per_element, 1,
                                     nb_quad_points_per_element, u_el.size());

  for_each_elements(nb_element, filter_elements, [&](auto && el) {
    auto & uq = *out_it;
    const auto & u = *u_it;
    auto N = Tensor3<Real>(shapes_it[el]);

    for (auto && q : arange(uq.size(2))) {
      auto uq_q = Matrix<Real>(uq(q));
      auto u_q = Matrix<Real>(u(q));
      auto N_q = Matrix<Real>(N(q));

      uq_q.mul<false, false>(N, u);
    }

    ++out_it;
    ++u_it;
  });
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeStructural<kind>::gradientOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_nablauq,
    UInt nb_dof, const GhostType & ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  auto itp_type = FEEngine::getInterpolationType(type);
  const auto & shapesd = shapes_derivatives(itp_type, ghost_type);

  auto nb_element = mesh.getNbElement(type, ghost_type);
  auto element_dimension = ElementClass<type>::getSpatialDimension();
  auto nb_quad_points_per_element = integration_points(type, ghost_type).cols();
  auto nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();

  Array<Real> u_el(0, nb_nodes_per_element * nb_dof);
  FEEngine::extractNodalToElementField<type>(mesh, in_u, u_el,
                                             ghost_type, filter_elements);

  auto nb_quad_points = nb_quad_points_per_element * u_el.size();
  out_nablauq.resize(nb_quad_points);

  auto out_it = out_nablauq.begin_reinterpret(element_dimension, 1, nb_quad_points_per_element,
                                              u_el.size());
  auto shapesd_it =
      shapesd.begin_reinterpret(element_dimension, nb_dof * nb_nodes_per_element,
                                nb_quad_points_per_element, nb_element);
  auto u_it = u_el.begin_reinterpret(nb_dof * nb_nodes_per_element, 1,
                                     nb_quad_points_per_element, u_el.size());

  for_each_elements(nb_element, filter_elements, [&](auto && el) {
    auto & nablau = *out_it;
    const auto & u = *u_it;
    auto B = Tensor3<Real>(shapesd_it[el]);

    for (auto && q : arange(nablau.size(2))) {
      auto nablau_q = Matrix<Real>(uq(q));
      auto u_q = Matrix<Real>(u(q));
      auto B_q = Matrix<Real>(N(q));

      nablau_q.mul<false, false>(B, u);
    }

    ++out_it;
    ++u_it;
  });


  AKANTU_DEBUG_OUT();
}

} // namespace akantu

#endif /* __AKANTU_SHAPE_STRUCTURAL_INLINE_IMPL_CC__ */
