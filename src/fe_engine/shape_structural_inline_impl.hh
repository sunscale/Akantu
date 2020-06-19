/**
 * @file   shape_structural_inline_impl.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Oct 11 2017
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  ShapeStructural inline implementation
 *
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
#include "mesh_iterators.hh"
#include "shape_structural.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SHAPE_STRUCTURAL_INLINE_IMPL_HH__
#define __AKANTU_SHAPE_STRUCTURAL_INLINE_IMPL_HH__

namespace akantu {

namespace {
  /// Extract nodal coordinates per elements
  template <ElementType type>
  std::unique_ptr<Array<Real>>
  getNodesPerElement(const Mesh & mesh, const Array<Real> & nodes,
                     const GhostType & ghost_type) {
    const auto dim = ElementClass<type>::getSpatialDimension();
    const auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    auto nodes_per_element =
        std::make_unique<Array<Real>>(0, dim * nb_nodes_per_element);
    FEEngine::extractNodalToElementField(mesh, nodes, *nodes_per_element, type,
                                         ghost_type);
    return nodes_per_element;
  }
} // namespace

template <ElementKind kind>
inline void ShapeStructural<kind>::initShapeFunctions(
    const Array<Real> & /* unused */, const Matrix<Real> & /* unused */,
    const ElementType & /* unused */, const GhostType & /* unused */) {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
#define INIT_SHAPE_FUNCTIONS(type)                                             \
  setIntegrationPointsByType<type>(integration_points, ghost_type);            \
  precomputeRotationMatrices<type>(nodes, ghost_type);                         \
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
void ShapeStructural<_ek_structural>::computeShapesOnIntegrationPoints(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    Array<Real> & shapes, const GhostType & ghost_type,
    const Array<UInt> & filter_elements) const {

  UInt nb_points = integration_points.cols();
  UInt nb_element = mesh.getConnectivity(type, ghost_type).size();

  shapes.resize(nb_element * nb_points);

  UInt ndof = ElementClass<type>::getNbDegreeOfFreedom();

#if !defined(AKANTU_NDEBUG)
  UInt size_of_shapes = ElementClass<type>::getShapeSize();
  AKANTU_DEBUG_ASSERT(shapes.getNbComponent() == size_of_shapes,
                      "The shapes array does not have the correct "
                          << "number of component");
#endif

  auto shapes_it = shapes.begin_reinterpret(
      ElementClass<type>::getNbNodesPerInterpolationElement(), ndof, nb_points,
      nb_element);

  auto shapes_begin = shapes_it;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }

  auto nodes_per_element = getNodesPerElement<type>(mesh, nodes, ghost_type);
  auto nodes_it = nodes_per_element->begin(mesh.getSpatialDimension(),
                                           Mesh::getNbNodesPerElement(type));
  auto nodes_begin = nodes_it;

  for (UInt elem = 0; elem < nb_element; ++elem) {
    if (filter_elements != empty_filter) {
      shapes_it = shapes_begin + filter_elements(elem);
      nodes_it = nodes_begin + filter_elements(elem);
    }

    Tensor3<Real> & N = *shapes_it;
    auto & real_coord = *nodes_it;
    ElementClass<type>::computeShapes(integration_points, real_coord, N);

    if (filter_elements == empty_filter) {
      ++shapes_it;
      ++nodes_it;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeStructural<kind>::precomputeRotationMatrices(
    const Array<Real> & nodes, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  const auto spatial_dimension = mesh.getSpatialDimension();
  const auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  const auto nb_element = mesh.getNbElement(type, ghost_type);
  const auto nb_dof = ElementClass<type>::getNbDegreeOfFreedom();

  if (not this->rotation_matrices.exists(type, ghost_type)) {
    this->rotation_matrices.alloc(0, nb_dof * nb_dof, type, ghost_type);
  }

  auto & rot_matrices = this->rotation_matrices(type, ghost_type);
  rot_matrices.resize(nb_element);

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type);

  bool has_extra_normal = mesh.hasData<Real>("extra_normal", type, ghost_type);
  Array<Real>::const_vector_iterator extra_normal;
  if (has_extra_normal)
    extra_normal = mesh.getData<Real>("extra_normal", type, ghost_type)
                       .begin(spatial_dimension);

  for (auto && tuple :
       zip(make_view(x_el, spatial_dimension, nb_nodes_per_element),
           make_view(rot_matrices, nb_dof, nb_dof))) {
    // compute shape derivatives
    auto & X = std::get<0>(tuple);
    auto & R = std::get<1>(tuple);

    if (has_extra_normal) {
      ElementClass<type>::computeRotationMatrix(R, X, *extra_normal);
      ++extra_normal;
    } else {
      ElementClass<type>::computeRotationMatrix(
          R, X, Vector<Real>(spatial_dimension));
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <ElementKind kind>
template <ElementType type>
void ShapeStructural<kind>::precomputeShapesOnIntegrationPoints(
    const Array<Real> & nodes, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  const auto & natural_coords = integration_points(type, ghost_type);
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_points = integration_points(type, ghost_type).cols();
  auto nb_element = mesh.getNbElement(type, ghost_type);
  auto nb_dof = ElementClass<type>::getNbDegreeOfFreedom();
  const auto dim = ElementClass<type>::getSpatialDimension();

  auto itp_type = FEEngine::getInterpolationType(type);
  if (not shapes.exists(itp_type, ghost_type)) {
    auto size_of_shapes = this->getShapeSize(type);
    this->shapes.alloc(0, size_of_shapes, itp_type, ghost_type);
  }

  auto & shapes_ = this->shapes(itp_type, ghost_type);
  shapes_.resize(nb_element * nb_points);

  auto nodes_per_element = getNodesPerElement<type>(mesh, nodes, ghost_type);

  for (auto && tuple :
       zip(make_view(shapes_, nb_dof, nb_dof * nb_nodes_per_element, nb_points),
           make_view(*nodes_per_element, dim, nb_nodes_per_element))) {
    auto & N = std::get<0>(tuple);
    auto & real_coord = std::get<1>(tuple);
    ElementClass<type>::computeShapes(natural_coords, real_coord, N);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeStructural<kind>::precomputeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  const auto & natural_coords = integration_points(type, ghost_type);
  const auto spatial_dimension = mesh.getSpatialDimension();
  const auto natural_spatial_dimension =
      ElementClass<type>::getNaturalSpaceDimension();
  const auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  const auto nb_points = natural_coords.cols();
  const auto nb_dof = ElementClass<type>::getNbDegreeOfFreedom();
  const auto nb_element = mesh.getNbElement(type, ghost_type);
  const auto nb_stress_components = ElementClass<type>::getNbStressComponents();

  auto itp_type = FEEngine::getInterpolationType(type);
  if (not this->shapes_derivatives.exists(itp_type, ghost_type)) {
    auto size_of_shapesd = this->getShapeDerivativesSize(type);
    this->shapes_derivatives.alloc(0, size_of_shapesd, itp_type, ghost_type);
  }

  auto & rot_matrices = this->rotation_matrices(type, ghost_type);

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type);

  auto & shapesd = this->shapes_derivatives(itp_type, ghost_type);
  shapesd.resize(nb_element * nb_points);

  for (auto && tuple :
       zip(make_view(x_el, spatial_dimension, nb_nodes_per_element),
           make_view(shapesd, nb_stress_components,
                     nb_nodes_per_element * nb_dof, nb_points),
           make_view(rot_matrices, nb_dof, nb_dof))) {
    // compute shape derivatives
    auto & X = std::get<0>(tuple);
    auto & B = std::get<1>(tuple);
    auto & RDOFs = std::get<2>(tuple);

    Tensor3<Real> dnds(natural_spatial_dimension,
                       ElementClass<type>::interpolation_property::dnds_columns,
                       B.size(2));
    ElementClass<type>::computeDNDS(natural_coords, X, dnds);

    Tensor3<Real> J(natural_spatial_dimension, natural_spatial_dimension,
                    natural_coords.cols());

    // Computing the coordinates of the element in the natural space
    auto R = RDOFs.block(0, 0, spatial_dimension, spatial_dimension);
    Matrix<Real> T(B.size(1), B.size(1), 0);

    for (UInt i = 0; i < nb_nodes_per_element; ++i) {
      T.block(RDOFs, i * RDOFs.rows(), i * RDOFs.rows());
    }

    // Rotate to local basis
    auto x =
        (R * X).block(0, 0, natural_spatial_dimension, nb_nodes_per_element);

    ElementClass<type>::computeJMat(natural_coords, x, J);
    ElementClass<type>::computeShapeDerivatives(J, dnds, T, B);
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
  FEEngine::extractNodalToElementField(mesh, in_u, u_el, type, ghost_type,
                                       filter_elements);

  auto nb_quad_points = nb_quad_points_per_element * u_el.size();
  out_uq.resize(nb_quad_points);

  auto out_it = out_uq.begin_reinterpret(nb_dof, 1, nb_quad_points_per_element,
                                         u_el.size());
  auto shapes_it =
      shapes_.begin_reinterpret(nb_dof, nb_dof * nb_nodes_per_element,
                                nb_quad_points_per_element, nb_element);
  auto u_it = u_el.begin_reinterpret(nb_dof * nb_nodes_per_element, 1,
                                     nb_quad_points_per_element, u_el.size());

  for_each_element(nb_element, filter_elements, [&](auto && el) {
    auto & uq = *out_it;
    const auto & u = *u_it;
    auto N = Tensor3<Real>(shapes_it[el]);

    for (auto && q : arange(uq.size(2))) {
      auto uq_q = Matrix<Real>(uq(q));
      auto u_q = Matrix<Real>(u(q));
      auto N_q = Matrix<Real>(N(q));

      uq_q.mul<false, false>(N_q, u_q);
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
    const Array<Real> & in_u, Array<Real> & out_nablauq, UInt nb_dof,
    const GhostType & ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  auto itp_type = FEEngine::getInterpolationType(type);
  const auto & shapesd = shapes_derivatives(itp_type, ghost_type);

  auto nb_element = mesh.getNbElement(type, ghost_type);
  auto element_dimension = ElementClass<type>::getSpatialDimension();
  auto nb_quad_points_per_element = integration_points(type, ghost_type).cols();
  auto nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();

  Array<Real> u_el(0, nb_nodes_per_element * nb_dof);
  FEEngine::extractNodalToElementField(mesh, in_u, u_el, type, ghost_type,
                                       filter_elements);

  auto nb_quad_points = nb_quad_points_per_element * u_el.size();
  out_nablauq.resize(nb_quad_points);

  auto out_it = out_nablauq.begin_reinterpret(
      element_dimension, 1, nb_quad_points_per_element, u_el.size());
  auto shapesd_it = shapesd.begin_reinterpret(
      element_dimension, nb_dof * nb_nodes_per_element,
      nb_quad_points_per_element, nb_element);
  auto u_it = u_el.begin_reinterpret(nb_dof * nb_nodes_per_element, 1,
                                     nb_quad_points_per_element, u_el.size());

  for_each_element(nb_element, filter_elements, [&](auto && el) {
    auto & nablau = *out_it;
    const auto & u = *u_it;
    auto B = Tensor3<Real>(shapesd_it[el]);

    for (auto && q : arange(nablau.size(2))) {
      auto nablau_q = Matrix<Real>(nablau(q));
      auto u_q = Matrix<Real>(u(q));
      auto B_q = Matrix<Real>(B(q));

      nablau_q.mul<false, false>(B_q, u_q);
    }

    ++out_it;
    ++u_it;
  });

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <>
template <ElementType type>
void ShapeStructural<_ek_structural>::computeBtD(
    const Array<Real> & Ds, Array<Real> & BtDs, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  auto itp_type = ElementClassProperty<type>::interpolation_type;

  auto nb_stress = ElementClass<type>::getNbStressComponents();
  auto nb_dof_per_element = ElementClass<type>::getNbDegreeOfFreedom() *
                            mesh.getNbNodesPerElement(type);

  const auto & shapes_derivatives =
      this->shapes_derivatives(itp_type, ghost_type);

  Array<Real> shapes_derivatives_filtered(0,
                                          shapes_derivatives.getNbComponent());
  auto && view = make_view(shapes_derivatives, nb_stress, nb_dof_per_element);
  auto B_it = view.begin();
  auto B_end = view.end();

  if (filter_elements != empty_filter) {
    FEEngine::filterElementalData(this->mesh, shapes_derivatives,
                                  shapes_derivatives_filtered, type, ghost_type,
                                  filter_elements);
    auto && view =
        make_view(shapes_derivatives_filtered, nb_stress, nb_dof_per_element);
    B_it = view.begin();
    B_end = view.end();
  }

  for (auto && values : zip(range(B_it, B_end), make_view(Ds, nb_stress),
                            make_view(BtDs, BtDs.getNbComponent()))) {
    const auto & B = std::get<0>(values);
    const auto & D = std::get<1>(values);
    auto & Bt_D = std::get<2>(values);
    Bt_D.template mul<true>(B, D);
  }
}

} // namespace akantu

#endif /* __AKANTU_SHAPE_STRUCTURAL_INLINE_IMPL_HH__ */
