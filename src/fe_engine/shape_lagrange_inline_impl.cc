/**
 * @file   shape_lagrange_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Oct 27 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  ShapeLagrange inline implementation
 *
 * @section LICENSE
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
#include "aka_iterators.hh"
#include "aka_voigthelper.hh"
#include "fe_engine.hh"
#include "shape_lagrange.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SHAPE_LAGRANGE_INLINE_IMPL_CC__
#define __AKANTU_SHAPE_LAGRANGE_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
#define INIT_SHAPE_FUNCTIONS(type)                                             \
  setIntegrationPointsByType<type>(integration_points, ghost_type);            \
  precomputeShapesOnIntegrationPoints<type>(nodes, ghost_type);                \
  if (ElementClass<type>::getNaturalSpaceDimension() ==                        \
          mesh.getSpatialDimension() ||                                        \
      kind != _ek_regular)                                                     \
    precomputeShapeDerivativesOnIntegrationPoints<type>(nodes, ghost_type);

template <ElementKind kind>
inline void ShapeLagrange<kind>::initShapeFunctions(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    const ElementType & type, const GhostType & ghost_type) {
  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(INIT_SHAPE_FUNCTIONS);
}

#undef INIT_SHAPE_FUNCTIONS

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
inline void ShapeLagrange<kind>::computeShapeDerivativesOnCPointsByElement(
    const Matrix<Real> & node_coords, const Matrix<Real> & natural_coords,
    Tensor3<Real> & shapesd) const {
  AKANTU_DEBUG_IN();

  // compute dnds
  Tensor3<Real> dnds(node_coords.rows(), node_coords.cols(),
                     natural_coords.cols());
  ElementClass<type>::computeDNDS(natural_coords, dnds);
  // compute jacobian
  Tensor3<Real> J(node_coords.rows(), natural_coords.rows(),
                  natural_coords.cols());
  ElementClass<type>::computeJMat(dnds, node_coords, J);

  // compute dndx
  ElementClass<type>::computeShapeDerivatives(J, dnds, shapesd);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::inverseMap(const Vector<Real> & real_coords,
                                     UInt elem, Vector<Real> & natural_coords,
                                     const GhostType & ghost_type) const {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  UInt * elem_val = mesh.getConnectivity(type, ghost_type).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(mesh.getNodes(), nodes_coord.storage(),
                                     elem_val + elem * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

  ElementClass<type>::inverseMap(real_coords, nodes_coord, natural_coords);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
bool ShapeLagrange<kind>::contains(const Vector<Real> & real_coords, UInt elem,
                                   const GhostType & ghost_type) const {

  UInt spatial_dimension = mesh.getSpatialDimension();
  Vector<Real> natural_coords(spatial_dimension);

  inverseMap<type>(real_coords, elem, natural_coords, ghost_type);
  return ElementClass<type>::contains(natural_coords);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::interpolate(const Vector<Real> & real_coords,
                                      UInt elem,
                                      const Matrix<Real> & nodal_values,
                                      Vector<Real> & interpolated,
                                      const GhostType & ghost_type) const {
  UInt nb_shapes = ElementClass<type>::getShapeSize();
  Vector<Real> shapes(nb_shapes);
  computeShapes<type>(real_coords, elem, shapes, ghost_type);
  ElementClass<type>::interpolate(nodal_values, shapes, interpolated);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::computeShapes(const Vector<Real> & real_coords,
                                        UInt elem, Vector<Real> & shapes,
                                        const GhostType & ghost_type) const {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  Vector<Real> natural_coords(spatial_dimension);

  inverseMap<type>(real_coords, elem, natural_coords, ghost_type);
  ElementClass<type>::computeShapes(natural_coords, shapes);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::computeShapeDerivatives(
    const Matrix<Real> & real_coords, UInt elem, Tensor3<Real> & shapesd,
    const GhostType & ghost_type) const {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_points = real_coords.cols();
  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  AKANTU_DEBUG_ASSERT(mesh.getSpatialDimension() == shapesd.size(0) &&
                          nb_nodes_per_element == shapesd.size(1),
                      "Shape size doesn't match");
  AKANTU_DEBUG_ASSERT(nb_points == shapesd.size(2),
                      "Number of points doesn't match shapes size");

  Matrix<Real> natural_coords(spatial_dimension, nb_points);

  // Creates the matrix of natural coordinates
  for (UInt i = 0; i < nb_points; i++) {
    Vector<Real> real_point = real_coords(i);
    Vector<Real> natural_point = natural_coords(i);

    inverseMap<type>(real_point, elem, natural_point, ghost_type);
  }

  UInt * elem_val = mesh.getConnectivity(type, ghost_type).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(mesh.getNodes(), nodes_coord.storage(),
                                     elem_val + elem * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

  computeShapeDerivativesOnCPointsByElement<type>(nodes_coord, natural_coords,
                                                  shapesd);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
ShapeLagrange<kind>::ShapeLagrange(const Mesh & mesh, const ID & id,
                                   const MemoryID & memory_id)
    : ShapeLagrangeBase(mesh, kind, id, memory_id) {}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::computeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    Array<Real> & shape_derivatives, const GhostType & ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  UInt nb_points = integration_points.cols();
  UInt nb_element = mesh.getConnectivity(type, ghost_type).size();

  UInt size_of_shapesd = ElementClass<type>::getShapeDerivativesSize();
  AKANTU_DEBUG_ASSERT(shape_derivatives.getNbComponent() == size_of_shapesd,
                      "The shapes_derivatives array does not have the correct "
                          << "number of component");
  shape_derivatives.resize(nb_element * nb_points);

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type,
                                       filter_elements);

  Real * shapesd_val = shape_derivatives.storage();
  Array<Real>::matrix_iterator x_it =
      x_el.begin(spatial_dimension, nb_nodes_per_element);

  if (filter_elements != empty_filter)
    nb_element = filter_elements.size();

  for (UInt elem = 0; elem < nb_element; ++elem, ++x_it) {
    if (filter_elements != empty_filter)
      shapesd_val = shape_derivatives.storage() +
                    filter_elements(elem) * size_of_shapesd * nb_points;

    Matrix<Real> & X = *x_it;
    Tensor3<Real> B(shapesd_val, spatial_dimension, nb_nodes_per_element,
                    nb_points);
    computeShapeDerivativesOnCPointsByElement<type>(X, integration_points, B);

    if (filter_elements == empty_filter)
      shapesd_val += size_of_shapesd * nb_points;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
void ShapeLagrange<kind>::computeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    Array<Real> & shape_derivatives, const ElementType & type,
    const GhostType & ghost_type, const Array<UInt> & filter_elements) const {
#define AKANTU_COMPUTE_SHAPES(type)                                            \
  computeShapeDerivativesOnIntegrationPoints<type>(                            \
      nodes, integration_points, shape_derivatives, ghost_type,                \
      filter_elements);

  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(AKANTU_COMPUTE_SHAPES);

#undef AKANTU_COMPUTE_SHAPES
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::precomputeShapesOnIntegrationPoints(
    const Array<Real> & nodes, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  Matrix<Real> & natural_coords = integration_points(type, ghost_type);
  UInt size_of_shapes = ElementClass<type>::getShapeSize();

  Array<Real> & shapes_tmp =
      shapes.alloc(0, size_of_shapes, itp_type, ghost_type);

  this->computeShapesOnIntegrationPoints<type>(nodes, natural_coords,
                                               shapes_tmp, ghost_type);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::precomputeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  Matrix<Real> & natural_coords = integration_points(type, ghost_type);
  UInt size_of_shapesd = ElementClass<type>::getShapeDerivativesSize();

  Array<Real> & shapes_derivatives_tmp =
      shapes_derivatives.alloc(0, size_of_shapesd, itp_type, ghost_type);

  this->computeShapeDerivativesOnIntegrationPoints<type>(
      nodes, natural_coords, shapes_derivatives_tmp, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::interpolateOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_uq, UInt nb_degree_of_freedom,
    const Array<Real> & shapes, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, in_u, u_el, type, ghost_type,
                                       filter_elements);

  this->interpolateElementalFieldOnIntegrationPoints<type>(
      u_el, out_uq, ghost_type, shapes, filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::interpolateOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_uq, UInt nb_degree_of_freedom,
    GhostType ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  AKANTU_DEBUG_ASSERT(shapes.exists(itp_type, ghost_type),
                      "No shapes for the type "
                          << shapes.printType(itp_type, ghost_type));

  this->interpolateOnIntegrationPoints<type>(in_u, out_uq, nb_degree_of_freedom,
                                             shapes(itp_type, ghost_type),
                                             ghost_type, filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::gradientOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_nablauq,
    UInt nb_degree_of_freedom, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  AKANTU_DEBUG_ASSERT(
      shapes_derivatives.exists(itp_type, ghost_type),
      "No shapes derivatives for the type "
          << shapes_derivatives.printType(itp_type, ghost_type));

  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, in_u, u_el, type, ghost_type,
                                       filter_elements);

  this->gradientElementalFieldOnIntegrationPoints<type>(
      u_el, out_nablauq, ghost_type, shapes_derivatives(itp_type, ghost_type),
      filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::computeBtD(
    const Array<Real> & Ds, Array<Real> & BtDs, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  auto itp_type = ElementClassProperty<type>::interpolation_type;
  const auto & shapes_derivatives =
      this->shapes_derivatives(itp_type, ghost_type);

  auto spatial_dimension = mesh.getSpatialDimension();
  auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);

  Array<Real> shapes_derivatives_filtered(0,
                                          shapes_derivatives.getNbComponent());
  auto && view =
      make_view(shapes_derivatives, spatial_dimension, nb_nodes_per_element);
  auto B_it = view.begin();
  auto B_end = view.end();

  if (filter_elements != empty_filter) {
    FEEngine::filterElementalData(this->mesh, shapes_derivatives,
                                  shapes_derivatives_filtered, type, ghost_type,
                                  filter_elements);
    auto && view = make_view(shapes_derivatives_filtered, spatial_dimension,
                             nb_nodes_per_element);
    B_it = view.begin();
    B_end = view.end();
  }

  for (auto && values :
       zip(range(B_it, B_end),
           make_view(Ds, Ds.getNbComponent() / spatial_dimension,
                     spatial_dimension),
           make_view(BtDs, BtDs.getNbComponent() / nb_nodes_per_element,
                     nb_nodes_per_element))) {
    const auto & B = std::get<0>(values);
    const auto & D = std::get<1>(values);
    auto & Bt_D = std::get<2>(values);
    // transposed due to the storage layout of B
    Bt_D.template mul<false, false>(D, B);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::computeBtDB(
    const Array<Real> & Ds, Array<Real> & BtDBs, UInt order_d,
    GhostType ghost_type, const Array<UInt> & filter_elements) const {
  auto itp_type = ElementClassProperty<type>::interpolation_type;
  const auto & shapes_derivatives =
      this->shapes_derivatives(itp_type, ghost_type);

  constexpr auto dim = ElementClass<type>::getSpatialDimension();
  auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);

  Array<Real> shapes_derivatives_filtered(0,
                                          shapes_derivatives.getNbComponent());
  auto && view = make_view(shapes_derivatives, dim, nb_nodes_per_element);
  auto B_it = view.begin();
  auto B_end = view.end();

  if (filter_elements != empty_filter) {
    FEEngine::filterElementalData(this->mesh, shapes_derivatives,
                                  shapes_derivatives_filtered, type, ghost_type,
                                  filter_elements);
    auto && view =
        make_view(shapes_derivatives_filtered, dim, nb_nodes_per_element);
    B_it = view.begin();
    B_end = view.end();
  }

  if (order_d == 4) {
    UInt tangent_size = VoigtHelper<dim>::size;
    Matrix<Real> B(tangent_size, dim * nb_nodes_per_element);
    Matrix<Real> Bt_D(dim * nb_nodes_per_element, tangent_size);

    for (auto && values :
         zip(range(B_it, B_end), make_view(Ds, tangent_size, tangent_size),
             make_view(BtDBs, dim * nb_nodes_per_element,
                       dim * nb_nodes_per_element))) {
      const auto & Bfull = std::get<0>(values);
      const auto & D = std::get<1>(values);
      auto & Bt_D_B = std::get<2>(values);

      VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(Bfull, B,
                                                         nb_nodes_per_element);
      Bt_D.template mul<true, false>(B, D);
      Bt_D_B.template mul<false, false>(Bt_D, B);
    }
  } else if (order_d == 2) {
    Matrix<Real> Bt_D(nb_nodes_per_element, dim);
    for (auto && values :
         zip(range(B_it, B_end), make_view(Ds, dim, dim),
             make_view(BtDBs, nb_nodes_per_element, nb_nodes_per_element))) {
      const auto & B = std::get<0>(values);
      const auto & D = std::get<1>(values);
      auto & Bt_D_B = std::get<2>(values);
      Bt_D.template mul<true, false>(B, D);
      Bt_D_B.template mul<false, false>(Bt_D, B);
    }
  }
}

template <>
template <>
inline void ShapeLagrange<_ek_regular>::computeBtDB<_point_1>(
    const Array<Real> & /*Ds*/, Array<Real> & /*BtDBs*/, UInt /*order_d*/,
    GhostType /*ghost_type*/, const Array<UInt> & /*filter_elements*/) const {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::computeNtb(
    const Array<Real> & bs, Array<Real> & Ntbs, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  Ntbs.resize(bs.size());

  UInt size_of_shapes = ElementClass<type>::getShapeSize();
  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  UInt nb_degree_of_freedom = bs.getNbComponent();

  Array<Real> shapes_filtered(0, size_of_shapes);
  auto && view = make_view(shapes(itp_type, ghost_type), 1, size_of_shapes);
  auto N_it = view.begin();
  auto N_end = view.end();

  if (filter_elements != empty_filter) {
    FEEngine::filterElementalData(this->mesh, shapes(itp_type, ghost_type),
                                  shapes_filtered, type, ghost_type,
                                  filter_elements);
    auto && view = make_view(shapes_filtered, 1, size_of_shapes);
    N_it = view.begin();
    N_end = view.end();
  }

  for (auto && values :
       zip(make_view(bs, nb_degree_of_freedom, 1), range(N_it, N_end),
           make_view(Ntbs, nb_degree_of_freedom, size_of_shapes))) {
    const auto & b = std::get<0>(values);
    const auto & N = std::get<1>(values);
    auto & Ntb = std::get<2>(values);

    Ntb.template mul<false, false>(b, N);
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

#endif /* __AKANTU_SHAPE_LAGRANGE_INLINE_IMPL_CC__ */
