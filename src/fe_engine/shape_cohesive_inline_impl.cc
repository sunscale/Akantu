/**
 * @file   shape_cohesive_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Feb 03 2012
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  ShapeCohesive inline implementation
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
#include "shape_cohesive.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SHAPE_COHESIVE_INLINE_IMPL_CC__
#define __AKANTU_SHAPE_COHESIVE_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline ShapeLagrange<_ek_cohesive>::ShapeLagrange(const Mesh & mesh,
                                                  const ID & id,
                                                  const MemoryID & memory_id)
    : ShapeLagrangeBase(mesh, _ek_cohesive, id, memory_id) {}

#define INIT_SHAPE_FUNCTIONS(type)                                             \
  setIntegrationPointsByType<type>(integration_points, ghost_type);            \
  precomputeShapesOnIntegrationPoints<type>(nodes, ghost_type);                \
  precomputeShapeDerivativesOnIntegrationPoints<type>(nodes, ghost_type);

/* -------------------------------------------------------------------------- */
inline void ShapeLagrange<_ek_cohesive>::initShapeFunctions(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    const ElementType & type, const GhostType & ghost_type) {
  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(INIT_SHAPE_FUNCTIONS);
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_cohesive>::computeShapeDerivativesOnIntegrationPoints(
    const Array<Real> &, const Matrix<Real> & integration_points,
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

  Real * shapesd_val = shape_derivatives.storage();

  if (filter_elements == empty_filter)
    nb_element = filter_elements.size();

  for (UInt elem = 0; elem < nb_element; ++elem) {
    Tensor3<Real> B(shapesd_val, spatial_dimension, nb_nodes_per_element,
                    nb_points);
    ElementClass<type>::computeDNDS(integration_points, B);

    if (filter_elements != empty_filter)
      shapesd_val += size_of_shapesd * nb_points;
    else {
      shapesd_val = shape_derivatives.storage() +
                    filter_elements(elem) * size_of_shapesd * nb_points;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void
ShapeLagrange<_ek_cohesive>::computeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    Array<Real> & shape_derivatives, const ElementType & type,
    const GhostType & ghost_type, const Array<UInt> & filter_elements) const {
#define AKANTU_COMPUTE_SHAPES(type)                                            \
  computeShapeDerivativesOnIntegrationPoints<type>(                            \
      nodes, integration_points, shape_derivatives, ghost_type,                \
      filter_elements);

  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(AKANTU_COMPUTE_SHAPES);

#undef AKANTU_COMPUTE_SHAPES
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_cohesive>::precomputeShapesOnIntegrationPoints(
    const Array<Real> & nodes, GhostType ghost_type) {
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
template <ElementType type>
void ShapeLagrange<_ek_cohesive>::precomputeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, GhostType ghost_type) {
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
template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::extractNodalToElementField(
    const Array<Real> & nodal_f, Array<Real> & elemental_f,
    const GhostType & ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_itp_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();
  UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();
  UInt nb_degree_of_freedom = nodal_f.getNbComponent();
  UInt nb_element = this->mesh.getNbElement(type, ghost_type);
  UInt * conn_val = this->mesh.getConnectivity(type, ghost_type).storage();

  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }

  elemental_f.resize(nb_element);

  Array<Real>::matrix_iterator u_it =
      elemental_f.begin(nb_degree_of_freedom, nb_nodes_per_itp_element);

  UInt * el_conn;
  ReduceFunction reduce_function;

  for (UInt el = 0; el < nb_element; ++el, ++u_it) {
    Matrix<Real> & u = *u_it;
    if (filter_elements != empty_filter)
      el_conn = conn_val + filter_elements(el) * nb_nodes_per_element;
    else
      el_conn = conn_val + el * nb_nodes_per_element;

    // compute the average/difference of the nodal field loaded from cohesive
    // element
    for (UInt n = 0; n < nb_nodes_per_itp_element; ++n) {
      UInt node_plus = *(el_conn + n);
      UInt node_minus = *(el_conn + n + nb_nodes_per_itp_element);
      for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
        Real u_plus = nodal_f(node_plus, d);
        Real u_minus = nodal_f(node_minus, d);
        u(d, n) = reduce_function(u_plus, u_minus);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::interpolateOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_uq, UInt nb_degree_of_freedom,
    GhostType ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

  AKANTU_DEBUG_ASSERT(this->shapes.exists(itp_type, ghost_type),
                      "No shapes for the type "
                          << this->shapes.printType(itp_type, ghost_type));

  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();
  Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  this->extractNodalToElementField<type, ReduceFunction>(in_u, u_el, ghost_type,
                                                         filter_elements);

  this->template interpolateElementalFieldOnIntegrationPoints<type>(
      u_el, out_uq, ghost_type, shapes(itp_type, ghost_type), filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::variationOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & nablauq, UInt nb_degree_of_freedom,
    GhostType ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

  AKANTU_DEBUG_ASSERT(
      this->shapes_derivatives.exists(itp_type, ghost_type),
      "No shapes for the type "
          << this->shapes_derivatives.printType(itp_type, ghost_type));

  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();
  Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  this->extractNodalToElementField<type, ReduceFunction>(in_u, u_el, ghost_type,
                                                         filter_elements);

  this->template gradientElementalFieldOnIntegrationPoints<type>(
      u_el, nablauq, ghost_type, shapes_derivatives(itp_type, ghost_type),
      filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::computeNormalsOnIntegrationPoints(
    const Array<Real> & u, Array<Real> & normals_u, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_element = this->mesh.getNbElement(type, ghost_type);
  UInt nb_points = this->integration_points(type, ghost_type).cols();
  UInt spatial_dimension = this->mesh.getSpatialDimension();

  if (filter_elements != empty_filter)
    nb_element = filter_elements.size();

  normals_u.resize(nb_points * nb_element);

  Array<Real> tangents_u(nb_element * nb_points,
                         (spatial_dimension * (spatial_dimension - 1)));

  this->template variationOnIntegrationPoints<type, ReduceFunction>(
      u, tangents_u, spatial_dimension, ghost_type, filter_elements);

  Array<Real>::vector_iterator normal = normals_u.begin(spatial_dimension);
  Array<Real>::vector_iterator normal_end = normals_u.end(spatial_dimension);

  Real * tangent = tangents_u.storage();

  if (spatial_dimension == 3)
    for (; normal != normal_end; ++normal) {
      Math::vectorProduct3(tangent, tangent + spatial_dimension,
                           normal->storage());

      (*normal) /= normal->norm();
      tangent += spatial_dimension * 2;
    }
  else if (spatial_dimension == 2)
    for (; normal != normal_end; ++normal) {
      Vector<Real> a1(tangent, spatial_dimension);

      (*normal)(0) = -a1(1);
      (*normal)(1) = a1(0);
      (*normal) /= normal->norm();

      tangent += spatial_dimension;
    }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
#endif /* __AKANTU_SHAPE_COHESIVE_INLINE_IMPL_CC__ */
