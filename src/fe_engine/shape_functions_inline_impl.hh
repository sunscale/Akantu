/**
 * @file   shape_functions_inline_impl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Oct 27 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  ShapeFunctions inline implementation
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
#include "fe_engine.hh"
#include "shape_functions.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SHAPE_FUNCTIONS_INLINE_IMPL_HH__
#define __AKANTU_SHAPE_FUNCTIONS_INLINE_IMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline const Array<Real> &
ShapeFunctions::getShapes(const ElementType & el_type,
                          const GhostType & ghost_type) const {
  return shapes(FEEngine::getInterpolationType(el_type), ghost_type);
}

/* -------------------------------------------------------------------------- */
inline const Array<Real> &
ShapeFunctions::getShapesDerivatives(const ElementType & el_type,
                                     const GhostType & ghost_type) const {
  return shapes_derivatives(FEEngine::getInterpolationType(el_type),
                            ghost_type);
}

/* -------------------------------------------------------------------------- */
inline UInt ShapeFunctions::getShapeSize(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt shape_size = 0;
#define GET_SHAPE_SIZE(type) shape_size = ElementClass<type>::getShapeSize()

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_SIZE); // ,

#undef GET_SHAPE_SIZE

  AKANTU_DEBUG_OUT();
  return shape_size;
}

/* -------------------------------------------------------------------------- */
inline UInt ShapeFunctions::getShapeDerivativesSize(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt shape_derivatives_size = 0;
#define GET_SHAPE_DERIVATIVES_SIZE(type)                                       \
  shape_derivatives_size = ElementClass<type>::getShapeDerivativesSize()

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_SIZE); // ,

#undef GET_SHAPE_DERIVATIVES_SIZE

  AKANTU_DEBUG_OUT();
  return shape_derivatives_size;
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeFunctions::setIntegrationPointsByType(const Matrix<Real> & points,
                                                const GhostType & ghost_type) {
  if (not this->integration_points.exists(type, ghost_type))
    this->integration_points(type, ghost_type).shallowCopy(points);
}

/* -------------------------------------------------------------------------- */
inline void
ShapeFunctions::buildInterpolationMatrix(const Matrix<Real> & coordinates,
                                         Matrix<Real> & coordMatrix,
                                         UInt integration_order) const {
  switch (integration_order) {
  case 1: {
    for (UInt i = 0; i < coordinates.cols(); ++i)
      coordMatrix(i, 0) = 1;
    break;
  }
  case 2: {
    UInt nb_quadrature_points = coordMatrix.cols();

    for (UInt i = 0; i < coordinates.cols(); ++i) {
      coordMatrix(i, 0) = 1;
      for (UInt j = 1; j < nb_quadrature_points; ++j)
        coordMatrix(i, j) = coordinates(j - 1, i);
    }
    break;
  }
  default: {
    AKANTU_TO_IMPLEMENT();
    break;
  }
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ShapeFunctions::buildElementalFieldInterpolationMatrix(
    const Matrix<Real> &, Matrix<Real> &, UInt) const {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
inline void ShapeFunctions::buildElementalFieldInterpolationMatrix<_segment_2>(
    const Matrix<Real> & coordinates, Matrix<Real> & coordMatrix,
    UInt integration_order) const {
  buildInterpolationMatrix(coordinates, coordMatrix, integration_order);
}

/* -------------------------------------------------------------------------- */
template <>
inline void ShapeFunctions::buildElementalFieldInterpolationMatrix<_segment_3>(
    const Matrix<Real> & coordinates, Matrix<Real> & coordMatrix,
    UInt integration_order) const {
  buildInterpolationMatrix(coordinates, coordMatrix, integration_order);
}

/* -------------------------------------------------------------------------- */
template <>
inline void ShapeFunctions::buildElementalFieldInterpolationMatrix<_triangle_3>(
    const Matrix<Real> & coordinates, Matrix<Real> & coordMatrix,
    UInt integration_order) const {
  buildInterpolationMatrix(coordinates, coordMatrix, integration_order);
}

/* -------------------------------------------------------------------------- */
template <>
inline void ShapeFunctions::buildElementalFieldInterpolationMatrix<_triangle_6>(
    const Matrix<Real> & coordinates, Matrix<Real> & coordMatrix,
    UInt integration_order) const {
  buildInterpolationMatrix(coordinates, coordMatrix, integration_order);
}

/* -------------------------------------------------------------------------- */
template <>
inline void
ShapeFunctions::buildElementalFieldInterpolationMatrix<_tetrahedron_4>(
    const Matrix<Real> & coordinates, Matrix<Real> & coordMatrix,
    UInt integration_order) const {
  buildInterpolationMatrix(coordinates, coordMatrix, integration_order);
}

/* -------------------------------------------------------------------------- */
template <>
inline void
ShapeFunctions::buildElementalFieldInterpolationMatrix<_tetrahedron_10>(
    const Matrix<Real> & coordinates, Matrix<Real> & coordMatrix,
    UInt integration_order) const {
  buildInterpolationMatrix(coordinates, coordMatrix, integration_order);
}

/**
 * @todo Write a more efficient interpolation for quadrangles by
 * dropping unnecessary quadrature points
 *
 */

/* -------------------------------------------------------------------------- */
template <>
inline void
ShapeFunctions::buildElementalFieldInterpolationMatrix<_quadrangle_4>(
    const Matrix<Real> & coordinates, Matrix<Real> & coordMatrix,
    UInt integration_order) const {

  if (integration_order !=
      ElementClassProperty<_quadrangle_4>::polynomial_degree) {
    AKANTU_TO_IMPLEMENT();
  } else {
    for (UInt i = 0; i < coordinates.cols(); ++i) {
      Real x = coordinates(0, i);
      Real y = coordinates(1, i);

      coordMatrix(i, 0) = 1;
      coordMatrix(i, 1) = x;
      coordMatrix(i, 2) = y;
      coordMatrix(i, 3) = x * y;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <>
inline void
ShapeFunctions::buildElementalFieldInterpolationMatrix<_quadrangle_8>(
    const Matrix<Real> & coordinates, Matrix<Real> & coordMatrix,
    UInt integration_order) const {

  if (integration_order !=
      ElementClassProperty<_quadrangle_8>::polynomial_degree) {
    AKANTU_TO_IMPLEMENT();
  } else {
    for (UInt i = 0; i < coordinates.cols(); ++i) {
      // UInt j = 0;
      Real x = coordinates(0, i);
      Real y = coordinates(1, i);

      coordMatrix(i, 0) = 1;
      coordMatrix(i, 1) = x;
      coordMatrix(i, 2) = y;
      coordMatrix(i, 3) = x * y;
      // for (UInt e = 0; e <= 2; ++e) {
      //   for (UInt n = 0; n <= 2; ++n) {
      //     coordMatrix(i, j) = std::pow(x, e) * std::pow(y, n);
      //     ++j;
      //   }
      // }
    }
  }
}
/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ShapeFunctions::interpolateElementalFieldFromIntegrationPoints(
    const Array<Real> & field,
    const Array<Real> & interpolation_points_coordinates_matrices,
    const Array<Real> & quad_points_coordinates_inv_matrices,
    ElementTypeMapArray<Real> & result, const GhostType & ghost_type,
    const Array<UInt> & element_filter) const {
  AKANTU_DEBUG_IN();

  auto nb_element = this->mesh.getNbElement(type, ghost_type);

  auto nb_quad_per_element =
      GaussIntegrationElement<type>::getNbQuadraturePoints();
  auto nb_interpolation_points_per_elem =
      interpolation_points_coordinates_matrices.getNbComponent() /
      nb_quad_per_element;

  if (not result.exists(type, ghost_type))
    result.alloc(nb_element * nb_interpolation_points_per_elem,
                 field.getNbComponent(), type, ghost_type);

  if (element_filter != empty_filter)
    nb_element = element_filter.size();

  Matrix<Real> coefficients(nb_quad_per_element, field.getNbComponent());

  auto & result_vec = result(type, ghost_type);

  auto field_it = field.begin_reinterpret(field.getNbComponent(),
                                          nb_quad_per_element, nb_element);

  auto interpolation_points_coordinates_it =
      interpolation_points_coordinates_matrices.begin(
          nb_interpolation_points_per_elem, nb_quad_per_element);

  auto result_begin = result_vec.begin_reinterpret(
      field.getNbComponent(), nb_interpolation_points_per_elem,
      result_vec.size() / nb_interpolation_points_per_elem);

  auto inv_quad_coord_it = quad_points_coordinates_inv_matrices.begin(
      nb_quad_per_element, nb_quad_per_element);

  /// loop over the elements of the current filter and element type
  for (UInt el = 0; el < nb_element; ++el, ++field_it, ++inv_quad_coord_it,
            ++interpolation_points_coordinates_it) {
    /**
     * matrix containing the inversion of the quadrature points'
     * coordinates
     */
    const auto & inv_quad_coord_matrix = *inv_quad_coord_it;

    /**
     * multiply it by the field values over quadrature points to get
     * the interpolation coefficients
     */
    coefficients.mul<false, true>(inv_quad_coord_matrix, *field_it);

    /// matrix containing the points' coordinates
    const auto & coord = *interpolation_points_coordinates_it;

    /// multiply the coordinates matrix by the coefficients matrix and store the
    /// result
    Matrix<Real> res(result_begin[element_filter(el)]);
    res.mul<true, true>(coefficients, coord);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ShapeFunctions::interpolateElementalFieldOnIntegrationPoints(
    const Array<Real> & u_el, Array<Real> & uq, const GhostType & ghost_type,
    const Array<Real> & shapes, const Array<UInt> & filter_elements) const {
  auto nb_element = mesh.getNbElement(type, ghost_type);
  auto nb_nodes_per_element = ElementClass<type>::getShapeSize();
  auto nb_points = shapes.size() / mesh.getNbElement(type, ghost_type);
  auto nb_degree_of_freedom = u_el.getNbComponent() / nb_nodes_per_element;

  Array<Real>::const_matrix_iterator N_it;
  Array<Real> * filtered_N = nullptr;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
    filtered_N = new Array<Real>(0, shapes.getNbComponent());
    FEEngine::filterElementalData(mesh, shapes, *filtered_N, type, ghost_type,
                                  filter_elements);
    N_it = filtered_N->begin_reinterpret(nb_nodes_per_element, nb_points,
                                         nb_element);
  } else {
    N_it =
        shapes.begin_reinterpret(nb_nodes_per_element, nb_points, nb_element);
  }

  uq.resize(nb_element * nb_points);

  auto u_it = u_el.begin(nb_degree_of_freedom, nb_nodes_per_element);
  auto inter_u_it =
      uq.begin_reinterpret(nb_degree_of_freedom, nb_points, nb_element);

  for (UInt el = 0; el < nb_element; ++el, ++N_it, ++u_it, ++inter_u_it) {
    const auto & u = *u_it;
    const auto & N = *N_it;
    auto & inter_u = *inter_u_it;

    inter_u.template mul<false, false>(u, N);
  }

  delete filtered_N;
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeFunctions::gradientElementalFieldOnIntegrationPoints(
    const Array<Real> & u_el, Array<Real> & out_nablauq,
    const GhostType & ghost_type, const Array<Real> & shapes_derivatives,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  auto nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();
  auto nb_points = integration_points(type, ghost_type).cols();
  auto element_dimension = ElementClass<type>::getNaturalSpaceDimension();
  auto nb_degree_of_freedom = u_el.getNbComponent() / nb_nodes_per_element;
  auto nb_element = mesh.getNbElement(type, ghost_type);

  Array<Real>::const_matrix_iterator B_it;

  Array<Real> * filtered_B = nullptr;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
    filtered_B = new Array<Real>(0, shapes_derivatives.getNbComponent());
    FEEngine::filterElementalData(mesh, shapes_derivatives, *filtered_B, type,
                                  ghost_type, filter_elements);
    B_it = filtered_B->begin(element_dimension, nb_nodes_per_element);
  } else {
    B_it = shapes_derivatives.begin(element_dimension, nb_nodes_per_element);
  }

  out_nablauq.resize(nb_element * nb_points);
  auto u_it = u_el.begin(nb_degree_of_freedom, nb_nodes_per_element);
  auto nabla_u_it = out_nablauq.begin(nb_degree_of_freedom, element_dimension);

  for (UInt el = 0; el < nb_element; ++el, ++u_it) {
    const auto & u = *u_it;
    for (UInt q = 0; q < nb_points; ++q, ++B_it, ++nabla_u_it) {
      const auto & B = *B_it;
      auto & nabla_u = *nabla_u_it;
      nabla_u.template mul<false, true>(u, B);
    }
  }

  delete filtered_B;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
} // namespace akantu

#endif /* __AKANTU_SHAPE_FUNCTIONS_INLINE_IMPL_HH__ */
