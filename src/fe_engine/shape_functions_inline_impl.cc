/**
 * @file   shape_functions_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Oct 27 2010
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  ShapeFunctions inline implementation
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
#include "fe_engine.hh"
#include "shape_functions.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SHAPE_FUNCTIONS_INLINE_IMPL_CC__
#define __AKANTU_SHAPE_FUNCTIONS_INLINE_IMPL_CC__

namespace akantu {

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
    AKANTU_DEBUG_TO_IMPLEMENT();
    break;
  }
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ShapeFunctions::buildElementalFieldInterpolationMatrix(
    const Matrix<Real> &, Matrix<Real> &, UInt) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
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
    AKANTU_DEBUG_TO_IMPLEMENT();
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
    AKANTU_DEBUG_TO_IMPLEMENT();
  } else {
    for (UInt i = 0; i < coordinates.cols(); ++i) {
      UInt j = 0;
      Real x = coordinates(0, i);
      Real y = coordinates(1, i);

      for (UInt e = 0; e <= 2; ++e) {
        for (UInt n = 0; n <= 2; ++n) {
          coordMatrix(i, j) = std::pow(x, e) * std::pow(y, n);
          ++j;
        }
      }
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

  UInt nb_element = this->mesh.getNbElement(type, ghost_type);

  UInt nb_quad_per_element =
      GaussIntegrationElement<type>::getNbQuadraturePoints();
  UInt nb_interpolation_points_per_elem =
      interpolation_points_coordinates_matrices.getNbComponent() /
      nb_quad_per_element;

  if (!result.exists(type, ghost_type))
    result.alloc(nb_element * nb_interpolation_points_per_elem,
                 field.getNbComponent(), type, ghost_type);

  if (element_filter != empty_filter)
    nb_element = element_filter.getSize();

  Matrix<Real> coefficients(nb_quad_per_element, field.getNbComponent());

  Array<Real> & result_vec = result(type, ghost_type);

  Array<Real>::const_matrix_iterator field_it = field.begin_reinterpret(
      field.getNbComponent(), nb_quad_per_element, nb_element);

  Array<Real>::const_matrix_iterator interpolation_points_coordinates_it =
      interpolation_points_coordinates_matrices.begin(
          nb_interpolation_points_per_elem, nb_quad_per_element);

  Array<Real>::matrix_iterator result_begin = result_vec.begin_reinterpret(
      field.getNbComponent(), nb_interpolation_points_per_elem,
      result_vec.getSize() / nb_interpolation_points_per_elem);

  Array<Real>::const_matrix_iterator inv_quad_coord_it =
      quad_points_coordinates_inv_matrices.begin(nb_quad_per_element,
                                                 nb_quad_per_element);

  /// loop over the elements of the current filter and element type
  for (UInt el = 0; el < nb_element; ++el, ++field_it, ++inv_quad_coord_it,
            ++interpolation_points_coordinates_it) {
    /**
     * matrix containing the inversion of the quadrature points'
     * coordinates
     */
    const Matrix<Real> & inv_quad_coord_matrix = *inv_quad_coord_it;

    /**
     * multiply it by the field values over quadrature points to get
     * the interpolation coefficients
     */
    coefficients.mul<false, true>(inv_quad_coord_matrix, *field_it);

    /// matrix containing the points' coordinates
    const Matrix<Real> & coord = *interpolation_points_coordinates_it;

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
  UInt nb_element;
  UInt nb_nodes_per_element = ElementClass<type>::getShapeSize();

  UInt nb_points = shapes.getSize() / mesh.getNbElement(type, ghost_type);
  UInt nb_degree_of_freedom = u_el.getNbComponent() / nb_nodes_per_element;

  Array<Real>::const_matrix_iterator N_it;
  Array<Real>::const_matrix_iterator u_it;
  Array<Real>::matrix_iterator inter_u_it;

  Array<Real> * filtered_N = NULL;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.getSize();
    filtered_N = new Array<Real>(0, shapes.getNbComponent());
    FEEngine::filterElementalData(mesh, shapes, *filtered_N, type, ghost_type,
                                  filter_elements);
    N_it = filtered_N->begin_reinterpret(nb_nodes_per_element, nb_points,
                                         nb_element);
  } else {
    nb_element = mesh.getNbElement(type, ghost_type);
    N_it =
        shapes.begin_reinterpret(nb_nodes_per_element, nb_points, nb_element);
  }

  uq.resize(nb_element * nb_points);

  u_it = u_el.begin(nb_degree_of_freedom, nb_nodes_per_element);
  inter_u_it =
      uq.begin_reinterpret(nb_degree_of_freedom, nb_points, nb_element);

  for (UInt el = 0; el < nb_element; ++el, ++N_it, ++u_it, ++inter_u_it) {
    const Matrix<Real> & u = *u_it;
    const Matrix<Real> & N = *N_it;
    Matrix<Real> & inter_u = *inter_u_it;

    inter_u.mul<false, false>(u, N);
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

  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();
  UInt nb_points = integration_points(type, ghost_type).cols();
  UInt element_dimension = ElementClass<type>::getNaturalSpaceDimension();
  UInt nb_degree_of_freedom = u_el.getNbComponent() / nb_nodes_per_element;

  Array<Real>::const_matrix_iterator B_it;
  Array<Real>::const_matrix_iterator u_it;
  Array<Real>::matrix_iterator nabla_u_it;

  UInt nb_element;
  Array<Real> * filtered_B = NULL;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.getSize();
    filtered_B = new Array<Real>(0, shapes_derivatives.getNbComponent());
    FEEngine::filterElementalData(mesh, shapes_derivatives, *filtered_B, type,
                                  ghost_type, filter_elements);
    B_it = filtered_B->begin(element_dimension, nb_nodes_per_element);
  } else {
    B_it = shapes_derivatives.begin(element_dimension, nb_nodes_per_element);
    nb_element = mesh.getNbElement(type, ghost_type);
  }

  out_nablauq.resize(nb_element * nb_points);

  u_it = u_el.begin(nb_degree_of_freedom, nb_nodes_per_element);
  nabla_u_it = out_nablauq.begin(nb_degree_of_freedom, element_dimension);

  for (UInt el = 0; el < nb_element; ++el, ++u_it) {
    const Matrix<Real> & u = *u_it;
    for (UInt q = 0; q < nb_points; ++q, ++B_it, ++nabla_u_it) {
      const Matrix<Real> & B = *B_it;
      Matrix<Real> & nabla_u = *nabla_u_it;

      nabla_u.mul<false, true>(u, B);
    }
  }

  delete filtered_B;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
} // namespace akantu

#endif /* __AKANTU_SHAPE_FUNCTIONS_INLINE_IMPL_CC__ */
