/**
 * @file   shape_functions.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 09 2017
 * @date last modification: Wed Oct 11 2017
 *
 * @brief  implementation of th shape functions interface
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
#include "shape_functions.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ShapeFunctions::ShapeFunctions(const Mesh & mesh, UInt spatial_dimension,
                               const ID & id, const MemoryID & memory_id)
    : Memory(id, memory_id), shapes("shapes_generic", id, memory_id),
      shapes_derivatives("shapes_derivatives_generic", id, memory_id),
      mesh(mesh), _spatial_dimension(spatial_dimension) {}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void
ShapeFunctions::initElementalFieldInterpolationFromIntegrationPoints(
    const Array<Real> & interpolation_points_coordinates,
    ElementTypeMapArray<Real> & interpolation_points_coordinates_matrices,
    ElementTypeMapArray<Real> & quad_points_coordinates_inv_matrices,
    const Array<Real> & quadrature_points_coordinates,
    GhostType ghost_type, const Array<UInt> & element_filter) const {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = this->mesh.getSpatialDimension();
  UInt nb_element = this->mesh.getNbElement(type, ghost_type);
  UInt nb_element_filter;

  if (element_filter == empty_filter) {
    nb_element_filter = nb_element;
  } else {
    nb_element_filter = element_filter.size();
  }

  auto nb_quad_per_element =
      GaussIntegrationElement<type>::getNbQuadraturePoints();
  auto nb_interpolation_points_per_elem =
      interpolation_points_coordinates.size() / nb_element;

  AKANTU_DEBUG_ASSERT(interpolation_points_coordinates.size() % nb_element == 0,
                      "Number of interpolation points should be a multiple of "
                      "total number of elements");

  if (not quad_points_coordinates_inv_matrices.exists(type, ghost_type)) {
    quad_points_coordinates_inv_matrices.alloc(
        nb_element_filter, nb_quad_per_element * nb_quad_per_element, type,
        ghost_type);
  } else {
    quad_points_coordinates_inv_matrices(type, ghost_type)
        .resize(nb_element_filter);
  }

  if (!interpolation_points_coordinates_matrices.exists(type, ghost_type)) {
    interpolation_points_coordinates_matrices.alloc(
        nb_element_filter,
        nb_interpolation_points_per_elem * nb_quad_per_element, type,
        ghost_type);
  } else {
    interpolation_points_coordinates_matrices(type, ghost_type)
        .resize(nb_element_filter);
  }

  Array<Real> & quad_inv_mat =
      quad_points_coordinates_inv_matrices(type, ghost_type);
  Array<Real> & interp_points_mat =
      interpolation_points_coordinates_matrices(type, ghost_type);

  Matrix<Real> quad_coord_matrix(nb_quad_per_element, nb_quad_per_element);

  Array<Real>::const_matrix_iterator quad_coords_it =
      quadrature_points_coordinates.begin_reinterpret(
          spatial_dimension, nb_quad_per_element, nb_element_filter);

  Array<Real>::const_matrix_iterator points_coords_begin =
      interpolation_points_coordinates.begin_reinterpret(
          spatial_dimension, nb_interpolation_points_per_elem, nb_element);

  Array<Real>::matrix_iterator inv_quad_coord_it =
      quad_inv_mat.begin(nb_quad_per_element, nb_quad_per_element);

  Array<Real>::matrix_iterator int_points_mat_it = interp_points_mat.begin(
      nb_interpolation_points_per_elem, nb_quad_per_element);

  /// loop over the elements of the current material and element type
  for (UInt el = 0; el < nb_element_filter;
       ++el, ++inv_quad_coord_it, ++int_points_mat_it, ++quad_coords_it) {
    /// matrix containing the quadrature points coordinates
    const Matrix<Real> & quad_coords = *quad_coords_it;
    /// matrix to store the matrix inversion result
    Matrix<Real> & inv_quad_coord_matrix = *inv_quad_coord_it;

    /// insert the quad coordinates in a matrix compatible with the
    /// interpolation
    buildElementalFieldInterpolationMatrix<type>(quad_coords,
                                                 quad_coord_matrix);

    /// invert the interpolation matrix
    inv_quad_coord_matrix.inverse(quad_coord_matrix);

    /// matrix containing the interpolation points coordinates
    const Matrix<Real> & points_coords =
        points_coords_begin[element_filter(el)];
    /// matrix to store the interpolation points coordinates
    /// compatible with these functions
    Matrix<Real> & inv_points_coord_matrix = *int_points_mat_it;

    /// insert the quad coordinates in a matrix compatible with the
    /// interpolation
    buildElementalFieldInterpolationMatrix<type>(points_coords,
                                                 inv_points_coord_matrix);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ShapeFunctions::initElementalFieldInterpolationFromIntegrationPoints(
    const ElementTypeMapArray<Real> & interpolation_points_coordinates,
    ElementTypeMapArray<Real> & interpolation_points_coordinates_matrices,
    ElementTypeMapArray<Real> & quad_points_coordinates_inv_matrices,
    const ElementTypeMapArray<Real> & quadrature_points_coordinates,
    const ElementTypeMapArray<UInt> * element_filter) const {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = this->mesh.getSpatialDimension();

  for (auto ghost_type : ghost_types) {
    auto types_iterable = mesh.elementTypes(spatial_dimension, ghost_type);
    if (element_filter != nullptr) {
      types_iterable =
          element_filter->elementTypes(spatial_dimension, ghost_type);
    }

    for (auto type : types_iterable) {
      UInt nb_element = mesh.getNbElement(type, ghost_type);
      if (nb_element == 0) {
        continue;
      }

      const Array<UInt> * elem_filter;
      if (element_filter != nullptr) {
        elem_filter = &((*element_filter)(type, ghost_type));
      } else {
        elem_filter = &(empty_filter);
      }

#define AKANTU_INIT_ELEMENTAL_FIELD_INTERPOLATION_FROM_C_POINTS(type)          \
  this->initElementalFieldInterpolationFromIntegrationPoints<type>(            \
      interpolation_points_coordinates(type, ghost_type),                      \
      interpolation_points_coordinates_matrices,                               \
      quad_points_coordinates_inv_matrices,                                    \
      quadrature_points_coordinates(type, ghost_type), ghost_type,             \
      *elem_filter)

      AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(
          AKANTU_INIT_ELEMENTAL_FIELD_INTERPOLATION_FROM_C_POINTS);
#undef AKANTU_INIT_ELEMENTAL_FIELD_INTERPOLATION_FROM_C_POINTS
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ShapeFunctions::interpolateElementalFieldFromIntegrationPoints(
    const ElementTypeMapArray<Real> & field,
    const ElementTypeMapArray<Real> & interpolation_points_coordinates_matrices,
    const ElementTypeMapArray<Real> & quad_points_coordinates_inv_matrices,
    ElementTypeMapArray<Real> & result, GhostType ghost_type,
    const ElementTypeMapArray<UInt> * element_filter) const {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = this->mesh.getSpatialDimension();

  auto types_iterable = mesh.elementTypes(spatial_dimension, ghost_type);
  if (element_filter != nullptr) {
    types_iterable =
        element_filter->elementTypes(spatial_dimension, ghost_type);
  }

  for (auto type : types_iterable) {
    UInt nb_element = mesh.getNbElement(type, ghost_type);
    if (nb_element == 0) {
      continue;
    }
    
    const Array<UInt> * elem_filter;
    if (element_filter != nullptr) {
      elem_filter = &((*element_filter)(type, ghost_type));
    } else {
      elem_filter = &(empty_filter);
    }

#define AKANTU_INTERPOLATE_ELEMENTAL_FIELD_FROM_C_POINTS(type)                 \
  interpolateElementalFieldFromIntegrationPoints<type>(                        \
      field(type, ghost_type),                                                 \
      interpolation_points_coordinates_matrices(type, ghost_type),             \
      quad_points_coordinates_inv_matrices(type, ghost_type), result,          \
      ghost_type, *elem_filter)

    AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(
        AKANTU_INTERPOLATE_ELEMENTAL_FIELD_FROM_C_POINTS);
#undef AKANTU_INTERPOLATE_ELEMENTAL_FIELD_FROM_C_POINTS
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
