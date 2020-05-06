/**
 * @file   geometrical_element_property.hh
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Thomas Menouillard <tmenouillard@stucky.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 29 2017
 *
 * @brief  Specialization of the geometrical types
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
#include "element_class.hh"
/* -------------------------------------------------------------------------- */
#include <array>
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace detail {
  template <typename properties> constexpr size_t sizeFacetConnectivity() {
    size_t s = 0;
    for (size_t n = 0; n < properties::nb_facet_types; ++n) {
      s += properties::nb_facets[n] * properties::nb_nodes_per_facet[n];
    }
    return s == 0 ? 1 : s;
  }
} // namespace detail

#if !defined(DOXYGEN)
template <> struct GeometricalElementProperty<_gt_not_defined> {
  static constexpr UInt spatial_dimension{0};
  static constexpr UInt nb_nodes_per_element{0};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{0}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{0}};
  static constexpr std::array<UInt, 1> facet_connectivity_vect{{0}};
};

template <> struct GeometricalElementProperty<_gt_point> {
  static constexpr UInt spatial_dimension{0};
  static constexpr UInt nb_nodes_per_element{1};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{1}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{1}};
  static constexpr std::array<UInt, 1> facet_connectivity_vect{{0}};
};

template <> struct GeometricalElementProperty<_gt_segment_2> {
  static constexpr UInt spatial_dimension{1};
  static constexpr UInt nb_nodes_per_element{2};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{2}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{1}};
  static constexpr std::array<UInt, 2> facet_connectivity_vect{{0, 1}};
};

template <> struct GeometricalElementProperty<_gt_segment_3> {
  static constexpr UInt spatial_dimension{1};
  static constexpr UInt nb_nodes_per_element{3};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{2}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{1}};
  // clang-format off
  static constexpr std::array<UInt, 2> facet_connectivity_vect{{0, 1}};
  // clang-format on
};

template <> struct GeometricalElementProperty<_gt_triangle_3> {
  static constexpr UInt spatial_dimension{2};
  static constexpr UInt nb_nodes_per_element{3};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{3}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{2}};
  // clang-format off
    static constexpr std::array<UInt, 6> facet_connectivity_vect{{
        0, 1, 2,
        1, 2, 0}};
  // clang-format on
};

template <> struct GeometricalElementProperty<_gt_triangle_6> {
  static constexpr UInt spatial_dimension{2};
  static constexpr UInt nb_nodes_per_element{6};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{3}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{3}};
  // clang-format off
    static constexpr std::array<UInt, 9> facet_connectivity_vect{{
        0, 1, 2,
        1, 2, 0,
        3, 4, 5}};
  // clang-format on
};

template <> struct GeometricalElementProperty<_gt_tetrahedron_4> {
  static constexpr UInt spatial_dimension{3};
  static constexpr UInt nb_nodes_per_element{4};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{4}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{3}};
  // clang-format off
    static constexpr std::array<UInt, 12> facet_connectivity_vect{{
        0, 1, 2, 0,
        2, 2, 0, 1,
        1, 3, 3, 3}};
  // clang-format on
};

template <> struct GeometricalElementProperty<_gt_tetrahedron_10> {
  static constexpr UInt spatial_dimension{3};
  static constexpr UInt nb_nodes_per_element{10};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{4}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{6}};
  // clang-format off
  static constexpr std::array<UInt, 6*4> facet_connectivity_vect{{
      0, 1, 2, 0,
      2, 2, 0, 1,
      1, 3, 3, 3,
      6, 5, 6, 4,
      5, 9, 7, 8,
      4, 8, 9, 7}};
  // clang-format on
};

template <> struct GeometricalElementProperty<_gt_quadrangle_4> {
  static constexpr UInt spatial_dimension{2};
  static constexpr UInt nb_nodes_per_element{4};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{4}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{2}};
  // clang-format off
  static constexpr std::array<UInt, 2*4> facet_connectivity_vect{{
      0, 1, 2, 3,
      1, 2, 3, 0}};
  // clang-format on
};

template <> struct GeometricalElementProperty<_gt_quadrangle_8> {
  static constexpr UInt spatial_dimension{2};
  static constexpr UInt nb_nodes_per_element{8};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{4}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{3}};
  // clang-format off
  static constexpr std::array<UInt, 4*3> facet_connectivity_vect{{
      0, 1, 2, 3,
      1, 2, 3, 0,
      4, 5, 6, 7}};
  // clang-format on
};

template <> struct GeometricalElementProperty<_gt_hexahedron_8> {
  static constexpr UInt spatial_dimension{3};
  static constexpr UInt nb_nodes_per_element{8};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{6}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{4}};
  // clang-format off
  static constexpr std::array<UInt, 4*6> facet_connectivity_vect{{
      0, 0, 1, 2, 3, 4,
      3, 1, 2, 3, 0, 5,
      2, 5, 6, 7, 4, 6,
      1, 4, 5, 6, 7, 7}};
  // clang-format on
};

template <> struct GeometricalElementProperty<_gt_hexahedron_20> {
  static constexpr UInt spatial_dimension{3};
  static constexpr UInt nb_nodes_per_element{20};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{6}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{8}};
  // clang-format off
  static constexpr std::array<UInt, 8*6> facet_connectivity_vect{{
      0,  1,  2,  3,  0,  4,
      1,  2,  3,  0,  3,  5,
      5,  6,  7,  4,  2,  6,
      4,  5,  6,  7,  1,  7,
      8,  9, 10, 11, 11, 16,
      13, 14, 15, 12, 10, 17,
      16, 17, 18, 19,  9, 18,
      12, 13, 14, 15,  8, 19}};
  // clang-format on
};

template <> struct GeometricalElementProperty<_gt_pentahedron_6> {
  static constexpr UInt spatial_dimension{3};
  static constexpr UInt nb_nodes_per_element{6};
  static constexpr UInt nb_facet_types{2};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{2, 3}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{3, 4}};
  // clang-format off
  static constexpr std::array<UInt, 3*2 + 4*3> facet_connectivity_vect{{
      // first type
      0, 3,
      2, 4,
      1, 5,
      // second type
      0, 0, 1,
      1, 3, 2,
      4, 5, 5,
      3, 2, 4}};
  // clang-format on
};

template <> struct GeometricalElementProperty<_gt_pentahedron_15> {
  static constexpr UInt spatial_dimension{3};
  static constexpr UInt nb_nodes_per_element{15};
  static constexpr UInt nb_facet_types{2};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{2, 3}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{6, 8}};
  // clang-format off
  static constexpr std::array<UInt, 6*2 + 8*3> facet_connectivity_vect{{
      // first type
      0, 3,
      2, 4,
      1, 5,
      8, 12,
      7, 13,
      6, 14,
      // second type
       0, 0, 1,
       1, 3, 2,
       4, 5, 5,
       3, 2, 4,
       6, 9, 7,
      10, 14, 11,
      12, 11, 13,
       9, 8, 10}};
  // clang-format on
};

#if defined(AKANTU_COHESIVE_ELEMENT)
/* -------------------------------------------------------------------------- */
template <> struct GeometricalElementProperty<_gt_cohesive_2d_4> {
  static constexpr UInt spatial_dimension{2};
  static constexpr UInt nb_nodes_per_element{4};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{2}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{2}};
  // clang-format off
  static constexpr std::array<UInt, 2 * 2> facet_connectivity_vect{{
      0, 2,
      1, 3}};
  // clang-format on
};

/* -------------------------------------------------------------------------- */
template <> struct GeometricalElementProperty<_gt_cohesive_2d_6> {
  static constexpr UInt spatial_dimension{2};
  static constexpr UInt nb_nodes_per_element{6};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{2}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{3}};
  // clang-format off
  static constexpr std::array<UInt, 3*2> facet_connectivity_vect{{
      0, 3,
      1, 4,
      2, 5}};
  // clang-format on
};

/* -------------------------------------------------------------------------- */
template <> struct GeometricalElementProperty<_gt_cohesive_1d_2> {
  static constexpr UInt spatial_dimension{1};
  static constexpr UInt nb_nodes_per_element{2};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{2}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{1}};
  // clang-format off
  static constexpr std::array<UInt, 2> facet_connectivity_vect{{0, 1}};
  // clang-format on
};

/* -------------------------------------------------------------------------- */
template <> struct GeometricalElementProperty<_gt_cohesive_3d_6> {
  static constexpr UInt spatial_dimension{3};
  static constexpr UInt nb_nodes_per_element{6};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{2}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{3}};
  // clang-format off
  static constexpr std::array<UInt, 3*2> facet_connectivity_vect{{
      0, 3,
      1, 4,
      2, 5}};
  // clang-format on
};

/* -------------------------------------------------------------------------- */
template <> struct GeometricalElementProperty<_gt_cohesive_3d_12> {
  static constexpr UInt spatial_dimension{3};
  static constexpr UInt nb_nodes_per_element{12};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{2}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{6}};
  // clang-format off
  static constexpr std::array<UInt, 6*2> facet_connectivity_vect{{
      0, 6,
      1, 7,
      2, 8,
      3, 9,
      4, 10,
      5, 11}};
  // clang-format on
};

/* -------------------------------------------------------------------------- */
template <> struct GeometricalElementProperty<_gt_cohesive_3d_8> {
  static constexpr UInt spatial_dimension{3};
  static constexpr UInt nb_nodes_per_element{8};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{2}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{4}};
  // clang-format off
  static constexpr std::array<UInt, 4*2> facet_connectivity_vect{{
      0, 4,
      1, 5,
      2, 6,
      3, 7}};
  // clang-format on
};

/* -------------------------------------------------------------------------- */
template <> struct GeometricalElementProperty<_gt_cohesive_3d_16> {
  static constexpr UInt spatial_dimension{3};
  static constexpr UInt nb_nodes_per_element{16};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{2}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{8}};
  // clang-format off
  static constexpr std::array<UInt, 8*2> facet_connectivity_vect{{
      0, 8,
      1, 9,
      2, 10,
      3, 11,
      4, 12,
      5, 13,
      6, 14,
      7, 15}};
  // clang-format on
};
#endif // AKANTU_COHESIVE_ELEMENT

/* -------------------------------------------------------------------------- */
template <> struct ElementClassExtraGeometryProperties<_not_defined> {
  static constexpr ElementType p1_type{_not_defined};
  static constexpr std::array<ElementType, 1> facet_type{{_not_defined}};
};

template <> struct ElementClassExtraGeometryProperties<_point_1> {
  static constexpr ElementType p1_type{_point_1};
  static constexpr std::array<ElementType, 1> facet_type{{_point_1}};
};

template <> struct ElementClassExtraGeometryProperties<_segment_2> {
  static constexpr ElementType p1_type{_segment_2};
  static constexpr std::array<ElementType, 1> facet_type{{_point_1}};
};

template <> struct ElementClassExtraGeometryProperties<_segment_3> {
  static constexpr ElementType p1_type{_segment_2};
  static constexpr std::array<ElementType, 1> facet_type{{_point_1}};
};

template <> struct ElementClassExtraGeometryProperties<_triangle_3> {
  static constexpr ElementType p1_type{_triangle_3};
  static constexpr std::array<ElementType, 1> facet_type{{_segment_2}};
};

template <> struct ElementClassExtraGeometryProperties<_triangle_6> {
  static constexpr ElementType p1_type{_triangle_3};
  static constexpr std::array<ElementType, 1> facet_type{{_segment_3}};
};

template <> struct ElementClassExtraGeometryProperties<_tetrahedron_4> {
  static constexpr ElementType p1_type{_tetrahedron_4};
  static constexpr std::array<ElementType, 1> facet_type{{_triangle_3}};
};

template <> struct ElementClassExtraGeometryProperties<_tetrahedron_10> {
  static constexpr ElementType p1_type{_tetrahedron_4};
  static constexpr std::array<ElementType, 1> facet_type{{_triangle_6}};
};

template <> struct ElementClassExtraGeometryProperties<_quadrangle_4> {
  static constexpr ElementType p1_type{_quadrangle_4};
  static constexpr std::array<ElementType, 1> facet_type{{_segment_2}};
};

template <> struct ElementClassExtraGeometryProperties<_quadrangle_8> {
  static constexpr ElementType p1_type{_quadrangle_4};
  static constexpr std::array<ElementType, 1> facet_type{{_segment_3}};
};

template <> struct ElementClassExtraGeometryProperties<_hexahedron_8> {
  static constexpr ElementType p1_type{_hexahedron_8};
  static constexpr std::array<ElementType, 1> facet_type{{_quadrangle_4}};
};

template <> struct ElementClassExtraGeometryProperties<_hexahedron_20> {
  static constexpr ElementType p1_type{_hexahedron_8};
  static constexpr std::array<ElementType, 1> facet_type{{_quadrangle_8}};
};

template <> struct ElementClassExtraGeometryProperties<_pentahedron_6> {
  static constexpr ElementType p1_type{_pentahedron_6};
  static constexpr std::array<ElementType, 2> facet_type{
      {_triangle_3, _quadrangle_4}};
};

template <> struct ElementClassExtraGeometryProperties<_pentahedron_15> {
  static constexpr ElementType p1_type{_pentahedron_6};
  static constexpr std::array<ElementType, 2> facet_type{
      {_triangle_6, _quadrangle_8}};
};

#if defined(AKANTU_COHESIVE_ELEMENT)
template <> struct ElementClassExtraGeometryProperties<_cohesive_2d_4> {
  static constexpr ElementType p1_type{_cohesive_2d_4};
  static constexpr std::array<ElementType, 1> facet_type{{_segment_2}};
};

template <> struct ElementClassExtraGeometryProperties<_cohesive_2d_6> {
  static constexpr ElementType p1_type{_cohesive_2d_4};
  static constexpr std::array<ElementType, 1> facet_type{{_segment_3}};
};

template <> struct ElementClassExtraGeometryProperties<_cohesive_1d_2> {
  static constexpr ElementType p1_type{_cohesive_1d_2};
  static constexpr std::array<ElementType, 1> facet_type{{_point_1}};
};

template <> struct ElementClassExtraGeometryProperties<_cohesive_3d_6> {
  static constexpr ElementType p1_type{_cohesive_3d_6};
  static constexpr std::array<ElementType, 1> facet_type{{_triangle_3}};
};

template <> struct ElementClassExtraGeometryProperties<_cohesive_3d_12> {
  static constexpr ElementType p1_type{_cohesive_3d_6};
  static constexpr std::array<ElementType, 1> facet_type{{_triangle_6}};
};

template <> struct ElementClassExtraGeometryProperties<_cohesive_3d_8> {
  static constexpr ElementType p1_type{_cohesive_3d_8};
  static constexpr std::array<ElementType, 1> facet_type{{_quadrangle_4}};
};

template <> struct ElementClassExtraGeometryProperties<_cohesive_3d_16> {
  static constexpr ElementType p1_type{_cohesive_3d_8};
  static constexpr std::array<ElementType, 1> facet_type{{_quadrangle_8}};
};
#endif // AKANTU_COHESIVE_ELEMENT

#endif // !defined(DOXYGEN)

} // namespace akantu
