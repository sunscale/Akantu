/**
 * @file   geometrical_data_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Jan 16 2013
 * @date last modification: Fri Jun 13 2014
 *
 * @brief  Specialization of the geometrical types
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
template<>
struct GeometricalElementData<_gt_point> {
  static const UInt spatial_dimension          = 0;
  static const UInt nb_nodes_per_element       = 1;
  static const GeometricalType p1_element_type = _gt_point;
  static const GeometricalType facet_type      = _gt_not_defined;
  static const UInt nb_facets                  = 0;
  static const UInt * facet_connectivity[]     = {};
};

/* -------------------------------------------------------------------------- */
template<>
struct GeometricalElementData<_gt_segment_2> {
  static const UInt spatial_dimension          = 1;
  static const UInt nb_nodes_per_element       = 2;
  static const GeometricalType p1_element_type = _gt_segment_2;
  static const UInt nb_facets                  = 2;
  static const GeometricalType facet_type      = _gt_point;
  static const UInt vec_facet_connectivity[]   = {0,
						  1};
  static const UInt * facet_connectivity[]     = {&vec_facet_connectivity[0],
						  &vec_facet_connectivity[1]};
};

/* -------------------------------------------------------------------------- */
template<>
struct GeometricalElementData<_gt_segment_3> {
  static const UInt spatial_dimension          = 1;
  static const UInt nb_nodes_per_element       = 3;
  static const GeometricalType p1_element_type = _gt_segment_2;
  static const UInt nb_facets                  = 2;
  static const GeometricalType facet_type      = _gt_point;
  static const UInt vec_facet_connectivity[]   = {0,
						  1};
  static const UInt * facet_connectivity[]     = {&vec_facet_connectivity[0],
						  &vec_facet_connectivity[1]};
};

/* -------------------------------------------------------------------------- */
template<>
struct GeometricalElementData<_gt_triangle_3> {
  static const UInt spatial_dimension          = 2;
  static const UInt nb_nodes_per_element       = 3;
  static const GeometricalType p1_element_type = _gt_triangle_3;
  static const UInt nb_facets                  = 3;
  static const GeometricalType facet_type      = _gt_segment_2;
  static const UInt vec_facet_connectivity[]   = {0, 1,
						  1, 2,
						  2, 0};
  static const UInt * facet_connectivity[]     = {&vec_facet_connectivity[0],
						  &vec_facet_connectivity[2],
						  &vec_facet_connectivity[4]};
};

/* -------------------------------------------------------------------------- */
template<>
struct GeometricalElementData<_gt_triangle_6> {
  static const UInt spatial_dimension          = 2;
  static const UInt nb_nodes_per_element       = 6;
  static const GeometricalType p1_element_type = _gt_triangle_3;
  static const UInt nb_facets                  = 3;
  static const GeometricalType facet_type      = _gt_segment_3;
  static const UInt vec_facet_connectivity[]   = {0, 1, 3,
						  1, 2, 4,
						  2, 0, 5};
  static const UInt * facet_connectivity[]     = {&vec_facet_connectivity[0],
						  &vec_facet_connectivity[3],
						  &vec_facet_connectivity[6]};
};

/* -------------------------------------------------------------------------- */
template<>
struct GeometricalElementData<_gt_tetrahedron_4> {
  static const UInt spatial_dimension          = 3;
  static const UInt nb_nodes_per_element       = 4;
  static const GeometricalType p1_element_type = _gt_tetrahedron_4;
  static const UInt nb_facets                  = 4;
  static const GeometricalType facet_type      = _gt_triangle_3;
  static const UInt vec_facet_connectivity[]   = {0, 2, 1,
						  1, 2, 3,
						  2, 0, 3,
						  0, 1, 3};
  static const UInt * facet_connectivity[]     = {&vec_facet_connectivity[0],
						  &vec_facet_connectivity[3],
						  &vec_facet_connectivity[6],
						  &vec_facet_connectivity[9]};
};

/* -------------------------------------------------------------------------- */
template<>
struct GeometricalElementData<_gt_tetrahedron_10> {
  static const UInt spatial_dimension          = 3;
  static const UInt nb_nodes_per_element       = 10;
  static const GeometricalType p1_element_type = _gt_tetrahedron_4;
  static const UInt nb_facets                  = 4;
  static const GeometricalType facet_type      = _gt_triangle_6;
  static const UInt vec_facet_connectivity[]   = {0, 2, 1, 6, 5, 4,
						  1, 2, 3, 5, 9, 8,
						  2, 0, 3, 6, 7, 9,
						  0, 1, 3, 4, 8, 7};
  static const UInt * facet_connectivity[]     = {&vec_facet_connectivity[0],
						  &vec_facet_connectivity[6],
						  &vec_facet_connectivity[12],
						  &vec_facet_connectivity[18]};
};

/* -------------------------------------------------------------------------- */
template<>
struct GeometricalElementData<_gt_quadrangle_4> {
  static const UInt spatial_dimension          = 2;
  static const UInt nb_nodes_per_element       = 4;
  static const GeometricalType p1_element_type = _gt_quadrangle_4;
  static const UInt nb_facets                  = 4;
  static const GeometricalType facet_type      = _gt_segment_2;
  static const UInt vec_facet_connectivity[]   = {0, 1,
						  1, 2,
						  2, 3,
						  3, 0};
  static const UInt * facet_connectivity[]     = {&vec_facet_connectivity[0],
						  &vec_facet_connectivity[2],
						  &vec_facet_connectivity[4],
						  &vec_facet_connectivity[6]};
};

/* -------------------------------------------------------------------------- */
template<>
struct GeometricalElementData<_gt_quadrangle_8> {
  static const UInt spatial_dimension          = 2;
  static const UInt nb_nodes_per_element       = 8;
  static const GeometricalType p1_element_type = _gt_quadrangle_4;
  static const UInt nb_facets                  = 4;
  static const GeometricalType facet_type      = _gt_segment_3;
  static const UInt vec_facet_connectivity[]   = {0, 1, 4,
						  1, 2, 5,
						  2, 3, 6,
						  3, 0, 7};
  static const UInt * facet_connectivity[]     = {vec_facet_connectivity + 0,
						  vec_facet_connectivity + 3,
						  vec_facet_connectivity + 6,
						  vec_facet_connectivity + 9};
};

/* -------------------------------------------------------------------------- */
template<>
struct GeometricalElementData<_gt_hexahedron_8> {
  static const UInt spatial_dimension          = 3;
  static const UInt nb_nodes_per_element       = 8;
  static const GeometricalType p1_element_type = _gt_hexahedron_8;
  static const UInt nb_facets                  = 6;
  static const GeometricalType facet_type      = _gt_quadrangle_4;
  static const UInt vec_facet_connectivity[]   = {0, 1, 2, 3,
						  0, 1, 5, 4,
						  1, 2, 6, 5,
						  2, 3, 7, 6,
						  3, 0, 4, 7,
						  4, 5, 6, 7};
  static const UInt * facet_connectivity[]     = {&vec_facet_connectivity[0],
						  &vec_facet_connectivity[4],
						  &vec_facet_connectivity[8],
						  &vec_facet_connectivity[12],
						  &vec_facet_connectivity[16],
						  &vec_facet_connectivity[20]};
};

