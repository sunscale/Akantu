/**
 * @file   geometrical_element.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Thomas Menouillard <tmenouillard@stucky.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Thu Jan 21 2016
 *
 * @brief  Specialization of the geometrical types
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* clang-format off */
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_not_defined>::spatial_dimension          = 0;
template<> UInt GeometricalElement<_gt_not_defined>::nb_nodes_per_element       = 0;
template<> UInt GeometricalElement<_gt_not_defined>::nb_facets[]                = { 0 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_point>::spatial_dimension            = 0;
template<> UInt GeometricalElement<_gt_point>::nb_nodes_per_element         = 1;
template<> UInt GeometricalElement<_gt_point>::nb_facets[]                  = { 1 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_segment_2>::spatial_dimension          = 1;
template<> UInt GeometricalElement<_gt_segment_2>::nb_nodes_per_element       = 2;
template<> UInt GeometricalElement<_gt_segment_2>::nb_facets[]                = { 2 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_segment_3>::spatial_dimension          = 1;
template<> UInt GeometricalElement<_gt_segment_3>::nb_nodes_per_element       = 3;
template<> UInt GeometricalElement<_gt_segment_3>::nb_facets[]                = { 2 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_triangle_3>::spatial_dimension          = 2;
template<> UInt GeometricalElement<_gt_triangle_3>::nb_nodes_per_element       = 3;
template<> UInt GeometricalElement<_gt_triangle_3>::nb_facets[]                = { 3 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_triangle_6>::spatial_dimension          = 2;
template<> UInt GeometricalElement<_gt_triangle_6>::nb_nodes_per_element       = 6;
template<> UInt GeometricalElement<_gt_triangle_6>::nb_facets[]                = { 3 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_tetrahedron_4>::spatial_dimension         = 3;
template<> UInt GeometricalElement<_gt_tetrahedron_4>::nb_nodes_per_element      = 4;
template<> UInt GeometricalElement<_gt_tetrahedron_4>::nb_facets[]               = { 4 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_tetrahedron_10>::spatial_dimension        = 3;
template<> UInt GeometricalElement<_gt_tetrahedron_10>::nb_nodes_per_element     = 10;
template<> UInt GeometricalElement<_gt_tetrahedron_10>::nb_facets[]              = { 4 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_quadrangle_4>::spatial_dimension          = 2;
template<> UInt GeometricalElement<_gt_quadrangle_4>::nb_nodes_per_element       = 4;
template<> UInt GeometricalElement<_gt_quadrangle_4>::nb_facets[]                = { 4 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_quadrangle_8>::spatial_dimension          = 2;
template<> UInt GeometricalElement<_gt_quadrangle_8>::nb_nodes_per_element       = 8;
template<> UInt GeometricalElement<_gt_quadrangle_8>::nb_facets[]                = { 4 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_hexahedron_8>::spatial_dimension          = 3;
template<> UInt GeometricalElement<_gt_hexahedron_8>::nb_nodes_per_element       = 8;
template<> UInt GeometricalElement<_gt_hexahedron_8>::nb_facets[]                = { 6 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_hexahedron_20>::spatial_dimension         = 3;
template<> UInt GeometricalElement<_gt_hexahedron_20>::nb_nodes_per_element      = 20;
template<> UInt GeometricalElement<_gt_hexahedron_20>::nb_facets[]               = { 6 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_pentahedron_6>::spatial_dimension         = 3;
template<> UInt GeometricalElement<_gt_pentahedron_6>::nb_nodes_per_element      = 6;
template<> UInt GeometricalElement<_gt_pentahedron_6>::nb_facets[]               = { 2, 3 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_pentahedron_15>::spatial_dimension        = 3;
template<> UInt GeometricalElement<_gt_pentahedron_15>::nb_nodes_per_element     = 15;
template<> UInt GeometricalElement<_gt_pentahedron_15>::nb_facets[]              = { 2, 3 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_not_defined>::nb_nodes_per_facet[]    = { 0 };
template<> UInt GeometricalElement<_gt_point>::nb_nodes_per_facet[]          = { 1 };
template<> UInt GeometricalElement<_gt_segment_2>::nb_nodes_per_facet[]      = { 1 };
template<> UInt GeometricalElement<_gt_segment_3>::nb_nodes_per_facet[]      = { 1 };
template<> UInt GeometricalElement<_gt_triangle_3>::nb_nodes_per_facet[]     = { 2 };
template<> UInt GeometricalElement<_gt_triangle_6>::nb_nodes_per_facet[]     = { 3 };
template<> UInt GeometricalElement<_gt_tetrahedron_4>::nb_nodes_per_facet[]  = { 3 };
template<> UInt GeometricalElement<_gt_tetrahedron_10>::nb_nodes_per_facet[] = { 6 };
template<> UInt GeometricalElement<_gt_quadrangle_4>::nb_nodes_per_facet[]   = { 2 };
template<> UInt GeometricalElement<_gt_quadrangle_8>::nb_nodes_per_facet[]   = { 3 };
template<> UInt GeometricalElement<_gt_hexahedron_8>::nb_nodes_per_facet[]   = { 4 };
template<> UInt GeometricalElement<_gt_hexahedron_20>::nb_nodes_per_facet[]  = { 8 };
template<> UInt GeometricalElement<_gt_pentahedron_6>::nb_nodes_per_facet[]  = { 3, 4 };
template<> UInt GeometricalElement<_gt_pentahedron_15>::nb_nodes_per_facet[] = { 6, 8 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_not_defined>::nb_facet_types    = 1;
template<> UInt GeometricalElement<_gt_point>::nb_facet_types          = 1;
template<> UInt GeometricalElement<_gt_segment_2>::nb_facet_types      = 1;
template<> UInt GeometricalElement<_gt_segment_3>::nb_facet_types      = 1;
template<> UInt GeometricalElement<_gt_triangle_3>::nb_facet_types     = 1;
template<> UInt GeometricalElement<_gt_triangle_6>::nb_facet_types     = 1;
template<> UInt GeometricalElement<_gt_tetrahedron_4>::nb_facet_types  = 1;
template<> UInt GeometricalElement<_gt_tetrahedron_10>::nb_facet_types = 1;
template<> UInt GeometricalElement<_gt_quadrangle_4>::nb_facet_types   = 1;
template<> UInt GeometricalElement<_gt_quadrangle_8>::nb_facet_types   = 1;
template<> UInt GeometricalElement<_gt_hexahedron_8>::nb_facet_types   = 1;
template<> UInt GeometricalElement<_gt_hexahedron_20>::nb_facet_types  = 1;
template<> UInt GeometricalElement<_gt_pentahedron_6>::nb_facet_types  = 2;
template<> UInt GeometricalElement<_gt_pentahedron_15>::nb_facet_types = 2;
/* !!! stored as a matrix nb_facets X nb_nodes_per_facet in COL MAJOR */
/* -------------------------------------------------------------------------- */
/*                                                                                  f1|f2|f3|f4|f5|f6 */
template<> UInt GeometricalElement<_gt_not_defined>::facet_connectivity_vect[]    = {0};
template<> UInt GeometricalElement<_gt_point>::facet_connectivity_vect[]          = {0};
template<> UInt GeometricalElement<_gt_segment_2>::facet_connectivity_vect[]      = {0, 1};
template<> UInt GeometricalElement<_gt_segment_3>::facet_connectivity_vect[]      = {0, 1};
template<> UInt GeometricalElement<_gt_triangle_3>::facet_connectivity_vect[]     = {0, 1, 2,
                                                                                     1, 2, 0};
template<> UInt GeometricalElement<_gt_triangle_6>::facet_connectivity_vect[]     = {0, 1, 2,
                                                                                     1, 2, 0,
                                                                                     3, 4, 5};
template<> UInt GeometricalElement<_gt_tetrahedron_4>::facet_connectivity_vect[]  = {0, 1, 2, 0,
                                                                                     2, 2, 0, 1,
                                                                                     1, 3, 3, 3};
template<> UInt GeometricalElement<_gt_tetrahedron_10>::facet_connectivity_vect[] = {0, 1, 2, 0,
                                                                                     2, 2, 0, 1,
                                                                                     1, 3, 3, 3,
                                                                                     6, 5, 6, 4,
                                                                                     5, 9, 7, 8,
                                                                                     4, 8, 9, 7};
template<> UInt GeometricalElement<_gt_quadrangle_4>::facet_connectivity_vect[]   = {0, 1, 2, 3,
                                                                                     1, 2, 3, 0};
template<> UInt GeometricalElement<_gt_quadrangle_8>::facet_connectivity_vect[]   = {0, 1, 2, 3,
                                                                                     1, 2, 3, 0,
                                                                                     4, 5, 6, 7};
template<> UInt GeometricalElement<_gt_hexahedron_8>::facet_connectivity_vect[]   = {0, 0, 1, 2, 3, 4,
                                                                                     3, 1, 2, 3, 0, 5,
                                                                                     2, 5, 6, 7, 4, 6,
                                                                                     1, 4, 5, 6, 7, 7};
template<> UInt GeometricalElement<_gt_hexahedron_20>::facet_connectivity_vect[]  = {0,   1,  2,  3,  0,  4,
                                                                                     1,   2,  3,  0,  3,  5,
                                                                                     5,   6,  7,  4,  2,  6,
                                                                                     4,   5,  6,  7,  1,  7,
                                                                                     8,   9, 10, 11, 11, 16,
                                                                                     13, 14, 15, 12, 10, 17,
                                                                                     16, 17, 18, 19,  9, 18,
                                                                                     12, 13, 14, 15,  8, 19};

template<> UInt GeometricalElement<_gt_pentahedron_6>::facet_connectivity_vect[]  = {// first type
                                                                                     0, 3,
                                                                                     2, 4,
                                                                                     1, 5,
                                                                                     // second type
                                                                                     0, 0, 1,
                                                                                     1, 3, 2,
                                                                                     4, 5, 5,
                                                                                     3, 2, 4};
template<> UInt GeometricalElement<_gt_pentahedron_15>::facet_connectivity_vect[] = {// first type
                                                                                     0,  3,
                                                                                     2,  4,
                                                                                     1,  5,
                                                                                     8,  12,
                                                                                     7,  13,
                                                                                     6,  14,
                                                                                     // second type
                                                                                     0,  0,  1,
                                                                                     1,  3,  2,
                                                                                     4,  5,  5,
                                                                                     3,  2,  4,
                                                                                     6,  9,  7,
                                                                                     10, 14, 11,
                                                                                     12, 11, 13,
                                                                                     9,  8,  10};

template<> UInt * GeometricalElement<_gt_not_defined>::facet_connectivity[]    = { &facet_connectivity_vect[0] };
template<> UInt * GeometricalElement<_gt_point>::facet_connectivity[]          = { &facet_connectivity_vect[0] };
template<> UInt * GeometricalElement<_gt_segment_2>::facet_connectivity[]      = { &facet_connectivity_vect[0] };
template<> UInt * GeometricalElement<_gt_segment_3>::facet_connectivity[]      = { &facet_connectivity_vect[0] };
template<> UInt * GeometricalElement<_gt_triangle_3>::facet_connectivity[]     = { &facet_connectivity_vect[0] };
template<> UInt * GeometricalElement<_gt_triangle_6>::facet_connectivity[]     = { &facet_connectivity_vect[0] };
template<> UInt * GeometricalElement<_gt_tetrahedron_4>::facet_connectivity[]  = { &facet_connectivity_vect[0] };
template<> UInt * GeometricalElement<_gt_tetrahedron_10>::facet_connectivity[] = { &facet_connectivity_vect[0] };
template<> UInt * GeometricalElement<_gt_quadrangle_4>::facet_connectivity[]   = { &facet_connectivity_vect[0] };
template<> UInt * GeometricalElement<_gt_quadrangle_8>::facet_connectivity[]   = { &facet_connectivity_vect[0] };
template<> UInt * GeometricalElement<_gt_hexahedron_8>::facet_connectivity[]   = { &facet_connectivity_vect[0] };
template<> UInt * GeometricalElement<_gt_hexahedron_20>::facet_connectivity[]  = { &facet_connectivity_vect[0] };
template<> UInt * GeometricalElement<_gt_pentahedron_6>::facet_connectivity[]  = { &facet_connectivity_vect[0], &facet_connectivity_vect[2*3] };
template<> UInt * GeometricalElement<_gt_pentahedron_15>::facet_connectivity[] = { &facet_connectivity_vect[0], &facet_connectivity_vect[2*6] };

/* clang-format on */
/* -------------------------------------------------------------------------- */

} // akantu
