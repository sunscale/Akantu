/**
 * @file   cohesive_element.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Sun Sep 26 2010
 * @date last modification: Mon Sep 14 2015
 *
 * @brief  CohesiveElement implementation
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
#include "aka_common.hh"
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <> UInt GeometricalElement<_gt_cohesive_2d_4>::spatial_dimension = 2;
template <> UInt GeometricalElement<_gt_cohesive_2d_4>::nb_nodes_per_element = 4;
template <> UInt GeometricalElement<_gt_cohesive_2d_4>::nb_facet_types = 1;
template <> UInt GeometricalElement<_gt_cohesive_2d_4>::nb_facets[] = {2};
template <>
UInt GeometricalElement<_gt_cohesive_2d_4>::nb_nodes_per_facet[] = {2};
template <>
UInt GeometricalElement<_gt_cohesive_2d_4>::facet_connectivity_vect[] = {0, 2,
                                                                         1, 3};
template <>
UInt * GeometricalElement<_gt_cohesive_2d_4>::facet_connectivity[] = {
    &facet_connectivity_vect[0]};

template <>
ElementType ElementClass<_cohesive_2d_4>::facet_type[] = {_segment_2};
template <> ElementType ElementClass<_cohesive_2d_4>::p1_type = _cohesive_2d_4;

/* -------------------------------------------------------------------------- */
template <> UInt GeometricalElement<_gt_cohesive_2d_6>::spatial_dimension = 2;
template <>
UInt GeometricalElement<_gt_cohesive_2d_6>::nb_nodes_per_element = 6;
template <> UInt GeometricalElement<_gt_cohesive_2d_6>::nb_facet_types = 1;
template <> UInt GeometricalElement<_gt_cohesive_2d_6>::nb_facets[] = {2};
template <>
UInt GeometricalElement<_gt_cohesive_2d_6>::nb_nodes_per_facet[] = {3};
template <>
UInt GeometricalElement<_gt_cohesive_2d_6>::facet_connectivity_vect[] = {
    0, 3, 1, 4, 2, 5};
template <>
UInt * GeometricalElement<_gt_cohesive_2d_6>::facet_connectivity[] = {
    &facet_connectivity_vect[0]};

template <>
ElementType ElementClass<_cohesive_2d_6>::facet_type[] = {_segment_3};
template <> ElementType ElementClass<_cohesive_2d_6>::p1_type = _cohesive_2d_4;

/* -------------------------------------------------------------------------- */
template <> UInt GeometricalElement<_gt_cohesive_1d_2>::spatial_dimension = 1;
template <>
UInt GeometricalElement<_gt_cohesive_1d_2>::nb_nodes_per_element = 2;
template <> UInt GeometricalElement<_gt_cohesive_1d_2>::nb_facet_types = 1;
template <> UInt GeometricalElement<_gt_cohesive_1d_2>::nb_facets[] = {2};
template <>
UInt GeometricalElement<_gt_cohesive_1d_2>::nb_nodes_per_facet[] = {1};
template <>
UInt GeometricalElement<_gt_cohesive_1d_2>::facet_connectivity_vect[] = {0, 1};
template <>
UInt * GeometricalElement<_gt_cohesive_1d_2>::facet_connectivity[] = {
    &facet_connectivity_vect[0]};

template <> ElementType ElementClass<_cohesive_1d_2>::facet_type[] = {_point_1};
template <> ElementType ElementClass<_cohesive_1d_2>::p1_type = _cohesive_1d_2;
/* -------------------------------------------------------------------------- */

template <> UInt GeometricalElement<_gt_cohesive_3d_6>::spatial_dimension = 3;
template <>
UInt GeometricalElement<_gt_cohesive_3d_6>::nb_nodes_per_element = 6;
template <> UInt GeometricalElement<_gt_cohesive_3d_6>::nb_facet_types = 1;
template <> UInt GeometricalElement<_gt_cohesive_3d_6>::nb_facets[] = {2};
template <>
UInt GeometricalElement<_gt_cohesive_3d_6>::nb_nodes_per_facet[] = {3};
template <>
UInt GeometricalElement<_gt_cohesive_3d_6>::facet_connectivity_vect[] = {
    0, 3, 1, 4, 2, 5};
template <>
UInt * GeometricalElement<_gt_cohesive_3d_6>::facet_connectivity[] = {
    &facet_connectivity_vect[0]};

template <>
ElementType ElementClass<_cohesive_3d_6>::facet_type[] = {_triangle_3};
template <> ElementType ElementClass<_cohesive_3d_6>::p1_type = _cohesive_3d_6;

/* -------------------------------------------------------------------------- */

template <> UInt GeometricalElement<_gt_cohesive_3d_12>::spatial_dimension = 3;
template <>
UInt GeometricalElement<_gt_cohesive_3d_12>::nb_nodes_per_element = 12;
template <> UInt GeometricalElement<_gt_cohesive_3d_12>::nb_facet_types = 1;
template <> UInt GeometricalElement<_gt_cohesive_3d_12>::nb_facets[] = {2};
template <>
UInt GeometricalElement<_gt_cohesive_3d_12>::nb_nodes_per_facet[] = {6};
template <>
UInt GeometricalElement<_gt_cohesive_3d_12>::facet_connectivity_vect[] = {
    0, 6, 1, 7, 2, 8, 3, 9, 4, 10, 5, 11};
template <>
UInt * GeometricalElement<_gt_cohesive_3d_12>::facet_connectivity[] = {
    &facet_connectivity_vect[0]};

template <>
ElementType ElementClass<_cohesive_3d_12>::facet_type[] = {_triangle_6};
template <> ElementType ElementClass<_cohesive_3d_12>::p1_type = _cohesive_3d_6;

/* -------------------------------------------------------------------------- */

template <> UInt GeometricalElement<_gt_cohesive_3d_8>::spatial_dimension = 3;
template <>
UInt GeometricalElement<_gt_cohesive_3d_8>::nb_nodes_per_element = 8;
template <> UInt GeometricalElement<_gt_cohesive_3d_8>::nb_facet_types = 1;
template <> UInt GeometricalElement<_gt_cohesive_3d_8>::nb_facets[] = {2};
template <>
UInt GeometricalElement<_gt_cohesive_3d_8>::nb_nodes_per_facet[] = {4};
template <>
UInt GeometricalElement<_gt_cohesive_3d_8>::facet_connectivity_vect[] = {
    0, 4, 1, 5, 2, 6, 3, 7};
template <>
UInt * GeometricalElement<_gt_cohesive_3d_8>::facet_connectivity[] = {
    &facet_connectivity_vect[0]};

template <>
ElementType ElementClass<_cohesive_3d_8>::facet_type[] = {_quadrangle_4};
template <> ElementType ElementClass<_cohesive_3d_8>::p1_type = _cohesive_3d_8;

/* -------------------------------------------------------------------------- */

template <> UInt GeometricalElement<_gt_cohesive_3d_16>::spatial_dimension = 3;
template <>
UInt GeometricalElement<_gt_cohesive_3d_16>::nb_nodes_per_element = 16;
template <> UInt GeometricalElement<_gt_cohesive_3d_16>::nb_facet_types = 1;
template <> UInt GeometricalElement<_gt_cohesive_3d_16>::nb_facets[] = {2};
template <>
UInt GeometricalElement<_gt_cohesive_3d_16>::nb_nodes_per_facet[] = {8};
template <>
UInt GeometricalElement<_gt_cohesive_3d_16>::facet_connectivity_vect[] = {
    0, 8, 1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15};
template <>
UInt * GeometricalElement<_gt_cohesive_3d_16>::facet_connectivity[] = {
    &facet_connectivity_vect[0]};

template <>
ElementType ElementClass<_cohesive_3d_16>::facet_type[] = {_quadrangle_8};
template <> ElementType ElementClass<_cohesive_3d_16>::p1_type = _cohesive_3d_8;

/* -------------------------------------------------------------------------- */

} // namespace akantu
