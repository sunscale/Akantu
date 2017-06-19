/**
 * @file   element_class.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Peter Spijker <peter.spijker@epfl.ch>
 *
 * @date creation: Tue Jul 20 2010
 * @date last modification: Thu Jul 16 2015
 *
 * @brief  Common part of element_classes
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
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

using std::sqrt;

namespace akantu {

/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_not_defined>::p1_type       = _not_defined;
template<> ElementType ElementClass<_not_defined>::facet_type[]  = {_not_defined};
/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_point_1>::p1_type           = _point_1;
template<> ElementType ElementClass<_point_1>::facet_type[]      = {_point_1};
/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_segment_2>::p1_type         = _segment_2;
template<> ElementType ElementClass<_segment_2>::facet_type[]    = {_point_1};
/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_segment_3>::p1_type         = _segment_2;
template<> ElementType ElementClass<_segment_3>::facet_type[]    = {_point_1};
/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_triangle_3>::p1_type        = _triangle_3;
template<> ElementType ElementClass<_triangle_3>::facet_type[]   = {_segment_2};
/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_triangle_6>::p1_type        = _triangle_3;
template<> ElementType ElementClass<_triangle_6>::facet_type[]   = {_segment_3};
/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_tetrahedron_4>::p1_type     = _tetrahedron_4;
template<> ElementType ElementClass<_tetrahedron_4>::facet_type[]= {_triangle_3};
/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_tetrahedron_10>::p1_type    = _tetrahedron_4;
template<> ElementType ElementClass<_tetrahedron_10>::facet_type[] = {_triangle_6};
/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_quadrangle_4>::p1_type      = _quadrangle_4;
template<> ElementType ElementClass<_quadrangle_4>::facet_type[] = {_segment_2};
/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_quadrangle_8>::p1_type      = _quadrangle_4;
template<> ElementType ElementClass<_quadrangle_8>::facet_type[] = { _segment_3 };
/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_hexahedron_8>::p1_type      = _hexahedron_8;
template<> ElementType ElementClass<_hexahedron_8>::facet_type[] = { _quadrangle_4 };
/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_hexahedron_20>::p1_type      = _hexahedron_8;
template<> ElementType ElementClass<_hexahedron_20>::facet_type[] = { _quadrangle_8 };
/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_pentahedron_6>::p1_type      = _pentahedron_6;
template<> ElementType ElementClass<_pentahedron_6>::facet_type[] = { _triangle_3, _quadrangle_4 };
/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_pentahedron_15>::p1_type      = _pentahedron_6;
template<> ElementType ElementClass<_pentahedron_15>::facet_type[] = { _triangle_6, _quadrangle_8 };
/* -------------------------------------------------------------------------- */
} // akantu
