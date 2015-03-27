/**
 * @file   aabb_primitive.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Mar 13 2015
 * @date last modification: Fri Mar 13 2015
 *
 * @brief  Macro classe (primitive) for AABB CGAL algos
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_AABB_PRIMITIVE_HH__
#define __AKANTU_AABB_PRIMITIVE_HH__

#include "aka_common.hh"
#include "triangle.hh"
#include "tetrahedron.hh"

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> Kernel;

#define AKANTU_AABB_CLASS(name) \
  class name##_primitive {      \
    typedef std::list< name<Kernel> >::iterator Iterator; \
                                                          \
  public:                                                 \
    typedef UInt Id;                                      \
    typedef Kernel::Point_3 Point;                        \
    typedef Kernel::name##_3 Datum;                       \
                                                          \
  public:                                                 \
    name##_primitive() : meshId(0), primitive() {}        \
    name##_primitive(Iterator it) : meshId(it->id()), primitive(*it) {} \
                                                                        \
  public:                                                               \
    const Datum & datum() const { return primitive; }                   \
    const Point & reference_point() const { return primitive.vertex(0); } \
    const Id & id() const { return meshId; }                              \
                                                                          \
  protected:                                                              \
    Id meshId;                                                            \
    name<Kernel> primitive;                                               \
                                                                          \
  }

AKANTU_AABB_CLASS(Triangle);
AKANTU_AABB_CLASS(Tetrahedron);

#undef AKANTU_AABB_CLASS

__END_AKANTU__

#endif // __AKANTU_AABB_PRIMITIVE_HH__
