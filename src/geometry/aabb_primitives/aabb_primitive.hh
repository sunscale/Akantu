/**
 * @file   aabb_primitive.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Clement Roux <clement.roux@epfl.ch>
 *
 * @date creation: Fri Mar 13 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Macro classe (primitive) for AABB CGAL algos
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_AABB_PRIMITIVE_HH_
#define AKANTU_AABB_PRIMITIVE_HH_

#include "aka_common.hh"
#include "line_arc.hh"
#include "tetrahedron.hh"
#include "triangle.hh"

#include "mesh_geom_common.hh"

namespace akantu {

/**
 * This macro defines a class that is used in the CGAL AABB tree algorithm.
 * All the `typedef`s and methods are required by the AABB module.
 *
 * The member variables are
 *  - the id of the element associated to the primitive
 *  - the geometric primitive of the element
 *
 *  @param name the name of the primitive type
 *  @param kernel the name of the kernel used
 */
#define AKANTU_AABB_CLASS(name, kernel)                                        \
  class name##_primitive {                                                     \
    using Iterator = std::list<name<kernel>>::iterator; /* NOLINT */           \
                                                                               \
  public:                                                                      \
    using Id = UInt;                                                           \
    using Point = kernel::Point_3;                                             \
    using Datum = kernel::name##_3;                                            \
                                                                               \
  public:                                                                      \
    name##_primitive() = default;                                              \
    name##_primitive(Iterator it) : meshId(it->id()), primitive(*it) {}        \
                                                                               \
  public:                                                                      \
    const Datum & datum() const { return primitive; }                          \
    Point reference_point() const;                                             \
    const Id & id() const { return meshId; }                                   \
                                                                               \
  protected:                                                                   \
    Id meshId{0};                                                              \
    name<kernel> primitive; /* NOLINT */                                       \
  }

// If the primitive is supported by CGAL::intersection() then the
// implementation process is really easy with this macro
AKANTU_AABB_CLASS(Triangle, cgal::Cartesian);
AKANTU_AABB_CLASS(Line_arc, cgal::Spherical);

#undef AKANTU_AABB_CLASS

} // namespace akantu

#endif // AKANTU_AABB_PRIMITIVE_HH_
