/**
 * @file   triangle_primitive.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Tue Mar 3 2015
 * @date last modification: Tue Mar 3 2015
 *
 * @brief  Triangle classe (primitive) for AABB CGAL algos
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

#ifndef __AKANTU_TRIANGLE_PRIMITIVE_HH__
#define __AKANTU_TRIANGLE_PRIMITIVE_HH__

#include "aka_common.hh"
#include "triangle.hh"

#include <CGAL/Cartesian.h>

__BEGIN_AKANTU__
  
typedef CGAL::Cartesian<Real> Kernel;

/* -------------------------------------------------------------------------- */

class Triangle_primitive {
  /// Type of storage::iterator
  typedef std::list< Triangle<Kernel> >::iterator Iterator;

/// Typedefs needed by CGAL
public:
  typedef UInt Id;
  typedef Kernel::Point_3 Point;
  typedef Kernel::Triangle_3 Datum;

public:
  /// Default constructor needed
  Triangle_primitive();

  /// Construct from iterator (*it is a Triangle<K>)
  Triangle_primitive(Iterator it);

/// Functions needed by CGAL
public:
  const Datum & datum() const;
  const Point & reference_point() const;
  const Id & id() const;

protected:
  Id meshId;
  Triangle<Kernel> triangle;
};

__END_AKANTU__

#endif // __AKANTU_TRIANGLE_PRIMITIVE_HH__
