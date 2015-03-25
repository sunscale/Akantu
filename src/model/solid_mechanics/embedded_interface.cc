/**
 * @file   embedded_interface.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Mon Mar 23 2015
 * @date last modification: Mon Mar 23 2015
 *
 * @brief Class used to represent embedded interfaces
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "aka_common.hh"
#include "embedded_interface.hh"

#include <CGAL/Cartesian.h>

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

EmbeddedInterface::EmbeddedInterface(UInt dim, const ID & id) :
  Parsable(_st_embedded_interface, id),
  points(2, dim)
{
  registerParam<std::string>("name", name, _pat_parsable | _pat_readable, "ID");
  registerParam<std::string>("material", mat_name, _pat_parsable | _pat_readable, "Material Name");
  registerParam< Matrix<Real> >("points", points, _pat_parsable | _pat_readable, "Interface Points");
}

EmbeddedInterface::~EmbeddedInterface()
{}

K::Segment_3 EmbeddedInterface::getPrimitive() const {
  if (points.cols() == 2) {
    K::Point_3 a(points(0, 0), points(0, 1), 0.);
    K::Point_3 b(points(1, 0), points(1, 1), 0.);
    return K::Segment_3(a, b);
  } else if (points.cols() == 3) {
    K::Point_3 a(points(0, 0), points(0, 1), points(0, 2));
    K::Point_3 b(points(1, 0), points(1, 1), points(1, 2));
    return K::Segment_3(a, b);
  } else {
    AKANTU_DEBUG_ERROR("Error in reading embedded interface : wrong number of components");
    return K::Segment_3(CGAL::ORIGIN, CGAL::ORIGIN);
  }
}

__END_AKANTU__
