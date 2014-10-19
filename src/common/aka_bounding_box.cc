/**
 * @file   aka_bounding_box.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Fri Mar 21 2014
 *
 * @brief  Implementation of the bounding box class
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

#include "aka_bounding_box.hh"
#include "aka_array.hh"

__BEGIN_AKANTU__


template <>
std::ostream& operator<< <1>(std::ostream& os, const BoundingBox<1>& bb) {
  
  os<<"Line["<<bb.min()<<","<<bb.max()<<"]";
  return os;
}


template <>
std::ostream& operator<< <2>(std::ostream& os, const BoundingBox<2>& bb) {
    
    os<<"Rectangle["<<bb.min()<<","<<bb.max()<<"]";
    return os;    
}

template <>
std::ostream& operator<< <3>(std::ostream& os, const BoundingBox<3>& bb) {
    
    os<<"Cuboid["<<bb.min()<<","<<bb.max()<<"]";
    return os;
}




/* -------------------------------------------------------------------------- */
template <int dim, class nodes_container>
BoundingBox<dim> createPointList(const nodes_container& nodes, const Array<Real>& coord) {
    
    //  AKANTU_DEBUG_ASSERT(nodes.getSize() != 0, "No nodes to create a bounding box with.");
    typedef typename nodes_container::const_iterator node_iterator;
    
    node_iterator it = nodes.begin();
    assert(it != nodes.end());
    
    BoundingBox<dim> bbox(Point<dim>(&coord(*it),coord.getNbComponent()));
    for (++it; it != nodes.end(); ++it) {
        Real * point_coord = &coord(*it);
        for (UInt d=0; d<coord.getNbComponent(); ++d) {
            ;
        }
        bbox += *it;
    }
    return bbox;
}

__END_AKANTU__
