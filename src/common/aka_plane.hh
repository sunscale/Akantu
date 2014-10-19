/**
 * @file   aka_plane.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Tue Jun 17 2014
 *
 * @brief  Geometry class representing planes
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

#ifndef __AKANTU_AKA_PLANE_HH__
#define __AKANTU_AKA_PLANE_HH__

#include "aka_point.hh"


__BEGIN_AKANTU__



//! Plane class
struct Plane {
  
  typedef Point<3> point_type;
  typedef point_type::value_type value_type;
  
  //! Parameter constructor creates a plane out of 3 points
  Plane(const point_type& a, const point_type& b, const point_type& c) {
    
    n_ = cross(b - a, c - a).normalize();
    d_ = n_*a;
  }
  
  //! Return the plane normal
  const point_type& normal() const
  { return n_; }
  
  //! Return the disctance to the origin
  value_type distance() const
  { return d_; }
  
  //! Standard output stream operator
  friend std::ostream& operator<<(std::ostream& os, const Plane& pi) {
    
    os<<"Plane[normal: "<<pi.normal()<<", origin distance: "<<pi.distance()<<"]"<<std::endl;
    return os;
  }
  
private:
  
  point_type n_;     //!< Plane unit normal
  value_type d_;     //!< Distance from the plane to the origin
  
};


__END_AKANTU__

#endif /* __AKANTU_AKA_PLANE_HH__ */
