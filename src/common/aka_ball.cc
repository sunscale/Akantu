/**
 * @file   aka_ball.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Tue May 07 2013
 * @date last modification: Tue May 07 2013
 *
 * @brief  bounding sphere classes
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

#include "aka_ball.hh"

__BEGIN_AKANTU__


template <>
std::ostream& Interval::print(std::ostream& os) const {
  os<<"Interval["<<c_<<", "<<r_<<"]";
  return os;
}


template <>
std::ostream& Circle::print(std::ostream& os) const {
  os<<"Disk["<<c_<<", "<<r_<<"]";
  return os;
}

template <>
std::ostream& Sphere::print(std::ostream& os) const {
  os<<"Sphere["<<c_<<", "<<r_<<"]";
  return os;
}



template <>
typename Ball<1>::value_type Ball<1>::measure() const
{ return r_; }

template <>
typename Ball<2>::value_type Ball<2>::measure() const
{ return pow(r_,2); }

template <>
typename Ball<3>::value_type Ball<3>::measure() const
{ return pow(r_,3); }



__END_AKANTU__
