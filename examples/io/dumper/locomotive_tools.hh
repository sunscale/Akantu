/**
 * @file   locomotive_tools.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Apr 13 2012
 * @date last modification: Mon Aug 17 2015
 *
 * @brief  interface for the common tools
 *
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
void applyRotation(const ::akantu::Vector<::akantu::Real> & center,
                   ::akantu::Real angle,
                   const ::akantu::Array<::akantu::Real> & nodes,
                   ::akantu::Array<::akantu::Real> & displacement,
                   const ::akantu::Array<::akantu::UInt> & node_group);

void fillColour(const ::akantu::Mesh & mesh,
                ::akantu::ElementTypeMapArray<::akantu::UInt> & colour);
