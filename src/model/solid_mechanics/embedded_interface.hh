/**
 * @file   embedded_interface.hh
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

#ifndef __AKANTU_EMBEDDED_INTERFACE_HH__
#define __AKANTU_EMBEDDED_INTERFACE_HH__

#include "aka_common.hh"
#include "parsable.hh"

#include <CGAL/Cartesian.h>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/**
 * @brief Parsable class for the embedded interfaces
 */
class EmbeddedInterface : public Parsable {

  typedef CGAL::Cartesian<Real> K;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  EmbeddedInterface(UInt dim, const ID & id = "");

  virtual ~EmbeddedInterface();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Constructs the CGAL primitive segment
  K::Segment_3 getPrimitive() const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Name, name, const std::string &);
  AKANTU_GET_MACRO(MaterialName, mat_name, const std::string &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Intetrface name
  std::string name;

  /// Material name
  std::string mat_name;

  /// Start and end points of interface
  Matrix<Real> points;

};

__END_AKANTU__

#endif // __AKANTU_EMBEDDED_INTERFACE_HH__
