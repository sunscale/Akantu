/**
 * @file   mesh_filter.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  the class representing the meshes
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef __AKANTU_MESH_FILTER_HH__
#define __AKANTU_MESH_FILTER_HH__

/* -------------------------------------------------------------------------- */
#include "element.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Filter Functors                                                            */
/* -------------------------------------------------------------------------- */

/// struct for the possible filter functors
struct FilterFunctor {
  enum Type { _node_filter_functor, _element_filter_functor };
};

/// class (functor) for the node filter
class NodeFilterFunctor : public FilterFunctor {
public:
  bool operator()(__attribute__((unused)) UInt node) { AKANTU_TO_IMPLEMENT(); }

public:
  static const Type type = _node_filter_functor;
};

/// class (functor) for the element filter
class ElementFilterFunctor : public FilterFunctor {
public:
  bool operator()(__attribute__((unused)) const Element & element) {
    AKANTU_TO_IMPLEMENT();
  }

public:
  static const Type type = _element_filter_functor;
};

} // namespace akantu

#endif /* __AKANTU_MESH_FILTER_HH__ */
