/**
 * @file   mesh_filter.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Mon Feb 10 2014
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  the class representing the meshes
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
#ifndef __AKANTU_MESH_FILTER_HH__
#define __AKANTU_MESH_FILTER_HH__

/* -------------------------------------------------------------------------- */
#include "element.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* Filter Fonctors                                                            */
/* -------------------------------------------------------------------------- */
struct FilterFunctor {
  enum Type {
    _node_filter_functor,
    _element_filter_functor
  };
};

class NodeFilterFunctor : public FilterFunctor {
public:
  bool operator()(UInt node) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
public:
  static const Type type = _node_filter_functor;
};

class ElementFilterFunctor : public FilterFunctor {
public:
  bool operator()(const Element & element) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
public:
  static const Type type = _element_filter_functor;
};

__END_AKANTU__

#endif /* __AKANTU_MESH_FILTER_HH__ */
