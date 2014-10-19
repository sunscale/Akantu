/**
 * @file   element_group.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Fri Sep 05 2014
 *
 * @brief  Stores information relevent to the notion of domain boundary and surfaces.
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
#include <sstream>
#include <algorithm>
#include <iterator>
#include "mesh.hh"
#include "group_manager_inline_impl.cc"
#include "dumpable_inline_impl.hh"
#include "aka_csr.hh"
#include "mesh_utils.hh"

#include "element_group.hh"
#if defined(AKANTU_USE_IOHELPER)
#  include "dumper_paraview.hh"
#endif
__BEGIN_AKANTU__
/* -------------------------------------------------------------------------- */
ElementGroup::ElementGroup(const std::string & group_name,
                           const Mesh & mesh,
                           NodeGroup & node_group,
                           UInt dimension,
                           const std::string & id,
                           const MemoryID & mem_id) :
  Memory(id, mem_id),
  mesh(mesh),
  name(group_name),
  elements("elements", id, mem_id),
  node_group(node_group),
  dimension(dimension) {
  AKANTU_DEBUG_IN();

#if defined(AKANTU_USE_IOHELPER)
  this->registerDumper<DumperParaview>("paraview_"  + group_name, group_name, true);
  this->addDumpFilteredMesh(mesh, elements, node_group.getNodes(), dimension);
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ElementGroup::empty() {
  elements.free();
}

/* -------------------------------------------------------------------------- */
void ElementGroup::append(const ElementGroup & other_group) {
  AKANTU_DEBUG_IN();

  node_group.append(other_group.node_group);

  /// loop on all element types in all dimensions
  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType ghost_type = *gt;

    type_iterator it = other_group.firstType(_all_dimensions, 
					     ghost_type,
					     _ek_not_defined);
    type_iterator last = other_group.lastType(_all_dimensions, 
					      ghost_type,
					      _ek_not_defined);

    for (; it != last; ++it) {
      ElementType type = *it;
      const Array<UInt> & other_elem_list = other_group.elements(type, ghost_type);
      UInt nb_other_elem = other_elem_list.getSize();

      Array<UInt> * elem_list;
      UInt nb_elem = 0;

      /// create current type if doesn't exists, otherwise get information
      if (elements.exists(type, ghost_type)) {
	elem_list = &elements(type, ghost_type);
	nb_elem = elem_list->getSize();
      }
      else {
	elem_list = &(elements.alloc(0, 1, type, ghost_type));
      }

      /// append new elements to current list
      elem_list->resize(nb_elem + nb_other_elem);
      std::copy(other_elem_list.begin(),
		other_elem_list.end(),
		elem_list->begin() + nb_elem);

      /// remove duplicates
      std::sort(elem_list->begin(), elem_list->end());
      Array<UInt>::iterator<> end = std::unique(elem_list->begin(), elem_list->end());
      elem_list->resize(end - elem_list->begin());
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ElementGroup::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "ElementGroup [" << std::endl;
  stream << space << " + name: "      << name << std::endl;
  stream << space << " + dimension: " << dimension << std::endl;
  elements.printself(stream, indent + 1);
  node_group.printself(stream, indent + 1);
  stream << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void ElementGroup::optimize() {
  // increasing the locality of data when iterating on the element of a group
  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;
    ElementList::type_iterator it = elements.firstType(_all_dimensions, ghost_type);
    ElementList::type_iterator last = elements.lastType(_all_dimensions, ghost_type);
    for (; it != last; ++it) {
      Array<UInt> & els = elements(*it, ghost_type);
      std::sort(els.begin(), els.end());

      Array<UInt>::iterator<> end = std::unique(els.begin(), els.end());
      els.resize(end - els.begin());
    }
  }

  node_group.optimize();
}

/* -------------------------------------------------------------------------- */




__END_AKANTU__

