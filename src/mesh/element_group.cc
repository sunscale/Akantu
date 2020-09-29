/**
 * @file   element_group.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Mon Jan 22 2018
 *
 * @brief  Stores information relevent to the notion of domain boundary and
 * surfaces.
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_csr.hh"
#include "dumpable.hh"
#include "dumpable_inline_impl.hh"
#include "group_manager.hh"
#include "group_manager_inline_impl.hh"
#include "mesh.hh"
#include "mesh_utils.hh"
#include <algorithm>
#include <iterator>
#include <sstream>

#include "element_group.hh"
#if defined(AKANTU_USE_IOHELPER)
#include "dumper_iohelper_paraview.hh"
#endif

namespace akantu {

/* -------------------------------------------------------------------------- */
ElementGroup::ElementGroup(const std::string & group_name, const Mesh & mesh,
                           NodeGroup & node_group, UInt dimension,
                           const std::string & id, const MemoryID & mem_id)
    : Memory(id, mem_id), mesh(mesh), name(group_name),
      elements("elements", id, mem_id), node_group(node_group),
      dimension(dimension) {
  AKANTU_DEBUG_IN();

#if defined(AKANTU_USE_IOHELPER)
  this->registerDumper<DumperParaview>("paraview_" + group_name, group_name,
                                       true);
  this->addDumpFilteredMesh(mesh, elements, node_group.getNodes(),
                            _all_dimensions);
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ElementGroup::ElementGroup(const ElementGroup & /*other*/) = default;

/* -------------------------------------------------------------------------- */
void ElementGroup::clear() { elements.free(); }

/* -------------------------------------------------------------------------- */
bool ElementGroup::empty() const { return elements.empty(); }

/* -------------------------------------------------------------------------- */
void ElementGroup::append(const ElementGroup & other_group) {
  AKANTU_DEBUG_IN();

  node_group.append(other_group.node_group);

  /// loop on all element types in all dimensions
  for (auto ghost_type : ghost_types) {
    for (auto type : other_group.elementTypes(_ghost_type = ghost_type,
                     _element_kind = _ek_not_defined)) {
      const Array<UInt> & other_elem_list =
          other_group.elements(type, ghost_type);
      UInt nb_other_elem = other_elem_list.size();

      Array<UInt> * elem_list;
      UInt nb_elem = 0;

      /// create current type if doesn't exists, otherwise get information
      if (elements.exists(type, ghost_type)) {
        elem_list = &elements(type, ghost_type);
        nb_elem = elem_list->size();
      } else {
        elem_list = &(elements.alloc(0, 1, type, ghost_type));
      }

      /// append new elements to current list
      elem_list->resize(nb_elem + nb_other_elem);
      std::copy(other_elem_list.begin(), other_elem_list.end(),
                elem_list->begin() + nb_elem);

      /// remove duplicates
      std::sort(elem_list->begin(), elem_list->end());
      Array<UInt>::iterator<> end =
          std::unique(elem_list->begin(), elem_list->end());
      elem_list->resize(end - elem_list->begin());
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ElementGroup::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT) {
    ;
  }

  stream << space << "ElementGroup [" << std::endl;
  stream << space << " + name: " << name << std::endl;
  stream << space << " + dimension: " << dimension << std::endl;
  elements.printself(stream, indent + 1);
  node_group.printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void ElementGroup::optimize() {
  // increasing the locality of data when iterating on the element of a group
  for (auto ghost_type : ghost_types) {
    for (auto type : elements.elementTypes(_ghost_type = ghost_type)) {
      Array<UInt> & els = elements(type, ghost_type);
      std::sort(els.begin(), els.end());

      Array<UInt>::iterator<> end = std::unique(els.begin(), els.end());
      els.resize(end - els.begin());
    }
  }

  node_group.optimize();
}

/* -------------------------------------------------------------------------- */
void ElementGroup::fillFromNodeGroup() {
  CSR<Element> node_to_elem;
  MeshUtils::buildNode2Elements(this->mesh, node_to_elem, this->dimension);

  std::set<Element> seen;

  Array<UInt>::const_iterator<> itn = this->node_group.begin();
  Array<UInt>::const_iterator<> endn = this->node_group.end();
  for (; itn != endn; ++itn) {
    CSR<Element>::iterator ite = node_to_elem.begin(*itn);
    CSR<Element>::iterator ende = node_to_elem.end(*itn);
    for (; ite != ende; ++ite) {
      const Element & elem = *ite;
      if (this->dimension != _all_dimensions &&
          this->dimension != Mesh::getSpatialDimension(elem.type)) {
        continue;
      }
      if (seen.find(elem) != seen.end()) {
        continue;
      }

      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(elem.type);
      Array<UInt>::const_iterator<Vector<UInt>> conn_it =
          this->mesh.getConnectivity(elem.type, elem.ghost_type)
              .begin(nb_nodes_per_element);
      const Vector<UInt> & conn = conn_it[elem.element];

      UInt count = 0;
      for (UInt n = 0; n < conn.size(); ++n) {
        count +=
            (this->node_group.getNodes().find(conn(n)) != UInt(-1) ? 1 : 0);
      }

      if (count == nb_nodes_per_element) {
        this->add(elem);
      }

      seen.insert(elem);
    }
  }

  this->optimize();
}

/* -------------------------------------------------------------------------- */
void ElementGroup::addDimension(UInt dimension) {
  this->dimension = std::max(dimension, this->dimension);
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
