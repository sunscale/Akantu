/**
 * @file   node_group.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Mon Jun 09 2014
 *
 * @brief  Node group definition
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
#include "aka_common.hh"
#include "aka_array.hh"
#include "aka_memory.hh"
#include "mesh_filter.hh"
/* -------------------------------------------------------------------------- */



#ifndef __AKANTU_NODE_GROUP_HH__
#define __AKANTU_NODE_GROUP_HH__

__BEGIN_AKANTU__

class NodeGroup : Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  NodeGroup(const std::string & name,
            const std::string & id = "node_group",
            const MemoryID & memory_id = 0);
  virtual ~NodeGroup();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  typedef Array<UInt>::const_iterator<UInt> const_node_iterator;

  /// empty the node group
  void empty();

  /// iterator to the beginning of the node group
  inline const_node_iterator begin() const;
  /// iterator to the end of the node group
  inline const_node_iterator end() const;

  /// add a node and give the local position through an iterator
  inline const_node_iterator add(UInt node, bool check_for_duplicate = true);

  /// remove duplicated nodes
  void optimize();

  /// append a group to current one
  void append(const NodeGroup & other_group);

  /// apply a filter on current node group
  template <typename T>
  void applyNodeFilter(T & filter);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO_NOT_CONST(Nodes, node_group, Array<UInt> &);
  AKANTU_GET_MACRO(Nodes, node_group, const Array<UInt> &);
  AKANTU_GET_MACRO(Name, name, const std::string &);

  /// give the number of nodes in the current group
  inline UInt getSize() const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// name of the group
  std::string name;

  /// list of nodes in the group
  Array<UInt> & node_group;
};

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const NodeGroup & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "node_group_inline_impl.cc"



#endif /* __AKANTU_NODE_GROUP_HH__ */
