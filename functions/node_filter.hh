/**
 * @file   node_filter.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Thu Jun 20 11:30:10 2013
 *
 * @brief  to filter nodes with functors
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef __AST_NODE_FILTER_HH__
#define __AST_NODE_FILTER_HH__

/* -------------------------------------------------------------------------- */
// akantu
#include "aka_common.hh"

// simtools
#include "ast_common.hh"

__BEGIN_SIMTOOLS__

using namespace akantu;

/* -------------------------------------------------------------------------- */
class GeometryFilter {

public:
  GeometryFilter(const Mesh & mesh,
		 UInt dir,
		 Real limit) :
    mesh(mesh), dir(dir), limit(limit) { 
    this->positions = &(mesh.getNodes());
  };
  ~GeometryFilter() {};

  bool operator() (UInt node) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };

protected:
  const Mesh & mesh;
  UInt dir;
  Real limit;
  const Array<Real> * positions;
};

/* -------------------------------------------------------------------------- */
class FilterPositionsGreaterThan : public GeometryFilter {
public:
  FilterPositionsGreaterThan(const Mesh & mesh,
			     UInt dir,
			     Real limit) : GeometryFilter(mesh, dir, limit) {};
  ~FilterPositionsGreaterThan() {};

  bool operator() (UInt node) {
    AKANTU_DEBUG_IN();
    bool to_filter = false;
    if((*this->positions)(node,this->dir) > this->limit)
      to_filter = true;
    AKANTU_DEBUG_OUT();
    return to_filter;
  };
};

/* -------------------------------------------------------------------------- */
class FilterPositionsLessThan : public GeometryFilter {
public:
  FilterPositionsLessThan(const Mesh & mesh,
			     UInt dir,
			     Real limit) : GeometryFilter(mesh, dir, limit) {};
  ~FilterPositionsLessThan() {};

  bool operator() (UInt node) {
    AKANTU_DEBUG_IN();
    bool to_filter = false;
    if((*this->positions)(node,this->dir) < this->limit)
      to_filter = true;
    AKANTU_DEBUG_OUT();
    return to_filter;
  };
};

/* -------------------------------------------------------------------------- */
template<class FilterType>
void applyNodeFilter(Array<UInt> & nodes, FilterType & filter) {

  Array<UInt>::iterator<> it = nodes.begin();

  for (; it != nodes.end(); ++it) {
    if (filter(*it)) {
      it = nodes.erase(it);
    }
  }
}; 

__END_SIMTOOLS__

#endif /* __AST_NODE_FILTER_HH__ */
