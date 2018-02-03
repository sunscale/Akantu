/**
 * @file   node_filter.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  to filter nodes with functors
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AST_NODE_FILTER_HH__
#define __AST_NODE_FILTER_HH__

/* -------------------------------------------------------------------------- */
// akantu
#include "aka_common.hh"
#include "mesh_filter.hh"


namespace akantu {

/* -------------------------------------------------------------------------- */
class GeometryFilter : public NodeFilterFunctor {

public:
  GeometryFilter(const Mesh & mesh,
		 UInt dir,
		 Real limit) :
    NodeFilterFunctor(),
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
    bool to_filter = true;
    if((*this->positions)(node,this->dir) > this->limit)
      to_filter = false;
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
    bool to_filter = true;
    if((*this->positions)(node,this->dir) < this->limit)
      to_filter = false;
    AKANTU_DEBUG_OUT();
    return to_filter;
  };
};

/* -------------------------------------------------------------------------- */
// this filter is erase because the convention of filter has changed!!
// filter == true -> keep node

// template<class FilterType>
// void applyNodeFilter(Array<UInt> & nodes, FilterType & filter) {

//   Array<UInt>::iterator<> it = nodes.begin();

//   for (; it != nodes.end(); ++it) {
//     if (filter(*it)) {
//       it = nodes.erase(it);
//     }
//   }
// }; 

} // namespace akantu

#endif /* __AST_NODE_FILTER_HH__ */
