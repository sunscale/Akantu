/**
 * @file   surface.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Wed Mar 13 2013
 *
 * @brief  contact surface classes
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

#ifndef __AKANTU_CSURFACE_HH__
#define __AKANTU_CSURFACE_HH__

#include <set>

#include "aka_common.hh"
#include "aka_bounding_box.hh"
#include "contact_common.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__


template <int d, class model_type>
BoundingBox<d> getBoundingBox(const Element& el, model_type& model) {
  
  typedef BoundingBox<d> bbox_type;
  typedef typename BoundingBox<d>::point_type point_type;
  
  Mesh &mesh = model.getMesh();
  const Array<Real> &x = model.getCurrentPosition();
  
  // for each element, construct the bounding box
  bbox_type elem_bb = bbox_type();
  
  UInt nb_nodes = mesh.getNbNodesPerElement(el.type);
  const Array<UInt> &conn = mesh.getConnectivity(el.type);
  
  for (UInt n = 0; n<nb_nodes; ++n) {
    UInt node = conn(el.element, n);
    elem_bb += point_type(&x(node));
  }
  return elem_bb;
}



//! This class represents a contact element surface
template <int d>
class CSurface {
  
  struct Element_set_comparator {
    bool operator()(const Element& el1, const Element& el2) const
    { return el1.element < el2.element; }
  };

public:
  
  // type definitions
  
  typedef BoundingBox<d> bbox_type;
  typedef typename bbox_type::point_type point_type;
  typedef SolidMechanicsModel model_type;
  typedef std::set<Element, Element_set_comparator> element_container;
  typedef typename element_container::const_iterator element_iterator;
  
  typedef std::set<UInt> node_container;
  typedef typename node_container::const_iterator node_iterator;

private:
  
  element_container elements_;         //!< Surface sides
  node_container nodes_;               //!< Surface nodes
  bbox_type bbox_;                     //!< Bounding box
  const model_type &model_;            //!< Reference to model

public:
    
  CSurface(const model_type &model) : elements_(), bbox_(), model_(model) {}
  
  CSurface& operator=(const CSurface& cs) {
    
    if (this != &cs) {
      new (this) CSurface(cs);
    }
    return *this;
  }
    
  void update_bounding_box() {
    
    bbox_ = bbox_type();
    const Array<Real> &x = model_.getCurrentPosition();
    
    // loop over all nodes to find out bounding box
    for (node_iterator it = nodes_.begin(); it != nodes_.end(); ++it)
      bbox_ += point_type(&x(*it));
  }
  
  bbox_type const& bounding_box() const
  { return bbox_; }
  
  //! Add an element to the surface
  void add_element(ElementType type, UInt id)
  { elements_.insert(Element(type, id)); }
  
  void add_node(UInt n)
  { nodes_.insert(n); }
  
  //! Return the number of elements in the surface
  UInt size() const
  { return elements_.size(); }
  
  void intersects(const bbox_type& bb, std::set<const Element*>& intersected) const {
        
    // loop over elements
    for (element_iterator it = elements_.begin(); it != elements_.end(); ++it) {
      
      // for each element, construct the bounding box
      bbox_type elem_bb = getBoundingBox<d>(*it, model_);

      // check for bounding box intersection
      if (bb & elem_bb)
        intersected.insert(&*it);
      
    } // loop over elements
  }
  
  bool in_surface(UInt np, const Element* sp) const {

    element_iterator elit = elements_.find(*sp);
    if (elit == elements_.end())
      return false;
    
    node_iterator nit = nodes_.find(np);
    if (nit == nodes_.end())
      return false;

    return true;
  }
    
  //! Enable std output
  friend std::ostream& operator<<(std::ostream& os, const CSurface& s) {
    
    os<<"    Contact surface: "<<endl;
    os<<"      bounding box: "<<s.bounding_box()<<endl;
    os<<"      "<< s.elements_.size()<<" surface elements:"<<endl;
    for (typename CSurface::element_iterator it = s.elements_.begin(); it != s.elements_.end(); ++it)
      os<<" "<<it->element;
    os<<"\n      "<<s.nodes_.size()<<" surface nodes:"<<endl;
    for (typename CSurface::node_iterator it = s.nodes_.begin(); it != s.nodes_.end(); ++it)
      os<<" "<<*it;
    return os<<endl;
  }
};


__END_AKANTU__

#endif /* __AKANTU_CSURFACE_HH__ */
