/**
 * @file   zone.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Tue May 13 2014
 *
 * @brief  contact zone classes
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

#ifndef __AKANTU_CZONE_HH__
#define __AKANTU_CZONE_HH__

#include <unordered_set>

#include "aka_common.hh"
#include "surface.hh"

__BEGIN_AKANTU__


//! This class represents a contact element surface
template <int d>
class CZone {
  
public:
  
  typedef CSurface<d> contact_surface;
  typedef typename std::list<const contact_surface*> surface_list;
  typedef typename surface_list::const_iterator surface_iterator;
  
  typedef typename CSurface<d>::point_type point_type;
  typedef typename CSurface<d>::bbox_type bbox_type;
  
//  typedef std::set<UInt> node_set;
  typedef std::unordered_set<UInt> node_set;
  typedef typename node_set::iterator node_iterator;
  typedef typename node_set::const_iterator const_node_iterator;
  typedef std::map<int, node_set> bucket_container;
  typedef typename bucket_container::iterator bucket_iterator;
  typedef typename bucket_container::const_iterator const_bucket_iterator;
  
  typedef std::set<const Element*> element_container;
  typedef typename element_container::const_iterator element_iterator;
  //    typedef std::map<UInt, element_container> node_element_map;
  //    typedef typename node_element_map::iterator node_element_iterator;
  
  typedef std::map<int, element_container> elementmap;
  typedef typename elementmap::iterator elementmap_iterator;
  typedef typename elementmap::const_iterator const_elementmap_iterator;
  
private:
  
  SolidMechanicsModel& model_;
  bbox_type bbox_;
  Real delta_[d];               //!< grid size length
  int size_[d];                //!< grid dimensions, number of buckets in each dimension
  int index_[d];                //!< indexes to convert to 1D equivalent array
  bucket_container buckets_;
  elementmap elements_;
  //    node_element_map node_element_map_;
  surface_list surfaces_;
  
public:
  
  CZone() {}
  
  // clear memory
  void clear() {
    buckets_.clear();
    //      node_element_map_.clear();
  }

  typename bucket_container::size_type buckets_size() const 
  { return buckets_.size(); }

  bucket_iterator buckets_begin()
  { return buckets_.begin(); }
  
  bucket_iterator buckets_end()
  { return buckets_.end(); }
  
  const_bucket_iterator buckets_begin() const
  { return buckets_.begin(); }
  
  const_bucket_iterator buckets_end() const
  { return buckets_.end(); }
  
  void erase_bucket(bucket_iterator bit)
  { buckets_.erase(bit); }
  
  CZone(SolidMechanicsModel& model, 
        const bbox_type& bb, 
        const element_container& elems,
        const contact_surface& s1,
        const contact_surface& s2)
  : model_(model), bbox_(bb), delta_(), size_(), index_() {
    
    assert(!elems.empty());
    
    surfaces_.push_back(&s1);
    surfaces_.push_back(&s2);
    
    // tolerance
    Real epsilon = 2*std::numeric_limits<Real>::epsilon();
    
    // loop over elements to get bucket increments
    typedef typename element_container::const_iterator element_iterator;
    for (element_iterator it = elems.begin(); it != elems.end(); ++it) {
      
      // get side bounding box
      bbox_type bb = getBoundingBox<d>(**it, model_);
      
      // process box to modify delta array
      const point_type &m = bb.min();
      const point_type &M = bb.max();
      for (int i=0; i<d; ++i)
        delta_[i] = std::max(delta_[i], (M[i] - m[i])/3. + epsilon);
    }
    
    // get the number of buckets in each direction
    const point_type& m = bbox_.min();
    const point_type& M = bbox_.max();
    for (int i=0; i<d; ++i) {
      assert(delta_[i] != 0.);
      size_[i] = static_cast<int>((M[i] - m[i]) / delta_[i])+1;
    }
    
    // get indexes to convert multi-dimensional indexes to a 1D array
    index_[0] = 1;
    for (int i=1; i<d; ++i)
      index_[i] = index_[i-1]*size_[i-1];

    // add elements to buckets
    add_elements(elems);
    
    std::set<UInt> nodes_inside;
    
    Mesh &mesh = model.getMesh();
    const Array<Real> &x = model.getCurrentPosition();
    
    // loop over elements
    for (element_iterator it = elems.begin(); it != elems.end(); ++it) {
      
      UInt nb_nodes = mesh.getNbNodesPerElement((*it)->type);
      const Array<UInt> &conn = mesh.getConnectivity((*it)->type);
      
      // loop over nodes
      for (UInt n = 0; n<nb_nodes; ++n) {
        
        UInt node = conn((*it)->element, n);
        point_type node_coord = point_type(&x(node));
        
        // if node is within bounding box
        if (bbox_ & node_coord) {
          
          nodes_inside.insert(node);
          
          // compute index into one-dimensional array
          int coord[d];
          for (int i=0; i<d; ++i)
            coord[i] = static_cast<int>((node_coord[i] - m[i])/delta_[i]);
          int idx = compute_index(coord);
          
          // add node to bucket
          buckets_[idx].insert(node);
          
          //            // add side to map
          //            node_element_map_[node].insert(*it);
        } // node within bounding box
      } // loop over nodes
    } // loop over elements
    
    // check nodes
    std::set<UInt> node_check;
    for (bucket_iterator it = buckets_.begin(); it != buckets_.end(); ++it)
      for (node_iterator nit = it->second.begin(); nit != it->second.end(); ++nit)
        node_check.insert(*nit);
    cout<<"*** INFO *** Node check after creating buckets passed: "<<(node_check.size() == nodes_inside.size())<<endl;
    assert(node_check.size() == nodes_inside.size());
    
#ifdef DEBUG_CONTACT
    cout<<"*** INFO *** A total of "<<nodes_inside.size()<<" nodes lie inside the contact zone."<<endl;
    cout<<"*** INFO *** A total of "<<buckets_.size()<<" node buckets were created."<<endl;
#endif

  }
  
  int compute_index(int *array) const {
    int idx = 0;
    for (int i=0; i<d; ++i)
      idx += index_[i]*array[i];
    return idx;
  }
  
  void decompute_index(int idx, int *array, Int2Type<2>) const {
    array[0] = idx % index_[1];
    array[1] = idx / index_[1];
  }
  
  void decompute_index(int idx, int *array, Int2Type<3>) const {
    array[0] = (idx % (index_[2])) % index_[1];
    array[1] = (idx % (index_[2])) / index_[1];
    array[2] = idx / index_[2];
  }
  
  
  void add_elements(const element_container& s) {
    
    const point_type& bbm = bbox_.min();
    
    // loop over elements
    for (element_iterator it = s.begin(); it != s.end(); ++it) {
      
      bbox_type elbb = getBoundingBox<d>(**it, model_);
      const point_type& elmin = elbb.min();
      const point_type& elmax = elbb.max();
      
      int min[d], max[d];
      for (int i=0; i<d; ++i) {
        min[i] = std::max(0, static_cast<int>((elmin[i] - bbm[i])/delta_[i]));
        max[i] = std::min(size_[i], static_cast<int>((elmax[i] - bbm[i])/delta_[i]) +1);
      }
      
      add_element(*it, min, max, Int2Type<d>());
      
    } // loop over elements
  }

  void add_element(const Element* el, int *min, int*max, Int2Type<2>) {
    
    for (int i=min[0]; i<max[0]; ++i) {
      for (int j=min[1]; j<max[1]; ++j) {
          int idxarray[d] = {i,j};
          int idx = compute_index(idxarray);
          elements_[idx].insert(el);
      }
    }
  }
  
  void add_element(const Element* el, int *min, int*max, Int2Type<3>) {
    
    for (int i=min[0]; i<max[0]; ++i) {
      for (int j=min[1]; j<max[1]; ++j) {
        for (int k=min[2]; k<max[2]; ++k) {
          int idxarray[d] = {i,j,k};
          int idx = compute_index(idxarray);
          elements_[idx].insert(el);
        }
      }
    }
  }
  
  element_container contiguous(int bucket) {
    
    element_container elements;
    collect(elements, bucket, Int2Type<d>());
    return elements;
  }
  
  bool in_surface(UInt np, const Element* sp) {
    
    for (surface_iterator it = surfaces_.begin(); it != surfaces_.end(); ++it) {
      
      const contact_surface& s = **it;
      if (s.in_surface(np, sp))
        return true;
    }
    return false;
  }
  
  bool in_element(UInt np, const Element* sp) {
    
    // loop over nodes of side
    Mesh& mesh = model_.getMesh();
    UInt nb_nodes = mesh.getNbNodesPerElement(sp->type);
    const Array<UInt> &conn = mesh.getConnectivity(sp->type);
    
    // loop over nodes
    for (UInt n = 0; n<nb_nodes; ++n) {
      UInt node = conn(sp->element, n);
      if (node == np)
        return true;
    }
    return false;
  }
  
  void set_union(int bucket, element_container& elements) const {
    
    // check if bucket exists
    const_elementmap_iterator it = elements_.find(bucket);
    
    // if bucket is not found, return since there are no elements to add
    if (it == elements_.end())
      return;
    
    // otherwise add elements
    for (element_iterator sit = it->second.begin(); sit != it->second.end(); ++sit)
      elements.insert(*sit);
  }
  
  
  
  template <class container_type>
  void collect(container_type &c, int idx, Int2Type<2>) const {
    
    const int &m = size_[0];
    const int &n = size_[1];
    
    int coord[d];
    decompute_index(idx, coord, Int2Type<2>());
    
    for (int i=coord[0]-1; i<=coord[0]+1; ++i)
      if (i >=0 && i<m)
        for (int j=coord[1]-1; j<=coord[1]+1; ++j)
          if (j >=0 && j<n) {
            int cont[d] = { i, j};
            set_union(compute_index(cont), c);
          }
  }
  
  template <class container_type>
  void collect(container_type &c, int idx, Int2Type<3>) const {
    
    const int &m = size_[0];
    const int &n = size_[1];
    const int &p = size_[2];
    
    int coord[d];
    decompute_index(idx, coord, Int2Type<3>());
    
    for (int i=coord[0]-1; i<=coord[0]+1; ++i)
      if (i >=0 && i<m)
        for (int j=coord[1]-1; j<=coord[1]+1; ++j)
          if (j >=0 && j<n)
            for (int k=coord[2]-1; k<=coord[2]+1; ++k)
              if (k >=0 && k<p) {
                int cont[d] = { i, j, k};
                set_union(compute_index(cont), c);
              }
  }
  
  //! Enable std output
  friend std::ostream& operator<<(std::ostream& os, const CZone& cz) {
    
    Mesh& mesh = cz.model_.getMesh();
    const Array<Real> &x = cz.model_.getCurrentPosition();
    
    int coord[d];
    
    os<<"    Contact zone: "<<endl;
    os<<"      origin: "<<cz.bbox_.min()<<endl;
    os<<"      bounding box: "<<cz.bbox_<<endl;
    os<<"      delta: (";
    for (int i=0; i<d-1; ++i)
      os<<cz.delta_[i]<<",";
    os<<cz.delta_[d-1]<<")"<<endl;
    os<<"      size: (";
    for (int i=0; i<d-1; ++i)
      os<<cz.size_[i]<<",";
    os<<cz.size_[d-1]<<")"<<endl;
    os<<"        number of node buckets: "<<cz.buckets_.size()<<endl;
    std::set<UInt> nodes;
    for (const_bucket_iterator it = cz.buckets_.begin(); it != cz.buckets_.end(); ++it) {
      os<<"          node bucket "<<it->first;
      
      cz.decompute_index(it->first, coord, Int2Type<d>());
      os<<" ("<<coord[0];
      for (int i=1; i<d; ++i)
        os<<","<<coord[i];
      os<<") contains "<<it->second.size()<<" nodes:";
      
      for (const_node_iterator nit = it->second.begin(); nit != it->second.end(); ++nit) {
        os<<" "<<*nit<<" ("<<point_type(&x(*nit))<<")";
        nodes.insert(*nit);
      }
      os<<endl;
    }
    os<<"        there are a total of "<<nodes.size()<<" nodes in the node buckets."<<endl;
    os<<"        number of element buckets: "<<cz.elements_.size()<<endl;
    
    for (const_elementmap_iterator it = cz.elements_.begin(); it != cz.elements_.end(); ++it) {
      os<<"          element bucket "<<it->first;
      
      
      cz.decompute_index(it->first, coord, Int2Type<d>());
      os<<" ("<<coord[0];
      for (int i=1; i<d; ++i)
        os<<","<<coord[i];
      os<<") contains "<<it->second.size()<<" elements:"<<endl;
      
      for (element_iterator sit = it->second.begin(); sit != it->second.end(); ++sit) {
        os<<"\t\t"<<**sit;
        
        UInt nb_nodes = mesh.getNbNodesPerElement((*sit)->type);
        const Array<UInt> &conn = mesh.getConnectivity((*sit)->type);
        
        os<<", nodes:";
        
        // loop over nodes
        for (UInt n = 0; n<nb_nodes; ++n)
          os<<" "<<conn((*sit)->element, n)<<" ("<<
          point_type(&x(conn((*sit)->element, n)))<<")";
        os<<endl;
      }
    }
    return os;
  }
};


__END_AKANTU__

#endif /* __AKANTU_CZONE_HH__ */
