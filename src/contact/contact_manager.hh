/**
 * @file   contact_manager.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Wed Mar 13 2013
 *
 * @brief  contact manager
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

#ifndef __AKANTU_CMANAGER_HH__
#define __AKANTU_CMANAGER_HH__

#include "aka_common.hh"
#include "surface.hh"
#include "zone.hh"
#include "element.hh"


__BEGIN_AKANTU__


template <int d>
class BBox_sorter {
  
  typedef CSurface<d> surface_type;
  
  struct Comparator {
    
    template <class surface_pointer>
    bool operator()(surface_pointer s1, surface_pointer s2) const {
      
      typedef typename surface_type::bbox_type bbox_type;
      
      auto b1 = s1->bounding_box();
      auto b2 = s2->bounding_box();
      
      // consider min point for sorting (chosen arbitrarily)
      return b1.min(d-1) < b2.min(d-1);
      // NOTE: point coordinates are 0-based, so calling d-1 with
      // dimension 2 actually gets the second element of the
      // coordinate array
    }
  };
  
public:
  
  template <class pointer>
  static std::pair<pointer, pointer> check(pointer p1, pointer p2) {
    typedef std::pair<pointer, pointer> contact_pair;
    return p1 < p2 ? contact_pair(p1,p2) : contact_pair(p2,p1);
  }
  
  template <class surface_container, class intersection_container>
  static void sort(surface_container& surfaces, intersection_container& intersections) {
    
    typedef typename surface_container::value_type surface_ptr;
    typedef typename surface_type::bbox_type bbox_type;
    
    typedef typename surface_container::iterator surface_iterator;
    typedef typename intersection_container::value_type contact_pair;
    
    // sort container
    surfaces.sort(Comparator());
    
    // loop over sorted surfaces
    surface_iterator it1 = surfaces.begin();
    surface_iterator it2 = ++surfaces.begin();
    for (; it2 != surfaces.end(); ++it2) {
      
      // check for intersection of bounding boxes
      auto b1 = (*it1)->bounding_box();
      auto b2 = (*it2)->bounding_box();
      
      // if bounding boxes intersect
      if (b1 & b2) {
        // add pair to set
        intersections.insert(check(*it1,*it2));
        ++it1;
        
      } else {
        
        // remove pair from set
        intersections.erase(check(*it1,*it2));
        
        // remove first surface from container, as it does not
        // intersect other surfaces
        surfaces.erase(it1++);
      }
    }
    
    if (surfaces.size() > 1)
      BBox_sorter<d-1>::sort(surfaces, intersections);
  }
  
};

//! Partial template specialization for zero dimension
/*! This class finishes the recursion on dimension by doing nothing in the
 * sort function.
 */
template <>
struct BBox_sorter<0> {
  
  template <class surface_container, class intersection_container>
  static void sort(surface_container&, intersection_container&) {}
};

template <int d>
class CManager {
  
  typedef SolidMechanicsModel model_type;
  typedef CSurface<d> contact_surface;
  typedef std::vector<contact_surface> surface_container;
  typedef typename surface_container::iterator surface_iterator;
  typedef typename surface_container::const_iterator const_surface_iterator;
  
  typedef CZone<d> contact_zone;
  typedef std::list<contact_zone*> zone_container;
  typedef typename zone_container::iterator zone_iterator;
  typedef typename zone_container::const_iterator const_zone_iterator;
  
  typedef CElement<d> contact_element;
  typedef std::list<contact_element> celement_container;
  typedef typename celement_container::const_iterator celement_iterator;
  
  SolidMechanicsModel &model_;
  surface_container surfaces_;
  zone_container zones_;
  celement_container celements_;
  Contact_type c_;
  
public:
  
  void clear() {
    
    zones_.clear();
    celements_.clear();
  }
  
  CManager(model_type& model, Contact_type c = No_self_contact_t) : model_(model), c_(c) {
    
    // get mesh from model
    Mesh &mesh = model.getMesh();
    
    // call update current position to be able to call later
    // the function to get current positions
    model.updateCurrentPosition();
    
    // obtain exterior surfaces
    MeshUtils::buildFacets(mesh);
    
    // assign surface ids
    MeshUtils::buildSurfaceID(mesh);
    
    // allocate memory for surfaces
    UInt nb_surfaces = mesh.getNbSurfaces();
    surfaces_.reserve(nb_surfaces);
    
    for (UInt i=0; i<nb_surfaces; ++i)
      surfaces_.push_back(contact_surface(model));
    
    // iterate over elements of lower dimension
    Mesh::type_iterator it = mesh.firstType(d-1);
    Mesh::type_iterator end = mesh.lastType(d-1);
    
    for(; it != end; ++it) {
      
      UInt nb_element = mesh.getNbElement(*it);
      UInt nb_nodes = mesh.getNbNodesPerElement(*it);
      const Array<UInt> &conn = mesh.getConnectivity(*it);
      
      Array<UInt> &surf_id = mesh.getSurfaceID(*it);
      
      // add elements to corresponding surface
      for(UInt e = 0; e < nb_element; ++e) {
        
        CSurface<d> &surface = surfaces_.at(surf_id(e));
        surface.add_element(*it,e);
        
        // add element nodes to surface
        for (UInt n = 0; n<nb_nodes; ++n)
          surface.add_node(conn(e, n));
      }
    }
  }
  
  ~CManager() {
    
    // delete over contact zones
    for (zone_iterator it = zones_.begin(); it != zones_.end(); ++it)
      delete *it;
  }
  
  void global_search() {
    
#ifdef DEBUG_CONTACT
    cout<<"__________________________________________________________"<<endl;
    cout<<"*** INFO *** Printing contact manager after global search."<<endl;
#endif
    
    typedef const contact_surface* surface_ptr;
    typedef std::list<surface_ptr> surface_list;
    typedef std::pair<surface_ptr, surface_ptr> contact_pair;
    typedef std::set<contact_pair> intersection_container;
    typedef typename intersection_container::iterator intersection_iterator;
    typedef typename contact_surface::bbox_type bbox_type;
    
    // create container of bounding boxes used for the sort and
    // the container to store intersections
    surface_list list;
    intersection_container intersections;
    
    // loop over surfaces to update bounding boxes
    for (surface_iterator it = surfaces_.begin(); it != surfaces_.end(); ++it) {
      
      const contact_surface &surface = *it;
      it->update_bounding_box();
      list.push_back(&surface);
    }
    
    // carry out sort in all dimensions to check for intersections
    BBox_sorter<d>::sort(list, intersections);
    
    if (!intersections.empty()) {
      
      // loop over intersections to find contact zones
      for (intersection_iterator it = intersections.begin();
           it != intersections.end(); ++it) {
        
        // get contact surfaces that intersect
        const contact_surface &surface1 = *it->first;
        const contact_surface &surface2 = *it->second;
        
        // find intersection bounding box
        bbox_type bb = surface1.bounding_box() && surface2.bounding_box();
        
        // find elements that intersect with the above bounding box
        std::set<const Element*> intersected_elements;
        surface1.intersects(bb, intersected_elements);
        
#if DEBUG_CONTACT
        UInt inter = intersected_elements.size();
        cout<<"Surface 1 contains "<<inter<<" elements in contact zone."<<endl;
#endif
        surface2.intersects(bb, intersected_elements);
#if DEBUG_CONTACT
        cout<<"Surface 2 contains "<<(intersected_elements.size() - inter)<<" elements in contact zone."<<endl;
#endif
        
        // if intersected_elements container is not empty, a contact zone can be created
        
        if (!intersected_elements.empty()) {
          
          // now that the intersection has been found, get contact
          // zone object using bounding box and intersecting elements
          zones_.push_back(new contact_zone(model_, bb,intersected_elements, surface1, surface2));
        }
      }
    }
#ifdef DEBUG_CONTACT
    cout<<*this;
#endif
  }
  
  void local_search() {
    
    typedef typename contact_zone::node_set node_set;
    typedef typename contact_zone::node_iterator node_iterator;
    typedef typename contact_zone::element_container element_set;
    typedef typename contact_zone::element_iterator element_iterator;
    
    // loop over contact zones
    for (zone_iterator it = zones_.begin(); it != zones_.end(); ++it) {
      
      contact_zone& cs = **it;
      
      // loop over buckets in current contact zone
      // incrementation of the iterator is done when removing the bucket taking
      // care of not invalidating the iterators
      for (typename contact_zone::bucket_iterator bit = cs.buckets_begin(); bit != cs.buckets_end();) {
        
        // get nodes of bucket plus those of contiguous buckets
        element_set contiguous = cs.contiguous(bit->first);
        const Element* closest = NULL;
        
        // flagged nodes in case that a node belongs to several buckets
        // (node in bucket boundaries)
        node_set flaggedNodes;
        
        // loop over nodes
        for (node_iterator nit1 = bit->second.begin(); nit1 != bit->second.end(); ++nit1) {
          
          auto np = *nit1;
          
          // node has already been considered for a contact element,
          // do nothing further with it
          node_iterator fit = flaggedNodes.find(np);
          if (fit != flaggedNodes.end())
            continue;
          
          Real m = std::numeric_limits<Real>::infinity();
          
          // loop over elements
          for (element_iterator sit = contiguous.begin(); sit != contiguous.end(); ++sit) {
            
            // checkf if no self-contact is allowed
            if (c_ == No_self_contact_t)
              if (cs.in_surface(np, *sit)) {
#ifdef DEBUG_CONTACT
                cout<<"*** INFO *** Node "<<np<<" and element "<<(*sit)->element<<" belong to the same surface"<<endl;
#endif
                continue;
              }
            
            // check if node pointed by np lies in the same segment
            if (cs.in_element(np, *sit)) {
#ifdef DEBUG_CONTACT
              cout<<"*** INFO *** Node "<<np<<" belongs to element "<<*sit<<endl;
#endif
              continue;
            }
            
            // compute distance from node to element using SQP
            Real length = distance<d>(np, *sit, model_);
            
            if (length < m) {
              m = length;
              closest = *sit;
            }
            
          } // loop over elements
          
          // if a close node was found that does not belong to the
          // same element, add contact element for resolution
          if (closest) {
            celements_.push_back(contact_element(np, closest, model_));
            flaggedNodes.insert(np);
          }
#ifdef DEBUG_CONTACT
          else
            cout<<"*** INFO *** No close element was found for node "<<np<<endl;
#endif
          
        } // loop over nodes in bucket
        
        // remove bucket for further search as it is not needed
        cs.erase_bucket(bit++);
        
      } // loop over buckests
      
#ifdef DEBUG_CONTACT
      // print contact elements
      print_celements(cout);
#endif
      
    } // loop over contact zones
  }
  
  
  void remove_penetrations() {
    
    // loop over contact elements
    for (celement_iterator eit = celements_.begin(); eit != celements_.end(); ++eit) {
      
      bool flag = eit->penetrates();
      
      if (flag) {
#ifdef DEBUG_CONTACT
        cout<<"*** WARNING *** Penetration occurs for element "<<*eit<<endl;
#endif
        
        cout<<"*** INFO *** Collision detected, aborting..."<<endl;
        exit(1);
        //        eit->remove_penetration();
      }
      
    }
    
    // clear non-permanent objects
    clear();
  }
  
  void print_celements(std::ostream& os) {
    os<<"  Contact elements: "<<celements_.size()<<endl;
    for (celement_iterator eit = celements_.begin(); eit != celements_.end(); ++eit)
      os<<*eit;
  }
  
  
  
  //! Enable std output
  friend std::ostream& operator<<(std::ostream& os, const CManager& cm) {
    
    os<<"Contact manager info: "<<endl;
    os<<"  Contact surfaces: "<<cm.surfaces_.size()<<endl;
    typename CManager::const_surface_iterator it = cm.surfaces_.begin();
    for (; it != cm.surfaces_.end(); ++it)
      os<<*it;
    os<<"  Contact zones: "<<cm.zones_.size()<<endl;
    for (typename CManager::const_zone_iterator it = cm.zones_.begin(); it != cm.zones_.end(); ++it)
      os<<**it;
    
    return os;
  }
};

__END_AKANTU__

#endif /* __AKANTU_CMANAGER_HH__ */
