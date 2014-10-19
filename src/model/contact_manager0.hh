/**
 * @file   contact_manager0.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Tue May 13 2014
 * @date last modification: Tue May 13 2014
 *
 * @brief  zeroth level of contact (simplest implementation)
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

#ifndef __AKANTU_CONTACT_MANAGER_0_HH__
#define __AKANTU_CONTACT_MANAGER_0_HH__

#include <unordered_set>


//#include "aka_common.hh"
#include "model_manager.hh"

#define DEBUG_MANAGER 1


__BEGIN_AKANTU__


template <class pair_type>
class PairComp {
  
public:
  bool operator()(pair_type const &a, pair_type const &b) {
    return a.first->first < b.first->first
    || (!(b.first->first < a.first->first) && a.second->first < b.second->first);
  }
};







template <class Bounding_policy, Discretization_type DT, Consider_acceleration accel = Consider_t, template <class> class Cost_policy = Cost_functor>
class Contact0_model_manager : public Model_manager<SolidMechanicsModel>, public Kinematic_traits<accel> {
  
  
  using Kinematic_traits<accel>::resolve_time;
  
  typedef Real time_type;
  typedef ctimer chronograph_type;
  
  // model type
  typedef SolidMechanicsModel model_type;
  typedef model_type* model_pointer;
  typedef model_type& model_reference;
  
  // geometric types
  typedef Bounding_policy volume_type;
  typedef typename volume_type::point_type point_type;
  typedef typename point_type::value_type value_type;
  typedef typename volume_type::aabb_type aabb_type;
  
  // element type
  typedef ModelElement<model_type> element_type;
  
  // Bounding volume hierarchy related types
  typedef Cost_policy<volume_type> cost_functor;
  typedef DataTree<volume_type, element_type, Cost_policy > tree_type;
  typedef typename tree_type::leaves_data_iterator leaves_data_iterator;
  typedef typename tree_type::leaf_iterator tree_leaf_iterator;
  typedef typename tree_type::iterator tree_iterator;
  typedef typename tree_type::const_iterator const_tree_iterator;
  typedef std::list<tree_type*> forest_container;
  typedef typename forest_container::iterator forest_iterator;
  typedef typename forest_container::const_iterator const_forest_iterator;
  
  // data structure for detailed check
  typedef std::pair<leaves_data_iterator, leaves_data_iterator> check_type;
  
  
  //  typedef std::unordered_set<check_type, PairComp2<check_type> > check_container;
  typedef std::set<check_type, PairComp<check_type> > check_container;
  typedef typename check_container::iterator check_iterator;
  
  // kinematic types
  typedef point_type velocity_type;
  typedef std::list<velocity_type> velocity_container;
  typedef typename velocity_container::iterator velocity_iterator;
  
  typedef unsigned long mask_size;
  
  // timer type
  typedef std::priority_queue<time_type, std::vector<time_type>, std::greater<time_type> > timer_type;
  
  // tuple type
  typedef std::tuple<time_type, tree_iterator, tree_iterator> tuple_type;
  
  struct Tuple_compare {
    bool operator()(const tuple_type& t1, const tuple_type& t2) const
    { return std::get<0>(t1) > std::get<0>(t2); }
  };
  
  typedef typename std::priority_queue<tuple_type, std::vector<tuple_type>, Tuple_compare> hierarchy_timer;
  
  
  //! Structure used to do a postorder update of tree hierarchies
  struct Updater {
    
    tree_type &t_;        //!< Reference to hierarchy
    
    Updater(tree_type& t) : t_(t) {}
    
    void operator()(tree_iterator it) {
      if (!it.is_leaf()) {
        volume_type& v = *it;
        volume_type& lv = *it.left();
        volume_type& rv = *it.right();
        v = lv + rv;
        assert(lv.last_time_ == rv.last_time_);
        v.last_time_ = lv.last_time_;
        v.velocity_ = 0.5 * (lv.velocity_ + rv.velocity_);
        v.acceleration_ = 0.5 * (lv.acceleration_ + rv.acceleration_);
      }
    }
  };
  
  struct Printer {
    void operator()(tree_iterator it)
    { cout<<*it<<", "; }
  };
  
  class Time_exception : public std::exception {
    virtual const char* what() const throw()
    { return "*** EXCEPTION *** Zero time increment."; }
  };
  
  struct Continuator : public std::exception {
    
    tuple_type best_;
    
    Continuator(const tuple_type& b) : best_(b) {}
    
    virtual const char* what() const throw()
    { return "*** EXCEPTION *** Continue."; }
  };
  
  template <bool flag>
  struct Bool2Type {
    enum { value = flag };
  };
  
private:
  
  forest_container forest_;              //!< Bounding volume hierarchies
  hierarchy_timer timer_;                //!< Priority queue of times
  mask_size masks_;                      //!< Variable used for static objects
  tree_iterator null_;
  time_type last_;                       //!< Keep time of last detection engine reset
  check_container detailed_;
  bool detect_;
  
  //! Data structure that holds contact data
  struct ContactData {
    
    typedef SolidMechanicsModel model_type;
    typedef ModelElement<model_type> element_type;
    
    typedef std::set<UInt> node_set;
    
    struct Comparator {
      template <class iterator>
      bool operator()(iterator const &a, iterator const &b)
      { return a->first < b->first; }
    };
    
    typedef std::set<leaves_data_iterator, Comparator> elem_set;
    
    leaves_data_iterator left_;
    leaves_data_iterator right_;
    
    mutable elem_set ssurface_;
    mutable elem_set msurface_;
    mutable node_set slaves_;
    
    mutable tree_type *stree_;
    mutable tree_type *mtree_;
    
    ContactData(leaves_data_iterator l, leaves_data_iterator r) : left_(l), right_(r), ssurface_(), msurface_(), slaves_(), stree_(nullptr), mtree_(nullptr) {}
    
    template <class neighbor_type>
    void initialize(neighbor_type& ln, neighbor_type& rn) const {
      
      stree_ = ln.tree_;
      mtree_ = rn.tree_;
      
      // insert slave and master elements
      for (auto it = ln.elems_.begin(); it != ln.elems_.end(); ++it)
        ssurface_.insert(*it);
      for (auto it = rn.elems_.begin(); it != rn.elems_.end(); ++it)
        msurface_.insert(*it);
    }
    
    
    
    bool operator<(const ContactData& cd) const {
      
      return this->left_->first < cd.left_->first
      || (!(cd.left_->first < this->left_->first) && this->right_->first < cd.right_->first);
    }
    
    
    friend std::ostream& operator<<(std::ostream& os, const ContactData& cd) {
      
      os<<"Contact data info:"<<endl;
      os<<"  "<<cd.slaves_.size()<<" slave nodes";
      //      for (std::set<UInt>::const_iterator it = cd.slaves_.begin(); it!=cd.slaves_.end(); ++it)
      //        os<<" "<<*it<<endl;
      //        os<<" "<<*it;
      return os;
    }
    
  };
  
  std::set<ContactData> contact_;
  typedef typename std::set<ContactData>::iterator contact_iterator;
  
public:
  
  //! Default constructor
  Contact0_model_manager()
  : Model_manager(), forest_(), timer_(), masks_(), null_(tree_iterator(nullptr)), last_(), detect_(true) {}
  
  //! Destructor
  ~Contact0_model_manager() {
    
    // delete trees
    for (forest_iterator it = forest_.begin(); it != forest_.end(); ++it)
      delete *it;
  }
  
  virtual void add_model(model_pointer m, Kinematic_type k = dynamic_object_t) {
    
    m->initializeUpdateResidualData();
    models_.push_back(m);
    
    if (models_.size() > 8*sizeof(mask_size)) {
      cout<<"*** ERROR *** Type used for masks is too small to handle all models."<<endl;
      cout<<"Aborting..."<<endl;
      exit(1);
    }
    
    // create tree
    tree_type* tp = construct_tree_bottom_up<tree_type, model_type, element_type>(*m);
    forest_.push_back(tp);
    
    
#ifdef DEBUG_MANAGER
    //    cout<<"tree "<<*tp<<endl;
    //    print_mathematica(*tp);
#endif
    
    // mask model as dynamic or static
    masks_ |= (k << (models_.size()-1));
  }
  
  
  void update_forest(time_type t) {
    
    int k = 0;
    for (forest_iterator fit = forest_.begin(); fit != forest_.end(); ++fit) {
      
      // check if the object is dynamic to update
      if (!(masks_ & (1 << k++)))
        continue;
      
      // loop over leaves
      for (leaves_data_iterator lit = (*fit)->leaves_data_begin();
           lit != (*fit)->leaves_data_end(); ++lit) {
        
        Real t_old = lit->first->last_time_;
        
        std::vector<const Real*> c = lit->second.coordinates();
        volume_type v = Volume_creator<volume_type>::create(c);
        
        // get positions
        const point_type& p0 = lit->first->center();
        const point_type& p1 = v.center();
        
        // get velocities
        const point_type& v0 = lit->first->velocity_;
        const point_type& v1 = v.velocity_;
        
        // new velocity and acceleration
        v.velocity_ = 1/(t - t_old) * (p1-p0);
        v.acceleration_ = 1/(t - t_old) * (v1-v0);
        
        v.last_time_ = t;
        
        // set new volume
        *lit->first = v;
      }
      
      tree_type &t = **fit;
      Updater u(t);
      postorder(t.root(),u);
    }
  }
  
  
  void reset() {
    
    // clear queue
    while (!timer_.empty())
      timer_.pop();
    
    timer_.push(std::make_tuple(last_, null_, null_));
    timer_.push(std::make_tuple(inf, null_, null_));
  }
  
  
  template <class queue_type>
  void print_queue(queue_type copy) {
    
    cout<<"Printing queue values:";
    while (!copy.empty()) {
      const tuple_type& tuple = copy.top();
      cout<<" "<<std::get<0>(tuple);
      copy.pop();
    }
    cout<<endl;
  }
  
  
  
  /*! \param t - Current elapsed time
   * \param Dt - Time step
   */
  void intersect(time_type t, time_type Dt) {
    
#ifdef DEBUG_MANAGER
    cout<<"t = "<<t<<endl;
//    cout<<"t = "<<t<<", Dt = "<<Dt<<", timer:";
    
//    hierarchy_timer copy = timer_;
//    while (!copy.empty()) {
//      const tuple_type& tuple = copy.top();
//      cout<<" "<<std::get<0>(tuple);
//      copy.pop();
//    }
//    cout<<endl;
    
#endif
    
    static time_type Dt1 = Dt;
    
    // reset if first enter the function
    if (t == last_ || t == Dt1) {
      Dt = Dt1;
      reset();
    }
    
    const tuple_type& tuple = timer_.top();
    time_type top = std::get<0>(tuple);
    
    // update hierarchies and get positions
    // note that the updating starts before the next intetersection check
    if (models_.size() > 1 && top <= t + 3*Dt1) {
      update_forest(t);
//#ifdef DEBUG_MANAGER
//      cout<<"Updating forest"<<endl;
//#endif
    }
    
    // check if detection is shut off
    if (t + Dt < top)
      return;
    
    if (detect_) {
      
      // get iterators from priority element
      tree_iterator it1 = std::get<1>(tuple);
      tree_iterator it2 = std::get<2>(tuple);
      
      // remove time from timer
      timer_.pop();
      
      // check for intersection becase:
      //   1. there are enough models
      //   2. intersection happens before the next increment
      if (models_.size() > 1 && top <= t + Dt1) {
        
        // if first step, add next time to timer and return because there
        // is not enough information for the computation of intersection times
        if (t <= last_ + (accel == Consider_t ? 2 : 1)*Dt) {
          // remove time from timer
          timer_.push(std::make_tuple(t, null_, null_));
#ifdef DEBUG_MANAGER
          cout<<"Early out"<<endl;
#endif
          return;
        }
        
        // check if iterators are null to compute O(n^2) collision times
        if (it1 == null_ || it2 == null_) {
          
          // do O(n^2) operation to obtain next time of intersections
          for (forest_iterator it1 = forest_.begin(); it1 != --forest_.end(); ++it1) {
            
            forest_iterator it2 = it1;
            for (++it2; it2 != forest_.end(); ++it2) {
              
              // get collision time
#ifdef DEBUG_MANAGER
              cout<<"Calling check_collision in O(n^2) branch"<<endl;
#endif
              time_type tstar = resolve_time((*it1)->root(), (*it2)->root());
              if (tstar != inf)
                timer_.push(std::make_tuple(t+tstar, (*it1)->root(), (*it2)->root()));
#ifdef DEBUG_MANAGER
              else
                cout<<"*** INFO *** Objects do not intersect:\n  "<<*(*it1)->root()<<"\n  "<<*(*it2)->root()<<endl;
#endif
            } // inner hierarchy loop
          } // outer hierarchy loop
          
        }
        
        // else use collision information previously computed (avoids O(n^2) operation above)
        else {
          
#ifdef DEBUG_MANAGER
          cout<<"Calling check_collision in non-O(n^2) branch"<<endl;
#endif
          
          // temporary queue for tree traversal
          hierarchy_timer pq;
          pq.push(std::make_tuple(top, it1, it2));
          
          try {
            
            // enter infinite loop
            while (!pq.empty()) {
              
#ifdef DEBUG_MANAGER
              cout<<"______________________________________________"<<endl;
              cout<<"queue empty? -> "<<pq.empty()<<endl;
              if (!pq.empty())
                print_queue(pq);
#endif
              
              const tuple_type& tuple = pq.top();
              
              time_type tstar = std::get<0>(tuple);
              tree_iterator left = std::get<1>(tuple);
              tree_iterator right = std::get<2>(tuple);
              
#ifdef DEBUG_MANAGER
              cout<<"Queue time "<<tstar<<", items: "<<*left<<", "<<*right<<endl;
#endif
              
              pq.pop();
              
              check_collision(t, left, right, pq);
              
            } // infinite loop
            
          }
          catch (Continuator& e) {
            
            cout<<"In continuator"<<endl;
            
            tuple_type& best = e.best_;
            time_type& best_time = std::get<0>(best);
            best_time += t;
            
            // clean hierarchy timer until best time
            while (std::get<0>(timer_.top()) <= best_time)
              timer_.pop();
            
            timer_.push(best);
            timer_.push(best);
            
          }
          
        }
        
      } // if statement on enough models
      // else do nothing as there are not enough models to carry out intersection
      // or the check engine is shut down until the next time in timer
    }
    
    // detailed check
    detailed_check(t, Dt);
  }
  
  
  struct Neighbor_finder {
    
    int count_;
    element_type *el_;
    tree_type *tree_;
    
    class Comparator {
      
    public:
      template <class iterator>
      bool operator()(iterator const &a, iterator const &b) {
        return a->first < b->first;
      }
    };
    
    std::set<leaves_data_iterator, Comparator > elems_;
    
    Neighbor_finder() : count_(), el_(nullptr), tree_(nullptr) {}
    
    UInt size() const
    { return elems_.size(); }
    
    template <class container>
    void add_slaves(container& c) {
      
      for (auto it = elems_.begin(); it != elems_.end(); ++it) {
        auto elem = (*it)->second;
        for (size_t i=0; i<elem.numNodes(); ++i)
          c.insert(elem.node(i));
      }
    }
    
    template <class container>
    void add_masters(container& c) {
      
      for (auto it = elems_.begin(); it != elems_.end(); ++it)
        c.insert(&(*it)->second);
    }
    
    bool operator()() { assert(el_ != nullptr); return elems_.size() == el_->numNodes()+1; }
    
    void operator()(tree_iterator it) {
      
      assert(tree_ != nullptr);
      assert(el_ != nullptr);
      auto eit = tree_->find_data(it);
      if (eit != tree_->leaves_data_end())
        if (el_->shareNodes(eit->second))
          elems_.insert(eit);
      ++count_;
    }
    
    friend std::ostream& operator<<(std::ostream& os, const Neighbor_finder& nf) {
      
      os<<"Neighbor finder info:"<<endl;
      os<<"  found neighbors in "<<nf.count_<<" evaluations"<<endl;
      
      if (nf.el_)
        os<<"  reference element "<<*nf.el_<<endl;
      
      cout<<"  neighbors: "<<nf.size()<<endl;
      for (auto n:nf.elems_)
        os<<'\t'<<n->second<<'\n';
      return os;
    }
    
  };
  
  
  void detailed_check(time_type t, time_type Dt) {
    
    constexpr int dim = point_type::dim();
    
    if (contact_.size() > 0) {
      
      // loop over contact structures
      for (contact_iterator cit = contact_.begin(); cit != contact_.end(); ++cit) {
        
        ContactData& c = const_cast<ContactData&>(*cit);
        
        typedef std::set<element_type> elem_list;
        typedef std::pair<point_type, vector_type> closest_type;
        typedef std::tuple<UInt, Real, element_type, closest_type, elem_list> test_type;
        
        std::map<UInt, test_type > map;
        std::set<UInt> checked;
        
        typename ContactData::elem_set snew, mnew;
        
        // loop over slave elements to determine slave nodes
        for (auto sit:c.ssurface_) {
          
          // get slave element
          auto sel = sit->second;
          auto sbb = sel.template boundingBox<dim>();
          
          // loop over master elements
          for (auto mit:c.msurface_) {
            
            auto mel = mit->second;
            
            // tighter check with AABBs
            auto mbb = mel.template boundingBox<dim>();
            if (!(sbb & mbb))
              continue;
            
            // loop over slave nodes
            for (UInt i=0; i<sel.numNodes(); ++i) {
              
              UInt s = sel.node(i);
              
              // treat first element as slave
              auto coord = sel.coordinates();
              
              // create point
              point_type p(coord[i]);
              
              if (!penetrates(p,mel)) {
                continue;
              } else {
                
                // find node in contact data structure
                auto fslave = c.slaves_.find(s);
                // if not found, add it to container and search for new contact elements
                if (fslave == c.slaves_.end()) {
                  
                  c.slaves_.insert(s);
                  
                  Neighbor_finder sc, mc;
                  sc.tree_ = c.stree_;
                  mc.tree_ = c.mtree_;
                  sc.el_ = &sit->second;
                  mc.el_ = &mit->second;
                  c.stree_->collect_neighbors(sit->first, sc);
                  c.mtree_->collect_neighbors(mit->first, mc);
                  
                  for (auto i:sc.elems_)
                    snew.insert(i);
                  for (auto i:mc.elems_)
                    mnew.insert(i);
                }
                
                // compute closest point
                closest_type r = closest_point_to_element(p, mel);
                
                Real nd = (p-r.first).sq_norm();
                
                auto elit = map.find(s);
                if (elit == map.end()) {
                  map[s] = test_type(i, nd, sel, r, elem_list());
                  
                  test_type &tuple = map[s];
                  std::get<4>(tuple).insert(mel);
                }
                else {
                  
                  auto tuple = map[s];
                  Real dist = std::get<1>(tuple);
                  element_type sel = std::get<2>(tuple);
                  
                  if (std::abs(nd - dist) < 1.0e-6) {
                    
                    test_type &tuple = map[s];
                    
                    // edit closest point
                    elem_list &els = std::get<4>(tuple);
                    els.insert(mel);
                    
                    if (els.size() == 2) {
                      
                      closest_type &p = std::get<3>(tuple);
                      
                      p = commonPonit<point_type>(els);
                    }
                    
                    
                  } else if (nd < dist) {
                    
                    map[s] = test_type(i, nd, sel, r, elem_list());
                    
                    test_type &tuple = map[s];
                    std::get<4>(tuple).insert(mel);
                    
                  }
                  
                }
                
              }
              
            } // loop over master elements
            
            
          } // loop over slave nodes
        } // loop over slave elements to determine slave nodes
        
        // add new elements if found
        if (!snew.empty() || !mnew.empty()) {
          for (auto s:snew)
            c.ssurface_.insert(s);
          for (auto m:mnew)
            c.msurface_.insert(m);
        }
        
        
        // loop over map to balance slave nodes with closest element
        for (auto it = map.begin(); it != map.end(); ++it) {
          
          UInt id = std::get<0>(it->second);
          element_type sel = std::get<2>(it->second);
          auto cp = std::get<3>(it->second);
          auto mel = std::get<4>(it->second);
          
          
          // process slave
          auto el = *mel.begin();
          balance<point_type>(Dt, id, cp, sel, const_cast<element_type&>(el));
          
        }
      } // loop over contact structures
    } // contact_ > 0
    
    // WORKING CODE FOR SINGLE PASS
    //    if (contact_.size() > 0) {
    //
    //      // loop over contact structures
    //      for (auto c: contact_) {
    //
    //        typedef std::tuple<UInt, Real, element_type, element_type> test_type;
    //
    //        std::map<UInt, test_type > map;
    //        std::set<UInt> checked;
    //
    //        typename ContactData::elem_set snew, mnew;
    //
    //        cout<<"------------------------------------------------------------------"<<endl;
    //        cout<<"slaves -> "<<c.ssurface_.size()<<endl;
    //        cout<<"masters -> "<<c.msurface_.size()<<endl;
    //
    //        // loop over slave elements to determine slave nodes
    //        for (auto sit:c.ssurface_) {
    //
    //          // get slave element
    //          auto sel = sit->second;
    //          auto sbb = sel.template boundingBox<dim>();
    //
    //          cout<<"sel -> "<<sel<<endl;
    //
    //          // loop over master elements
    //          for (auto mit:c.msurface_) {
    //
    //            auto mel = mit->second;
    //            cout<<"mel -> "<<mel<<endl;
    //
    //            // tighter check with AABBs
    //            auto mbb = mel.template boundingBox<dim>();
    //            if (!(sbb & mbb)) {
    //              cout<<"AABBS"<<endl;
    //              continue;
    //            }
    //
    //            // loop over slave nodes
    //            for (UInt i=0; i<sel.numNodes(); ++i) {
    //
    //              UInt s = sel.node(i);
    //
    //              cout<<"***NODE "<<s<<endl;
    //              cout<<"***MBB -> "<<mbb<<endl;
    //              cout<<"***MN -> "<<mel.normal()<<endl;
    //
    //              // treat first element as slave
    //              auto coord = sel.coordinates();
    //
    //              // create point
    //              point_type p(coord[i]);
    //
    //              if (!penetrates(p,mel)) {
    //                cout<<">SLAVE GETOUT"<<endl;
    //                continue;
    //              } else {
    //
    //
    //                // find node in contact data structure
    //                auto fslave = c.slaves_.find(s);
    //                // if not found, add it to container and search for new contact elements
    //                if (fslave == c.slaves_.end()) {
    //                  c.slaves_.insert(s);
    //
    //                  cout<<"type of sit -> "<<typeid(sit).name()<<endl;
    //
    //                  Neighbor_finder sc, mc;
    //                  sc.tree_ = c.stree_;
    //                  mc.tree_ = c.mtree_;
    //                  sc.el_ = &sit->second;
    //                  mc.el_ = &mit->second;
    //                  c.stree_->collect_neighbors(sit->first, sc);
    //                  c.mtree_->collect_neighbors(mit->first, mc);
    //
    //                  for (auto i:sc.elems_)
    //                    snew.insert(i);
    //                  for (auto i:mc.elems_)
    //                    mnew.insert(i);
    //
    //                  for (auto i:sc.elems_)
    //                    cout<<i->second;
    //                  for (auto i:mc.elems_)
    //                    cout<<i->second;
    //                }
    //
    //                cout<<">SLAVE "<<s<<": "<<p<<", with id "<<i<<", belonging to "<<sel<<" PENETRATES"<<endl;
    //                cout<<"master element "<<mel<<", master bb -> "<<mbb<<endl;
    //
    //                // compute closest point
    //                std::pair<point_type, vector_type> r = closest_point_to_element(p, mel);
    //
    //                cout<<"r -> "<<r.first<<endl;
    //
    //                Real nd = (p-r.first).sq_norm();
    //
    //                auto elit = map.find(s);
    //                if (elit == map.end()) {
    //                  cout<<"NO ELEMENT IN MAP"<<endl;
    //                  map[s] = test_type(i, nd, sel, mel);
    //                }
    //                else {
    //                  cout<<"ELEMENT IN MAP"<<endl;
    //                  cout<<"stored info in map: ";
    //
    //                  auto tuple = map[s];
    //                  UInt id = std::get<0>(tuple);
    //                  Real dist = std::get<1>(tuple);
    //                  element_type sel = std::get<2>(tuple);
    //                  element_type mel = std::get<3>(tuple);
    //
    //                  cout<<"id "<<id<<", dist "<<dist<<"slave element "<<sel<<", master el "<<mel<<endl;
    //
    //                  cout<<"comparing new distance "<<nd<<" with stored value "<<dist<<endl;
    //                  if (nd < dist) {
    //                    cout<<"new distance smaller, inserting new element"<<endl;
    //                    map[s] = test_type(i, nd, sel, mel);
    //                  } else
    //                    cout<<"new distance is bigger!"<<endl;
    //                }
    //
    //              }
    //
    //            } // loop over master elements
    //
    //
    //          } // loop over slave nodes
    //        } // loop over slave elements to determine slave nodes
    //
    //        // add new elements if found
    //        if (!snew.empty() || !mnew.empty()) {
    //          for (auto s:snew)
    //            c.ssurface_.insert(s);
    //          for (auto m:mnew)
    //            c.msurface_.insert(m);
    //        }
    //
    //
    //        // loop over map to balance slave nodes with closest element
    //        for (auto it = map.begin(); it != map.end(); ++it) {
    //
    //          UInt id = std::get<0>(it->second);
    //          element_type sel = std::get<2>(it->second);
    //          element_type mel = std::get<3>(it->second);
    //
    //          cout<<"balance node with id -> "<<id<<" in slave element "<<sel<<endl;
    //          cout<<"master el -> "<<mel<<endl;
    //
    //          // process slave
    //          balance<point_type>(Dt, id, sel, mel);
    //
    //        }
    //      } // loop over contact structures
    //    } // contact_ > 0
  }
  
  
  template <class contact_type, class bbox_type, class slave_container>
  void solve_contact(time_type Dt, contact_type& c1, contact_type& c2, const bbox_type& bb, slave_container& slaves) {
    
    // treat first element as slave
    auto coord = c1.coordinates();
    
    int slave = 0;
    for (const double* c:coord) {
      
      // create point
      point_type p(c);
      
      // if point lies outside of collision zone, continue
      if (!(bb & p)) {
        ++slave;
        continue;
      }
      
      UInt id = c1.node(slave);
      auto it = slaves.find(id);
      if (it != slaves.end())
        continue;
      slaves.insert(id);
      
      // else find closest distance from p to contacting element c2
      std::pair<point_type, vector_type> r = closest_point_to_element(p, c2);
      const point_type& q = r.first;
      const vector_type& n = r.second;
      
      // get distance from current position
      Real delta = sqrt((q-p).sq_norm());
      
      auto mass = c1.getMass(slave)[0];
      
      // compute force at slave node
      vector_type N = 2 * delta * mass / pow(Dt,2.) * n;
      
      // compute forces in master element balancing linear and angular momenta
      balance(Dt, slave, r, N, c1, c2);
      
      ++slave;
    }
    
  }
  
  
private:
  
  
  
  template <class iterator>
  void traverse_right(iterator it1, iterator it2, hierarchy_timer& pq) {
    
#ifdef DEBUG_MANAGER
    cout<<"      traversing right hierarchy"<<endl;
#endif
    
    iterator lit = it2.left();
    iterator rit = it2.right();
    assert(lit != null_);
    assert(rit != null_);
    
    time_type tstar1 = resolve_time(it1, lit);
    if (tstar1 != inf) {
#ifdef DEBUG_MANAGER
      cout<<"    Adding queue time "<<tstar1<<endl;
#endif
      pq.push(std::make_tuple(tstar1, it1, lit));
    }
    time_type tstar2 = resolve_time(it1, rit);
    if (tstar2 != inf) {
#ifdef DEBUG_MANAGER
      cout<<"    Adding queue time "<<tstar2<<endl;
#endif
      pq.push(std::make_tuple(tstar2, it1, rit));
    }
  }
  
  template <class iterator>
  void traverse_left(iterator it1, iterator it2, hierarchy_timer& pq) {
    
#ifdef DEBUG_MANAGER
    cout<<"      traversing left hierarchy"<<endl;
#endif
    
    iterator lit = it1.left();
    iterator rit = it1.right();
    assert(lit != null_);
    assert(rit != null_);
    
    time_type tstar1 = resolve_time(lit, it2);
    if (tstar1 != inf) {
#ifdef DEBUG_MANAGER
      cout<<"    Adding queue time "<<tstar1<<endl;
#endif
      pq.push(std::make_tuple(tstar1, lit, it2));
    }
    time_type tstar2 = resolve_time(rit, it2);
    if (tstar2 != inf) {
#ifdef DEBUG_MANAGER
      cout<<"    Adding queue time "<<tstar2<<endl;
#endif
      pq.push(std::make_tuple(tstar2, rit, it2));
    }
  }
  
  
  template <class iterator>
  void check_collision(time_type t, iterator it1, iterator it2, hierarchy_timer& pq) {
    
    // if volumes are leaves, change the timer and time step
    if (it1.is_leaf() && it2.is_leaf()) {
      
      cout<<"*** FOUND LEAVES ***"<<endl;
      
      // case where objects intersect
      if (*it1 & *it2) {
        
        cout<<"*** BOUNDING SPHERE INTERSECTION ***"<<endl;
        
        // add leaves to carry out penetration tests
        leaves_data_iterator lit1(nullptr), lit2(nullptr);
        
        Neighbor_finder lc, rc;
        
        for (forest_iterator it = forest_.begin(); it != forest_.end(); ++it) {
          
          leaves_data_iterator lit = (*it)->find_data(it1);
          if (lit != (*it)->leaves_data_end()) {
            lit1 = lit;
            
            lc.tree_ = *it;
            lc.el_ = &lit->second;
            (*it)->collect_neighbors(lit->first, lc);
          }
          
          lit = (*it)->find_data(it2);
          if (lit != (*it)->leaves_data_end()) {
            lit2 = lit;
            
            rc.tree_ = *it;
            rc.el_ = &lit->second;
            (*it)->collect_neighbors(lit->first, rc);
          }
        }
        
        assert (lit1 != leaves_data_iterator(nullptr));
        assert (lit2 != leaves_data_iterator(nullptr));
        
        // check for new ContactData
        auto ins = contact_.insert(ContactData(lit1, lit2));
        
        if (ins.second) {
          auto cd = ins.first;
          cd->initialize(lc,rc);
          cout<<"*** INFO *** Adding contact data."<<endl;
        } else {
          cout<<"*** WARNING *** Contact data not inserted."<<endl;
        }
        
        // add to data structure for detailed check
        detailed_.insert(std::make_pair(lit1, lit2));
        
        detect_ = false;
        
        cout<<"***DETAILED SIZE -> "<<detailed_.size()<<endl;
        
      }
      else
        cout<<"*** NO INTERSECTION BETWEEN BOUNDING SPHERES ***"<<endl;
      
      
      time_type tstar = resolve_time(it1, it2);
      
#ifdef DEBUG_MANAGER
      
      if (tstar == inf)
        cout<<"    Leaves found that DO NOT collide"<<endl;
      else {
        cout<<"    Leaves found that collide at time "<<(tstar)<<endl;
      }
#endif
      
      if (tstar == inf) {
        cout<<"    Leaves found that DO NOT collide"<<endl;
      }
      
      throw Continuator(std::make_tuple(tstar, it1, it2));
    }
    
    // found left leaf, traverse right hierarchy
    else if (it1.is_leaf() && !it2.is_leaf()) {
      
#ifdef DEBUG_MANAGER
      cout<<"      s1 is leaf"<<endl;
#endif
      traverse_right(it1, it2, pq);
    }
    
    // found right leaf, traverse left hierarchy
    else if (!it1.is_leaf() && it2.is_leaf()) {
      
#ifdef DEBUG_MANAGER
      cout<<"      s2 is leaf"<<endl;
#endif
      traverse_left(it1, it2, pq);
    }
    
    // else non-leaf case found, check volume sizes
    else {
      
      value_type m1 = it1->measure();
      value_type m2 = it2->measure();
      
      // volumes are equal to numerical error, traverse both hierarchies
      if (equal(m1, m2)) {
        
#ifdef DEBUG_MANAGER
        cout<<"      "<<m1<<" == "<<m2<<endl;
#endif
        traverse_right(it1, it2, pq);
        traverse_left(it1, it2, pq);
      }
      // left volume is bigger, traverse right hierarchy
      else if (m1 > m2) {
        
#ifdef DEBUG_MANAGER
        cout<<"      "<<m1<<" > "<<m2<<endl;
#endif
        traverse_left(it1, it2, pq);
      }
      // right volume is bigger, traverse left hierarchy
      else if (m1 < m2) {
        
#ifdef DEBUG_MANAGER
        cout<<"      "<<m1<<" < "<<m2<<endl;
#endif
        traverse_right(it1, it2, pq);
      }
    } // non-leaf case
    
  }
  
  
  friend std::ostream& operator<<(std::ostream& os, const Contact0_model_manager& mm) {
    
    os<<"Contact model manager info:"<<endl;
    os<<"  models: "<<mm.models_.size()<<endl;
    size_t i=0;
    const_forest_iterator tit = mm.forest_.begin();
    for (const_model_iterator it = mm.models_.begin(); it != mm.models_.end(); ++it) {
      os<<"\tmodel "<<++i<<" memory address: "<<*it<<endl;
      os<<"\tmodel: "<<**it<<endl;
      os<<"\ttree: ";
      print_mathematica(**tit++);
    }
    return os;
  }
};

__END_AKANTU__

#endif /* __AKANTU_CONTACT_MANAGER_0_HH__ */
