/**
 * @file   model_manager.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Mon Jan 07 2013
 * @date last modification: Fri Sep 05 2014
 *
 * @brief  higher order object that deals with collections of models
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

#ifndef __AKANTU_MODEL_MANAGER_HH__
#define __AKANTU_MODEL_MANAGER_HH__

#include "aka_common.hh"
#include "mesh.hh"
#include "model.hh"
#include "aka_tree.hh"
#include "solid_mechanics_model.hh"
#include "aka_plane.hh"
#include "aka_geometry.hh"
#include "solid_mechanics_model_element.hh"
#include "aka_timer.hh"

#include <cstring>
#include <queue>
#include <unordered_set>
#include <array/expr.hpp>



#define DEBUG_MANAGER 1


__BEGIN_AKANTU__


typedef array::vector_type<Real> vector_type;
typedef array::matrix_type<Real> matrix_type;
using array::transpose;


enum Discretization_type { Node_to_node, Node_to_segment };


template <class Model_policy>
class Model_manager {
  
public:
  
  typedef Model_policy model_type;
  typedef model_type* model_pointer;
  typedef model_type& model_reference;
  
  typedef std::list<model_type*> model_container;
  typedef typename model_container::iterator model_iterator;
  typedef typename model_container::const_iterator const_model_iterator;
  
protected:
  
  model_container models_;            //!< Models
  
  
public:
  
  //! Default constructor
  Model_manager() : models_() {}
  
  virtual void add_model(model_reference m)
  { models_.push_back(&m); }
  
  virtual void add_model(model_pointer m)
  { models_.push_back(m); }
  
  model_iterator models_begin()
  { return models_.begin(); }

  model_iterator models_end()
  { return models_.end(); }

  friend std::ostream& operator<<(std::ostream& os, const Model_manager& mm) {
    
    os<<"Model manager info:"<<endl;
    os<<"  models: "<<mm.models_.size()<<endl;
    size_t i=0;
    for (const_model_iterator it = mm.models_.begin(); it != mm.models_.end(); ++it) {
      os<<"\tmodel "<<++i<<" memory address: "<<*it<<endl;
      os<<"\tmodel "<<**it<<endl;
    }
    
    return os;
  }
  
};



enum Kinematic_type { static_object_t, dynamic_object_t};
enum Consider_acceleration { Consider_t, Neglect_t };


template <Consider_acceleration>
class Kinematic_traits;



template<>
class Kinematic_traits<Consider_t> {
  
protected:
  
  typedef Real time_type;
  
  // check collision considering acceleration
  // the equation to solve is
  //
  //   alpha t^4 + beta t^3 + gamma t^2 + delta t + epsilon = 0
  //
  // where alpha = (a.a)/4, beta = (a.v), gamma = a.s + v.v, delta = 2(s.v),
  // epsilon = s.s - r^2, a = a1-a2, v = v1-v2, s = c1-c2, r = r1-r2
  //
  template <class iterator>
  time_type resolve_time(iterator it1, iterator it2) {
    
    typedef typename iterator::value_type volume_type;
    typedef typename volume_type::point_type point_type;
    
    const point_type& v1 = it1->velocity_;
    const point_type& v2 = it2->velocity_;
    
    point_type a1 = it1->acceleration_;
    point_type a2 = it2->acceleration_;
    
    Real r = it1->radius() + it2->radius();
    point_type s = it2->center() - it1->center();
    point_type v = v2 - v1;
    point_type a = a2 - a1;
    
    // get coefficients
    Real alpha = (a*a)/4;
    Real beta = a*v;
    Real gamma = a*s + v*v;
    Real delta = 2*(s*v);
    Real epsilon = s*s - r*r;
    
    // obtain roots from quartic equation by calling the function such that
    // the coefficient for the quartic term is equal to one
    std::vector<Real> x(4, inf);
    
    
    uint32_t roots = solve_quartic(beta/alpha, gamma/alpha, delta/alpha, epsilon/alpha,
                                   &x[0], &x[1], &x[2], &x[3]);
    
    Real tmin = inf;
    
    // if there are roots, take the first one as an indication of the collision
    if (roots > 0) {
      
      for (size_t i=0; i<x.size(); ++i)
        tmin = std::min(tmin, x[i]);
      
      if (tmin > 0) {
#ifdef DEBUG_MANAGER
        cout<<"        Solve quartic with coefficients ";
        cout<<(beta/alpha)<<", "<<(gamma/alpha)<<", "<<(delta/alpha)<<", "<<(epsilon/alpha)<<endl;
        cout<<"          Approximate collision time -> tmin = "<<tmin<<" = "<<(tmin)<<endl;
        cout<<"          Roots:";
        for (size_t i=0; i<x.size(); ++i)
          cout<<" "<<x[i];
        cout<<endl;
#endif
        for (size_t i=0; i<x.size(); ++i)
          tmin = std::min(tmin, x[i]);
        
      }
    } // if roots
    return tmin > 0 ? tmin : inf;
  }
  
private:
  
  /*! \brief Solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d = 0.
   *
   *  Solves the quartic equation. Returns the number of roots
   *  found. Roots are filled in ascending order.
   *  Code taken from the Ion Beam Simulator, which is distrubuted under
   *  the terms of the GNU General Public License as published by the Free
   *  Software Foundation; either version 2 of the License, or (at your
   *  option) any later version.
   */
  uint32_t solve_quartic(Real a, Real b, Real c, Real d,
                         Real *x0, Real *x1, Real *x2, Real *x3 );
  
};



template<>
class Kinematic_traits<Neglect_t> {
  
protected:
  
  typedef Real time_type;
  
  // check collision neglecting acceleration
  template <class iterator>
  time_type resolve_time(iterator it1, iterator it2) {
    
    typedef typename iterator::value_type volume_type;
    typedef typename volume_type::point_type point_type;
    typedef typename point_type::value_type value_type;
    
    const volume_type& s1 = *it1;
    const volume_type& s2 = *it2;
    
    const point_type& v1 = s1.velocity_;
    const point_type& v2 = s2.velocity_;
    
    // vector between spheress
    point_type s = s2.center() - s1.center();
    // relative motion of s1 with respect to stationary s0
    point_type v = v2 - v1;
    // sum of radii
    value_type r = s1.radius() + s2.radius();
    value_type c = s*s - r*r;
    
#ifdef DEBUG_MANAGER
    cout<<"Checking collisiong between:"<<endl;
    cout<<"  "<<s1<<", velocity: "<<v1<<endl;
    cout<<"  "<<s2<<", velocity: "<<v2<<endl;
    cout<<"  Relative velocity: "<<v<<endl;
#endif
    
    value_type epsilon = 1e-8;
    
    // already intersecting
    if (c < -epsilon) {
      
#ifdef DEBUG_MANAGER
      cout<<"  Intersecting bounding volumes"<<endl;
#endif
      // should not get to this point
      return time_type();
    }
    
    value_type a = v*v;
    // if spheres not moving relative to each other
    if (a < epsilon) {
#ifdef DEBUG_MANAGER
      cout<<"    Objects not moving relative to each other"<<endl;
      cout<<"    a: "<<a<<endl;
#endif
      return inf;
    }
    
    value_type b = v*s;
    // if spheres not moving towards each other
    if (b >= 0.) {
#ifdef DEBUG_MANAGER
      cout<<"    Objects not moving towards each other"<<endl;
      cout<<"    b: "<<b<<endl;
#endif
      return inf;
    }
    
    value_type d = b*b - a*c;
    // if no real-valued root (d < 0), spheres do not intersect
    // otherwise add time to timer
    if (d >= 0.) {
      time_type ts = (-b - sqrt(d))/a;
      if (ts > -epsilon) {
#ifdef DEBUG_MANAGER
        cout<<"    Objects intersects at time "<<(ts)<<endl;
#endif
        return ts;
      }
#ifdef DEBUG_MANAGER
      else {
        cout<<"    ts negative: "<<ts<<endl;
      }
#endif
      
    }
#ifdef DEBUG_MANAGER
    else {
      cout<<"    Objects do not intersect"<<endl;
      cout<<"    discriminant: "<<d<<endl;
    }
#endif
    
    return inf;
  }
};



template <class VolumeType, class DataPolicy, template <class> class CostPolicy = Cost_functor>
class DataTree :
public Tree<VolumeType, CostPolicy> {
  
public:
  
  typedef DataPolicy data_type;
  typedef VolumeType volume_type;
  typedef CostPolicy<volume_type> cost_functor;
  typedef Tree<volume_type, CostPolicy> tree_type;
  typedef typename tree_type::iterator tree_iterator;
  typedef typename tree_type::const_iterator const_tree_iterator;
  
  // leaf information
  typedef std::map<tree_iterator, data_type> leaves_data;
  typedef typename leaves_data::iterator leaves_data_iterator;
  
  bool add_data(tree_iterator it, data_type& data) {
    auto i = data_.insert(std::make_pair(it, data));
    return i.second;
  }
  
  leaves_data_iterator leaves_data_begin()
  { return data_.begin(); }
  
  leaves_data_iterator leaves_data_end()
  { return data_.end(); }
  
  leaves_data_iterator find_data(tree_iterator it)
  { return data_.find(it); }
  
private:
  
  leaves_data data_;
  
};

template <class U, class T>
std::pair<U, T> minmax(const U& u, const T& t) {
  return u < t ? std::make_pair(u,t) : std::make_pair(t,u);
}

template <Discretization_type, class, class>
class ContactElement;
  

template <class point_type, class element_type>
class ContactElement<Node_to_node, point_type, element_type> {
    
public:
    
  typedef Real time_type;
  typedef typename vector_type::value_type value_type;
  typedef std::tuple<time_type, point_type> impact_tuple;
  typedef typename element_type::model_type model_type;
  
  struct Comparator {
    
    bool operator()(const ContactElement *c1, const ContactElement *c2) const {
      
      auto p1 = minmax(c1->el1_, c1->el2_);
      auto p2 = minmax(c2->el1_, c2->el2_);
      
      if (p1.first < p2.first || p1.second < p2.second)
        return true;
      
      auto i1 = minmax(c1->id1_, c1->id2_);
      auto i2 = minmax(c2->id1_, c2->id2_);
      
      if (i1.first < i2.first || i1.second < i2.second)
        return true;
      
      return false;
    }
  };
  
  typedef Comparator comparator_type;
  
  //! Parameter constructor
  template <class parameter_type>
  ContactElement(const parameter_type& p) :
  id1_(std::get<0>(p)), id2_(std::get<1>(p)),
  el1_(std::get<2>(p)), el2_(std::get<3>(p)),
  impact_(std::get<4>(p)), m1_(), m2_(), linked_(), release_(true) {}
  
  
  //! Resolve impact
  /*! At the moment of impact, obtain velocities after impact from
   * contacting nodes, and join masses (nodes behave as one)
   */
  void resolve_impact(time_type t, time_type Dt) {
    
    // if at moment of impact
    if (equal(t, std::get<0>(impact_)) && !linked_) {
      
      
#ifdef DEBUG_MANAGER
      cout<<"Linking nodes"<<endl;
#endif
      
      // get models
      model_type &model1 = el1_->model();
      model_type &model2 = el2_->model();
      
      // get global node ids
      UInt n1 = el1_->node(id1_);
      UInt n2 = el2_->node(id2_);
      
      // get references to mass and velocity vectors
      Array<Real> &mass1 = model1.getMass();
      Array<Real> &mass2 = model2.getMass();
      Array<Real> &velocity1 = model1.getVelocity();
      Array<Real> &velocity2 = model2.getVelocity();
      
      // get references to masses and velocities involved in the collision
      Real &m1 = mass1(n1);
      Real &m2 = mass2(n2);
      value_type &v1 = velocity1(n1);
      value_type &v2 = velocity2(n2);
      
      // set correct velocity
      Real vc = (m1*v1 + m2*v2) / (m1+m2);
      velocity1(n1) = vc;
      velocity2(n2) = vc;
      
      // save mass values and add masses to treat them as a single node
      m1_ = m1;
      m2_ = m2;
      mass1(n1) += m2_;
      mass2(n2) += m1_;
      
      // set link flag
      linked_ = true;
    }
  }

  //! Treat linked nodes
  /*! During the time that the nodes are in contact, treat them as a single
   * node of joint mass and synchronize residual values
   */
  bool resolve(time_type t, time_type Dt) {
    
    // get models
    model_type &model1 = el1_->model();
    model_type &model2 = el2_->model();
    
    // get global node ids
    UInt n1 = el1_->node(id1_);
    UInt n2 = el2_->node(id2_);
    
    // get references to residual vectors
    Array<Real> & r1 = model1.getResidual();
    Array<Real> & r2 = model2.getResidual();
    
    // check condition for delinking of nodes
    vector_type vec = el2_->barycenter() - el1_->barycenter();
    
    if (std::signbit(vec[0]) != std::signbit(r1(n1))) {

#ifdef DEBUG_MANAGER
      cout<<"Unlinking nodes"<<endl;
#endif

      if (release_) {
        
        // get masses
        Array<Real> &mass1 = model1.getMass();
        Array<Real> &mass2 = model2.getMass();
                
        mass1(id1_) = m1_;
        mass2(id2_) = m2_;
        
        release_ = false;
        return true;
      }
    }
    
    // synchronize if necessary
    if (linked_ && release_) {
      Real tmp = r1(n1);
      
      r1(n1) += r2(n2);
      r2(n2) += tmp;
    }
    return false;
  }
  
  friend std::ostream& operator<<(std::ostream& os, const ContactElement& ce) {
    cout<<"Contact element info:\n  linking elements ("<<ce.el1_<<" - "<<ce.el2_<<")"<<endl;
    cout<<"  colliding nodes ("<<ce.id1_<<" - "<<ce.id2_<<")"<<endl;
    cout<<"  first impact at location "<<std::get<1>(ce.impact_)<<" at time "<<std::get<0>(ce.impact_)<<endl;
    if (ce.m1_ > 0. && ce.m2_ > 0)
      cout<<"  saved masses: ("<<ce.m1_<<","<<ce.m2_<<")"<<endl;
    if (ce.linked_)
      cout<<"  linked state (treating contacting nodes as a single node)"<<endl;
    
    return os;
  }
  
private:
  size_t id1_, id2_;                    //!< Ids of element nodes involved in the collision
  element_type *el1_, *el2_;            //!< Pointers to elements involved in teh collision
  impact_tuple impact_;                 //!< Impact information, tuple containing the time and point of contact
  
  value_type m1_, m2_;                  //!< Mass values stored after linking of the nodes
  bool linked_, release_;               //!< Flags used to specify the state of the linking
};


template <class Bounding_policy, Discretization_type DT, Consider_acceleration accel = Consider_t, template <class> class Cost_policy = Cost_functor>
class Contact_model_manager : public Model_manager<SolidMechanicsModel>, public Kinematic_traits<accel> {
  
  
  using Kinematic_traits<accel>::resolve_time;
  
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
  typedef ContactElement<DT, point_type, element_type> contact_element_type;
  typedef std::set<contact_element_type*, typename contact_element_type::comparator_type > contact_element_container;
  typedef typename contact_element_container::iterator contact_element_iterator;
  typedef typename contact_element_type::time_type time_type;
  typedef typename contact_element_type::impact_tuple impact_tuple;
    
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
  
  struct Contactor : public std::exception {
    
    impact_tuple data_;
    
    Contactor(const impact_tuple& c) : data_(c) {}
    
    virtual const char* what() const throw()
    { return "*** EXCEPTION *** Contact."; }
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
  contact_element_container celems_;     //!< Contact elements
  
public:
  
  //! Default constructor
  Contact_model_manager()
  : Model_manager(), forest_(), timer_(), masks_(), null_(tree_iterator(nullptr)), last_() {}
  
  //! Destructor
  ~Contact_model_manager() {

    // delete trees
    for (forest_iterator it = forest_.begin(); it != forest_.end(); ++it)
      delete *it;
    
    // delete contact element
    for (contact_element_iterator it = celems_.begin(); it != celems_.end(); ++it)
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
    cout<<"tree "<<*tp<<endl;
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
      
  void resolve(time_type t, time_type Dt) {
    
    // loop over contact elements
    contact_element_iterator elit = celems_.begin();
    while (elit != celems_.end()) {

      cout<<**elit<<endl;
      if ((*elit)->resolve(t, Dt)) {
        celems_.erase(elit++);
 
        last_ = t + Dt;
        
      } else
        ++elit;
    }
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
  void intersect(time_type t, time_type& Dt) {
    
#ifdef DEBUG_MANAGER
    cout<<"t = "<<t<<", Dt = "<<Dt<<", timer:";
    
    hierarchy_timer copy = timer_;
    while (!copy.empty()) {
      const tuple_type& tuple = copy.top();
      cout<<" "<<std::get<0>(tuple);
      copy.pop();
    }
    cout<<endl;
    
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
#ifdef DEBUG_MANAGER
      cout<<"Updating forest"<<endl;
#endif
    }
    
    // check if detection is shut off
    if (t + Dt < top)
      return;
    
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
          while (true) {
            
#ifdef DEBUG_MANAGER
            cout<<"______________________________________________"<<endl;
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
          
          tuple_type& best = e.best_;
          time_type& best_time = std::get<0>(best);
          Dt = best_time;
          best_time += t;
          
          // clean hierarchy timer until best time
          while (std::get<0>(timer_.top()) <= best_time)
            timer_.pop();
          
          timer_.push(best);
          timer_.push(best);
          
          try {
            // set time step if required
            cout<<"calling set_time_set in Continuator"<<endl;

            set_time_step(t, Dt, Dt1);
          } catch (Time_exception& e) {
            cout<<"catching  inner Time_exception"<<endl;
            
          }
          
          
        }
        catch (Contactor& c) {
          
          // get collision impact
          Dt = std::get<0>(c.data_);
          
          try {
            
            cout<<"calling set_time_set in Contactor"<<endl;
            
            // set time step if required
            set_time_step(t, Dt, Dt1);
            //            last_ = t + Dt;
            
          } catch (Time_exception& e) {
            // if too small a time step, do nothing, next time step
            // will carry out DCR
            
            cout<<"catching outter Time_exception"<<endl;
            
            //            last_ = t + Dt;
            
            Dt = Dt1;
          }
          cout<<"last -> "<<last_<<endl;
        }
      }
      
    } // if statement on enough models
    // else do nothing as there are not enough models to carry out intersection
    // or the check engine is shut down until the next time in timer
    
    
    // resolve collision if necessary
    for (contact_element_iterator elit = celems_.begin(); elit != celems_.end(); ++elit)
      (*elit)->resolve_impact(t, Dt);

    
  }
  
  
private:
  
  void set_time_step(time_type t, time_type& Dt, time_type Dt1) {
    
    if (Dt > Dt1)
      Dt = Dt1;
    
    if (Dt < 1e-10) {
      cout<<"*** INFO *** New time step is too small. Throwing exception..."<<endl;
      throw Time_exception();
    }
#ifdef DEBUG_MANAGER
    cout<<"    Leaves found that collide at time "<<(t + Dt)<<endl;
    cout<<"    Setting new time step: "<<Dt<<endl;
#endif
    for (model_iterator mit = models_.begin(); mit != models_.end(); ++mit)
      (*mit)->setTimeStep(Dt);
  }
  
  // check time of collision between points
  impact_tuple check_points(time_type t, leaves_data_iterator it1, leaves_data_iterator it2) {
    
    const volume_type& v1 = *it1->first;
    const volume_type& v2 = *it2->first;
    
    std::vector<const Real*> c1 = it1->second.coordinates();
    std::vector<const Real*> c2 = it2->second.coordinates();
    
    // compute intersection
    volume_type vint = v1 && v2;
    
    // get indices of colliding nodes
    size_t ii=0, jj=0;
    for (size_t i=1; i<c1.size(); ++i) {
      
      if (vint & point_type(c1[i]))
        ii = i;
      if (vint & point_type(c2[i]))
        jj = i;
    }
    
    impact_tuple impact =
    moving_point_against_point(point_type(c1[ii]), point_type(c2[jj]), it1->first->velocity_, it2->first->velocity_);
    
    time_type &timpact = std::get<0>(impact);
    timpact += t;
    
    celems_.insert(new contact_element_type(std::make_tuple(ii, jj, &it1->second, &it2->second, impact)));
    
    return impact;
    
  }
  
  // check time of collision between 2D segments
  impact_tuple check_2D_sides(time_type t, leaves_data_iterator it1, leaves_data_iterator it2) {
    
    typedef Point<3> test_point;

      
    const volume_type& v1 = *it1->first;
    const volume_type& v2 = *it2->first;
    
    std::vector<const Real*> c1 = it1->second.coordinates();
    std::vector<const Real*> c2 = it2->second.coordinates();
    
    // create plane from second container node
    // THIS WON'T WORK, SIDES HAVE ONLY TWO NODES, CONTINUE DEVELOPING FROM HERE
    assert(c2.size() == 3);
    
    // form plane from three points
    test_point o,p,q;
    
    for (size_t j=0; j<point_type::dim(); ++j) {
      o[j] = c2[0][j];
      p[j] = c2[1][j];
      q[j] = c2[2][j];
    }
    
//    point_type o(c2[0]);
//    point_type p(c2[1]);
//    point_type q(c2[2]);
    
    Plane pi(o,p,q);
    
    // loop over the nodes of the container
    for (size_t i=0; i<c1.size(); ++i) {

      test_point v;
      point_type v2d = v1.velocity_ -  v2.velocity_;

      // create 3D point, used to check for plane intersection
      test_point x;
      for (size_t j=0; j<point_type::dim(); ++j) {
        x[j] = c1[i][j];
        v[j] = v2d[j];
      }
      
      cout<<"x -> "<<x<<endl;
      
      cout<<"v -> "<<v<<endl;
      
      
      std::tuple<time_type, test_point> impact =
      moving_point_against_plane(x, v, pi);

      
      
    }

    exit(1);
    
////    impact_tuple impact =
////    moving_point_against_point(point_type(c1[ii]), point_type(c2[jj]), it1->first->velocity_, it2->first->velocity_);
////    
////    time_type &timpact = std::get<0>(impact);
////    timpact += t;
////    
////    celems_.insert(new contact_element_type(std::make_tuple(ii, jj, &it1->second, &it2->second, impact)));
////    
//    return impact;
    
  }
  
    
  // check time of collision between triangles
  impact_tuple check_triangles(leaves_data_iterator it1, leaves_data_iterator it2) {
    
    std::vector<const Real*> c1 = it1->second.coordinates();
    std::vector<const Real*> c2 = it2->second.coordinates();
    
    assert(c1.size() == 3);
    assert(c1.size() == c2.size());
    impact_tuple min = std::make_tuple(inf,point_type());
    
    for (size_t i=0; i<c1.size(); ++i) {
      
      point_type r(c1[i]);
      
      // form plane from three points
      point_type o(c2[0]);
      point_type p(c2[1]);
      point_type q(c2[2]);
      
      Plane pi(o,p,q);
      
      // relative velocity
      point_type v = it1->first->velocity_ - it2->first->velocity_;
      
      // intersect point r with velocity v with plane pi
      // the function returns collision time and point of contact
      impact_tuple impact = moving_point_against_plane(r, v, pi);
      
      // make sure intersection point lies within the triangle
      if (is_point_in_triangle(std::get<1>(impact), o, p, q))
        if (std::get<0>(impact) < std::get<0>(min))
          min = impact;
    }
    return min;
  }
  
  impact_tuple fine_collision_time(time_type t, leaves_data_iterator it1, leaves_data_iterator it2, Int2Type<1>) {
    
    // check nodes of first segment against second segment
    return check_points(t, it1, it2);
  }
  
  impact_tuple fine_collision_time(time_type t, leaves_data_iterator it1, leaves_data_iterator it2, Int2Type<2>) {

    // check sides of the colliding elements
    return check_2D_sides(t, it1, it2);
  }
  
  impact_tuple fine_collision_time(time_type t, leaves_data_iterator it1, leaves_data_iterator it2, Int2Type<3>) {
    
    // check nodes of first triangle against second triangle
    impact_tuple impact1 = check_triangles(it1, it2);
    impact_tuple impact2 = check_triangles(it2, it1);
    
    // check nodes of second triangle against first triangle
    return std::get<0>(impact1) < std::get<0>(impact2) ? impact1 : impact2;
  }
  
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
        
        for (forest_iterator it = forest_.begin(); it != forest_.end(); ++it) {
          
          leaves_data_iterator lit = (*it)->find_data(it1);
          if (lit != (*it)->leaves_data_end())
            lit1 = lit;
          
          lit = (*it)->find_data(it2);
          if (lit != (*it)->leaves_data_end())
            lit2 = lit;
          
        }
        
        assert (lit1 != leaves_data_iterator(nullptr));
        assert (lit2 != leaves_data_iterator(nullptr));
        
        // determine collision time at the lowest level of detection
        impact_tuple impact = fine_collision_time(t, lit1, lit2, Int2Type<volume_type::dim()>());
        
#ifdef DEBUG_MANAGER
        cout<<"  Fine intersection time: "<<std::get<0>(impact)<<endl;
#endif
        throw Contactor(impact);
      }
      else
      cout<<"*** NO INTERSECTION BETWEEN BOUNDING SPHERES ***"<<endl;

      
      time_type tstar = resolve_time(it1, it2);
      
#ifdef DEBUG_MANAGER
      
      if (tstar == inf)
        cout<<"    Leaves found that DO NOT collide"<<endl;
      else {
        cout<<"    Leaves found that collide at time "<<(tstar)<<endl;
        cout<<"    Modifying timer and time step..."<<endl;
      }
#endif
      
      if (tstar == inf) {
        cout<<"    Leaves found that DO NOT collide"<<endl;
        return;
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
  
  
  friend std::ostream& operator<<(std::ostream& os, const Contact_model_manager& mm) {
    
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

#endif /* __AKANTU_MODEL_MANAGER_HH__ */
