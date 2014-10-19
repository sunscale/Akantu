/**
 * @file   aka_optimize.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  Objects that can be used to carry out optimization
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

#ifndef __AKANTU_OPTIMIZE_HH__
#define __AKANTU_OPTIMIZE_HH__

#include "aka_config.hh"
#include "aka_common.hh"
#include "aka_point.hh"
#include "solid_mechanics_model.hh"

#include <iostream>
#include <nlopt.hpp>
#include <array/expr.hpp>


//#define DEBUG_OPTIMIZE 1


__BEGIN_AKANTU__

using std::cout;
using std::endl;

typedef array::Array<1,Real> vector_type;
typedef array::Array<2,Real> matrix_type;

std::ostream& operator<<(std::ostream&, nlopt::result);



enum Optimizator_type { Min_t, Max_t };


//! Class used for optimization
/*! This class is a convenience object that inherits from nlopt::opt and carries
 * some routines that are common to a nonlinear optimization. The objects sets 
 * the optimization algorithm as nlopt::LD_SLSQP, a sequential quadratic programming
 * (SQP) algorithm for nonlinearly constrained gradient-based optimization
 * (supporting both inequality and equality constraints), based on the 
 * implementation by Dieter Kraft.
 */
class Optimizator : public nlopt::opt {
  
  typedef std::vector<Real> point_type;
  typedef nlopt::opt base_type;
  
  point_type& x_;     //!< Initial guess for the optimization
  Real val_;          //!< Function value
  
public:
  
  //! Parameter constructor that takes an initial guess and a functor
  template <class functor_type>
  Optimizator(point_type& x0,
              functor_type& fn,
              Optimizator_type t = Min_t,
              nlopt::algorithm alg = nlopt::LD_SLSQP) :
  base_type(alg, x0.size()), x_(x0) {
    
    if (t == Min_t)
      this->set_min_objective(functor_type::wrap, &fn);
    else
      this->set_max_objective(functor_type::wrap, &fn);
    
    this->set_xtol_rel(1e-4);
  }
  
  //! Carry out the optimization and print result
  Real result() {
    
    optimize(x_, val_);
    
    cout<<"Optimum value found at location";
    for (size_t i=0; i<x_.size(); ++i)
      cout<<" "<<x_[i];
    cout<<"\nFunction value: "<<val_;
    
    return val_;
  }
  
};

//! Traits class used as a base class for the Distance_minimizator class template
template <ElementType>
struct Distance_minimizator_traits;


//! Partial template specialization for a segment
template <>
struct Distance_minimizator_traits<_segment_2> {
  
  //! Set lower and upper bounds for the master coordinate
  static void set_bounds(nlopt::opt &opt) {
    
    std::vector<Real> lb(1,-1), ub(1,1);
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
  }
  
  //! Select the start point between the center or the ends of the segmet
  template <class object_type>
  static void start(object_type &obj) {
    
    Real min = std::numeric_limits<Real>::infinity();
    std::vector<Real> xstart = { -1., 0., 1. }; // check center and extremes of element
    int idx = -1;
    for (size_t i=0; i<xstart.size(); ++i) {
      obj.xi_[0] = xstart[i];
      std::vector<Real> grad; // empty vector
      Real new_dist = obj(obj.xi_, grad);
      if (new_dist < min) {
        min = new_dist;
        idx = i;
      }
    }
    obj.xi_[0] = xstart[idx];
  }
};


//! Constraint function for a triangle
static Real constrain_triangle_3(const std::vector<Real> &xi, std::vector<Real> &grad, void *data) {
  if (!grad.empty()) {
    grad[0] = 1.;
    grad[1] = 1.;
  }
  return (xi[0] + xi[1] - 1.);
}


//! Helper class used for the definition of triangle partial template specializations
struct Triangle_minimizator_traits {
  
  //! Set lower and upper bounds for the master coordinate
  static void set_bounds(nlopt::opt &opt) {
    
    std::vector<Real> lb(2, Real());
    opt.set_lower_bounds(lb);
    opt.add_inequality_constraint(constrain_triangle_3, NULL, 1e-4);
  }
  
  //! Select the start point between the center or the vertices of the triangle
  template <class object_type>
  static void start(object_type &obj) {
    
    Real min = std::numeric_limits<Real>::infinity();
    Real xstart[4][2] = { {0.,0.}, {1.,0.}, {0.,1.}, {1./3.,1./3.} }; // check center and corners of element
    int idx = -1;
    for (int i=0; i<4; ++i) {
      obj.xi_[0] = xstart[i][0];
      obj.xi_[1] = xstart[i][1];
      std::vector<Real> grad; // empty vector
      Real new_dist = obj(obj.xi_, grad);
      if (new_dist < min) {
        min = new_dist;
        idx = i;
      }
    }
    obj.xi_[0] = xstart[idx][0];
    obj.xi_[1] = xstart[idx][1];
  }
};


//! Partial template specialization for a 3-node triangle
template <>
struct Distance_minimizator_traits<_triangle_3> : public Triangle_minimizator_traits {};


//! Partial template specialization for a 6-node triangle
template <>
struct Distance_minimizator_traits<_triangle_6> : public Triangle_minimizator_traits {};



/*! The Distance_minimizator class template can be used to obtain the closest point
 * to a a finite element using the NLopt optimization library.
 * \tparam d - The dimension of the problem
 * \tparam element_policy - The element type to which the closest point is sought
 * The class inherits from Distance_minimizator_traits to take care of functionality
 * specific to elements of certain type. The code can be used for elements of type 
 * _segment_2, _triangle_3, and _triangle_6. The optimization stage is done during the
 * constructor by calling the function constructor_common. The closest point is then
 * obtained by calling the function point.
 */
template <int d, ElementType element_policy>
class Distance_minimizator : public Distance_minimizator_traits<element_policy> {
  
  friend class Triangle_minimizator_traits;
  friend class Distance_minimizator_traits<element_policy>;
  
  const UInt nb_nodes = ElementClass<element_policy>::getNbNodesPerElement();
  
  typedef Distance_minimizator_traits<element_policy> traits_type;
  
  typedef Point<d> point_type;
  
  nlopt::opt opt_;         //!< Optimizator reference
  std::vector<Real> xi_;   //!< Master coordinate closest to point
  vector_type p_;          //!< Point to which the distance is minimized
  matrix_type XX_;         //!< Triangle coordinates
  UInt counter_;           //!< Optimization iteration counter
  Real fmin_;              //!< Minimum distance value
  
  //! Common function called during construction to carry out the minimization
  void constructor_common() {
    
    traits_type::set_bounds(opt_);
    
    opt_.set_min_objective(wrap, this);
    opt_.set_ftol_abs(1e-4);
    
    // compute start point
    traits_type::start(*this);
    
    // optimize
#ifdef DEBUG_OPTIMIZE
    nlopt::result result = opt_.optimize(xi_, fmin_);
    if (result > 0)
    cout<<"Optimium found in "<<counter_<<" iterations: "<<fmin_<<endl;
    cout<<"Point at master coordinate "<<xi_[0]<<": "<<point()<<endl;
    cout<<result<<endl;
#else
    opt_.optimize(xi_, fmin_);
#endif
  }
  
  public:
  
  //! Parameter constructor

  /*! This parameter constructor takes the point to which the minimum distance is
   * sought, and a container of points points that are the coordinates of the finite
   * element
   * \param r - Point coordinates
   * \param pts - Container of triangle points
   */
  template <class point_type, class point_container>
  Distance_minimizator(const point_type& p, const point_container& pts)
  : opt_(nlopt::LD_SLSQP, d-1), xi_(d-1), p_(d), XX_(nb_nodes,d), counter_() {
    
    // get triangle and point coordinates
    for (UInt i=0; i<d; ++i) {
      p_[i] = p[i];
      for (UInt j=0; j<nb_nodes; ++j)
      XX_(j,i) = pts[j][i];
    }
    
    // common constructor operation
    constructor_common();
  }
  
  //! Parameter constructor
  /*! This parameter constructor uses a pointer to the coordinates of the point, 
   * an Element pointer, and the SolidMechanicsModel.
   * \param r - Point coordinates
   * \param el - Finite element to which the distance is minimized
   * \param model - Solid mechanics model
   */
  Distance_minimizator(const Real *r, const Element *el, SolidMechanicsModel &model)
  : opt_(nlopt::LD_SLSQP, d-1), xi_(d-1), p_(d), XX_(nb_nodes,d), counter_() {
    
    Mesh& mesh = model.getMesh();
    const Array<Real> &X = model.getCurrentPosition();
    const Array<UInt> &conn = mesh.getConnectivity(el->type);
    for (UInt i=0; i<nb_nodes; ++i) {
      p_(i) = r[i];
      for (UInt j=0; j<d; ++j)
        XX_(i,j) = X(conn(el->element,i),j);
    }
    
    constructor_common();
  }
  
  
  vector_type operator()(const std::vector<Real> &xi)
  {
    vector_type N(nb_nodes);
    vector_type xi2(xi.size(), const_cast<Real*>(&xi[0]));
    ElementClass<element_policy>::computeShapes(xi2, N);
    return transpose(XX_)*N;
  }
  
  //! Evaluation of the function and its gradient
  Real operator()(const std::vector<Real> &xi, std::vector<Real> &grad)
  {
    // increment function evaluation counter
    ++counter_;
    
    vector_type x = (*this)(xi);
    vector_type diff = x-p_;
    
    if (!grad.empty()) {
      
      // compute shape function derivatives
      matrix_type DN(d-1,nb_nodes);
      vector_type xi2(xi.size(), const_cast<Real*>(&xi[0]));
      
      ElementClass<element_policy>::computeDNDS(xi2, DN);
      DN = transpose(DN);
      
      // compute jacobian
      matrix_type J = transpose(XX_)*DN;
      
      // compute function gradient
      vector_type gradF = transpose(J) * diff;
      for (UInt i=0; i<gradF.size(); ++i)
      grad[i] = gradF[i];
      
    }
    // return function value
    return 0.5 * transpose(diff)*diff;
  }
  
  //! Return point at current master coordinate
  point_type point() {
    vector_type x = (*this)(xi_);
    point_type p;
    for (UInt i=0; i<x.size(); ++i)
      p[i] = x[i];
    return p;
  }
  
  //! Return the number of function evaluations
  UInt iterations() const
  { return counter_; }
  
  //! Return the master coordinate
  const std::vector<Real>& master_coordinates()
  {  return xi_; }
  
  //! Function that is used to have this class working with nlopt
  static double wrap(const std::vector<double> &x, std::vector<double> &grad, void *data) {
    return (*reinterpret_cast<Distance_minimizator<d,element_policy>*>(data))(x, grad); }
  
};


__END_AKANTU__


#endif /* __AKANTU_OPTIMIZE_HH__ */


