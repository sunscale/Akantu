/**
 * @file   solid_mechanics_model_element.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue May 07 2013
 * @date last modification: Tue May 13 2014
 *
 * @brief  elements for solid mechanics models
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


#ifndef AKANTU_SOLID_MECHANICS_MODEL_ELEMENT_HH
#define AKANTU_SOLID_MECHANICS_MODEL_ELEMENT_HH

#include <array/expr.hpp>

#include "mesh.hh"
#include "aka_error.hh"
#include "solid_mechanics_model.hh"
#include "aka_bounding_box.hh"


__BEGIN_AKANTU__


typedef array::vector_type<Real> vector_type;
typedef array::matrix_type<Real> matrix_type;




template <>
class ModelElement<SolidMechanicsModel> : public Element {
  
public:
  
  
  typedef Element element_base;
  
  SolidMechanicsModel *model_;       //!< Reference to model
  UInt *connectivity_;               //!< Ponter to connectivity array
  
  
  typedef SolidMechanicsModel model_type;
  
  ModelElement() : element_base(), model_(nullptr), connectivity_() {}
  
  ModelElement(SolidMechanicsModel& m, Element &el)
  : element_base(el), model_(&m) {
    connectivity_ = &m.getMesh().getConnectivity(this->type, this->ghost_type)(this->element);
  }
  
  ModelElement(SolidMechanicsModel& m, ElementType type, UInt id, GhostType gt = _not_ghost)
  : element_base(type, id, gt), model_(&m) {
    connectivity_ = &m.getMesh().getConnectivity(type, gt)(id);
  }

  ModelElement(const ModelElement& p) : element_base(static_cast<element_base>(p)), model_(p.model_), connectivity_(p.connectivity_) {}
  
  
  model_type& model() { return *model_; }
  
  UInt numNodes() const
  { return Mesh::getNbNodesPerElement(type); }
  
  template <class element_type>
  bool shareNodes(element_type& el) {
    
    for (UInt i=0; i<Mesh::getNbNodesPerElement(type); ++i)
      for (UInt j=0; j<Mesh::getNbNodesPerElement(el.type); ++j)
        if (connectivity_[i] == el.connectivity_[j])
          return true;
    return false;
  }
  
  template <class element_type>
  UInt shareNode(element_type& el) {
    
    for (UInt i=0; i<Mesh::getNbNodesPerElement(type); ++i)
      for (UInt j=0; j<Mesh::getNbNodesPerElement(el.type); ++j)
        if (connectivity_[i] == el.connectivity_[j])
          return connectivity_[i];
    return -1;
  }
  
  
  UInt& node(UInt n) {
    AKANTU_DEBUG_ASSERT(n < Mesh::getNbNodesPerElement(type),
                        "Node "<<n<<" is larger than element number of nodes: "<<Mesh::getNbNodesPerElement(type));
    return connectivity_[n];
  }
  
  UInt node(UInt n) const {
    AKANTU_DEBUG_ASSERT(n < Mesh::getNbNodesPerElement(type),
                        "Node "<<n<<" is larger than element number of nodes: "<<Mesh::getNbNodesPerElement(type));
    return connectivity_[n];
  }
  
  
  // vector of pointers to nodes' first coordinates
  std::vector<const Real*> coordinates() {
    
    UInt nb_nodes = Mesh::getNbNodesPerElement(this->type);
    const Array<Real> &position = model_->getCurrentPosition();
    std::vector<const Real*> coord(nb_nodes);
    for (size_t i=0; i<nb_nodes; ++i)
      coord[i] = &position(connectivity_[i]);
    return coord;
  }
  
  // barycenter
  vector_type barycenter() const {
    
    typedef typename vector_type::value_type value_type;
    
    UInt nb_nodes = Mesh::getNbNodesPerElement(this->type);
    const Array<Real> &position = model_->getCurrentPosition();
    vector_type sum(model_->getSpatialDimension());
    for (size_t i=0; i<nb_nodes; ++i) {
      Real * p = const_cast<Real*>(&position(connectivity_[i]));
      sum += vector_type(model_->getSpatialDimension(), p);
    }
    
    return (1./static_cast<value_type>(nb_nodes)) * sum;
  }
  
  //! Returns the closest point to an element and the element normal
  vector_type normal() {
    
    vector_type n;
    auto coord = coordinates();
    
    switch (type) {
      case _segment_2:
      {
        // create points from segment
        Point<2> x(coord[0]);
        Point<2> y(coord[1]);        
        Point<2> t = y-x;
        n = vector_type(2);
        // normal2 normalizes the normal so that its magnitude is 1
        Math::normal2(&t[0], &n[0]);
        
        break;
      }
      case _triangle_3:
      {
        Point<3> x(coord[0]);
        Point<3> y(coord[1]);
        Point<3> z(coord[2]);
        Point<3> t1 = y-x;
        Point<3> t2 = z-x;
        n = vector_type(3);
        Math::vectorProduct3(&t1[0], &t2[0], &n[0]);
        Math::normalize3(&n[0]);
        
        break;
      }
      default:
        AKANTU_DEBUG_ERROR("No element type found");
    }
    return n;
  }
  
  
  
  template <int dim>
  BoundingBox<dim> boundingBox()  {
    
    typedef typename BoundingBox<dim>::point_type point_type;
    
    assert(dim == model_->getSpatialDimension());
    BoundingBox<dim> bb;
    UInt nb_nodes = Mesh::getNbNodesPerElement(this->type);
    const Array<Real> &position = model_->getCurrentPosition();
    for (size_t i=0; i<nb_nodes; ++i) {
      point_type p(&position(connectivity_[i]));
      bb += p;
    }
    return bb;
  }
  
  template <int dim>
  Point<dim> point(UInt nid) {
    
    AKANTU_DEBUG_ASSERT(dim == model_->getSpatialDimension(),
                        "Point and model dimensions do not match");
    UInt n = node(nid);
    const Array<Real> &position = model_->getCurrentPosition();
    Real * p = const_cast<Real*>(&position(n));
    return Point<dim>(p);
  }

  
  // mass
  vector_type getMass(UInt nid) {
    UInt n = node(nid);
    Array<Real> &mass = model_->getMass();
    return vector_type(model_->getSpatialDimension(), &mass(n));
  }
  
  // mass for const objects
  const vector_type getMass(UInt nid) const {
    UInt n = node(nid);
    Array<Real> &mass = model_->getMass();
    return vector_type(model_->getSpatialDimension(), &mass(n));
  }
  
  // initial coordinates
  vector_type getInitialCoordinates(UInt nid) {
    UInt n = node(nid);
    Array<Real> &coord = model_->getMesh().getNodes();
    return vector_type(model_->getSpatialDimension(), &coord(n));
  }
  
  const vector_type getInitialCoordinates(UInt nid) const {
    UInt n = node(nid);
    Array<Real> &coord = model_->getMesh().getNodes();
    return vector_type(model_->getSpatialDimension(), &coord(n));
  }
  
  
  // displacement
  vector_type getDisplacement(UInt nid) {
    UInt n = node(nid);
    Array<Real> &displacement = model_->getDisplacement();
    return vector_type(model_->getSpatialDimension(), &displacement(n));
  }
  
  // displacement for const objects
  const vector_type getDisplacement(UInt nid) const {
    UInt n = node(nid);
    Array<Real> &displacement = model_->getDisplacement();
    return vector_type(model_->getSpatialDimension(), &displacement(n));
  }
  
  // velocity
  vector_type getVelocity(UInt nid) {
    UInt n = node(nid);
    Array<Real> &velocity = model_->getVelocity();
    return vector_type(model_->getSpatialDimension(), &velocity(n));
  }
  
  // velocity for const objects
  const vector_type getVelocity(UInt nid) const {
    UInt n = node(nid);
    Array<Real> &velocity = model_->getVelocity();
    return vector_type(model_->getSpatialDimension(), &velocity(n));
  }
  
  // acceleration
  vector_type getAcceleration(UInt nid) {
    UInt n = node(nid);
    Array<Real> &acceleration = model_->getAcceleration();
    return vector_type(model_->getSpatialDimension(), &acceleration(n));
  }
  
  // acceleration for const objects
  const vector_type getAcceleration(UInt nid) const {
    UInt n = node(nid);
    Array<Real> &acceleration = model_->getAcceleration();
    return vector_type(model_->getSpatialDimension(), &acceleration(n));
  }
  
  // position (location + displacement)
  vector_type getCurrentPosition(UInt nid) {
    UInt n = node(nid);
    const Array<Real> &position = model_->getCurrentPosition();
    Real * p = const_cast<Real*>(&position(n));
    return vector_type(model_->getSpatialDimension(), p);
  }
  
  // position (location + displacement) for const objects
  const vector_type getCurrentPosition(UInt nid) const {
    UInt n = node(nid);
    const Array<Real> &position = model_->getCurrentPosition();
    Real * p = const_cast<Real*>(&position(n));
    return vector_type(model_->getSpatialDimension(), p);
  }
  
  // residual
  vector_type getResidual(UInt nid) {
    UInt n = node(nid);
    const Array<Real> & residual = model_->getResidual();
    Real * p = const_cast<Real*>(&residual(n));
    return vector_type(model_->getSpatialDimension(), p);
  }
  
  // residual for const objects
  const vector_type getResidual(UInt nid) const {
    UInt n = node(nid);
    const Array<Real> & residual = model_->getResidual();
    Real * p = const_cast<Real*>(&residual(n));
    return vector_type(model_->getSpatialDimension(), p);
  }
  
  // momentum
  vector_type getMomentum(UInt nid) {
    UInt n = node(nid);
    Array<Real> &velocity = model_->getVelocity();
    Array<Real> &mass = model_->getMass();
    vector_type p(model_->getSpatialDimension());
    for (size_t i=0; i<p.size(); ++i)
      p[i] = mass(n,i)*velocity(n,i);
    return p;
  }
  
  // momentum for const objects
  const vector_type getMomentum(UInt nid) const {
    return const_cast<ModelElement&>(*this).getMomentum(nid);
  }
  
};


//! Returns the closest point to an element and the element normal
template <class point_type, class element_type>
std::pair<point_type, vector_type> closest_point_to_element(const point_type& p, element_type& el) {
  
  point_type r;
  vector_type n;
  
  switch (el.type) {
    case _segment_2:
    {
      // create points from segment
      auto coord = el.coordinates();
      assert (coord.size() == 2);
      point_type x(coord[0]);
      point_type y(coord[1]);
      
      r = closest_point_to_segment(p, x, y);
      
      point_type t = y-x;
      n = vector_type(2);
      Math::normal2(&t[0], &n[0]);
      
      break;
    }
      
    default:
      AKANTU_DEBUG_ERROR("No element type found");
  }
  
  return std::make_pair(r, n);
}





template <class point_type, class pair_type>
bool balance(Real Dt, int id, const pair_type& r, ModelElement<SolidMechanicsModel>& sel, ModelElement<SolidMechanicsModel>& mel) {
  
  
//  constexpr int dim = point_type::dim();
  
  // treat first element as slave
  auto coord = sel.coordinates();
  
  // create point
  point_type p(coord[id]);
  
  // else find closest distance from p to contacting element c2
  //  std::pair<point_type, vector_type> r = closest_point_to_element(p, mel);
  const point_type& q = r.first;
  const vector_type& n = r.second;
  
  // get distance from current position
  Real delta = sqrt((q-p).sq_norm());
  
  auto mass = sel.getMass(id)[0];
  
  // compute force at slave node
  vector_type N = 2 * delta * mass / pow(Dt,2.) * n;
    
  // update residual and velocity for slave
  vector_type sr = sel.getResidual(id);
  vector_type m = sel.getMass(id);
  vector_type v = sel.getVelocity(id);
  vector_type a = sel.getAcceleration(id);
  for (size_t i=0; i<N.size(); ++i) {
    //    sr[i] = N[i];
    v[i] += N[i]/m[i] * Dt;
    //    a[i] -= N[i]/m[i];
  }
  
  // set location of slave node
  auto xs = sel.getDisplacement(id);
  auto oc = sel.getInitialCoordinates(id);
  for (size_t i=0; i<xs.size(); ++i)
    xs[i] = q[i] - oc[i];
  
  // balance with slave forces
  switch (mel.type) {
    case _segment_2:
    {
      // create points from segment
      auto coord = mel.coordinates();
      assert (coord.size() == 2);
      point_type X1(coord[0]);
      point_type X2(coord[1]);
      
      // get weights for distribution of loads
      Real alpha = (q - X1).sq_norm() / (X2 - X1).sq_norm();
      
      //get vectors
      vector_type R1 = mel.getResidual(0);
      vector_type R2 = mel.getResidual(1);
      vector_type V1 = mel.getVelocity(0);
      vector_type V2 = mel.getVelocity(1);
      //      vector_type A1 = mel.getAcceleration(0);
      //      vector_type A2 = mel.getAcceleration(1);
      vector_type M1 = mel.getMass(0);
      vector_type M2 = mel.getMass(1);
      
      for (size_t i=0; i<N.size(); ++i) {
        Real r1 = (alpha - 1)*N[i];
        Real r2 = -alpha * N[i];
        //        R1[i] = r1;
        //        R2[i] = r2;
        V1[i] += r1/M1[i] * Dt;
        V2[i] += r2/M2[i] * Dt;
        //        A1[i] -= R1[i]/M1[i];
        //        A2[i] -= R2[i]/M2[i];
      }
      
      break;
    }
      
    default:
      AKANTU_DEBUG_ERROR("No element type found");
  }
  
  
  return true;
  
}



template <class point_type, class element_container>
std::pair<point_type, vector_type> commonPonit(const element_container& els) {
  
  typedef std::pair<point_type, vector_type> pair_type;
  
  // balance with slave forces
  switch (els.size()) {
    case 2:
    {
      
      auto it = els.begin();
      auto el1 = *it++;
      auto el2 = *it;
      
      // get normals
      auto n1 = el1.normal();
      auto n2 = el2.normal();
      
      // average normal
      vector_type n = n1 + n2;
      n *= (1/n.norm());
      
      // get common point
      assert(el1.shareNodes(el2));
      
      UInt c = el1.shareNode(el2);
      
      assert(el1.shareNode(el2) == el2.shareNode(el1));
      
      // set location of slave node
      Array<Real> &displacement = el1.model().getDisplacement();
      Array<Real> &coord = el1.model().getMesh().getNodes();
      
      point_type p;
      for (int i=0; i<point_type::dim(); ++i)
        p[i] = coord(c,i) + displacement(c,i);
      
      return std::make_pair(p,n);
      
      break;
    }
      
    default:
      AKANTU_DEBUG_ERROR("Invalid size");
  }
  
  
  return pair_type();
}



template <class point_type, class element_type>
bool penetrates(point_type& r, element_type& el) {
  
  // balance with slave forces
  switch (el.type) {
    case _segment_2:
    {
      
      // create points from master
      auto coord = el.coordinates();
      assert (coord.size() == 2);
      point_type p(coord[0]);
      point_type q(coord[1]);
      
      return left_turn(p,q,r) > 0 && has_projection(r, p, q);
      break;
    }
      
    default:
      AKANTU_DEBUG_ERROR("No element type found");
  }
  
  // should never get here
  assert(false);
  return false;
}

__END_AKANTU__


#endif /* AKANTU_SOLID_MECHANICS_MODEL_ELEMENT_HH */
