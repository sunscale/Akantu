/**
 * @file   element.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Wed Mar 13 2013
 *
 * @brief  contact element classes
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

#ifndef __AKANTU_CELEMENT_HH__
#define __AKANTU_CELEMENT_HH__

#include <iostream>

#include "aka_common.hh"
#include "aka_point.hh"
#include "solid_mechanics_model.hh"
#include "mesh.hh"
#include "contact_common.hh"
#if AKANTU_OPTIMIZATION
#  include "aka_optimize.hh"
#endif

__BEGIN_AKANTU__

using std::cout;
using std::endl;


// compute signed measure, used to find out if a node penetrates or not
template <int d>
bool check_penetration(UInt node, const Element* el, SolidMechanicsModel& model);

template <int d>
Point<d> minimize_distance(UInt node, const Element* el, SolidMechanicsModel& model);


template <int d>
Real distance(UInt node, const Element* el, SolidMechanicsModel& model) {

  typedef Point<d> point_type;
  
  const Array<Real> &x = model.getCurrentPosition();

  // get point of node
  point_type o(&x(node));

  // compute closest location within the element to master node
  point_type p = minimize_distance<d>(node, el, model);

  return sqrt((o-p).sq_norm());
}
  

template <int d>
class CElement {
  
  typedef Point<d> point_type;
  typedef const Element* element_pointer;
  
  UInt master_;                  //!< Slave node
  element_pointer element_;      //!< Master surface element
  SolidMechanicsModel& model_;
  
public:
  
  CElement(UInt n, element_pointer el, SolidMechanicsModel& model) : master_(n), element_(el), model_(model) {}
  
  //! Checks if there is penetration between the master node and the element
  bool penetrates() const { 
#ifdef DEBUG_CONTACT
    cout<<"*** INFO *** Checking penetration for"<<*this<<endl;
#endif
    return check_penetration<d>(master_, element_, model_);
  }
  
  void remove_penetration() const {
    
#ifdef DEBUG_CONTACT
    cout<<"*** INFO *** Removing penetration for master node "<<master_<<endl;
#endif
    
    cout<<"before calling minimize distance"<<endl;
    
    // compute closest location within the element to master node
    point_type o = minimize_distance<d>(master_, element_, model_);
    
    cout<<"after calling minimize"<<endl;
    cout<<"point obtained -> "<<o<<endl;
    
    // modify master node coordinates
    Array<Real> &X = model_.getDisplacement();
    const Array<Real> &coord = model_.getMesh().getNodes();
    for (UInt i=0; i<d; ++i)
      X(master_, i) = o[i] - coord(master_, i);
    
    cout<<"after modifying master node coords"<<endl;
  }
  
  
  //! Enable std output
  friend std::ostream& operator<<(std::ostream& os, const CElement& cel) {
    
    Mesh& mesh = cel.model_.getMesh();
    const Array<Real> &x = cel.model_.getCurrentPosition();
    
    UInt nb_nodes = mesh.getNbNodesPerElement((cel.element_)->type);
    const Array<UInt> &conn = mesh.getConnectivity((cel.element_)->type);
    
    os<<"    Contact element: "<<endl;
    os<<"      slave node: "<<cel.master_<<", coordinates: "<<point_type(&x(cel.master_))<<endl;
    os<<"      master surface element: "<<*cel.element_<<", with extreme nodes";
    
    for (int i=0; i<d; ++i)
      os<<" - "<<point_type(&x(conn((cel.element_)->element, i)));
    os<<endl;
    return os;
  }
};



__END_AKANTU__

#endif /* __AKANTU_CELEMENT_HH__ */
