/**
 * @file   discretization.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  contact discretization classes
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

#ifndef __AKANTU_DISCRETIZATION_HH__
#define __AKANTU_DISCRETIZATION_HH__

//#include <cassert>
//#include <stdio.h>
//#include <stdlib.h>

#include <map>


#include "aka_common.hh"
//#include "aka_types.hh"
//#include "contact_common.hh"
//#include "solid_mechanics_model.hh"
#include "solid_mechanics_model_element.hh"
//#include "friction.hh"

using std::endl;
using std::cout;

__BEGIN_AKANTU__



class Discretization_base {
  //! Support the visitor design pattern
//  virtual void accept(Contact_discretization_visitor&) = 0;
};


template <Discretization_type d>
class Contact_discretization;

//template <>
//class Contact_discretization<Node_to_node> : public Discretization_base {
//  
//  template <int> friend class Contact_scheme;
//  
//private:
//  
//  // pair structure
//  struct N2N_pair {
//    
//    typedef UInt node_id;
//    
//    node_id master_, slave_;
//    bool contact_, stick_;
//    vector_type n_;
//  };
//  
//  
//};



//template <>
//class Contact_discretization<Node_to_segment> : public Discretization_base {
//
//  
//  typedef SolidMechanicsModel model_type;
//  typedef ModelElement <model_type> element_type;
//  
//  typedef std::map <UInt, element_type> slave_master_map;
//  typedef std::map <UInt, Real> real_map;
//  
//  typedef typename real_map::iterator real_iterator;
//  typedef std::map <UInt, Real> gap_map;
//  
//  slave_master_map sm_;
//  real_map areas_, gaps_;
//  
//public:
//  
//  //! Add slave
//	void addSlave(UInt s)
//  { sm_[s] = element_type(); }
//  
//	//! Add slave-master pair
//	void addPair(UInt s, element_type el)
//  { sm_[s] = el; }
//
//  //! Add area to a slave node
//	void addArea(UInt n, Real a)
//  { if (a != 0.) areas_[n] = a; }
//
//  
////  //! Member function to support the visitor design pattern.
////  void accept(Contact_discretization_visitor& guest)
////  { guest.GenericVisit(*this); }
//
//};




__END_AKANTU__

#endif /* __AKANTU_DISCRETIZATION_HH__ */
