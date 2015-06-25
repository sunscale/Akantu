/**
 * @file   dumper_igfem_element_iterator.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Helper class to return sub element information
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#ifndef __AKANTU_IGFEM_HELPEER_HH__
#define __AKANTU_IGFEM_HELPER_HH__
/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

struct IGFEMHelper {

  /// get the number of nodes for a given sub-element
  static UInt getNbNodesPerSubElement(const ElementType & type, const UInt sub_element) {
  UInt nb_nodes_per_sub_element = 0;
#define GET_NB_NODES_PER_SUB_ELEMENT(type)			\
  switch (sub_element) {					\
  case 0:							\
    nb_nodes_per_sub_element = ElementClass<ElementClassProperty<type>::sub_element_type_1>::getNbNodesPerInterpolationElement(); \
  case 1:							\
    nb_nodes_per_sub_element = ElementClass<ElementClassProperty<type>::sub_element_type_2>::getNbNodesPerInterpolationElement(); \
  }			
  AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(GET_NB_NODES_PER_SUB_ELEMENT);
#undef GET_NB_NODES_PER_SUB_ELEMENT
  return nb_nodes_per_sub_element;
  }

  /// get the connectivity for a given sub-element
  static UInt * getSubElementConnectivity(const ElementType & type, const UInt sub_element) {
  UInt * sub_element_connectivity = NULL;
#define GET_SUB_ELEMENT_CONNECTIVITY(type)				\
  sub_element_connectivity = ElementClass<type>::getSubElementConnectivity(sub_element);	
  AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(GET_SUB_ELEMENT_CONNECTIVITY);
#undef GET_SUB_ELEMENT_CONNECTIVITY
  return sub_element_connectivity;
  }

  /// get the sub-element type
  static ElementType getSubElementType(const ElementType & type, const UInt sub_element) {
    ElementType sub_type = _not_defined;
#define GET_SUB_ELEMENT_TYPE(type)						\
    switch (sub_element) {						\
    case 0:								\
      sub_type = ElementClassProperty<type>::sub_element_type_1;	\
    case 1:								\
      sub_type = ElementClassProperty<type>::sub_element_type_2;	\
    }
    AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(GET_SUB_ELEMENT_TYPE);    
#undef GET_SUB_ELEMENT_TYPE   
    return sub_type;
  }

};

__END_AKANTU__
/* -------------------------------------------------------------------------- */



#endif /* __AKANTU_IGFEM_HELPER_HH__ */
