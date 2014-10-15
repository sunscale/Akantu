/**
 * @file   element_class_igfem.hh
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu May 22 13:50:46 2014
 *
 * @brief  Specialization for interface-enriched finite elements
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "element_class.hh"

#ifndef __AKANTU_ELEMENT_CLASS_IGFEM_HH__
#define __AKANTU_ELEMENT_CLASS_IGFEM_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
#define AKANTU_DEFINE_IGFEM_ELEMENT_CLASS_PROPERTY_(elem_type,		\
						    regular_el_type,	\
						    is_sub)		\
  static const GeometricalType geometrical_type = ElementClassProperty<regular_el_type>::geometrical_type; \
  static const InterpolationType interpolation_type = ElementClassProperty<regular_el_type>::interpolation_type; \
  static const ElementType regular_element_type = regular_el_type;	\
    static const ElementKind element_kind = _ek_igfem;			\
    static const UInt spatial_dimension = ElementClassProperty<regular_el_type>::spatial_dimension; \
    static const GaussIntergrationType gauss_integration_type = ElementClassProperty<regular_el_type>::gauss_integration_type; \
    static const UInt minimal_integration_order = ElementClassProperty<regular_el_type>::minimal_integration_order; \
    static const bool is_subelement = is_sub;


/// define ElementClassProperty for parent elements
#define AKANTU_DEFINE_IGFEM_PARENT_ELEMENT_CLASS_PROPERTY(elem_type,	\
							  regular_el_type) \
  template<>								\
  struct ElementClassProperty<elem_type> {				\
  AKANTU_DEFINE_IGFEM_ELEMENT_CLASS_PROPERTY_(elem_type,		\
					      regular_el_type,		\
					      false)			\
  }

/// define ElementClassProperty for subelements
#define AKANTU_DEFINE_IGFEM_SUBELEMENT_CLASS_PROPERTY(elem_type,	\
						      regular_el_type,	\
						      parent_el_type)	\
  template<>								\
  struct ElementClassProperty<elem_type> {				\
  AKANTU_DEFINE_IGFEM_ELEMENT_CLASS_PROPERTY_(elem_type,		\
					      regular_el_type,		\
					      true)			\
  static const ElementType parent_element_type = parent_el_type;	\
  }

/* -------------------------------------------------------------------------- */
/* ElementClass                                                               */
/* -------------------------------------------------------------------------- */
template<ElementType element_type>
class ElementClass<element_type, _ek_igfem> :
  public ElementClass<ElementClassProperty<element_type>::regular_element_type> {
protected:
  typedef ElementClassProperty<element_type> element_property;
  typedef ElementClass<element_property::regular_element_type> regular_element_class; 
  typedef typename regular_element_class::geometrical_element geometrical_element;
  typedef typename regular_element_class::interpolation_element interpolation_element;

  typedef typename interpolation_element::interpolation_property interpolation_property;
};

/* -------------------------------------------------------------------------- */

__END_AKANTU__

#endif /* __AKANTU_ELEMENT_CLASS_IGFEM_HH__ */
