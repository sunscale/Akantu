/**
 * @file   internal_field_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  Material internal properties
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
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_INTERNAL_FIELD_TMPL_HH__
#define __AKANTU_INTERNAL_FIELD_TMPL_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<typename T>
InternalField<T>::InternalField(const ID & id, Material & material) :
  ElementTypeMapArray<T>(id, material.getID(), material.getMemoryID()),
  material(material),
  fem(material.getModel().getFEEngine()),
  element_filter(material.getElementFilter()),
  default_value(T()),
  spatial_dimension(material.getModel().getSpatialDimension()),
  element_kind(_ek_regular),
  nb_component(0),
  is_init(false),
  previous_values(NULL) {
}

/* -------------------------------------------------------------------------- */
template<typename T>
InternalField<T>::InternalField(const ID & id, Material & material, FEEngine & fem,
				const ElementTypeMapArray<UInt> & element_filter) :
  ElementTypeMapArray<T>(id, material.getID(), material.getMemoryID()),
  material(material),
  fem(fem),
  element_filter(element_filter),
  default_value(T()),
  spatial_dimension(material.getModel().getSpatialDimension()),
  element_kind(_ek_regular),
  nb_component(0),
  is_init(false),
  previous_values(NULL) {
}

/* -------------------------------------------------------------------------- */
template<typename T>
InternalField<T>::InternalField(const ID & id, const InternalField<T> & other) :
  ElementTypeMapArray<T>(id, other.material.getID(), other.material.getMemoryID()),
  material(other.material),
  fem(other.fem),
  element_filter(other.element_filter),
  default_value(other.default_value),
  spatial_dimension(other.spatial_dimension),
  element_kind(other.element_kind),
  nb_component(other.nb_component),
  is_init(false),
  previous_values(NULL) {

  AKANTU_DEBUG_ASSERT(other.is_init, "Cannot create a copy of a non initialized field");
  this->internalInitialize(this->nb_component);
}


/* -------------------------------------------------------------------------- */
template<typename T>
InternalField<T>::~InternalField() {
  if(this->is_init) {
    this->material.unregisterInternal(*this);
  }

  delete previous_values;
}

/* -------------------------------------------------------------------------- */
template<typename T>
void InternalField<T>::initialize(UInt nb_component) {
  internalInitialize(nb_component);
}

/* -------------------------------------------------------------------------- */
template<typename T>
void InternalField<T>::initializeHistory() {
  if(!previous_values)
    previous_values = new InternalField<T>("previous_" + this->getID(), *this);
}

/* -------------------------------------------------------------------------- */
template<typename T>
void InternalField<T>::resize() {
  if(!this->is_init) return;

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    typename ElementTypeMapArray<UInt>::type_iterator it
      = element_filter.firstType(spatial_dimension, gt, element_kind);
    typename ElementTypeMapArray<UInt>::type_iterator end
      = element_filter.lastType(spatial_dimension, gt, element_kind);
    for(; it != end; ++it) {
      UInt nb_element = element_filter(*it, gt).getSize();

      UInt nb_quadrature_points = fem.getNbQuadraturePoints(*it, gt);
      UInt new_size = nb_element * nb_quadrature_points;

      UInt old_size = 0;
      Array<T> * vect = NULL;

      if(this->exists(*it, gt)) {
	vect = &(this->operator()(*it, gt));
	old_size = vect->getSize();
	vect->resize(new_size);
      } else {
	vect = &(this->alloc(nb_element * nb_quadrature_points, nb_component, *it, gt));
      }

      this->setArrayValues(vect->storage() + old_size * vect->getNbComponent(),
			   vect->storage() + new_size * vect->getNbComponent());
    }
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
void InternalField<T>::setDefaultValue(const T & value) {
  this->default_value = value;
  this->reset();
}


/* -------------------------------------------------------------------------- */
template<typename T>
void InternalField<T>::reset() {
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    typename ElementTypeMapArray<T>::type_iterator it  = this->firstType(spatial_dimension, gt, element_kind);
    typename ElementTypeMapArray<T>::type_iterator end = this->lastType(spatial_dimension, gt, element_kind);
    for(; it != end; ++it) {
      Array<T> & vect = this->operator()(*it, gt);
      vect.clear();
      this->setArrayValues(vect.storage(),
			   vect.storage() + vect.getSize() * vect.getNbComponent());
    }
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
void InternalField<T>::internalInitialize(UInt nb_component) {
  if(!this->is_init) {
    this->nb_component = nb_component;
    for(UInt g = _not_ghost; g <= _ghost; ++g) {
      GhostType gt = (GhostType) g;
      typename ElementTypeMapArray<UInt>::type_iterator it
	= element_filter.firstType(spatial_dimension, gt, element_kind);
      typename ElementTypeMapArray<UInt>::type_iterator end
	= element_filter.lastType(spatial_dimension, gt, element_kind);

      for(; it != end; ++it) {
	UInt nb_element = element_filter(*it, gt).getSize();
	UInt nb_quadrature_points = fem.getNbQuadraturePoints(*it, gt);
	if(this->exists(*it, gt))
	  this->operator()(*it, gt).resize(nb_element * nb_quadrature_points);
	else
	  this->alloc(nb_element * nb_quadrature_points, nb_component, *it, gt);
      }
    }

    this->material.registerInternal(*this);
    this->is_init = true;
  }
  this->reset();

  if(previous_values) previous_values->internalInitialize(nb_component);
}

/* -------------------------------------------------------------------------- */
template<typename T>
void InternalField<T>::setArrayValues(T * begin, T * end) {
  for(; begin < end; ++begin) *begin = default_value;
}

/* -------------------------------------------------------------------------- */
template<typename T>
void InternalField<T>::saveCurrentValues() {
  AKANTU_DEBUG_ASSERT(previous_values != NULL,
		      "The history of the internal " << this->getID()
		      << " has not been activated");

  if(!is_init) return;

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    typename ElementTypeMapArray<T>::type_iterator it  = this->firstType(spatial_dimension, gt, element_kind);
    typename ElementTypeMapArray<T>::type_iterator end = this->lastType(spatial_dimension, gt, element_kind);
    for(; it != end; ++it) {
      this->previous_values->operator()(*it, gt).copy(this->operator()(*it, gt));
    }
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
void InternalField<T>::removeQuadraturePoints(const ElementTypeMapArray<UInt> & new_numbering) {
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    ElementTypeMapArray<UInt>::type_iterator it  = new_numbering.firstType(_all_dimensions, gt, _ek_not_defined);
    ElementTypeMapArray<UInt>::type_iterator end = new_numbering.lastType(_all_dimensions, gt, _ek_not_defined);
    for (; it != end; ++it) {
      ElementType type = *it;
      if(this->exists(type, gt)){
	const Array<UInt> & renumbering = new_numbering(type, gt);

	Array<T> & vect = this->operator()(type, gt);

	UInt nb_quad_per_elem = fem.getNbQuadraturePoints(type, gt);
	UInt nb_component = vect.getNbComponent();

	Array<T> tmp(renumbering.getSize()*nb_quad_per_elem, nb_component);

	AKANTU_DEBUG_ASSERT(tmp.getSize() == vect.getSize(), "Something strange append some mater was created from nowhere!!");

	AKANTU_DEBUG_ASSERT(tmp.getSize() == vect.getSize(), "Something strange append some mater was created or disappeared in "<< vect.getID() << "("<< vect.getSize() <<"!=" << tmp.getSize() <<") ""!!");

	UInt new_size = 0;
	for (UInt i = 0; i < renumbering.getSize(); ++i) {
	  UInt new_i = renumbering(i);
	  if(new_i != UInt(-1)) {
	    memcpy(tmp.storage() + new_i * nb_component * nb_quad_per_elem,
		   vect.storage() + i * nb_component * nb_quad_per_elem,
		   nb_component * nb_quad_per_elem * sizeof(T));
	    ++new_size;
	  }
	}
	tmp.resize(new_size * nb_quad_per_elem);
	vect.copy(tmp);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
void InternalField<T>::printself(std::ostream & stream, unsigned int indent) const {
  stream << "InternalField [ " << this->getID();
#if !defined(AKANTU_NDEBUG)
  if(AKANTU_DEBUG_TEST(dblDump)) {
    stream << std::endl;
    InternalField<T>::printself(stream, indent + 3);
  } else {
#endif
    stream << " {"
	   << this->getData(_not_ghost).size() << " types - "
	   << this->getData(_ghost).size() << " ghost types"
	   << "}";
#if !defined(AKANTU_NDEBUG)
  }
#endif
  stream << " ]";
}

/* -------------------------------------------------------------------------- */
template<>
inline void ParsableParamTyped< InternalField<Real> >::parseParam(const ParserParameter & in_param) {
  ParsableParam::parseParam(in_param);
  Real r = in_param;
  param.setDefaultValue(r);
}

__END_AKANTU__

#endif /* __AKANTU_INTERNAL_FIELD_TMPL_HH__ */
