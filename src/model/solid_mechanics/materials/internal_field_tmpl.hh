/**
 * @file   internal_field_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Material internal properties
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_INTERNAL_FIELD_TMPL_HH_
#define AKANTU_INTERNAL_FIELD_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
InternalFieldTmpl<Material, T>::InternalFieldTmpl(const ID & id, Material & material)
    : ElementTypeMapArray<T>(id, material.getID(), material.getMemoryID()),
      material(material), fem(&(material.getModel().getFEEngine())),
      element_filter(material.getElementFilter()),
      spatial_dimension(material.getModel().getSpatialDimension()) {}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
InternalFieldTmpl<Material, T>::InternalFieldTmpl(
    const ID & id, Material & material, FEEngine & fem,
    const ElementTypeMapArray<UInt> & element_filter)
    : ElementTypeMapArray<T>(id, material.getID(), material.getMemoryID()),
      material(material), fem(&fem), element_filter(element_filter),
      spatial_dimension(material.getSpatialDimension()) {}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
InternalFieldTmpl<Material, T>::InternalFieldTmpl(
    const ID & id, Material & material, UInt dim, FEEngine & fem,
    const ElementTypeMapArray<UInt> & element_filter)
    : ElementTypeMapArray<T>(id, material.getID(), material.getMemoryID()),
      material(material), fem(&fem), element_filter(element_filter),
      spatial_dimension(dim) {}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
InternalFieldTmpl<Material, T>::InternalFieldTmpl(const ID & id,
						  const InternalFieldTmpl<Material, T> & other)
    : ElementTypeMapArray<T>(id, other.material.getID(),
                             other.material.getMemoryID()),
      material(other.material), fem(other.fem),
      element_filter(other.element_filter), default_value(other.default_value),
      spatial_dimension(other.spatial_dimension),
      element_kind(other.element_kind), nb_component(other.nb_component) {

  AKANTU_DEBUG_ASSERT(other.is_init,
                      "Cannot create a copy of a non initialized field");
  this->internalInitialize(this->nb_component);
}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
InternalFieldTmpl<Material, T>::~InternalFieldTmpl() {
  if (this->is_init) {
    this->material.unregisterInternal(*this);
  }
}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
void InternalFieldTmpl<Material, T>::setFEEngine(FEEngine & fe_engine) {
  this->fem = &fe_engine;
}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
void InternalFieldTmpl<Material, T>::setElementKind(ElementKind element_kind) {
  this->element_kind = element_kind;
}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
void InternalFieldTmpl<Material, T>::initialize(UInt nb_component) {
  internalInitialize(nb_component);
}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
void InternalFieldTmpl<Material, T>::initializeHistory() {
  if (!previous_values)
    previous_values =
      std::make_unique<InternalFieldTmpl<Material, T>>("previous_" + this->getID(), *this);
}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
void InternalFieldTmpl<Material, T>::resize() {
  if (!this->is_init){
    return;
  }

  for (auto ghost : ghost_types) {
    for (const auto & type : this->filterTypes(ghost)) {
      UInt nb_element = this->element_filter(type, ghost).size();

      UInt nb_quadrature_points =
          this->fem->getNbIntegrationPoints(type, ghost);
      UInt new_size = nb_element * nb_quadrature_points;

      UInt old_size = 0;
      Array<T> * vect = nullptr;

      if (this->exists(type, ghost)) {
        vect = &(this->operator()(type, ghost));
        old_size = vect->size();
        vect->resize(new_size);
      } else {
        vect = &(this->alloc(nb_element * nb_quadrature_points, nb_component,
                             type, ghost));
      }

      this->setArrayValues(vect->storage() + old_size * vect->getNbComponent(),
                           vect->storage() + new_size * vect->getNbComponent());
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
void InternalFieldTmpl<Material, T>::setDefaultValue(const T & value) {
  this->default_value = value;
  this->reset();
}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
void InternalFieldTmpl<Material, T>::reset() {
  for (auto ghost_type : ghost_types){
    for (const auto & type : this->elementTypes(ghost_type)) {
      Array<T> & vect = (*this)(type, ghost_type);
      //vect.zero();
      this->setArrayValues(
          vect.storage(), vect.storage() + vect.size() * vect.getNbComponent());
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
void InternalFieldTmpl<Material, T>::internalInitialize(UInt nb_component) {
  if (!this->is_init) {
    this->nb_component = nb_component;

    for (auto ghost : ghost_types) {
      for (const auto & type : this->filterTypes(ghost)) {
        UInt nb_element = this->element_filter(type, ghost).size();
        UInt nb_quadrature_points =
            this->fem->getNbIntegrationPoints(type, ghost);
        if (this->exists(type, ghost)) {
          this->operator()(type, ghost)
              .resize(nb_element * nb_quadrature_points);
        } else {
          this->alloc(nb_element * nb_quadrature_points, nb_component, type,
                      ghost);
        }
      }
    }

    this->material.registerInternal(*this);
    this->is_init = true;
  }
  this->reset();

  if (this->previous_values) {
    this->previous_values->internalInitialize(nb_component);
  }
}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
void InternalFieldTmpl<Material, T>::setArrayValues(T * begin, T * end) {
  for (; begin < end; ++begin){
    *begin = this->default_value;
  }
}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
void InternalFieldTmpl<Material, T>::saveCurrentValues() {
  AKANTU_DEBUG_ASSERT(this->previous_values != nullptr,
                      "The history of the internal "
                          << this->getID() << " has not been activated");

  if (not this->is_init) {
    return;
  }

  for (auto ghost_type : ghost_types) {
    for (const auto & type : this->elementTypes(ghost_type)) {
      (*this->previous_values)(type, ghost_type)
          .copy((*this)(type, ghost_type));
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
void InternalFieldTmpl<Material, T>::restorePreviousValues() {
  AKANTU_DEBUG_ASSERT(this->previous_values != nullptr,
                      "The history of the internal "
                          << this->getID() << " has not been activated");

  if (not this->is_init) {
    return;
  }

  for (auto ghost_type : ghost_types) {
    for (const auto & type : this->elementTypes(ghost_type)) {
      (*this)(type, ghost_type)
          .copy((*this->previous_values)(type, ghost_type));
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
void InternalFieldTmpl<Material, T>::removeIntegrationPoints(
    const ElementTypeMapArray<UInt> & new_numbering) {
  for (auto ghost_type : ghost_types) {
    for (auto type : new_numbering.elementTypes(_all_dimensions, ghost_type,
                                                _ek_not_defined)) {
      if (not this->exists(type, ghost_type)) {
        continue;
      }

      Array<T> & vect = (*this)(type, ghost_type);
      if (vect.empty()) {
        continue;
      }

      const Array<UInt> & renumbering = new_numbering(type, ghost_type);

      UInt nb_quad_per_elem = fem->getNbIntegrationPoints(type, ghost_type);
      UInt nb_component = vect.getNbComponent();

      Array<T> tmp(renumbering.size() * nb_quad_per_elem, nb_component);

      AKANTU_DEBUG_ASSERT(
          tmp.size() == vect.size(),
          "Something strange append some mater was created from nowhere!!");

      AKANTU_DEBUG_ASSERT(
          tmp.size() == vect.size(),
          "Something strange append some mater was created or disappeared in "
              << vect.getID() << "(" << vect.size() << "!=" << tmp.size()
              << ") "
                 "!!");

      UInt new_size = 0;
      for (UInt i = 0; i < renumbering.size(); ++i) {
        UInt new_i = renumbering(i);
        if (new_i != UInt(-1)) {
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

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
void InternalFieldTmpl<Material, T>::printself(std::ostream & stream,
                                 int indent [[gnu::unused]]) const {
  stream << "InternalField [ " << this->getID();
#if !defined(AKANTU_NDEBUG)
  if (AKANTU_DEBUG_TEST(dblDump)) {
    stream << std::endl;
    ElementTypeMapArray<T>::printself(stream, indent + 3);
  } else {
#endif
    stream << " {" << this->getData(_not_ghost).size() << " types - "
           << this->getData(_ghost).size() << " ghost types"
           << "}";
#if !defined(AKANTU_NDEBUG)
  }
#endif
  stream << " ]";
}

/* -------------------------------------------------------------------------- */
template <>
inline void
ParameterTyped<InternalField<Real>>::setAuto(const ParserParameter & in_param) {
  Parameter::setAuto(in_param);
  Real r = in_param;
  param.setDefaultValue(r);
}

/* -------------------------------------------------------------------------- */
template <class Material, typename T>
inline InternalFieldTmpl<Material, T>::operator T() const {
  return default_value;
}

} // namespace akantu

#endif /* AKANTU_INTERNAL_FIELD_TMPL_HH_ */
