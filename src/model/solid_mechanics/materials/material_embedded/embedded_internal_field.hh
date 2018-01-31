/**
 * @file   embedded_internal_field.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Jun 30 2015
 *
 * @brief  Embedded Material internal properties
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "aka_common.hh"
#include "internal_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_EMBEDDED_INTERNAL_FIELD_HH__
#define __AKANTU_EMBEDDED_INTERNAL_FIELD_HH__

namespace akantu {

class Material;
class FEEngine;

/// This class is used for MaterialReinforcement internal fields
template<typename T>
class EmbeddedInternalField : public InternalField<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// Constructor
  EmbeddedInternalField(const ID & id, Material & material):
    InternalField<T>(id,
                     material,
                     material.getModel().getFEEngine("EmbeddedInterfaceFEEngine"),
                     material.getElementFilter()) {
    this->spatial_dimension = 1;
  }

  /// Copy constructor
  EmbeddedInternalField(const ID & id, const EmbeddedInternalField & other):
    InternalField<T>(id, other) {
    this->spatial_dimension = 1;
  }

  void operator=(const EmbeddedInternalField & other) {
    InternalField<T>::operator=(other);
    this->spatial_dimension = 1;
  }
};

/// Method used to initialise the embedded internal fields from material file
template <>
inline void ParameterTyped<EmbeddedInternalField<Real>>::setAuto(
    const ParserParameter & in_param) {
  Parameter::setAuto(in_param);
  Real r = in_param;
  param.setDefaultValue(r);
}

} // akantu

#endif // __AKANTU_EMBEDDED_INTERNAL_FIELD_HH__
