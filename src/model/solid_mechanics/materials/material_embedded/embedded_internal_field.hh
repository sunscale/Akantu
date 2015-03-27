/**
 * @file   embedded_internal_field.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Mar 19 2015
 * @date last modification: Thu Mar 19 2015
 *
 * @brief  Embedded Material internal properties
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "internal_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_EMBEDDED_INTERNAL_FIELD_HH__
#define __AKANTU_EMBEDDED_INTERNAL_FIELD_HH__

__BEGIN_AKANTU__

class Material;
class FEEngine;

template<typename T>
class EmbeddedInternalField : public InternalField<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  EmbeddedInternalField(const ID & id, Material & material):
    InternalField<T>(id,
                     material,
                     material.getModel().getFEEngine("EmbeddedInterfaceFEEngine"),
                     material.getElementFilter())
  {
    this->spatial_dimension = 1;
  }

  EmbeddedInternalField(const ID & id, const EmbeddedInternalField & other):
    InternalField<T>(id, other)
  {
    this->spatial_dimension = 1;
  }

  void operator=(const EmbeddedInternalField & other) {
    InternalField<T>::operator=(other);
    this->spatial_dimension = 1;
  }

};

__END_AKANTU__

#endif // __AKANTU_EMBEDDED_INTERNAL_FIELD_HH__
