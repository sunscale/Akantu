/**
 * @file   dumper_elemental_field.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Dec 03 2017
 *
 * @brief  description of elemental fields
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_DUMPER_ELEMENTAL_FIELD_HH_
#define AKANTU_DUMPER_ELEMENTAL_FIELD_HH_
/* -------------------------------------------------------------------------- */
#include "communicator.hh"
#include "dumper_field.hh"
#include "dumper_generic_elemental_field.hh"
#ifdef AKANTU_IGFEM
#include "dumper_igfem_elemental_field.hh"
#endif
/* -------------------------------------------------------------------------- */
namespace akantu {
namespace dumpers {
/* -------------------------------------------------------------------------- */

template <typename T, template <class> class ret = Vector,
          bool filtered = false>
class ElementalField
    : public GenericElementalField<SingleType<T, ret, filtered>,
                                   elemental_field_iterator> {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using types = SingleType<T, ret, filtered>;
  using field_type = typename types::field_type;
  using iterator = elemental_field_iterator<types>;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ElementalField(const field_type & field,
                 UInt spatial_dimension = _all_dimensions,
                 GhostType ghost_type = _not_ghost,
                 ElementKind element_kind = _ek_not_defined)
      : GenericElementalField<types, elemental_field_iterator>(
            field, spatial_dimension, ghost_type, element_kind) {}
};

/* -------------------------------------------------------------------------- */

} // namespace dumpers
} // namespace akantu

#endif /* AKANTU_DUMPER_ELEMENTAL_FIELD_HH_ */
