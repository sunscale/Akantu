/**
 * @file   dumper_igfem_elemental_field.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  description of IGFEM elemental fields
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#ifndef AKANTU_DUMPER_IGFEM_ELEMENTAL_FIELD_HH_
#define AKANTU_DUMPER_IGFEM_ELEMENTAL_FIELD_HH_
/* -------------------------------------------------------------------------- */
#include "dumper_field.hh"
#include "dumper_igfem_generic_elemental_field.hh"
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
namespace dumpers {
/* -------------------------------------------------------------------------- */

template <typename T, template <class> class ret = Vector,
          bool filtered = false>
class IGFEMElementalField
    : public IGFEMGenericElementalField<SingleType<T, ret, filtered>,
                                        igfem_elemental_field_iterator> {

public:
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

  typedef SingleType<T, ret, filtered> types;
  typedef typename types::field_type field_type;
  typedef elemental_field_iterator<types> iterator;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  IGFEMElementalField(const field_type & field,
                      UInt spatial_dimension = _all_dimensions,
                      GhostType ghost_type = _not_ghost,
                      ElementKind element_kind = _ek_igfem)
      : IGFEMGenericElementalField<types, igfem_elemental_field_iterator>(
            field, spatial_dimension, ghost_type, element_kind) {}
};

} // namespace dumpers
} // namespace akantu

#endif /* AKANTU_DUMPER_IGFEM_ELEMENTAL_FIELD_HH_ */
