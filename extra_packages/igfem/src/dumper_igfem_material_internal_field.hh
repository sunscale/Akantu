/**
 * @file   dumper_igfem_material_internal_field.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  description of  IGFEM material internal field
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#ifndef AKANTU_DUMPER_IGFEM_MATERIAL_INTERNAL_FIELD_HH_
#define AKANTU_DUMPER_IGFEM_MATERIAL_INTERNAL_FIELD_HH_
/* -------------------------------------------------------------------------- */
#include "dumper_igfem_quadrature_points_field.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
namespace dumpers {
/* -------------------------------------------------------------------------- */

template <typename T, bool filtered = false>
class IGFEMInternalMaterialField
    : public IGFEMGenericElementalField<SingleType<T, Vector, filtered>,
                                        igfem_quadrature_point_iterator> {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

public:
  typedef SingleType<T, Vector, filtered> types;
  typedef IGFEMGenericElementalField<types, igfem_quadrature_point_iterator>
      parent;
  typedef typename types::field_type field_type;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  IGFEMInternalMaterialField(const field_type & field,
                             UInt spatial_dimension = _all_dimensions,
                             GhostType ghost_type = _not_ghost,
                             ElementKind kind = _ek_igfem)
      : parent(field, spatial_dimension, ghost_type, kind) {}
};

} // namespace dumpers
} // namespace akantu

#endif /* AKANTU_DUMPER_IGFEM_MATERIAL_INTERNAL_FIELD_HH_ */
