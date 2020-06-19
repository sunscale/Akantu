/**
 * @file   material_igfem_saw_tooth_damage_inline_impl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief Implementation of inline functions of the squentially linear
 * saw-tooth damage model for IGFEM elements
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
inline void
MaterialIGFEMSawToothDamage<spatial_dimension>::computeDamageAndStressOnQuad(
    Matrix<Real> & sigma, Real & dam) {
  sigma *= 1 - dam;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
UInt MaterialIGFEMSawToothDamage<spatial_dimension>::updateDamage(
    UInt quad_index, const Real eq_stress, const ElementType & el_type,
    const GhostType & ghost_type) {
  AKANTU_DEBUG_ASSERT(prescribed_dam > 0.,
                      "Your prescribed damage must be greater than zero");

  Array<Real> & dam = this->damage(el_type, ghost_type);
  Real & dam_on_quad = dam(quad_index);

  /// check if damage occurs
  if (equivalent_stress(el_type, ghost_type)(quad_index) >=
      (1 - dam_tolerance) * norm_max_equivalent_stress) {
    /// damage the entire sub-element -> get the element index
    UInt el_index =
        quad_index / this->element_filter(el_type, ghost_type).getSize();
    UInt nb_quads = this->fem->getNbIntegrationPoints(el_type, ghost_type);
    UInt start_idx = el_index * nb_quads;
    Array<UInt> & sub_mat = this->sub_material(el_type, ghost_type);
    UInt damaged_quads = 0;
    if (dam_on_quad < dam_threshold) {
      for (UInt q = 0; q < nb_quads; ++q, ++start_idx) {
        if (sub_mat(start_idx)) {
          dam(start_idx) += prescribed_dam;
          damaged_quads += 1;
        }
      }
    } else {
      for (UInt q = 0; q < nb_quads; ++q, ++start_idx) {
        if (sub_mat(start_idx)) {
          dam(start_idx) += max_damage;
          damaged_quads += 1;
        }
      }
    }
    return damaged_quads;
  }

  return 0;
}

/* -------------------------------------------------------------------------- */
