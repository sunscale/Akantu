/**
 * @file   material_vreepeerlings_non_local.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 *
 * @brief  MaterialVreePeerlings header for non-local damage
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material_damage_non_local.hh"
#include "material_vreepeerlings.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_VREEPEERLINGS_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_VREEPEERLINGS_NON_LOCAL_HH__

namespace akantu {

/**
 * Material VreePeerlings Non local
 *
 * parameters in the material files :
 */
template <UInt spatial_dimension,
          template <UInt> class MatParent = MaterialElastic>
class MaterialVreePeerlingsNonLocal
    : public MaterialDamageNonLocal<
          spatial_dimension,
          MaterialVreePeerlings<spatial_dimension, MatParent>> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef MaterialVreePeerlings<spatial_dimension, MatParent> Parent;
  typedef MaterialDamageNonLocal<spatial_dimension, Parent>
      MaterialVreePeerlingsNonLocalParent;

  MaterialVreePeerlingsNonLocal(SolidMechanicsModel & model,
                                const ID & id = "");

  virtual ~MaterialVreePeerlingsNonLocal(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial();

  /// constitutive law for all element of a type
  // void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// constitutive law
  virtual void computeNonLocalStress(ElementType el_type,
                                     GhostType ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// non local version of equivalent strain
  InternalField<Real> equi_strain_non_local;

  /// non local version of equivalent strain rate
  InternalField<Real> equi_strain_rate_non_local;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_vreepeerlings_non_local_inline_impl.hh"

} // namespace akantu

#endif /* __AKANTU_MATERIAL_VREEPEERLINGS_NON_LOCAL_HH__ */
