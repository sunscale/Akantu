/**
 * @file   material_brittle_non_local.hh
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 *
 * @brief  MaterialBrittleNonLocal header for non-local damage
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material_brittle.hh"
#include "material_damage_non_local.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_BRITTLE_NON_LOCAL_HH_
#define AKANTU_MATERIAL_BRITTLE_NON_LOCAL_HH_

namespace akantu {

/**
 * Material Brittle Non local
 *
 * parameters in the material files :
 */
template <UInt spatial_dimension>
class MaterialBrittleNonLocal
    : public MaterialDamageNonLocal<spatial_dimension,
                                    MaterialBrittle<spatial_dimension>> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef MaterialDamageNonLocal<spatial_dimension,
                                 MaterialBrittle<spatial_dimension>>
      MaterialBrittleNonLocalParent;
  MaterialBrittleNonLocal(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialBrittleNonLocal(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial();

protected:
  /// constitutive law
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  void computeNonLocalStress(ElementType type,
                             GhostType ghost_type = _not_ghost);

  /// associate the non-local variables of the material to their neighborhoods
  virtual void nonLocalVariableToNeighborhood();

private:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Sigma_max, Sigma_max, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  InternalField<Real> Sigma_max;
  InternalField<Real> Sigma_maxnl;
  InternalField<Real> Sigma_fracture;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_brittle_non_local_inline_impl.hh"

} // namespace akantu

#endif /* AKANTU_MATERIAL_BRITTLE_NON_LOCAL_HH_ */
