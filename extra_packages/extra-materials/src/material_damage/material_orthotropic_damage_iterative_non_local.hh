/**
 * @file   material_orthotropic_damage_iterative_non_local.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  MaterialOrthotropicDamageIterativeNonLocal header for non-local
 * damage
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material_damage_non_local.hh"
#include "material_orthotropic_damage_iterative.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_ITERATIVE_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_ITERATIVE_NON_LOCAL_HH__

namespace akantu {

/**
 * Material Damage Iterative Non local
 *
 * parameters in the material files :
 */
template <UInt spatial_dimension>
class MaterialOrthotropicDamageIterativeNonLocal
    : public MaterialDamageNonLocal<
          spatial_dimension,
          MaterialOrthotropicDamageIterative<spatial_dimension>> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef MaterialDamageNonLocal<
      spatial_dimension, MaterialOrthotropicDamageIterative<spatial_dimension>>
      MaterialOrthotropicDamageIterativeNonLocalParent;
  MaterialOrthotropicDamageIterativeNonLocal(SolidMechanicsModel & model,
                                             const ID & id = "");

  virtual ~MaterialOrthotropicDamageIterativeNonLocal(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial();

protected:
  void computeStress(ElementType type, GhostType ghost_type);

  void computeNonLocalStress(ElementType type,
                             GhostType ghost_type = _not_ghost);

  /// associate the non-local variables of the material to their neighborhoods
  virtual void nonLocalVariableToNeighborhood();

private:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  InternalField<Real> grad_u_nl;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_orthotropic_damage_iterative_non_local_inline_impl.hh"

} // namespace akantu

#endif /* __AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_ITERATIVE_NON_LOCAL_HH__ */
