/**
 * @file   material_orthtropic_damage_iterative.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date   Sun Mar  8 12:54:30 2015
 *
 * @brief Specialization of the class material orthotropic damage to
 * damage only one gauss point at a time and propagate damage in a
 * linear way. Max principal stress criterion is used as a failure
 * criterion.
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
#include "material_orthotropic_damage.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_ITERATIVE_HH__
#define __AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_ITERATIVE_HH__

namespace akantu {

/**
 * Material damage iterative
 *
 * parameters in the material files :
 *   - Sc
 */
template <UInt spatial_dimension>
class MaterialOrthotropicDamageIterative
    : public MaterialOrthotropicDamage<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialOrthotropicDamageIterative(SolidMechanicsModel & model,
                                     const ID & id = "");

  virtual ~MaterialOrthotropicDamageIterative(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  ///  virtual void updateInternalParameters();

  virtual void computeAllStresses(GhostType ghost_type = _not_ghost);

  /// update internal field damage
  UInt updateDamage();

  /// update energies after damage has been updated
  virtual void updateEnergiesAfterDamage(ElementType el_type);

protected:
  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type,
                             GhostType ghost_type = _not_ghost);

  /// compute the equivalent stress on each Gauss point (i.e. the max prinicpal
  /// stress) and normalize it by the tensile strength
  void computeNormalizedEquivalentStress(const Array<Real> & grad_u,
                                         ElementType el_type,
                                         GhostType ghost_type = _not_ghost);

  /// find max normalized equivalent stress
  void findMaxNormalizedEquivalentStress(ElementType el_type,
                                         GhostType ghost_type = _not_ghost);

  inline void computeDamageAndStressOnQuad(Matrix<Real> & sigma,
                                           Matrix<Real> & one_minus_D,
                                           Matrix<Real> & root_one_minus_D,
                                           Matrix<Real> & damage,
                                           Matrix<Real> & first_term,
                                           Matrix<Real> & third_term);
  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get max normalized equivalent stress
  AKANTU_GET_MACRO(NormMaxEquivalentStress, norm_max_equivalent_stress, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// resistance to damage
  RandomInternalField<Real> Sc;

  /// internal field to store equivalent stress on each Gauss point
  InternalField<Real> equivalent_stress;

  /// internal field to store the direction of the current damage frame
  InternalField<Real> stress_dir;

  /// damage increment
  Real prescribed_dam;

  /// maximum equivalent stress
  Real norm_max_equivalent_stress;

  /// define damage threshold at which damage will be set to 1
  Real dam_threshold;

  /// quadrature point with highest equivalent Stress
  IntegrationPoint q_max;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_orthotropic_damage_iterative_inline_impl.hh"

} // namespace akantu

#endif /* __AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_ITERATIVE_HH__ */
