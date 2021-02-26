/**
 * @file   material_vreepeerlings.hh
 *
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 *
 * @brief  Specialization of the material class for the VreePeerlings material
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material_damage.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_VREEPEERLINGS_HH_
#define AKANTU_MATERIAL_VREEPEERLINGS_HH_

namespace akantu {

/**
 * Material vreepeerlings
 *
 * parameters in the material files :
 *   - Kapaoi  : (default: 0.0001) Initial threshold (of the equivalent strain)
 * >= Crate
 *   - Kapac  : (default: 0.0002) Final threshold (of the equivalent strain)
 *   - Arate  : (default: 1.) Fitting parameter (must be close to 1 to do tend
 * to 0 the stress in the damaged element)
 *   - Brate   : (default: 1.) This parameter determines the rate at which the
 * damage grows
 *   - Crate   : (default: 0.0001) This parameter determines the rate at which
 * the damage grows
 *   - Kct    : (default: 1.) Ratio between compressive and tensile strength
 */
template <UInt spatial_dimension,
          template <UInt> class MatParent = MaterialElastic>
class MaterialVreePeerlings
    : public MaterialDamage<spatial_dimension, MatParent> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef MaterialDamage<spatial_dimension, MatParent>
      MaterialVreePeerlingsParent;
  MaterialVreePeerlings(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialVreePeerlings(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial();

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(Matrix<Real> & F, Matrix<Real> & sigma,
                                  Real & dam, Real & Equistrain_rate,
                                  Real & Equistrain, Real & Kapaq, Real dt,
                                  Matrix<Real> & strain_rate_vrpgls,
                                  Real & FullDam_ValStrain,
                                  Real & FullDam_ValStrain_rate,
                                  Real & Nb_damage);

  inline void computeDamageAndStressOnQuad(Matrix<Real> & sigma, Real & dam,
                                           Real & Equistrain_rate,
                                           Real & Equistrain, Real & Kapaq,
                                           Real dt, Real & FullDam_Valstrain,
                                           Real & FullDam_Valstrain_rate,
                                           Real & Nb_damage);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Initial threshold (of the equivalent strain) (used in the initial step)
  Real Kapaoi;

  /// Final threshold (of the equivalent strain) (used in the initial step)
  Real Kapac;

  /// This parameter determines the rate at which the damage grows
  Real Arate;

  /// This parameter determines the rate at which the damage grows
  Real Brate;

  /// This parameter determines the rate at which the damage grows
  Real Crate;

  /// Ratio between compressive and tensile strength
  Real Kct;

  /// Randomness on Kapaoi
  Real Kapao_randomness;

  /// Kapa vector which contains the initial damage threshold
  RandomInternalField<Real> Kapa;

  /// Strain rate tensor to compute the rate dependent damage law
  InternalField<Real> strain_rate_vreepeerlings;

  /// Value of the equivalent strain when damage = 1
  InternalField<Real> Full_dam_value_strain;

  /// Value of the equivalent strain rate when damage = 1
  InternalField<Real> Full_dam_value_strain_rate;

  /// Count the number of times that the material is damaged to damage = 0 until
  /// damage = 1
  InternalField<Real> Number_damage;

  /// Equivalent strain used to compute the damage evolution
  InternalField<Real> equi_strain;

  /// Equivalent strain rate used to compute the damage evolution
  InternalField<Real> equi_strain_rate;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_vreepeerlings_inline_impl.hh"
#include "material_vreepeerlings_tmpl.hh"

} // namespace akantu

#endif /* AKANTU_MATERIAL_VREEPEERLINGS_HH_ */
