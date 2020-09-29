/**
 * @file   material_damage_linear.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 *
 *
 * @brief  Material isotropic elastic + linear softening
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

#ifndef AKANTU_MATERIAL_DAMAGE_LINEAR_HH_
#define AKANTU_MATERIAL_DAMAGE_LINEAR_HH_

namespace akantu {

/**
 * Material liner damage
 *
 * parameters in the material files :
 *   - Sigc : (default: 1e5)
 *   - Gc  : (default: 2)
 */
template <UInt spatial_dimension>
class MaterialDamageLinear : public MaterialDamage<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialDamageLinear(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialDamageLinear(){};

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
                                  Real & damage, Real & K);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// kind of toughness
  Real Gc;

  /// critical stress
  Real Sigc;

  /// damage internal variable
  InternalField<Real> K;

  Real Epsmin, Epsmax;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_damage_linear_inline_impl.hh"

} // namespace akantu

#endif /* AKANTU_MATERIAL_DAMAGE_LINEAR_HH_ */
