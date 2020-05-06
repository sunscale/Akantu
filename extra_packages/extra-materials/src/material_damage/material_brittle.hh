/**
 * @file   material_brittle.hh
 *
 * @author Aranda Ruiz Josue <josue.arandaruiz@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 *
 * @brief  Brittle damage law
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
#include "material_damage.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_BRITTLE_HH__
#define __AKANTU_MATERIAL_BRITTLE_HH__

namespace akantu {

/**
 * Material brittle
 *
 * parameters in the material files :
 *   - S_0      : Critical stress at low strain rate (default: 157e6)
 *   - E_0      : Low strain rate threshold (default: 27e3)
 *   - A,B,C,D  : Fitting parameters for the critical stress at high strain
 * rates
 *                (default: 1.622e-11, -1.3274e-6, 3.6544e-2, -181.38)
 */
template <UInt spatial_dimension>
class MaterialBrittle : public MaterialDamage<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialBrittle(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial();

  virtual void updateInternalParameters();

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(Matrix<Real> & grad_u, Matrix<Real> & grad_v,
                                  Matrix<Real> & sigma, Real & dam,
                                  Real & sigma_equivalent,
                                  Real & fracture_stress);

  inline void computeDamageAndStressOnQuad(Matrix<Real> & sigma, Real & dam,
                                           Real & sigma_c,
                                           Real & fracture_stress);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// strain rate arrays ordered by element types
  InternalField<Real> strain_rate_brittle;

  // polynome constants for critical stress value
  Real A;
  Real B;
  Real C;
  Real D;

  // minimum strain rate
  Real E_0;

  // Critical stress at low strain rates
  Real S_0;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_brittle_inline_impl.hh"

} // namespace akantu

#endif /* __AKANTU_MATERIAL_brittle_HH__ */
