/**
 * @file   material_viscoplastic.hh
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 *
 *
 * @brief  Specialization of the material class for
 * MaterialLinearIsotropicHardening to include viscous effects (small
 * deformation)
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_voigthelper.hh"
#include "material_plastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_VISCOPLASTIC_HH_
#define AKANTU_MATERIAL_VISCOPLASTIC_HH_

namespace akantu {

/**
 * Material plastic isotropic
 *
 * parameters in the material files :
 *   - h : Hardening parameter (default: 0)
 *   - sigmay : Yield stress
 *   - rate : Rate sensitivity
 *   - edot0 : Reference strain rate
 *
 *   - ts: Time step
 */

template <UInt spatial_dimension>
class MaterialViscoPlastic : public MaterialPlastic<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialViscoPlastic(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type,
                             GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost);

protected:
  inline void
  computeStressOnQuad(const Matrix<Real> & grad_u,
                      const Matrix<Real> & previous_grad_u,
                      Matrix<Real> & sigma, const Matrix<Real> & previous_sigma,
                      Matrix<Real> & inelastic_strain,
                      const Matrix<Real> & previous_inelastic_strain,
                      Real & iso_hardening) const;

  inline void computeTangentModuliOnQuad(
      Matrix<Real> & tangent, const Matrix<Real> & grad_u,
      const Matrix<Real> & previous_grad_u, const Matrix<Real> & sigma_tensor,
      const Matrix<Real> & previous_sigma_tensor,
      const Real & iso_hardening) const;
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// Rate sensitivity component (rate)
  Real rate;

  /// Reference strain rate (edot0)
  Real edot0;

  /// Time step (ts)
  Real ts;
};
/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "material_viscoplastic_inline_impl.hh"

} // namespace akantu

#endif /* AKANTU_MATERIAL_VISCOPLASTIC_HH_ */
