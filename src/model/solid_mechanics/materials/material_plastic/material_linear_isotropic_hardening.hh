/**
 * @file   material_linear_isotropic_hardening.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 *
 * @date creation: Thu Oct 03 2013
 * @date last modification: Mon Apr 07 2014
 *
 * @brief  Specialization of the material class for isotropic finite deformation linear hardening plasticity
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_voigthelper.hh"
#include "material_plastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_LINEAR_ISOTROPIC_HARDENING_HH__
#define __AKANTU_MATERIAL_LINEAR_ISOTROPIC_HARDENING_HH__

__BEGIN_AKANTU__

/**
 * Material plastic with a linear evolution of the yielding stress
 */
template <UInt spatial_dimension>
class MaterialLinearIsotropicHardening : public MaterialPlastic<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialLinearIsotropicHardening(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(const ElementType & el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost);

protected:
  inline void computeStressOnQuad(const Matrix<Real> & grad_u,
                                  const Matrix<Real> & previous_grad_u,
                                  Matrix<Real> & sigma,
                                  const Matrix<Real> & previous_sigma,
                                  Matrix<Real> & inelas_strain,
                                  const Matrix<Real> & previous_inelas_strain,
                                  Real & iso_hardening,
                                  const Real & previous_iso_hardening,
                                  const Real & sigma_th,
                                  const Real & previous_sigma_th);

  inline void computeTangentModuliOnQuad(Matrix<Real> & tangent,
                                         const Matrix<Real> & grad_u,
                                         const Matrix<Real> & previous_grad_u,
                                         const Matrix<Real> & sigma_tensor,
                                         const Matrix<Real> & previous_sigma_tensor,
                                         const Real & iso_hardening) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

};
/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "material_linear_isotropic_hardening_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_LINEAR_ISOTROPIC_HARDENING_HH__ */
