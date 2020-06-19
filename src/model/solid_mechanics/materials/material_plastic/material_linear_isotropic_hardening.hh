/**
 * @file   material_linear_isotropic_hardening.hh
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Benjamin Paccaud <benjamin.paccaud@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sat Dec 02 2017
 *
 * @brief  Specialization of the material class for isotropic finite deformation
 * linear hardening plasticity
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
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

namespace akantu {

/**
 * Material plastic with a linear evolution of the yielding stress
 */
template <UInt spatial_dimension>
class MaterialLinearIsotropicHardening
    : public MaterialPlastic<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialLinearIsotropicHardening(SolidMechanicsModel & model,
                                   const ID & id = "");
  MaterialLinearIsotropicHardening(SolidMechanicsModel & model, UInt dim,
                                   const Mesh & mesh, FEEngine & fe_engine,
                                   const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(const ElementType & el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

protected:
  /// Infinitesimal deformations
  inline void computeStressOnQuad(
      const Matrix<Real> & grad_u, const Matrix<Real> & previous_grad_u,
      Matrix<Real> & sigma, const Matrix<Real> & previous_sigma,
      Matrix<Real> & inelas_strain, const Matrix<Real> & previous_inelas_strain,
      Real & iso_hardening, const Real & previous_iso_hardening,
      const Real & sigma_th, const Real & previous_sigma_th);
  /// Finite deformations
  inline void computeStressOnQuad(
      const Matrix<Real> & grad_u, const Matrix<Real> & previous_grad_u,
      Matrix<Real> & sigma, const Matrix<Real> & previous_sigma,
      Matrix<Real> & inelas_strain, const Matrix<Real> & previous_inelas_strain,
      Real & iso_hardening, const Real & previous_iso_hardening,
      const Real & sigma_th, const Real & previous_sigma_th,
      const Matrix<Real> & F_tensor);

  inline void computeTangentModuliOnQuad(
      Matrix<Real> & tangent, const Matrix<Real> & grad_u,
      const Matrix<Real> & previous_grad_u, const Matrix<Real> & sigma_tensor,
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

} // namespace akantu

#include "material_linear_isotropic_hardening_inline_impl.hh"

#endif /* __AKANTU_MATERIAL_LINEAR_ISOTROPIC_HARDENING_HH__ */
