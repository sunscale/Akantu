/**
 * @file   material_plastic.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Thu Dec 07 2017
 *
 * @brief  Common interface for plastic materials
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
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_PLASTIC_HH_
#define AKANTU_MATERIAL_PLASTIC_HH_

namespace akantu {

/**
 * Parent class for the plastic constitutive laws
 * parameters in the material files :
 *   - h : Hardening parameter (default: 0)
 *   - sigmay : Yield stress
 */
template <UInt dim> class MaterialPlastic : public MaterialElastic<dim> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialPlastic(SolidMechanicsModel & model, const ID & id = "");
  MaterialPlastic(SolidMechanicsModel & model, UInt a_dim, const Mesh & mesh,
                  FEEngine & fe_engine, const ID & id = "");

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /**
   * @brief Return potential or plastic energy
   *
   * Plastic dissipated energy is integrated over time.
   */
  Real getEnergy(const std::string & type) override;

  /// Update the plastic energy for the current timestep
  void updateEnergies(ElementType el_type) override;

  /// Compute the true potential energy (w/ elastic strain)
  void computePotentialEnergy(ElementType el_type) override;

protected:
  /// compute the stress and inelastic strain for the quadrature point
  inline void computeStressAndInelasticStrainOnQuad(
      const Matrix<Real> & grad_u, const Matrix<Real> & previous_grad_u,
      Matrix<Real> & sigma, const Matrix<Real> & previous_sigma,
      Matrix<Real> & inelastic_strain,
      const Matrix<Real> & previous_inelastic_strain,
      const Matrix<Real> & delta_inelastic_strain) const;

  inline void computeStressAndInelasticStrainOnQuad(
      const Matrix<Real> & delta_grad_u, Matrix<Real> & sigma,
      const Matrix<Real> & previous_sigma, Matrix<Real> & inelastic_strain,
      const Matrix<Real> & previous_inelastic_strain,
      const Matrix<Real> & delta_inelastic_strain) const;

  /// Get the integrated plastic energy for the time step
  Real getPlasticEnergy();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Yield stresss
  Real sigma_y;

  /// hardening modulus
  Real h;

  /// isotropic hardening, r
  InternalField<Real> iso_hardening;

  /// inelastic strain arrays ordered by element types (inelastic deformation)
  InternalField<Real> inelastic_strain;

  /// Plastic energy
  InternalField<Real> plastic_energy;

  /// @todo : add a coefficient beta that will multiply the plastic energy
  /// increment
  /// to compute the energy converted to heat

  /// Plastic energy increment
  InternalField<Real> d_plastic_energy;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

} // namespace akantu

#include "material_plastic_inline_impl.hh"
#endif /* AKANTU_MATERIAL_PLASTIC_HH_ */
