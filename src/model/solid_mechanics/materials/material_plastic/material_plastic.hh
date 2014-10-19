/**
 * @file   material_plastic.hh
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date creation: Mon Apr 07 2014
 * @date last modification: Mon Apr 07 2014
 *
 * @brief  Common interface for plastic materials
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
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */


#ifndef __AKANTU_MATERIAL_PLASTIC_HH__
#define __AKANTU_MATERIAL_PLASTIC_HH__

__BEGIN_AKANTU__

/**
 * Parent class for the plastic constitutive laws
 * parameters in the material files :
 *   - h : Hardening parameter (default: 0)
 *   - sigmay : Yield stress
 */
template<UInt dim>
class MaterialPlastic : public MaterialElastic<dim> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialPlastic(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual Real getEnergy(std::string type);

  /// Compute the plastic energy
  virtual void updateEnergies(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// Compute the true potential energy
  virtual void computePotentialEnergy(ElementType el_type, GhostType ghost_type);

protected:
  inline void computeStressAndInelasticStrainOnQuad(const Matrix<Real> & grad_u,
                                                    const Matrix<Real> & previous_grad_u,
                                                    Matrix<Real> & sigma,
                                                    const Matrix<Real> & previous_sigma,
                                                    Matrix<Real> & inelas_strain,
                                                    const Matrix<Real> & previous_inelas_strain,
                                                    const Matrix<Real> & delta_inelastic_strain) const;

  inline void computeStressAndInelasticStrainOnQuad(const Matrix<Real> & delta_grad_u,
                                                    Matrix<Real> & sigma,
                                                    const Matrix<Real> & previous_sigma,
                                                    Matrix<Real> & inelas_strain,
                                                    const Matrix<Real> & previous_inelas_strain,
                                                    const Matrix<Real> & delta_inelastic_strain) const;


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

  /// @todo : add a coefficient beta that will multiply the plastic energy increment
  /// to compute the energy converted to heat

  /// Plastic energy increment
  InternalField<Real> d_plastic_energy;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_plastic_inline_impl.cc"


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_PLASTIC_HH__ */
