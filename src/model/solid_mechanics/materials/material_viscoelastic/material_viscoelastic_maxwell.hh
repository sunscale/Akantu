/**
 * @file   material_viscoelastic_maxwell.hh
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Tue May 08 2018
 * @date last modification: Tue May 08 2018
 *
 * @brief  Material Visco-elastic, based on Maxwell chain,
 * see
 * [] R. de Borst and A.H. van den Boogaard "Finite-element modeling of
 * deformation and cracking in early-age concrete", J.Eng.Mech., 1994
 * as well as
 * [] Manual of DIANA FEA Theory manual v.10.2 Section 37.6
 *
 * @section LICENSE
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
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_VISCOELASTIC_MAXWELL_HH__
#define __AKANTU_MATERIAL_VISCOELASTIC_MAXWELL_HH__

namespace akantu {

/**
 * Material Viscoelastic based on Maxwell chain
 *
 *
 * @verbatim

              E_0
      ------|\/\/\|-------
      |                  |
   ---|                  |---
      |                  |
      ----|\/\/\|--[|-----
      |   E_v1  \Eta_1|
   ---|                  |---
      |                  |
      ----|\/\/\|--[|-----
      |   E_v2 \Eta_2 |
   ---|                  |---
      |                  |
      ----|\/\/\|--[|----
          E_vN \Eta_N

 @endverbatim
 *
 * keyword : viscoelastic_maxwell
 *
 * parameters in the material files :
 *   - N   : number of Maxwell elements
 *   - Einf  : one spring element stiffness
 *   - Ev1 : stiffness of the 1st viscous element
 *   - Eta1: viscosity of the 1st Maxwell element
 *   ...
 *   - Ev<N> : stiffness of the Nst viscous element
 *   - Eta<N>: viscosity of the Nst Maxwell element
 */

template <UInt spatial_dimension>
class MaterialViscoelasticMaxwell : public MaterialElastic<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialViscoelasticMaxwell(SolidMechanicsModel & model, const ID & id = "");
  ~MaterialViscoelasticMaxwell() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material computed parameter
  void initMaterial() override;

  /// recompute the lame coefficient and effective tangent moduli
  void updateInternalParameters() override;

  /// update internal variable on a converged Newton
  void afterSolveStep(bool converged) override;

  /// update internal variable based on previous and current strain values
  void updateIntVariables();

  /// update the internal variable sigma_v on quadrature point
  void updateIntVarOnQuad(const Matrix<Real> & grad_u,
                          const Matrix<Real> & previous_grad_u,
                          Tensor3<Real> & sigma_v, Tensor3<Real> & epsilon_v);

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(const ElementType & el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

  /// save previous stress and strain values into "previous" arrays
  void savePreviousState() override;

  /// change flag of updateIntVar to true
  void forceUpdateVariable();

  /// change flag of updateIntVar to false
  void forceNotUpdateVariable();

  /// compute the elastic potential energy
  void computePotentialEnergy(ElementType el_type) override;

protected:
  void computePotentialEnergyOnQuad(const Matrix<Real> & grad_u, Real & epot,
                                    Tensor3<Real> & sigma_v,
                                    Tensor3<Real> & epsilon_v);

  /// update the dissipated energy, is called after the stress have been
  /// computed
  void updateDissipatedEnergy(ElementType el_type);

  void updateDissipatedEnergyOnQuad(const Matrix<Real> & grad_u,
                                    const Matrix<Real> & previous_grad_u,
                                    const Matrix<Real> & sigma,
                                    const Matrix<Real> & previous_sigma,
                                    Real & dis_energy, Real & mech_work,
                                    const Real & pot_energy);

  /// compute stresses on a quadrature point
  void computeStressOnQuad(const Matrix<Real> & grad_u,
                           const Matrix<Real> & previous_grad_u,
                           Matrix<Real> & sigma, Tensor3<Real> & sigma_v,
                           const Real & sigma_th);

  /// compute tangent moduli on a quadrature point
  void computeTangentModuliOnQuad(Matrix<Real> & tangent);

  bool hasStiffnessMatrixChanged() override {
    Real dt = this->model.getTimeStep();

    return ((this->previous_dt == dt)
                ? (!(this->previous_dt == dt)) * (this->was_stiffness_assembled)
                : (!(this->previous_dt == dt)));
    //  return (!(this->previous_dt == dt));
  }

  MatrixType getTangentType() override {
    return _symmetric;
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// give the dissipated energy
  Real getDissipatedEnergy() const;
  Real getDissipatedEnergy(ElementType type, UInt index) const;

  /// get the potential energy
  Real getPotentialEnergy() const;
  Real getPotentialEnergy(ElementType type, UInt index) const;

  /// get the potential energy
  Real getMechanicalWork() const;
  Real getMechanicalWork(ElementType type, UInt index) const;

  /// get the energy using an energy type string for the time step
  Real getEnergy(const std::string & type) override;
  Real getEnergy(const std::string & energy_id, ElementType type,
                 UInt index) override;
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  using voigt_h = VoigtHelper<spatial_dimension>;

  /// Vectors of viscosity, viscous elastic modulus, one spring element elastic
  /// modulus
  Vector<Real> Eta;
  Vector<Real> Ev;
  Real Einf;

  /// time step from previous solveStep
  Real previous_dt;

  /// Stiffness matrix template
  Matrix<Real> C;
  /// Compliance matrix template
  Matrix<Real> D;

  /// Internal variable: viscous_stress
  InternalField<Real> sigma_v;

  /// Internal variable: spring strain in Maxwell element
  InternalField<Real> epsilon_v;

  /// Dissipated energy
  InternalField<Real> dissipated_energy;

  /// Mechanical work
  InternalField<Real> mechanical_work;

  /// Update internal variable after solve step or not
  bool update_variable_flag;
};

} // namespace akantu

#endif /* __AKANTU_MATERIAL_VISCOELASTIC_MAXWELL_HH__ */
