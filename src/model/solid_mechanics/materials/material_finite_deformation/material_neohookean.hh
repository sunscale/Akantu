/**
 * @file   material_neohookean.hh
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 29 2017
 *
 * @brief  Material isotropic elastic
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
#include "material.hh"
#include "plane_stress_toolbox.hh"

/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_NEOHOOKEAN_HH_
#define AKANTU_MATERIAL_NEOHOOKEAN_HH_

namespace akantu {

/**
 * Material elastic isotropic
 *
 * parameters in the material files :
 *   - rho : density (default: 0)
 *   - E   : Young's modulus (default: 0)
 *   - nu  : Poisson's ratio (default: 1/2)
 *   - Plane_Stress : if 0: plane strain, else: plane stress (default: 0)
 */
template <UInt spatial_dimension>
class MaterialNeohookean : public PlaneStressToolbox<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialNeohookean(SolidMechanicsModel & model, const ID & id = "");

  ~MaterialNeohookean() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material computed parameter
  void initMaterial() override;

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  /// Computation of the cauchy stress for plane strain materials
  void
  computeCauchyStressPlaneStress(ElementType el_type,
                                 GhostType ghost_type = _not_ghost) override;

  /// Non linear computation of the third direction strain in 2D plane stress
  /// case
  void computeThirdAxisDeformation(ElementType el_type,
                                   GhostType ghost_type = _not_ghost) override;

  /// compute the elastic potential energy
  void computePotentialEnergy(ElementType el_type) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

  /// compute the p-wave speed in the material
  Real getPushWaveSpeed(const Element & element) const override;

  /// compute the s-wave speed in the material
  Real getShearWaveSpeed(const Element & element) const override;

  MatrixType getTangentType() override {
    return _symmetric;
  }

protected:
  /// constitutive law for a given quadrature point
  inline void computePiolaKirchhoffOnQuad(const Matrix<Real> & E,
                                          Matrix<Real> & S);

  /// constitutive law for a given quadrature point (first piola)
  inline void computeFirstPiolaKirchhoffOnQuad(const Matrix<Real> & grad_u,
                                               const Matrix<Real> & S,
                                               Matrix<Real> & P);

  /// constitutive law for a given quadrature point
  inline void computeDeltaStressOnQuad(const Matrix<Real> & grad_u,
                                       const Matrix<Real> & grad_delta_u,
                                       Matrix<Real> & delta_S);

  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(Matrix<Real> & grad_u, Matrix<Real> & S,
                                  const Real & C33 = 1.0);

  /// constitutive law for a given quadrature point
  inline void computeThirdAxisDeformationOnQuad(Matrix<Real> & grad_u,
                                                Real & c33_value);

  /// constitutive law for a given quadrature point
  // inline void updateStressOnQuad(const Matrix<Real> & sigma,
  //                              Matrix<Real> & cauchy_sigma);

  /// compute the potential energy for a quadrature point
  inline void computePotentialEnergyOnQuad(const Matrix<Real> & grad_u,
                                           Real & epot);

  /// compute the tangent stiffness matrix for an element
  void computeTangentModuliOnQuad(Matrix<Real> & tangent, Matrix<Real> & grad_u,
                                  const Real & C33 = 1.0);

  /// recompute the lame coefficient if E or nu changes
  void updateInternalParameters() override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the young modulus
  Real E;

  /// Poisson coefficient
  Real nu;

  /// First Lamé coefficient
  Real lambda;

  /// Second Lamé coefficient (shear modulus)
  Real mu;

  /// Bulk modulus
  Real kpa;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "material_neohookean_inline_impl.hh"

#endif /* AKANTU_MATERIAL_NEOHOOKEAN_HH_ */
