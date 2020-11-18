/**
 * @file   material_elastic.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Nov 17 2017
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
#include "material_thermal.hh"
#include "plane_stress_toolbox.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_ELASTIC_HH_
#define AKANTU_MATERIAL_ELASTIC_HH_

namespace akantu {

/**
 * Material elastic isotropic
 *
 * parameters in the material files :
 *   - E   : Young's modulus (default: 0)
 *   - nu  : Poisson's ratio (default: 1/2)
 *   - Plane_Stress : if 0: plane strain, else: plane stress (default: 0)
 */
template <UInt spatial_dimension>
class MaterialElastic
    : public PlaneStressToolbox<spatial_dimension,
                                MaterialThermal<spatial_dimension>> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:
  using Parent =
      PlaneStressToolbox<spatial_dimension, MaterialThermal<spatial_dimension>>;

public:
  MaterialElastic(SolidMechanicsModel & model, const ID & id = "");
  MaterialElastic(SolidMechanicsModel & model, UInt dim, const Mesh & mesh,
                  FEEngine & fe_engine, const ID & id = "");

  ~MaterialElastic() override = default;

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial() override;

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

  /// compute the elastic potential energy
  void computePotentialEnergy(ElementType el_type) override;

  void
  computePotentialEnergyByElement(ElementType type, UInt index,
                                  Vector<Real> & epot_on_quad_points) override;

  /// compute the p-wave speed in the material
  Real getPushWaveSpeed(const Element & element) const override;

  /// compute the s-wave speed in the material
  Real getShearWaveSpeed(const Element & element) const override;

protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(const Matrix<Real> & grad_u,
                                  Matrix<Real> & sigma,
                                  Real sigma_th = 0) const;

  /// compute the tangent stiffness matrix for an element
  inline void computeTangentModuliOnQuad(Matrix<Real> & tangent) const;

  /// recompute the lame coefficient if E or nu changes
  void updateInternalParameters() override;

  static inline void computePotentialEnergyOnQuad(const Matrix<Real> & grad_u,
                                                  const Matrix<Real> & sigma,
                                                  Real & epot);

  bool hasStiffnessMatrixChanged() override {
    return (not was_stiffness_assembled);
  }

  MatrixType getTangentType() override {
    return _symmetric;
  }
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get first Lame constant
  AKANTU_GET_MACRO(Lambda, lambda, Real);

  /// get second Lame constant
  AKANTU_GET_MACRO(Mu, mu, Real);

  /// get bulk modulus
  AKANTU_GET_MACRO(Kappa, kpa, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// First Lamé coefficient
  Real lambda;

  /// Second Lamé coefficient (shear modulus)
  Real mu;

  /// Bulk modulus
  Real kpa;

  /// defines if the stiffness was computed
  bool was_stiffness_assembled;
};

} // namespace akantu

#include "material_elastic_inline_impl.hh"

#endif /* AKANTU_MATERIAL_ELASTIC_HH_ */
