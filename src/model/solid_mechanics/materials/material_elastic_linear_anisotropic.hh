/**
 * @file   material_elastic_linear_anisotropic.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Enrico Milanese <enrico.milanese@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Feb 16 2018
 *
 * @brief  Orthotropic elastic material
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
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_HH_
#define AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_HH_

namespace akantu {

/**
 * General linear anisotropic elastic material
 * The only constraint on the elastic tensor is that it can be represented
 * as a symmetric 6x6 matrix (3D) or 3x3 (2D).
 *
 * parameters in the material files :
 *   - rho  : density (default: 0)
 *   - C_ij  : entry on the stiffness
 */
template <UInt Dim> class MaterialElasticLinearAnisotropic : public Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialElasticLinearAnisotropic(SolidMechanicsModel & model,
                                   const ID & id = "", bool symmetric = true);

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

  void updateInternalParameters() override;

  bool hasStiffnessMatrixChanged() override {
    return (not was_stiffness_assembled);
  }

  MatrixType getTangentType() override {
    return _symmetric;
  }

protected:
  // compute C from Cprime
  void rotateCprime();

  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(const Matrix<Real> & grad_u,
                                  Matrix<Real> & sigma) const;

  /// tangent matrix for a given quadrature point
  inline void computeTangentModuliOnQuad(Matrix<Real> & tangent) const;

  inline void computePotentialEnergyOnQuad(const Matrix<Real> & grad_u,
                                           const Matrix<Real> & sigma,
                                           Real & epot);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// compute max wave celerity
  Real getCelerity(const Element & element) const override;

  AKANTU_GET_MACRO(VoigtStiffness, C, Matrix<Real>);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  using voigt_h = VoigtHelper<Dim>;

  /// direction matrix and vectors
  std::vector<std::unique_ptr<Vector<Real>>> dir_vecs;

  Matrix<Real> rot_mat;
  /// Elastic stiffness tensor in material frame and full vectorised notation
  Matrix<Real> Cprime;
  /// Elastic stiffness tensor in voigt notation
  Matrix<Real> C;
  /// eigenvalues of stiffness tensor
  Vector<Real> eigC;

  bool symmetric;

  /// viscous proportion
  Real alpha;

  /// defines if the stiffness was computed
  bool was_stiffness_assembled;
};
} // namespace akantu

#include "material_elastic_linear_anisotropic_inline_impl.hh"

#endif /* AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_HH_ */
