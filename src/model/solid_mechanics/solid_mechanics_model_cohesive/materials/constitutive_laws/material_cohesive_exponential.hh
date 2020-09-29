/**
 * @file   material_cohesive_exponential.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Dec 03 2017
 *
 * @brief  Exponential irreversible cohesive law of mixed mode loading
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
#include "material_cohesive.hh"

/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_COHESIVE_EXPONENTIAL_HH_
#define AKANTU_MATERIAL_COHESIVE_EXPONENTIAL_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * Cohesive material Exponential damage
 *
 * parameters in the material files :
 *   - sigma_c   : critical stress sigma_c  (default: 0)
 *   - beta      : weighting parameter for sliding and normal opening (default:
 * 0)
 *   - delta_c   : critical opening (default: 0)
 */
template <UInt spatial_dimension>
class MaterialCohesiveExponential : public MaterialCohesive {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialCohesiveExponential(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// Initialization
  void initMaterial() override;

  /// constitutive law
  void computeTraction(const Array<Real> & normal, ElementType el_type,
                       GhostType ghost_type = _not_ghost) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentTraction(ElementType el_type,
                              Array<Real> & tangent_matrix,
                              const Array<Real> & normal,
                              GhostType ghost_type = _not_ghost) override;

private:
  void computeCoupledTraction(Vector<Real> & tract, const Vector<Real> & normal,
                              Real delta, const Vector<Real> & opening,
                              Real & delta_max_new, Real delta_max);

  void computeCompressiveTraction(Vector<Real> & tract,
                                  const Vector<Real> & normal, Real delta_n,
                                  const Vector<Real> & opening);

  void computeCoupledTangent(Matrix<Real> & tangent,
                             const Vector<Real> & normal, Real delta,
                             const Vector<Real> & opening, Real delta_max_new);

  void computeCompressivePenalty(Matrix<Real> & tangent,
                                 const Vector<Real> & normal, Real delta_n);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// beta parameter
  Real beta;

  /// contact penalty = initial slope ?
  bool exp_penalty;

  /// Ratio of contact tangent over the initial exponential tangent
  Real contact_tangent;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

// #include "material_cohesive_exponential_inline_impl.hh"

} // namespace akantu

#endif /* AKANTU_MATERIAL_COHESIVE_EXPONENTIAL_HH_ */
