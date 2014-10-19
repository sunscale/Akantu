/**
 * @file   material_cohesive_linear.hh
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Tue Jul 29 2014
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
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

#include "material_cohesive.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_COHESIVE_LINEAR_HH__
#define __AKANTU_MATERIAL_COHESIVE_LINEAR_HH__

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/**
 * Cohesive material linear damage for extrinsic case
 *
 * parameters in the material files :
 *   - sigma_c   : critical stress sigma_c  (default: 0)
 *   - beta      : weighting parameter for sliding and normal opening (default: 0)
 *   - G_cI      : fracture energy for mode I (default: 0)
 *   - G_cII     : fracture energy for mode II (default: 0)
 *   - penalty   : stiffness in compression to prevent penetration
 */
template<UInt spatial_dimension>
class MaterialCohesiveLinear : public MaterialCohesive {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialCohesiveLinear(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// check stress for cohesive elements' insertion
  virtual void checkInsertion();

protected:

  /// constitutive law
  void computeTraction(const Array<Real> & normal,
		       ElementType el_type,
		       GhostType ghost_type = _not_ghost);

  /// compute stress norms on quadrature points for each facet for stress check
  virtual void computeStressNorms(const Array<Real> & facet_stress,
				  Array<Real> & stress_check,
				  Array<Real> & normal_stress,
				  ElementType type_facet);

  /// compute effective stress norm for insertion check
  inline void computeEffectiveNorm(const Matrix<Real> & stress,
				   const Vector<Real> & normal,
				   const Vector<Real> & tangent,
				   Vector<Real> & normal_stress,
				   Real & effective_norm);

  /**
   * Scale insertion traction sigma_c according to the volume of the
   * two elements surrounding a facet
   */
  void scaleInsertionTraction();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// get sigma_c_eff
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(InsertionTraction, sigma_c_eff, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// beta parameter
  Real beta;

  /// beta square inverse to compute effective norm
  Real beta2_inv;

  /// mode I fracture energy
  Real G_cI;

  /// mode II fracture energy
  Real G_cII;

  /// kappa parameter
  Real kappa;

  /// constitutive law scalar to compute delta
  Real beta2_kappa2;

  /// constitutive law scalar to compute traction
  Real beta2_kappa;

  /// penalty coefficient
  Real penalty;

  /// reference volume used to scale sigma_c
  Real volume_s;

  /// weibull exponent used to scale sigma_c
  Real m_s;

  /// critical effective stress
  RandomInternalField<Real, CohesiveInternalField> sigma_c_eff;

  /// critical displacement
  CohesiveInternalField<Real> delta_c;

  /// stress at insertion
  CohesiveInternalField<Real> insertion_stress;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "material_cohesive_linear_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_COHESIVE_LINEAR_HH__ */
