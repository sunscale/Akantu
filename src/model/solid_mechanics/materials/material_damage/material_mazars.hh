/**
 * @file   material_mazars.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 *
 * @date creation: Wed Apr 06 2011
 * @date last modification: Fri Dec 20 2013
 *
 * @brief  Material Following the Mazars law for damage evolution
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
#include "material_damage.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_MAZARS_HH__
#define __AKANTU_MATERIAL_MAZARS_HH__

__BEGIN_AKANTU__

/**
 * Material Mazars
 *
 * parameters in the material files :
 *   - rho  : density (default: 0)
 *   - E    : Young's modulus (default: 0)
 *   - nu   : Poisson's ratio (default: 1/2)
 *   - K0   : Damage threshold
 *   - At   : Parameter damage traction 1
 *   - Bt   : Parameter damage traction 2
 *   - Ac   : Parameter damage compression 1
 *   - Bc   : Parameter damage compression 2
 *   - beta : Parameter for shear
 */
template<UInt spatial_dimension>
class MaterialMazars : public MaterialDamage<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialMazars(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialMazars() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(const Matrix<Real> & grad_u,
				  Matrix<Real> & sigma,
				  Real & damage,
				  Real & Ehat);

  inline void computeDamageAndStressOnQuad(const Matrix<Real> & grad_u,
					   Matrix<Real> & sigma,
					   Real & damage,
					   Real & Ehat);

  inline void computeDamageOnQuad(const Real & epsilon_equ,
				  const Matrix<Real> & sigma,
				  const Vector<Real> & epsilon_princ,
				  Real & dam);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// damage threshold
  RandomInternalField<Real>  K0;
  ///parameter damage traction 1
  Real At ;
  ///parameter damage traction 2
  Real Bt ;
  ///parameter damage compression 1
  Real Ac ;
  ///parameter damage compression 2
  Real Bc ;
  ///parameter for shear
  Real beta ;

  /// specify the variable to average false = ehat, true = damage (only valid for non local version)
  bool damage_in_compute_stress;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_mazars_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_MAZARS_HH__ */
