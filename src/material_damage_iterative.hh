/**
 * @file   material_damage_iterative.hh
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Fri Feb 14 14:34:37 2014
 *
 * @brief  Specialization of the class material damage to damage only one gauss 
 * point at a time and propagate damage in a linear way. Max principal stress 
 * criterion is used as a failure criterion.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_MATERIAL_DAMAGE_ITERATIVE_HH__
#define __AKANTU_MATERIAL_DAMAGE_ITERATIVE_HH__

__BEGIN_AKANTU__

/**
 * Material damage iterative
 *
 * parameters in the material files :
 *   - Sc  
 */
template<UInt spatial_dimension>
class MaterialDamageIterative : public MaterialDamage<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialDamageIterative(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialDamageIterative() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  ///  virtual void updateInternalParameters();

  virtual void computeAllStresses(GhostType ghost_type = _not_ghost);

  /// update internal field damage
  UInt updateDamage();

  /// update energies after damage has been updated
  virtual void updateEnergiesAfterDamage(ElementType el_type, GhostType ghost_typ);
  
  virtual void onBeginningSolveStep(const AnalysisMethod & method) { };

  virtual void onEndSolveStep(const AnalysisMethod & method) { };
protected:
  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  ///compute the equivalent stress on each Gauss point (i.e. the max prinicpal stress) and normalize it by the tensile strength
  void computeNormalizedEquivalentStress(const Array<Real> & grad_u,
                                         ElementType el_type, GhostType ghost_type = _not_ghost);

  /// find max normalized equivalent stress
  void findMaxNormalizedEquivalentStress(ElementType el_type, GhostType ghost_type = _not_ghost);


  inline void computeDamageAndStressOnQuad(Matrix<Real> & sigma,
                                           Real & dam);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get max normalized equivalent stress
  AKANTU_GET_MACRO(NormMaxEquivalentStress, norm_max_equivalent_stress, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// resistance to damage
  RandomInternalField<Real> Sc;

  /// internal field to store equivalent stress on each Gauss point
  InternalField<Real> equivalent_stress;

  /// damage increment
  Real prescribed_dam;

  /// maximum equivalent stress
  Real norm_max_equivalent_stress;

  /// deviation from max stress at which Gauss point will still get damaged 
  Real dam_tolerance;

  /// define damage threshold at which damage will be set to 1 
  Real dam_threshold; 

  /// maximum damage value
  Real max_damage;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_damage_iterative_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_DAMAGE_ITERATIVE_HH__ */
