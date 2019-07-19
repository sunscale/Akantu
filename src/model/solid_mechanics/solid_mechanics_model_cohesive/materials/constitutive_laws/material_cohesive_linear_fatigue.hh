/**
 * @file   material_cohesive_linear_fatigue.hh
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Dec 03 2017
 *
 * @brief  Linear irreversible cohesive law with dissipative
 * unloading-reloading cycles
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

#include "material_cohesive_linear.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_COHESIVE_LINEAR_FATIGUE_HH__
#define __AKANTU_MATERIAL_COHESIVE_LINEAR_FATIGUE_HH__

/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * Linear irreversible cohesive law with dissipative
 * unloading-reloading cycles
 *
 * This law uses two different stiffnesses during unloading and
 * reloading. The implementation is based on the article entitled "A
 * cohesive model for fatigue crack growth" by Nguyen, Repetto, Ortiz
 * and Radovitzky (2001). This law is identical to the
 * MaterialCohesiveLinear one except for the unloading-reloading
 * phase.
 *
 * input parameter:
 *
 * - delta_f : it must be greater than delta_c and it is inversely
 *      proportional to the dissipation in the unloading-reloading
 *      cycles (default: delta_c)
 */

template <UInt spatial_dimension>
class MaterialCohesiveLinearFatigue
    : public MaterialCohesiveLinear<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialCohesiveLinearFatigue(SolidMechanicsModel & model,
                                const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material parameters
  void initMaterial() override;

protected:
  /// constitutive law
  void computeTraction(const Array<Real> & normal, ElementType el_type,
                       GhostType ghost_type = _not_ghost) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the switches
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Switches, switches, UInt);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// delta_f parameter
  Real delta_f;

  /// variable saying if delta_f is equal to delta_max for each
  /// element when the traction is computed
  bool progressive_delta_f;

  /// count the opening/closing switches per element
  bool count_switches;

  /// delta of the previous step
  CohesiveInternalField<Real> delta_prec;

  /// stiffness for reloading
  CohesiveInternalField<Real> K_plus;

  /// stiffness for unloading
  CohesiveInternalField<Real> K_minus;

  /// 1D traction in the cohesive law
  CohesiveInternalField<Real> T_1d;

  /// Number of opening/closing switches
  CohesiveInternalField<UInt> switches;

  /// delta increment of the previous time step
  CohesiveInternalField<Real> delta_dot_prec;

  /// has the element passed to normal regime (not in fatigue anymore)
  CohesiveInternalField<bool> normal_regime;

  /// ratio indicating until what point fatigue is applied in the cohesive law
  Real fatigue_ratio;
};

} // namespace akantu

#endif /* __AKANTU_MATERIAL_COHESIVE_LINEAR_FATIGUE_HH__ */
