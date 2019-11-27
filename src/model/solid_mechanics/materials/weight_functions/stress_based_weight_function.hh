/**
 * @file   stress_based_weight_function.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Mon Aug 24 2015
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Removed damaged weight function for non local materials
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "base_weight_function.hh"
/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_STRESS_BASED_WEIGHT_FUNCTION_HH__
#define __AKANTU_STRESS_BASED_WEIGHT_FUNCTION_HH__

namespace akantu {
/* -------------------------------------------------------------------------- */
/* Stress Based Weight */
/* -------------------------------------------------------------------------- */
/// based on based on Giry et al.: Stress-based nonlocal damage model,
/// IJSS, 48, 2011
class StressBasedWeightFunction : public BaseWeightFunction {
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  StressBasedWeightFunction(NonLocalManager & manager);

  /* --------------------------------------------------------------------------
   */
  /* Base Weight Function inherited methods */
  /* --------------------------------------------------------------------------
   */
  void init() override;

  inline void updateInternals() override;

  void updatePrincipalStress(GhostType ghost_type);

  inline void updateQuadraturePointsCoordinates(
      ElementTypeMapArray<Real> & quadrature_points_coordinates);

  inline Real operator()(Real r, const IntegrationPoint & q1,
                         const IntegrationPoint & q2);

  /// computation of ellipsoid
  inline Real computeRhoSquare(Real r, Vector<Real> & eigs,
                               Matrix<Real> & eigenvects, Vector<Real> & x_s);

protected:
  inline void setInternal();

private:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  /// tensile strength
  Real ft;

  /// prinicipal stresses
  ElementTypeMapReal * stress_diag;

  /// for preselection of types (optimization)
  ElementTypeMapReal * selected_stress_diag;

  /// principal directions
  ElementTypeMapReal * stress_base;

  /// lenght intrinisic to the material
  ElementTypeMapReal * characteristic_size;
};

#if defined(AKANTU_INCLUDE_INLINE_IMPL)
#include "stress_based_weight_function_inline_impl.cc"
#endif

} // namespace akantu

#endif /* __AKANTU_STRESS_BASED_WEIGHT_FUNCTION_HH__ */
