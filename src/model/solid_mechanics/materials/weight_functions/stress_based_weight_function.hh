/**
 * @file   stress_based_weight_function.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Fri Apr 13 2012
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  Removed damaged weight function for non local materials
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
#include "base_weight_function.hh"
/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_STRESS_BASED_WEIGHT_FUNCTION_HH__
#define __AKANTU_STRESS_BASED_WEIGHT_FUNCTION_HH__

__BEGIN_AKANTU__
/* -------------------------------------------------------------------------- */
/* Stress Based Weight                                                         */
/* -------------------------------------------------------------------------- */
/// based on based on Giry et al.: Stress-based nonlocal damage model,
/// IJSS, 48, 2011
template<UInt spatial_dimension>
class StressBasedWeightFunction : public BaseWeightFunction<spatial_dimension> {
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  StressBasedWeightFunction(Material & material);

  /* -------------------------------------------------------------------------- */
  /* Base Weight Function inherited methods                                     */
  /* -------------------------------------------------------------------------- */
  void init();

  virtual inline void updateInternals(__attribute__((unused)) const ElementTypeMapArray<Real> & quadrature_points_coordinates);

  void updatePrincipalStress(GhostType ghost_type);

  inline void updateQuadraturePointsCoordinates(ElementTypeMapArray<Real> & quadrature_points_coordinates);

  inline void selectType(ElementType type1, GhostType ghost_type1,
                         ElementType type2, GhostType ghost_type2);

  inline Real operator()(Real r,
			 const QuadraturePoint & q1,
			 const QuadraturePoint & q2);

  /// computation of ellipsoid
  inline Real computeRhoSquare(Real r,
                               Vector<Real> & eigs,
                               Matrix<Real> & eigenvects,
                               Vector<Real> & x_s);

private:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  /// tensile strength
  Real ft;

  /// prinicipal stresses
  InternalField<Real> stress_diag;

  /// for preselection of types (optimization)
  Array<Real> * selected_stress_diag;

  /// principal directions
  InternalField<Real> stress_base;

  /// for preselection of types (optimization)
  Array<Real> * selected_stress_base;

  //  InternalField<Real> quadrature_points_coordinates;
  /// for preselection of types (optimization)
  Array<Real> * selected_position_1;
  Array<Real> * selected_position_2;

  /// lenght intrinisic to the material
  InternalField<Real> characteristic_size;

  /// for preselection of types (optimization)
  Array<Real> * selected_characteristic_size;
};

#include "stress_based_weight_function_tmpl.hh"
#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "stress_based_weight_function_inline_impl.cc"
#endif

__END_AKANTU__

#endif /* __AKANTU_STRESS_BASED_WEIGHT_FUNCTION_HH__ */
