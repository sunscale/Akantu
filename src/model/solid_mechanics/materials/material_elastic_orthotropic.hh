/**
 * @file   material_elastic_orthotropic.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Enrico Milanese <enrico.milanese@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Feb 16 2018
 *
 * @brief  Orthotropic elastic material
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material_elastic_linear_anisotropic.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_ELASTIC_ORTHOTROPIC_HH__
#define __AKANTU_MATERIAL_ELASTIC_ORTHOTROPIC_HH__

namespace akantu {

/**
 * Orthotropic elastic material
 *
 * parameters in the material files :
 *   - n1   : direction of x-axis in material base, normalisation not necessary
 * (default: {1, 0, 0})
 *   - n2   : direction of y-axis in material base, normalisation not necessary
 * (default: {0, 1, 0})
 *   - n3   : direction of z-axis in material base, normalisation not necessary
 * (default: {0, 0, 1})
 *   - rho  : density (default: 0)
 *   - E1   : Young's modulus along n1 (default: 0)
 *   - E2   : Young's modulus along n2 (default: 0)
 *   - E3   : Young's modulus along n3 (default: 0)
 *   - nu12 : Poisson's ratio along 12 (default: 0)
 *   - nu13 : Poisson's ratio along 13 (default: 0)
 *   - nu23 : Poisson's ratio along 23 (default: 0)
 *   - G12  : Shear modulus along 12 (default: 0)
 *   - G13  : Shear modulus along 13 (default: 0)
 *   - G23  : Shear modulus along 23 (default: 0)
 */

template <UInt Dim>
class MaterialElasticOrthotropic
    : public MaterialElasticLinearAnisotropic<Dim> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialElasticOrthotropic(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial() override;

  void updateInternalParameters() override;

  void
  computePotentialEnergyByElement(ElementType type, UInt index,
                                  Vector<Real> & epot_on_quad_points) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(E1, E1, Real);
  AKANTU_GET_MACRO(E2, E2, Real);
  AKANTU_GET_MACRO(E3, E3, Real);
  AKANTU_GET_MACRO(Nu12, nu12, Real);
  AKANTU_GET_MACRO(Nu13, nu13, Real);
  AKANTU_GET_MACRO(Nu23, nu23, Real);
  AKANTU_GET_MACRO(G12, G12, Real);
  AKANTU_GET_MACRO(G13, G13, Real);
  AKANTU_GET_MACRO(G23, G23, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the n1 young modulus
  Real E1{0.};

  /// the n2 young modulus
  Real E2{0.};

  /// the n3 young modulus
  Real E3{0.};

  /// 12 Poisson coefficient
  Real nu12{0.};

  /// 13 Poisson coefficient
  Real nu13{0.};

  /// 23 Poisson coefficient
  Real nu23{0.};

  /// 12 shear modulus
  Real G12{0.};

  /// 13 shear modulus
  Real G13{0.};

  /// 23 shear modulus
  Real G23{0.};
};

} // namespace akantu

#endif /* __AKANTU_MATERIAL_ELASTIC_ORTHOTROPIC_HH__ */
