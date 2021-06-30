/**
 * @file   local_material_damage.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Mon Sep 11 2017
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
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_LOCAL_MATERIAL_DAMAGE_HH_
#define AKANTU_LOCAL_MATERIAL_DAMAGE_HH_

namespace akantu {

class LocalMaterialDamage : public Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  LocalMaterialDamage(SolidMechanicsModel & model, const ID & id = "");

  virtual ~LocalMaterialDamage(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial();

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(Matrix<Real> & grad_u, Matrix<Real> & sigma,
                                  Real & damage);

  /// compute tangent stiffness
  virtual void computeTangentStiffness(__attribute__((unused))
                                       ElementType el_type,
                                       __attribute__((unused))
                                       Array<Real> & tangent_matrix,
                                       __attribute__((unused))
                                       GhostType ghost_type = _not_ghost){};

  /// compute the potential energy for all elements
  void computePotentialEnergy(ElementType el_type);

  /// compute the potential energy for on element
  inline void computePotentialEnergyOnQuad(Matrix<Real> & grad_u,
                                           Matrix<Real> & sigma, Real & epot);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// compute the celerity of wave in the material
  inline Real getCelerity(const Element & element) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Damage, damage, Real);

private:
  /// the young modulus
  Real E;

  /// Poisson coefficient
  Real nu;

  /// First Lamé coefficient
  Real lambda;

  /// Second Lamé coefficient (shear modulus)
  Real mu;

  /// resistance to damage
  Real Yd;

  /// damage threshold
  Real Sd;

  /// Bulk modulus
  Real kpa;

  /// damage internal variable
  InternalField<Real> damage;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "local_material_damage_inline_impl.hh"

} // namespace akantu

#endif /* AKANTU_LOCAL_MATERIAL_DAMAGE_HH_ */
