/**
 * @file   material_orthotropic_damage.hh
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Sun Mar  8 12:49:56 2015
 *
 * @brief  Material for orthotropic damage
 *
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
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_HH__
#define __AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_HH__

namespace akantu {
template <UInt spatial_dimension,
          template <UInt> class Parent = MaterialElastic>
class MaterialOrthotropicDamage : public Parent<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialOrthotropicDamage(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialOrthotropicDamage(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void initMaterial();

  /// compute the tangent stiffness matrix for an element type
  virtual void computeTangentModuli(const ElementType & el_type,
                                    Array<Real> & tangent_matrix,
                                    GhostType ghost_type = _not_ghost);

protected:
  /// update the dissipated energy, must be called after the stress have been
  /// computed
  void updateEnergies(ElementType el_type) override;

  /// compute the tangent stiffness matrix for a given quadrature point
  inline void computeTangentModuliOnQuad(
      Matrix<Real> & tangent, const Matrix<Real> C, const Matrix<Real> & dam,
      const Matrix<Real> & dam_directions, Matrix<Real> & O_prime,
      Matrix<Real> & S_prime, Matrix<Real> & O, Matrix<Real> & S,
      Matrix<Real> & rotation_tmp);

  inline void computeDamageAndStressOnQuad(Matrix<Real> & sigma,
                                           Matrix<Real> & one_minus_D,
                                           Matrix<Real> & root_one_minus_D,
                                           Matrix<Real> & damage,
                                           Matrix<Real> & first_term,
                                           Matrix<Real> & third_term);

  /// rotate a Matrix of size dim*dim into the coordinate system of the FE
  /// computation
  inline void rotateIntoComputationFrame(const Matrix<Real> & to_rotate,
                                         Matrix<Real> & rotated,
                                         const Matrix<Real> & damage_directions,
                                         Matrix<Real> & rotation_tmp);

  /// rotate a Matrix of size dim*dim into the coordinate system of the damage
  inline void rotateIntoNewFrame(const Matrix<Real> & to_rotate,
                                 Matrix<Real> & rotated,
                                 const Matrix<Real> & damage_directions,
                                 Matrix<Real> & rotation_tmp);

  /// compute (1-D)
  inline void computeOneMinusD(Matrix<Real> & one_minus_D,
                               const Matrix<Real> & damage);

  /// compute (1-D)^(1/2)
  inline void computeSqrtOneMinusD(const Matrix<Real> & one_minus_D,
                                   Matrix<Real> & sqrt_one_minus_D);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// give the dissipated energy for the time step
  Real getDissipatedEnergy() const;

  // virtual Real getEnergy(std::string type);
  // virtual Real getEnergy(std::string energy_id, ElementType type, UInt index)
  // {
  //   return Parent<spatial_dimension>::getEnergy(energy_id, type, index);
  // };

  AKANTU_GET_MACRO_NOT_CONST(Damage, damage, ElementTypeMapArray<Real> &);
  AKANTU_GET_MACRO(Damage, damage, const ElementTypeMapArray<Real> &);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Damage, damage, Real)

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// damage internal variable
  InternalField<Real> damage;

  /// dissipated energy
  InternalField<Real> dissipated_energy;

  /// contain the current value of @f$ \int_0^{\epsilon}\sigma(\omega)d\omega
  /// @f$ the dissipated energy
  InternalField<Real> int_sigma;

  /// direction vectors for the damage frame
  InternalField<Real> damage_dir_vecs;

  Real eta;

  /// maximum damage value
  Real max_damage;
};

} // namespace akantu

#include "material_orthotropic_damage_tmpl.hh"

#endif /* __AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_HH__ */
