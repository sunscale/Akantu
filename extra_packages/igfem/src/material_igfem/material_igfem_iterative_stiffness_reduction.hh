/**
 * @file   material_igfem_iterative_stiffness_reduction.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Mar  9 19:44:22 2016
 *
 * @brief  Material for iterative stiffness reduction by contant factor
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
#include "material_igfem_saw_tooth_damage.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_IGFEM_ITERATIVE_STIFFNESS_REDUCTION_HH__
#define __AKANTU_MATERIAL_IGFEM_ITERATIVE_STIFFNESS_REDUCTION_HH__

namespace akantu {

/**
 * Material damage iterative
 *
 * parameters in the material files :
 *   - Gfx
 *   - h
 *   - Sc
 */
/// Proposed by Rots and Invernizzi, 2004: Regularized sequentially linear
// saw-tooth softening model (section 4.2)
template <UInt spatial_dimension>
class MaterialIGFEMIterativeStiffnessReduction
    : public MaterialIGFEMSawToothDamage<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialIGFEMIterativeStiffnessReduction(SolidMechanicsModel & model,
                                           const ID & id = "");

  virtual ~MaterialIGFEMIterativeStiffnessReduction(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// set the material parameters
  virtual void initMaterial();

  /// compute the equivalent stress on each Gauss point (i.e. the max prinicpal
  /// stress) and normalize it by the tensile stiffness
  virtual void
  computeNormalizedEquivalentStress(const Array<Real> & grad_u,
                                    ElementType el_type,
                                    GhostType ghost_type = _not_ghost);

  /// update internal field damage
  virtual UInt updateDamage();

  virtual void onElementsAdded(const Array<Element> & element_list,
                               const NewElementsEvent & event);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the ultimate strain
  IGFEMInternalField<Real> eps_u;

  /// the reduction
  IGFEMInternalField<UInt> reduction_step;

  /// the tangent of the tensile stress-strain softening
  IGFEMInternalField<Real> D;

  /// fracture energy
  Real Gf;

  /// crack band width for normalization of fracture energy
  Real crack_band_width;

  /// the number of total reductions steps until complete failure
  UInt max_reductions;

  /// the reduction constant (denoated by a in the paper of rots)
  Real reduction_constant;
};

} // namespace akantu

#endif /* __AKANTU_MATERIAL_IGFEM_ITERATIVE_STIFFNESS_REDUCTION_HH__ */
