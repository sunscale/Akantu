/**
 * @file   material_reinforcement_template.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Wed Mar 25 2015
 * @date last modification: Mon Jun 01 2015
 *
 * @brief  Reinforcement material templated with constitutive law
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

#ifndef __AKANTU_MATERIAL_REINFORCEMENT_TEMPLATE_HH__
#define __AKANTU_MATERIAL_REINFORCEMENT_TEMPLATE_HH__

/* -------------------------------------------------------------------------- */

#include "aka_common.hh"
#include "material_reinforcement.hh"
#include "material_elastic.hh"
#include "material_linear_isotropic_hardening.hh"

namespace akantu {

/**
 * @brief Implementation of MaterialReinforcement with 1D constitutive law
 * @see MaterialReinforcement, MaterialElastic
 *
 * This class is a reinforcement featuring a constitutive law.
 * <strong>Be careful !</strong> Because of multiple inheritance, this class
 * forms a diamond.
 */
template<UInt dim, class ConstLaw = MaterialElastic<1> >
class MaterialReinforcementTemplate : public MaterialReinforcement<dim>,
                                      public ConstLaw {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// Constructor
  MaterialReinforcementTemplate(SolidMechanicsModel & a_model, const ID & id = "");

  /// Destructor
  virtual ~MaterialReinforcementTemplate();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Initialises the material
  void initMaterial();

  /// Compute the stiffness parameter for elements of a type
  virtual void computeTangentModuli(const ElementType & type,
                                    Array<Real> & tangent,
                                    GhostType ghost_type);

  /// Computes stress used by constitutive law
  virtual void computeStress(ElementType type, GhostType ghost_type);

  /// Computes gradu to be used by the constitutive law
  virtual void computeGradU(const ElementType & type, GhostType ghost_type);

  /// Compute the potential energy of the reinforcement
  virtual void computePotentialEnergy(ElementType type, GhostType ghost_type = _not_ghost);

  /// Get energy in reinforcement (currently limited to potential)
  virtual Real getEnergy(std::string id);

protected:
  /**
   * @brief Compute interface gradu from bulk gradu
   * \f[
   *  \varepsilon_s = C \varepsilon_c
   * \f]
   */
  inline void computeInterfaceGradUOnQuad(const Matrix<Real> & full_gradu,
                                          Real & gradu,
                                          const Matrix<Real> & C);

};

#include "material_reinforcement_template_tmpl.hh"

} // akantu


#endif // __AKANTU_MATERIAL_REINFORCEMENT_TEMPLATE_HH__
