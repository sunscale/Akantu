/**
 * @file   material_reinforcement_template.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Mon Mar 16 2015
 * @date last modification: Mon Mar 16 2015
 *
 * @brief  Reinforcement material templated with constitutive law
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_MATERIAL_REINFORCEMENT_TEMPLATE_HH__
#define __AKANTU_MATERIAL_REINFORCEMENT_TEMPLATE_HH__

/* -------------------------------------------------------------------------- */

#include "aka_common.hh"
#include "material_reinforcement.hh"
#include "material_elastic.hh"

__BEGIN_AKANTU__

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

  /// TODO determine the significance and how to compute the energies
  virtual Real getEnergy(std::string id) { return 0.; }

protected:
  /// Compute interface gradu from bulk gradu
  inline void computeInterfaceGradUOnQuad(const Matrix<Real> & full_gradu,
                                          Real & gradu,
                                          const Matrix<Real> & C);

};

#include "material_reinforcement_template_inline_impl.cc"

__END_AKANTU__

#endif // __AKANTU_MATERIAL_REINFORCEMENT_TEMPLATE_HH__
