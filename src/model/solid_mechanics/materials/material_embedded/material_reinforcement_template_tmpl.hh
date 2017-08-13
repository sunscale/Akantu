/**
 * @file   material_reinforcement_template_tmpl.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Wed Mar 25 2015
 * @date last modification: Thu Oct 15 2015
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

// /!\ no namespace here !

/* -------------------------------------------------------------------------- */

template<UInt dim, class ConstLaw>
MaterialReinforcementTemplate<dim, ConstLaw>::MaterialReinforcementTemplate(SolidMechanicsModel & a_model,
                                                                            const ID & id):
  Material(a_model, 1,
           dynamic_cast<EmbeddedInterfaceModel &>(a_model).getInterfaceMesh(),
           a_model.getFEEngine("EmbeddedInterfaceFEEngine"), id),
  MaterialReinforcement<dim>(a_model, 1,
                             dynamic_cast<EmbeddedInterfaceModel &>(a_model).getInterfaceMesh(),
                             a_model.getFEEngine("EmbeddedInterfaceFEEngine"), id),
  ConstLaw(a_model, 1,
           dynamic_cast<EmbeddedInterfaceModel &>(a_model).getInterfaceMesh(),
           a_model.getFEEngine("EmbeddedInterfaceFEEngine"), id)
{}

/* -------------------------------------------------------------------------- */

template<UInt dim, class ConstLaw>
MaterialReinforcementTemplate<dim, ConstLaw>::~MaterialReinforcementTemplate()
{}

/* -------------------------------------------------------------------------- */

template<UInt dim, class ConstLaw>
void MaterialReinforcementTemplate<dim, ConstLaw>::initMaterial() {
  MaterialReinforcement<dim>::initMaterial();
  ConstLaw::initMaterial();

  // Needed for plasticity law
  this->ConstLaw::nu = 0.5;
  this->ConstLaw::updateInternalParameters();
}

/* -------------------------------------------------------------------------- */

template<UInt dim, class ConstLaw>
void MaterialReinforcementTemplate<dim, ConstLaw>::computeGradU(const ElementType & el_type,
                                                                GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  
  MaterialReinforcement<dim>::computeGradU(el_type, ghost_type);

  const UInt voigt_size = Material::getTangentStiffnessVoigtSize(dim);

  Array<Real>::matrix_iterator full_gradu_it =
    this->MaterialReinforcement<dim>::gradu_embedded(el_type, ghost_type).begin(dim, dim);
  Array<Real>::matrix_iterator full_gradu_end =
    this->MaterialReinforcement<dim>::gradu_embedded(el_type, ghost_type).end(dim, dim);

  Array<Real>::scalar_iterator gradu_it =
    this->ConstLaw::gradu(el_type, ghost_type).begin();

  Array<Real>::matrix_iterator cosines_it = 
    this->directing_cosines(el_type, ghost_type).begin(voigt_size, voigt_size);


  for (; full_gradu_it != full_gradu_end ; ++full_gradu_it,
                                          ++gradu_it,
                                          ++cosines_it) {
    Matrix<Real> & full_gradu = *full_gradu_it;
    Real & gradu = *gradu_it;
    Matrix<Real> & C = *cosines_it;
    
    computeInterfaceGradUOnQuad(full_gradu, gradu, C);
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template<UInt dim, class ConstLaw>
void MaterialReinforcementTemplate<dim, ConstLaw>::computeInterfaceGradUOnQuad(const Matrix<Real> & full_gradu,
                                                                               Real & gradu,
                                                                               const Matrix<Real> & C) {
  const UInt voigt_size = Material::getTangentStiffnessVoigtSize(dim);

  Matrix<Real> epsilon(dim, dim);
  Vector<Real> e_voigt(voigt_size);
  Vector<Real> e_interface_voigt(voigt_size);

  epsilon = 0.5 * (full_gradu + full_gradu.transpose());
  MaterialReinforcement<dim>::strainTensorToVoigtVector(epsilon, e_voigt);
  e_interface_voigt.mul<false>(C, e_voigt);

  gradu = e_interface_voigt(0);
}
/* -------------------------------------------------------------------------- */

template<UInt dim, class ConstLaw>
void MaterialReinforcementTemplate<dim, ConstLaw>::computeTangentModuli(const ElementType & el_type,
                                                                        Array<Real> & tangent,
                                                                        GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(tangent.getNbComponent() == 1, "Reinforcements only work in 1D");

  ConstLaw::computeTangentModuli(el_type, tangent, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template<UInt dim, class ConstLaw>
void MaterialReinforcementTemplate<dim, ConstLaw>::computeStress(ElementType type,
                                                                 GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  
  ConstLaw::computeStress(type, ghost_type);

  Array<Real>::matrix_iterator full_sigma_it =
    this->MaterialReinforcement<dim>::stress_embedded(type, ghost_type).begin(dim, dim);
  Array<Real>::matrix_iterator full_sigma_end =
    this->MaterialReinforcement<dim>::stress_embedded(type, ghost_type).end(dim, dim);
  Array<Real>::scalar_iterator sigma_it =
    this->ConstLaw::stress(type, ghost_type).begin();
  Array<Real>::scalar_iterator pre_stress_it =
    this->MaterialReinforcement<dim>::pre_stress(type, ghost_type).begin();

  for (; full_sigma_it != full_sigma_end ; ++full_sigma_it, ++sigma_it, ++pre_stress_it) {
    Matrix<Real> & sigma = *full_sigma_it;

    sigma(0, 0) = *sigma_it + *pre_stress_it;
  }
  
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */

template<UInt dim, class ConstLaw>
Real MaterialReinforcementTemplate<dim, ConstLaw>::getEnergy(std::string id) {
  return MaterialReinforcement<dim>::getEnergy(id);
}

/* -------------------------------------------------------------------------- */

template <UInt dim, class ConstLaw>
void MaterialReinforcementTemplate<dim, ConstLaw>::computePotentialEnergy(ElementType type,
                                                                          GhostType ghost_type) {
  const UInt nb_elements = this->element_filter(type, ghost_type).size();
  const UInt nb_quad = this->model->getFEEngine("EmbeddedInterfaceFEEngine").getNbIntegrationPoints(type);
  this->ConstLaw::potential_energy.alloc(nb_quad * nb_elements, 1, type, ghost_type, 0.);

  ConstLaw::computePotentialEnergy(type, ghost_type);
}

