/**
 * @file   material_igfem_elastic.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Specializaton of material class for the igfem elastic material
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */
/* -------------------------------------------------------------------------- */
#include "material_igfem_elastic.hh"
#include "material_elastic.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
MaterialIGFEMElastic<dim>::MaterialIGFEMElastic(SolidMechanicsModel & model,
                                                const ID & id)
    : Material(model, id), Parent(model, id), lambda("lambda", *this),
      mu("mu", *this), kpa("kappa", *this) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim> void MaterialIGFEMElastic<dim>::initialize() {
  this->lambda.initialize(1);
  this->mu.initialize(1);
  this->kpa.initialize(1);
}

/* -------------------------------------------------------------------------- */
template <UInt dim> void MaterialIGFEMElastic<dim>::initMaterial() {
  AKANTU_DEBUG_IN();

  Parent::initMaterial();

  /// insert the sub_material names into the map
  this->sub_material_names[0] = this->name_sub_mat_1;
  this->sub_material_names[1] = this->name_sub_mat_2;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::updateElasticInternals(
    const Array<Element> & element_list) {

  /// compute the Lamé constants for both sub-materials
  Vector<Real> lambda_per_sub_mat(this->nb_sub_materials);
  Vector<Real> mu_per_sub_mat(this->nb_sub_materials);
  Vector<Real> kpa_per_sub_mat(this->nb_sub_materials);
  for (UInt i = 0; i < this->nb_sub_materials; ++i) {
    ID mat_name = this->sub_material_names[i];
    const MaterialElastic<spatial_dimension> & mat =
        dynamic_cast<MaterialElastic<spatial_dimension> &>(
            this->model->getMaterial(mat_name));
    lambda_per_sub_mat(i) = mat.getLambda();
    mu_per_sub_mat(i) = mat.getMu();
    kpa_per_sub_mat(i) = mat.getKappa();
  }

  for (ghost_type_t::iterator g = ghost_type_t::begin();
       g != ghost_type_t::end(); ++g) {
    GhostType ghost_type = *g;
    /// loop over all types in the material
    typedef ElementTypeMapArray<UInt>::type_iterator iterator;
    iterator it = this->element_filter.firstType(spatial_dimension, ghost_type,
                                                 _ek_igfem);
    iterator last_type =
        this->element_filter.lastType(spatial_dimension, ghost_type, _ek_igfem);
    /// loop over all types in the filter
    for (; it != last_type; ++it) {
      ElementType el_type = *it;
      if (el_type == _igfem_triangle_4)
        this->template setSubMaterial<_igfem_triangle_4>(element_list,
                                                         ghost_type);
      else if (el_type == _igfem_triangle_5)
        this->template setSubMaterial<_igfem_triangle_5>(element_list,
                                                         ghost_type);
      else
        AKANTU_ERROR("There is currently no other IGFEM type implemented");

      UInt nb_element = this->element_filter(el_type, ghost_type).getSize();
      UInt nb_quads = this->fem->getNbIntegrationPoints(el_type);
      /// get pointer to internals for given type
      Real * lambda_ptr = this->lambda(el_type, ghost_type).storage();
      Real * mu_ptr = this->mu(el_type, ghost_type).storage();
      Real * kpa_ptr = this->kpa(el_type, ghost_type).storage();
      UInt * sub_mat_ptr = this->sub_material(el_type, ghost_type).storage();
      for (UInt q = 0; q < nb_element * nb_quads;
           ++q, ++lambda_ptr, ++mu_ptr, ++kpa_ptr, ++sub_mat_ptr) {
        UInt index = *sub_mat_ptr;
        *lambda_ptr = lambda_per_sub_mat(index);
        *mu_ptr = mu_per_sub_mat(index);
        *kpa_ptr = kpa_per_sub_mat(index);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Parent::computeStress(el_type, ghost_type);

  if (!this->finite_deformation) {
    /// get pointer to internals
    Real * lambda_ptr = this->lambda(el_type, ghost_type).storage();
    Real * mu_ptr = this->mu(el_type, ghost_type).storage();
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
    this->computeStressOnQuad(grad_u, sigma, *lambda_ptr, *mu_ptr);
    ++lambda_ptr;
    ++mu_ptr;
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  } else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::computeTangentModuli(
    __attribute__((unused)) ElementType el_type,
    Array<Real> & tangent_matrix,
    __attribute__((unused)) GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// get pointer to internals
  Real * lambda_ptr = this->lambda(el_type, ghost_type).storage();
  Real * mu_ptr = this->mu(el_type, ghost_type).storage();
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  this->computeTangentModuliOnQuad(tangent, *lambda_ptr, *mu_ptr);
  ++lambda_ptr;
  ++mu_ptr;
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::computePotentialEnergy(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // MaterialThermal<spatial_dimension>::computePotentialEnergy(el_type,
  // ghost_type);

  // if(ghost_type != _not_ghost) return;
  // Array<Real>::scalar_iterator epot = this->potential_energy(el_type,
  // ghost_type).begin();

  // if (!this->finite_deformation) {
  //   MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  //   this->computePotentialEnergyOnQuad(grad_u, sigma, *epot);
  //   ++epot;

  //   MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  // } else {
  //   Matrix<Real> E(spatial_dimension, spatial_dimension);

  //   MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  //   this->template gradUToGreenStrain<spatial_dimension>(grad_u, E);

  //   this->computePotentialEnergyOnQuad(E, sigma, *epot);
  //   ++epot;

  //   MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  // }

  AKANTU_DEBUG_TO_IMPLEMENT();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::computePotentialEnergyByElement(
    ElementType type, UInt index, Vector<Real> & epot_on_quad_points) {
  // Array<Real>::matrix_iterator gradu_it =
  //   this->gradu(type).begin(spatial_dimension,
  //                            spatial_dimension);
  // Array<Real>::matrix_iterator gradu_end =
  //   this->gradu(type).begin(spatial_dimension,
  //                            spatial_dimension);
  // Array<Real>::matrix_iterator stress_it =
  //   this->stress(type).begin(spatial_dimension,
  //                            spatial_dimension);

  // if (this->finite_deformation)
  //   stress_it = this->piola_kirchhoff_2(type).begin(spatial_dimension,
  //                                                        spatial_dimension);

  // UInt nb_quadrature_points =
  // this->model->getFEEngine().getNbQuadraturePoints(type);

  // gradu_it  += index*nb_quadrature_points;
  // gradu_end += (index+1)*nb_quadrature_points;
  // stress_it  += index*nb_quadrature_points;

  // Real * epot_quad = epot_on_quad_points.storage();

  // Matrix<Real> grad_u(spatial_dimension, spatial_dimension);

  // for(;gradu_it != gradu_end; ++gradu_it, ++stress_it, ++epot_quad) {

  //   if (this->finite_deformation)
  //     this->template gradUToGreenStrain<spatial_dimension>(*gradu_it,
  //     grad_u);
  //   else
  //     grad_u.copy(*gradu_it);

  //   this->computePotentialEnergyOnQuad(grad_u, *stress_it, *epot_quad);
  // }
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::onElementsAdded(
    const Array<Element> & element_list, const NewElementsEvent & event) {

  Parent::onElementsAdded(element_list, event);
  updateElasticInternals(element_list);
};

INSTANTIATE_MATERIAL(MaterialIGFEMElastic);

} // namespace akantu
