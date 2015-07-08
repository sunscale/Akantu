/**
 * @file   material_igfem_elastic.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Specializaton of material class for the igfem elastic material
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */
/* -------------------------------------------------------------------------- */
#include "material_igfem_elastic.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt dim>
MaterialIGFEMElastic<dim>::MaterialIGFEMElastic(SolidMechanicsModel & model, const ID & id)  :
  Material(model, id),
  Parent(model, id),
  E(this->nb_sub_materials),
  nu(this->nb_sub_materials),
  lambda("lambda", *this),
  mu("mu", *this),
  kpa("kappa", *this) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialIGFEMElastic<dim>::initialize() {
  this->registerParam("E"      , E      , _pat_parsable | _pat_modifiable, "Young's modulus"        );
  this->registerParam("nu"     , nu     , _pat_parsable | _pat_modifiable, "Poisson's ratio"        );
  ///  this->registerParam("Plane_Stress", plane_stress, false, _pat_parsmod, "Is plane stress");

  this->lambda           .initialize(1);
  this->mu           .initialize(1);
  this->kpa           .initialize(1);
}


/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialIGFEMElastic<dim>::initMaterial() {
  AKANTU_DEBUG_IN();

  Parent::initMaterial();

  if (dim == 1) this->nu.clear(); /// set the Poisson ratios to zero

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::updateElasticInternals(GhostType ghost_type) {

  SolidMechanicsModelIGFEM * igfem_model = static_cast<SolidMechanicsModelIGFEM*>(this->model);
  const Mesh & mesh = this->model->getMesh();
  Array<Real> nodes_coordinates(mesh.getNodes(), true);
  Array<Real>::const_vector_iterator nodes_it = nodes_coordinates.begin(spatial_dimension);

  /// compute the Lamé constants for both sub-materials
  Vector<Real> lambda_per_sub_mat(this->nb_sub_materials);
  Vector<Real> mu_per_sub_mat(this->nb_sub_materials);
  Vector<Real> kpa_per_sub_mat(this->nb_sub_materials);
  
  this->updateElasticConstants(lambda_per_sub_mat, mu_per_sub_mat, kpa_per_sub_mat);


  /// loop over all types in the material
  typedef ElementTypeMapArray<UInt>:: type_iterator iterator;
  iterator it = this->element_filter.firstType(spatial_dimension, ghost_type, _ek_igfem);
  iterator last_type = this->element_filter.lastType(spatial_dimension, ghost_type, _ek_igfem);
  UInt index_1 = 0;
  UInt index_2 = 0;
  /// loop over all types in the filter
  for(; it != last_type; ++it) {
    ElementType el_type = *it;
    UInt nb_nodes_per_el = mesh.getNbNodesPerElement(el_type);
    UInt nb_parent_nodes = IGFEMHelper::getNbParentNodes(el_type);
    Vector<bool> is_inside(nb_parent_nodes);
    const Array<UInt> & connectivity = mesh.getConnectivity(el_type, ghost_type);
    Array<UInt>::const_vector_iterator connec_it = connectivity.begin(nb_nodes_per_el);

    /// get the number of quadrature points for the two sub-elements
    UInt quads_1 = IGFEMHelper::getNbQuadraturePoints(el_type, 0);
    UInt quads_2 = IGFEMHelper::getNbQuadraturePoints(el_type, 1);
  
    /// get pointer to internals for given type
    Real * lambda_ptr = this->lambda(el_type, ghost_type).storage();
    Real * mu_ptr = this->mu(el_type, ghost_type).storage();
    Real * kpa_ptr = this->kpa(el_type, ghost_type).storage();

    /// loop all elements for the given type
    const Array<UInt> & filter   = this->element_filter(el_type,ghost_type);
    UInt nb_elements = filter.getSize();
    for (UInt e = 0; e < nb_elements; ++e, ++connec_it) {
      for (UInt i = 0; i < nb_parent_nodes; ++i) {
	Vector<Real> node = nodes_it[(*connec_it)(i)];
	is_inside(i) = igfem_model->isInside(node, this->name_sub_mat_1);
      }

      UInt orientation = IGFEMHelper::getElementOrientation(el_type, is_inside);

      if (orientation) {
	index_1 = 0;
	index_2 = 1;
      }
      else {
	index_1 = 1;
	index_2 = 0;	
      }

      for (UInt q = 0; q < quads_1; ++q, ++lambda_ptr, ++mu_ptr, ++kpa_ptr) {
	*lambda_ptr = lambda_per_sub_mat(index_1);
	*mu_ptr = mu_per_sub_mat(index_1);
	*kpa_ptr = kpa_per_sub_mat(index_1); 
      }
      for (UInt q = 0; q < quads_2; ++q, ++lambda_ptr, ++mu_ptr, ++kpa_ptr) {
	*lambda_ptr = lambda_per_sub_mat(index_2);
	*mu_ptr = mu_per_sub_mat(index_2);
	*kpa_ptr = kpa_per_sub_mat(index_2); 
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialIGFEMElastic<dim>::updateElasticConstants(Vector<Real> & lambda_vec, Vector<Real> & mu_vec, Vector<Real> & kpa_vec) {

  for (UInt i = 0; i < this->nb_sub_materials; ++i) {
    lambda_vec(i) = this->nu(i) * this->E(i) / ((1 + this->nu(i)) * (1 - 2*this->nu(i)));
    mu_vec(i)     = this->E(i) / (2 * (1 + this->nu(i)));

    kpa_vec(i)    = lambda_vec(i) + 2./3. * mu_vec(i);
  }
}

/* -------------------------------------------------------------------------- */
template<>
void MaterialIGFEMElastic<2>::updateElasticConstants(Vector<Real> & lambda_vec, Vector<Real> & mu_vec, Vector<Real> & kpa_vec) {

  for (UInt i = 0; i < this->nb_sub_materials; ++i) {
    lambda_vec(i) = this->nu(i) * this->E(i) / ((1 + this->nu(1)) * (1 - 2*this->nu(i)));
    mu_vec(i)     = this->E(i) / (2 * (1 + this->nu(i)));

    if(this->plane_stress) lambda_vec(i) = this->nu(i) * this->E(i) / ((1 + this->nu(i))*(1 - this->nu(i)));

    kpa_vec(i)    = lambda_vec(i) + 2./3. * mu_vec(i);
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
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
template<UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::computeTangentModuli(__attribute__((unused)) const ElementType & el_type,
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
template<UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::computePotentialEnergy(ElementType el_type,
                                                                GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // MaterialThermal<spatial_dimension>::computePotentialEnergy(el_type, ghost_type);

  // if(ghost_type != _not_ghost) return;
  // Array<Real>::scalar_iterator epot = this->potential_energy(el_type, ghost_type).begin();

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
template<UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::computePotentialEnergyByElement(ElementType type, UInt index,
                                                                         Vector<Real> & epot_on_quad_points) {
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

  // UInt nb_quadrature_points = this->model->getFEEngine().getNbQuadraturePoints(type);

  // gradu_it  += index*nb_quadrature_points;
  // gradu_end += (index+1)*nb_quadrature_points;
  // stress_it  += index*nb_quadrature_points;

  // Real * epot_quad = epot_on_quad_points.storage();

  // Matrix<Real> grad_u(spatial_dimension, spatial_dimension);

  // for(;gradu_it != gradu_end; ++gradu_it, ++stress_it, ++epot_quad) {

  //   if (this->finite_deformation)
  //     this->template gradUToGreenStrain<spatial_dimension>(*gradu_it, grad_u);
  //   else
  //     grad_u.copy(*gradu_it);

  //   this->computePotentialEnergyOnQuad(grad_u, *stress_it, *epot_quad);
  // }
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::onElementsAdded(const Array<Element> & element_list,
							      const NewElementsEvent & event) {
  
  Parent::onElementsAdded(element_list, event);
  updateElasticInternals(_not_ghost);

};

INSTANTIATE_MATERIAL(MaterialIGFEMElastic);

__END_AKANTU__
