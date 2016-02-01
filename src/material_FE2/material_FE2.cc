/**
 * @file   material_FE2.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @brief Material for multi-scale simulations. It stores an
 * underlying RVE on each integration point of the material.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_FE2.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialFE2<spatial_dimension>::MaterialFE2(SolidMechanicsModel & model,
					    const ID & id)  :
  Material(model, id),
  C("material_stiffness", *this)
 {
  AKANTU_DEBUG_IN();

  this->C.initialize(voigt_h::size * voigt_h::size);
  this->initialize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialFE2<spatial_dimension>::~MaterialFE2() {
  for (UInt i = 0; i < RVEs.size(); ++i) {
    delete RVEs[i];
  }
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialFE2<dim>::initialize() {
  this->registerParam("element_type"      ,el_type, _triangle_3             ,  _pat_parsable | _pat_modifiable, "element type in RVE mesh" );
  this->registerParam("mesh_file"          ,mesh_file                , _pat_parsable | _pat_modifiable, "the mesh file for the RVE");
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::initMaterial() {

  Material::initMaterial();

  /// compute the number of integration points in this material and resize the RVE vector
  UInt nb_integration_points = this->element_filter(this->el_type, _not_ghost).getSize() 
    * this->fem->getNbIntegrationPoints(this->el_type);
  RVEs.resize(nb_integration_points);

  /// create a SolidMechanicsModel on each integration point of the material
  std::vector<SolidMechanicsModelRVE *>::iterator RVE_it = RVEs.begin();

  Array<Real>::matrix_iterator C_it =
    this->C(this->el_type).begin(voigt_h::size, voigt_h::size);
  for (UInt i = 1; i < nb_integration_points+1; ++RVE_it, ++i, ++C_it) {
    Mesh mesh(spatial_dimension, "RVE_mesh", i);
    mesh.read(mesh_file);
    *RVE_it = new SolidMechanicsModelRVE(mesh, true, _all_dimensions, "SMM_RVE", i); 
    (*RVE_it)->initFull();
    /// compute intial stiffness of the RVE
    (*RVE_it)->homogenizeStiffness(*C_it);
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::computeStress(ElementType el_type,
						   GhostType ghost_type) {
  // Wikipedia convention:
  // 2*eps_ij (i!=j) = voigt_eps_I
  // http://en.wikipedia.org/wiki/Voigt_notation
  AKANTU_DEBUG_IN();

  Array<Real>::const_matrix_iterator C_it = this->C(el_type, ghost_type).begin(voigt_h::size, 
									       voigt_h::size);

  // create vectors to store stress and strain in Voigt notation
  // for efficient computation of stress
  Vector<Real> voigt_strain(voigt_h::size);
  Vector<Real> voigt_stress(voigt_h::size);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  const Matrix<Real> & C_mat = *C_it;

  /// copy strains in Voigt notation
  for(UInt I = 0; I < voigt_h::size; ++I) {
    /// copy stress in
    Real voigt_factor = voigt_h::factors[I];
    UInt i = voigt_h::vec[I][0];
    UInt j = voigt_h::vec[I][1];

    voigt_strain(I) = voigt_factor * (grad_u(i, j) + grad_u(j, i)) / 2.;
  }

  // compute stresses in Voigt notation
  voigt_stress.mul<false>(C_mat, voigt_strain);

  /// copy stresses back in full vectorised notation
  for(UInt I = 0; I < voigt_h::size; ++I) {
    UInt i = voigt_h::vec[I][0];
    UInt j = voigt_h::vec[I][1];
    sigma(i, j) = sigma(j, i) = voigt_stress(I);
  }
  
  ++C_it;
  
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::computeTangentModuli(const ElementType & el_type,
							  Array<Real> & tangent_matrix,
							  GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<Real>::const_matrix_iterator C_it = this->C(el_type, ghost_type).begin(voigt_h::size, 
									       voigt_h::size);

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  tangent.copy(*C_it);
  ++C_it;
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::advanceASR(const Matrix<Real> & prestrain) {

  std::vector<SolidMechanicsModelRVE *>::iterator RVE_it = RVEs.begin();
  std::vector<SolidMechanicsModelRVE *>::iterator RVE_end = RVEs.end();

  Array<Real>::matrix_iterator C_it =
    this->C(this->el_type).begin(voigt_h::size, voigt_h::size);
  Array<Real>::matrix_iterator gradu_it =
    this->gradu(this->el_type).begin(spatial_dimension, spatial_dimension);
  Array<Real>::matrix_iterator eigen_gradu_it =
    this->eigengradu(this->el_type).begin(spatial_dimension, spatial_dimension);

  for (; RVE_it != RVE_end; ++RVE_it, ++C_it, ++gradu_it, ++eigen_gradu_it) {
    /// apply boundary conditions based on the current macroscopic displ. gradient
    (*RVE_it)->applyBoundaryConditions(*gradu_it);

    /// advance the ASR in every RVE
    (*RVE_it)->advanceASR(prestrain);

    /// compute the new effective stiffness of the RVE
    (*RVE_it)->homogenizeStiffness(*C_it);

    /// compute the average eigen_grad_u
    (*RVE_it)->homogenizeEigenGradU(*eigen_gradu_it);
  }
}

INSTANTIATE_MATERIAL(MaterialFE2);


__END_AKANTU__
