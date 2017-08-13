/**
 * @file   material_reinforcement.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Wed Mar 25 2015
 * @date last modification: Tue Dec 08 2015
 *
 * @brief  Reinforcement material
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

#include "aka_common.hh"
#include "aka_voigthelper.hh"
#include "material_reinforcement.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
MaterialReinforcement<dim>::MaterialReinforcement(SolidMechanicsModel & model,
                                                  UInt spatial_dimension,
                                                  const Mesh & mesh,
                                                  FEEngine & fe_engine,
                                                  const ID & id)
    : Material(model, dim, mesh, fe_engine,
               id), // /!\ dim, not spatial_dimension !
      model(NULL),
      stress_embedded("stress_embedded", *this, 1, fe_engine,
                      this->element_filter),
      gradu_embedded("gradu_embedded", *this, 1, fe_engine,
                     this->element_filter),
      directing_cosines("directing_cosines", *this, 1, fe_engine,
                        this->element_filter),
      pre_stress("pre_stress", *this, 1, fe_engine, this->element_filter),
      area(1.0), shape_derivatives() {
  AKANTU_DEBUG_IN();
  this->initialize(model);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialReinforcement<dim>::initialize(SolidMechanicsModel & a_model) {
  AKANTU_DEBUG_IN();
  this->model = dynamic_cast<EmbeddedInterfaceModel *>(&a_model);
  AKANTU_DEBUG_ASSERT(this->model != NULL,
                      "MaterialReinforcement needs an EmbeddedInterfaceModel");

  this->registerParam("area", area, _pat_parsable | _pat_modifiable,
                      "Reinforcement cross-sectional area");
  this->registerParam("pre_stress", pre_stress, _pat_parsable | _pat_modifiable,
                      "Uniform pre-stress");

  this->unregisterInternal(this->stress);

  // Fool the AvgHomogenizingFunctor
  // stress.initialize(dim * dim);

  // Reallocate the element filter
  this->element_filter.free();
  this->model->getInterfaceMesh().initElementTypeMapArray(this->element_filter,
                                                          1, 1, false, _ek_regular);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <UInt dim> MaterialReinforcement<dim>::~MaterialReinforcement() {
  AKANTU_DEBUG_IN();

  ElementTypeMap<ElementTypeMapArray<Real> *>::type_iterator
    it = shape_derivatives.firstType(),
    end = shape_derivatives.lastType();

  for (; it != end; ++it) {
    delete shape_derivatives(*it, _not_ghost);
    delete shape_derivatives(*it, _ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <UInt dim> void MaterialReinforcement<dim>::initMaterial() {
  Material::initMaterial();

  stress_embedded.initialize(dim * dim);
  gradu_embedded.initialize(dim * dim);
  pre_stress.initialize(1);

  /// We initialise the stuff that is not going to change during the simulation
  this->allocBackgroundShapeDerivatives();
  this->initBackgroundShapeDerivatives();
  this->initDirectingCosines();
}

/* -------------------------------------------------------------------------- */

/**
 * Background shape derivatives need to be stored per background element
 * types but also per embedded element type, which is why they are stored
 * in an ElementTypeMap<ElementTypeMapArray<Real> *>. The outer ElementTypeMap
 * refers to the embedded types, and the inner refers to the background types.
 */
template <UInt dim>
void MaterialReinforcement<dim>::allocBackgroundShapeDerivatives() {
  AKANTU_DEBUG_IN();

  Mesh & interface_mesh = model->getInterfaceMesh();
  Mesh & mesh = model->getMesh();

  ghost_type_t::iterator int_ghost_it = ghost_type_t::begin();

  // Loop over interface ghosts
  for (; int_ghost_it != ghost_type_t::end(); ++int_ghost_it) {
    Mesh::type_iterator interface_type_it =
        interface_mesh.firstType(1, *int_ghost_it);
    Mesh::type_iterator interface_type_end =
        interface_mesh.lastType(1, *int_ghost_it);

    // Loop over interface types
    for (; interface_type_it != interface_type_end; ++interface_type_it) {
      Mesh::type_iterator background_type_it =
          mesh.firstType(dim, *int_ghost_it);
      Mesh::type_iterator background_type_end =
          mesh.lastType(dim, *int_ghost_it);

      // Loop over background types
      for (; background_type_it != background_type_end ; ++background_type_it) {
        const ElementType & int_type = *interface_type_it;
        const ElementType & back_type = *background_type_it;
        const GhostType & int_ghost = *int_ghost_it;

        std::string shaped_id = "embedded_shape_derivatives";

        if (int_ghost == _ghost) shaped_id += ":ghost";

        ElementTypeMapArray<Real> * shaped_etma = new ElementTypeMapArray<Real>(shaped_id, this->name);

        UInt nb_points = Mesh::getNbNodesPerElement(back_type);
        UInt nb_quad_points = model->getFEEngine("EmbeddedInterfaceFEEngine").getNbIntegrationPoints(int_type);
        UInt nb_elements = element_filter(int_type, int_ghost).size();

        // Alloc the background ElementTypeMapArray
        shaped_etma->alloc(nb_elements * nb_quad_points,
                           dim * nb_points,
                           back_type);

        // Insert the background ElementTypeMapArray in the interface
        // ElementTypeMap
        shape_derivatives(shaped_etma, int_type, int_ghost);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <UInt dim>
void MaterialReinforcement<dim>::initBackgroundShapeDerivatives() {
  AKANTU_DEBUG_IN();

  Mesh & mesh = model->getMesh();

  Mesh::type_iterator type_it = mesh.firstType(dim, _not_ghost);
  Mesh::type_iterator type_end = mesh.lastType(dim, _not_ghost);

  for (; type_it != type_end; ++type_it) {
    computeBackgroundShapeDerivatives(*type_it, _not_ghost);
    // computeBackgroundShapeDerivatives(*type_it, _ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <UInt dim> void MaterialReinforcement<dim>::initDirectingCosines() {
  AKANTU_DEBUG_IN();

  Mesh & mesh = model->getInterfaceMesh();

  Mesh::type_iterator type_it = mesh.firstType(1, _not_ghost);
  Mesh::type_iterator type_end = mesh.lastType(1, _not_ghost);

  const UInt voigt_size = getTangentStiffnessVoigtSize(dim);
  directing_cosines.initialize(voigt_size * voigt_size);

  for (; type_it != type_end; ++type_it) {
    computeDirectingCosines(*type_it, _not_ghost);
    // computeDirectingCosines(*type_it, _ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <UInt dim>
void MaterialReinforcement<dim>::assembleStiffnessMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & interface_mesh = model->getInterfaceMesh();

  Mesh::type_iterator type_it = interface_mesh.firstType(1, _not_ghost);
  Mesh::type_iterator type_end = interface_mesh.lastType(1, _not_ghost);

  for (; type_it != type_end; ++type_it) {
    assembleStiffnessMatrix(*type_it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <UInt dim>
void MaterialReinforcement<dim>::assembleInternalForces(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & interface_mesh = model->getInterfaceMesh();

  Mesh::type_iterator type_it = interface_mesh.firstType(1, _not_ghost);
  Mesh::type_iterator type_end = interface_mesh.lastType(1, _not_ghost);

  for (; type_it != type_end; ++type_it) {
    this->assembleInternalForces(*type_it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialReinforcement<dim>::computeGradU(const ElementType & type,
                                              GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<UInt> & elem_filter = element_filter(type, ghost_type);
  UInt nb_element = elem_filter.size();
  UInt nb_quad_points = model->getFEEngine("EmbeddedInterfaceFEEngine")
                            .getNbIntegrationPoints(type);

  Array<Real> & gradu_vec = gradu_embedded(type, ghost_type);

  Mesh::type_iterator back_it = model->getMesh().firstType(dim, ghost_type);
  Mesh::type_iterator back_end = model->getMesh().lastType(dim, ghost_type);

  for (; back_it != back_end; ++back_it) {
    UInt nodes_per_background_e = Mesh::getNbNodesPerElement(*back_it);

    Array<Real> & shapesd =
        shape_derivatives(type, ghost_type)->operator()(*back_it, ghost_type);

    Array<UInt> * background_filter =
        new Array<UInt>(nb_element, 1, "background_filter");
    filterInterfaceBackgroundElements(*background_filter, *back_it, type,
                                      ghost_type);

    Array<Real> * disp_per_element = new Array<Real>(0, dim * nodes_per_background_e, "disp_elem");
    FEEngine::extractNodalToElementField(model->getMesh(),
                                         model->getDisplacement(),
                                         *disp_per_element,
                                         *back_it, ghost_type, *background_filter);

    Array<Real>::matrix_iterator disp_it =
        disp_per_element->begin(dim, nodes_per_background_e);
    Array<Real>::matrix_iterator disp_end =
        disp_per_element->end(dim, nodes_per_background_e);

    Array<Real>::matrix_iterator shapes_it =
        shapesd.begin(dim, nodes_per_background_e);
    Array<Real>::matrix_iterator grad_u_it = gradu_vec.begin(dim, dim);

    for (; disp_it != disp_end; ++disp_it) {
      for (UInt i = 0; i < nb_quad_points; i++, ++shapes_it, ++grad_u_it) {
        Matrix<Real> & B = *shapes_it;
        Matrix<Real> & du = *grad_u_it;
        Matrix<Real> & u = *disp_it;

        du.mul<false, true>(u, B);
      }
    }

    delete background_filter;
    delete disp_per_element;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialReinforcement<dim>::computeAllStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it = model->getInterfaceMesh().firstType();
  Mesh::type_iterator last_type = model->getInterfaceMesh().lastType();

  for (; it != last_type; ++it) {
    computeGradU(*it, ghost_type);
    computeStress(*it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialReinforcement<dim>::assembleInternalForces(
    const ElementType & type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = model->getMesh();

  Mesh::type_iterator type_it = mesh.firstType(dim, ghost_type);
  Mesh::type_iterator type_end = mesh.lastType(dim, ghost_type);

  for (; type_it != type_end; ++type_it) {
    assembleInternalForcesInterface(type, *type_it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

/**
 * Computes and assemble the residual. Residual in reinforcement is computed as
 * :
 * \f[
 * \vec{r} = A_s \int_S{\mathbf{B}^T\mathbf{C}^T \vec{\sigma_s}\,\mathrm{d}s}
 * \f]
 */
template <UInt dim>
void MaterialReinforcement<dim>::assembleInternalForcesInterface(
    const ElementType & interface_type, const ElementType & background_type,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt voigt_size = getTangentStiffnessVoigtSize(dim);

  Array<Real> & residual = const_cast<Array<Real> &>(model->getResidual());

  FEEngine & interface_engine = model->getFEEngine("EmbeddedInterfaceFEEngine");
  FEEngine & background_engine = model->getFEEngine();

  Array<UInt> & elem_filter = element_filter(interface_type, ghost_type);

  UInt nodes_per_background_e = Mesh::getNbNodesPerElement(background_type);
  UInt nb_quadrature_points =
      interface_engine.getNbIntegrationPoints(interface_type, ghost_type);
  UInt nb_element = elem_filter.size();

  UInt back_dof = dim * nodes_per_background_e;

  Array<Real> & shapesd = shape_derivatives(interface_type, ghost_type)
                              ->
                              operator()(background_type, ghost_type);

  Array<Real> * integrant = new Array<Real>(nb_quadrature_points * nb_element,
                                            back_dof,
                                            "integrant");

  Array<Real>::vector_iterator integrant_it = integrant->begin(back_dof);
  Array<Real>::vector_iterator integrant_end = integrant->end(back_dof);

  Array<Real>::matrix_iterator B_it =
      shapesd.begin(dim, nodes_per_background_e);
  Array<Real>::matrix_iterator C_it =
      directing_cosines(interface_type, ghost_type)
          .begin(voigt_size, voigt_size);
  Array<Real>::matrix_iterator sigma_it =
      stress_embedded(interface_type, ghost_type).begin(dim, dim);

  Vector<Real> sigma(voigt_size);
  Matrix<Real> Bvoigt(voigt_size, back_dof);
  Vector<Real> Ct_sigma(voigt_size);

  for (; integrant_it != integrant_end ; ++integrant_it,
         ++B_it,
         ++C_it,
         ++sigma_it) {
    VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(*B_it, Bvoigt, nodes_per_background_e);
    Matrix<Real> & C = *C_it;
    Vector<Real> & BtCt_sigma = *integrant_it;

    stressTensorToVoigtVector(*sigma_it, sigma);

    Ct_sigma.mul<true>(C, sigma);
    BtCt_sigma.mul<true>(Bvoigt, Ct_sigma);
    BtCt_sigma *= area;
  }

  Array<Real> * residual_interface =
      new Array<Real>(nb_element, back_dof, "residual_interface");
  interface_engine.integrate(*integrant, *residual_interface, back_dof,
                             interface_type, ghost_type, elem_filter);

  delete integrant;

  Array<UInt> * background_filter =
      new Array<UInt>(nb_element, 1, "background_filter");

  filterInterfaceBackgroundElements(*background_filter, background_type,
                                    interface_type, ghost_type);
  background_engine.assembleArray(
      *residual_interface, residual,
      model->getDOFSynchronizer().getLocalDOFEquationNumbers(), dim,
      background_type, ghost_type, *background_filter, -1.0);
  delete residual_interface;
  delete background_filter;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialReinforcement<dim>::filterInterfaceBackgroundElements(
    Array<UInt> & filter, const ElementType & type,
    const ElementType & interface_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  filter.resize(0);
  filter.clear();

  Array<Element> & elements =
      model->getInterfaceAssociatedElements(interface_type, ghost_type);
  Array<UInt> & elem_filter = element_filter(interface_type, ghost_type);

  Array<UInt>::scalar_iterator filter_it = elem_filter.begin(),
                               filter_end = elem_filter.end();

  for (; filter_it != filter_end; ++filter_it) {
    Element & elem = elements(*filter_it);
    if (elem.type == type)
      filter.push_back(elem.element);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <UInt dim>
void MaterialReinforcement<dim>::computeDirectingCosines(
    const ElementType & type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & interface_mesh = this->model->getInterfaceMesh();

  const UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  const UInt steel_dof = dim * nb_nodes_per_element;
  const UInt voigt_size = getTangentStiffnessVoigtSize(dim);
  const UInt nb_quad_points = model->getFEEngine("EmbeddedInterfaceFEEngine")
                                  .getNbIntegrationPoints(type, ghost_type);

  Array<Real> node_coordinates(this->element_filter(type, ghost_type).size(),
                               steel_dof);

  this->model->getFEEngine().template extractNodalToElementField<Real>(interface_mesh,
                                                                       interface_mesh.getNodes(),
                                                                       node_coordinates,
                                                                       type,
                                                                       ghost_type,
                                                                       this->element_filter(type, ghost_type));

  Array<Real>::matrix_iterator directing_cosines_it =
      directing_cosines(type, ghost_type).begin(voigt_size, voigt_size);

  Array<Real>::matrix_iterator node_coordinates_it =
      node_coordinates.begin(dim, nb_nodes_per_element);
  Array<Real>::matrix_iterator node_coordinates_end =
      node_coordinates.end(dim, nb_nodes_per_element);

  for (; node_coordinates_it != node_coordinates_end; ++node_coordinates_it) {
    for (UInt i = 0; i < nb_quad_points; i++, ++directing_cosines_it) {
      Matrix<Real> & nodes = *node_coordinates_it;
      Matrix<Real> & cosines = *directing_cosines_it;

      computeDirectingCosinesOnQuad(nodes, cosines);
    }
  }

  // Mauro: the directing_cosines internal is defined on the quadrature points of each element

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <UInt dim>
void MaterialReinforcement<dim>::assembleStiffnessMatrix(
    const ElementType & type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = model->getMesh();

  Mesh::type_iterator type_it = mesh.firstType(dim, ghost_type);
  Mesh::type_iterator type_end = mesh.lastType(dim, ghost_type);

  for (; type_it != type_end; ++type_it) {
    assembleStiffnessMatrixInterface(type, *type_it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Computes the reinforcement stiffness matrix (Gomes & Awruch, 2001)
 * \f[
 * \mathbf{K}_e = \sum_{i=1}^R{A_i\int_{S_i}{\mathbf{B}^T
 * \mathbf{C}_i^T \mathbf{D}_{s, i} \mathbf{C}_i \mathbf{B}\,\mathrm{d}s}}
 * \f]
 */
template <UInt dim>
void MaterialReinforcement<dim>::assembleStiffnessMatrixInterface(
    const ElementType & interface_type, const ElementType & background_type,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt voigt_size = getTangentStiffnessVoigtSize(dim);

  SparseMatrix & K = const_cast<SparseMatrix &>(model->getStiffnessMatrix());

  FEEngine & background_engine = model->getFEEngine();
  FEEngine & interface_engine = model->getFEEngine("EmbeddedInterfaceFEEngine");

  Array<UInt> & elem_filter = element_filter(interface_type, ghost_type);
  Array<Real> & grad_u = gradu_embedded(interface_type, ghost_type);

  UInt nb_element = elem_filter.size();
  UInt nodes_per_background_e = Mesh::getNbNodesPerElement(background_type);
  UInt nb_quadrature_points =
      interface_engine.getNbIntegrationPoints(interface_type, ghost_type);

  UInt back_dof = dim * nodes_per_background_e;

  UInt integrant_size = back_dof;

  grad_u.resize(nb_quadrature_points * nb_element);

  Array<Real> * tangent_moduli = new Array<Real>(
      nb_element * nb_quadrature_points, 1, "interface_tangent_moduli");
  tangent_moduli->clear();
  computeTangentModuli(interface_type, *tangent_moduli, ghost_type);

  Array<Real> & shapesd = shape_derivatives(interface_type, ghost_type)
                              ->
                              operator()(background_type, ghost_type);

  Array<Real> * integrant = new Array<Real>(nb_element * nb_quadrature_points,
                                            integrant_size * integrant_size,
                                            "B^t*C^t*D*C*B");
  integrant->clear();

  /// Temporary matrices for integrant product
  Matrix<Real> Bvoigt(voigt_size, back_dof);
  Matrix<Real> DC(voigt_size, voigt_size);
  Matrix<Real> DCB(voigt_size, back_dof);
  Matrix<Real> CtDCB(voigt_size, back_dof);

  Array<Real>::scalar_iterator D_it = tangent_moduli->begin();
  Array<Real>::scalar_iterator D_end = tangent_moduli->end();

  Array<Real>::matrix_iterator C_it =
      directing_cosines(interface_type, ghost_type)
          .begin(voigt_size, voigt_size);
  Array<Real>::matrix_iterator B_it =
      shapesd.begin(dim, nodes_per_background_e);
  Array<Real>::matrix_iterator integrant_it =
      integrant->begin(integrant_size, integrant_size);

  for (; D_it != D_end; ++D_it, ++C_it, ++B_it, ++integrant_it) {
    Real & D = *D_it;
    Matrix<Real> & C = *C_it;
    Matrix<Real> & B = *B_it;
    Matrix<Real> & BtCtDCB = *integrant_it;

    VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(B, Bvoigt,
                                                       nodes_per_background_e);

    DC.clear();
    DC(0, 0) = D * area;
    DC *= C;
    DCB.mul<false, false>(DC, Bvoigt);
    CtDCB.mul<true, false>(C, DCB);
    BtCtDCB.mul<true, false>(Bvoigt, CtDCB);
  }

  delete tangent_moduli;

  Array<Real> * K_interface = new Array<Real>(
      nb_element, integrant_size * integrant_size, "K_interface");
  interface_engine.integrate(*integrant, *K_interface,
                             integrant_size * integrant_size, interface_type,
                             ghost_type, elem_filter);

  delete integrant;

  // Mauro: Here K_interface contains the local stiffness matrices,
  // directing_cosines contains the information about the orientation
  // of the reinforcements, any rotation of the local stiffness matrix
  // can be done here

  Array<UInt> * background_filter = new Array<UInt>(nb_element, 1, "background_filter");

  filterInterfaceBackgroundElements(*background_filter, background_type,
                                    interface_type, ghost_type);

  background_engine.assembleMatrix(*K_interface, K, dim, background_type,
                                   ghost_type, *background_filter);

  delete K_interface;
  delete background_filter;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

/// In this function, type and ghost_type refer to background elements
template <UInt dim>
void MaterialReinforcement<dim>::computeBackgroundShapeDerivatives(
    const ElementType & type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & interface_mesh = model->getInterfaceMesh();

  FEEngine & engine = model->getFEEngine();
  FEEngine & interface_engine = model->getFEEngine("EmbeddedInterfaceFEEngine");

  Mesh::type_iterator interface_type = interface_mesh.firstType();
  Mesh::type_iterator interface_last = interface_mesh.lastType();

  for (; interface_type != interface_last; ++interface_type) {
    Array<UInt> & filter = element_filter(*interface_type, ghost_type);

    const UInt nb_elements = filter.size();
    const UInt nb_nodes = Mesh::getNbNodesPerElement(type);
    const UInt nb_quad_per_element =
        interface_engine.getNbIntegrationPoints(*interface_type);

    Array<Real> quad_pos(nb_quad_per_element * nb_elements, dim,
                         "interface_quad_points");
    quad_pos.resize(nb_quad_per_element * nb_elements);
    interface_engine.interpolateOnIntegrationPoints(
        interface_mesh.getNodes(), quad_pos, dim, *interface_type, ghost_type,
        filter);

    Array<Real> & background_shapesd =
        shape_derivatives(*interface_type, ghost_type)
            ->
            operator()(type, ghost_type);
    background_shapesd.clear();

    Array<UInt> * background_elements = new Array<UInt>(
        nb_elements, 1, "computeBackgroundShapeDerivatives:background_filter");
    filterInterfaceBackgroundElements(*background_elements, type,
                                      *interface_type, ghost_type);

    Array<UInt>::scalar_iterator back_it = background_elements->begin(),
                                 back_end = background_elements->end();

    Array<Real>::matrix_iterator shapesd_it =
        background_shapesd.begin(dim, nb_nodes);

    Array<Real>::vector_iterator quad_pos_it = quad_pos.begin(dim);

    for (; back_it != back_end; ++back_it) {
      for (UInt i = 0; i < nb_quad_per_element;
           i++, ++shapesd_it, ++quad_pos_it)
        engine.computeShapeDerivatives(*quad_pos_it, *back_it, type,
                                       *shapesd_it, ghost_type);
    }

    delete background_elements;
  }

  AKANTU_DEBUG_OUT();
}

template <UInt dim> Real MaterialReinforcement<dim>::getEnergy(std::string id) {
  AKANTU_DEBUG_IN();
  if (id == "potential") {
    Real epot = 0.;

    computePotentialEnergyByElements();

    Mesh::type_iterator it = element_filter.firstType(spatial_dimension),
                        end = element_filter.lastType(spatial_dimension);

    for (; it != end; ++it) {
      FEEngine & interface_engine =
          model->getFEEngine("EmbeddedInterfaceFEEngine");
      epot += interface_engine.integrate(potential_energy(*it, _not_ghost), *it,
                                         _not_ghost,
                                         element_filter(*it, _not_ghost));
      epot *= area;
    }

    return epot;
  }

  AKANTU_DEBUG_OUT();
  return 0;
}

// /* --------------------------------------------------------------------------
// */
// template<UInt dim>
// ElementTypeMap<UInt> MaterialReinforcement<dim>::getInternalDataPerElem(const
// ID & field_name,
//                                                                         const
//                                                                         ElementKind
//                                                                         &
//                                                                         kind,
//                                                                         const
//                                                                         ID &
//                                                                         fe_engine_id)
//                                                                         const
//                                                                         {
//   return Material::getInternalDataPerElem(field_name, kind,
//   "EmbeddedInterfaceFEEngine");
// }

// /* --------------------------------------------------------------------------
// */
// // Author is Guillaume Anciaux, see material.cc
// template<UInt dim>
// void MaterialReinforcement<dim>::flattenInternal(const std::string &
// field_id,
//                                                  ElementTypeMapArray<Real> &
//                                                  internal_flat,
//                                                  const GhostType ghost_type,
//                                                  ElementKind element_kind)
//                                                  const {
//   AKANTU_DEBUG_IN();

//   if (field_id == "stress_embedded" || field_id == "inelastic_strain") {
//     Material::flattenInternalIntern(field_id,
//                                     internal_flat,
//                                     1,
//                                     ghost_type,
//                                     _ek_not_defined,
//                                     &(this->element_filter),
//                                     &(this->model->getInterfaceMesh()));
//   }

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(MaterialReinforcement);

/* -------------------------------------------------------------------------- */

} // akantu
