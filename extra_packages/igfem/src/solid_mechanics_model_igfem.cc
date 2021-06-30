/**
 * @file   solid_mechanics_model_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  solid mechanics model for IGFEM analysis
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_igfem.hh"
#include "dumpable_inline_impl.hh"
#include "group_manager_inline_impl.hh"
#include "igfem_helper.hh"
#include "material_igfem.hh"
#ifdef AKANTU_USE_IOHELPER
#include "dumper_igfem_element_partition.hh"
#include "dumper_igfem_elemental_field.hh"
#include "dumper_igfem_material_internal_field.hh"
#include "dumper_material_padders.hh"
#include "dumper_paraview.hh"
#endif

/* -------------------------------------------------------------------------- */

namespace akantu {

const SolidMechanicsModelIGFEMOptions
    default_solid_mechanics_model_igfem_options(_static, false);

SolidMechanicsModelIGFEM::SolidMechanicsModelIGFEM(Mesh & mesh, UInt dim,
                                                   const ID & id)
    : SolidMechanicsModel(mesh, dim, id), IGFEMEnrichment(mesh),
      global_ids_updater(NULL) {
  AKANTU_DEBUG_IN();

  delete material_selector;
  material_selector = new DefaultMaterialIGFEMSelector(*this);

  this->registerEventHandler(*this);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("igfem elements", id);
  this->mesh.addDumpMeshToDumper("igfem elements", mesh, spatial_dimension,
                                 _not_ghost, _ek_igfem);
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModelIGFEM::~SolidMechanicsModelIGFEM() {
  AKANTU_DEBUG_IN();
  if (global_ids_updater)
    delete global_ids_updater;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::initFull(const ModelOptions & options) {
  AKANTU_DEBUG_IN();

  /// intialize the IGFEM enrichment
  this->initialize();

  SolidMechanicsModel::initFull(options);

  // set the initial condition to 0
  real_force->clear();
  real_displacement->clear();
  real_residual->clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Allocate all the needed vectors. By  default their are not necessarily set to
 * 0
 *
 */
void SolidMechanicsModelIGFEM::initArrays() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = mesh.getNbNodes();
  std::stringstream sstr_rdisp;
  sstr_rdisp << id << ":real_displacement";
  std::stringstream sstr_rforc;
  sstr_rforc << id << ":real_force";
  std::stringstream sstr_rresi;
  sstr_rresi << id << ":real_residual";

  real_displacement = &(alloc<Real>(sstr_rdisp.str(), nb_nodes,
                                    spatial_dimension, REAL_INIT_VALUE));
  real_force = &(alloc<Real>(sstr_rforc.str(), nb_nodes, spatial_dimension,
                             REAL_INIT_VALUE));
  real_residual = &(alloc<Real>(sstr_rresi.str(), nb_nodes, spatial_dimension,
                                REAL_INIT_VALUE));

  SolidMechanicsModel::initArrays();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::initParallel(MeshPartition * partition,
                                            DataAccessor * data_accessor) {
  SolidMechanicsModel::initParallel(partition, data_accessor);
  this->intersector_sphere.setDistributedSynchronizer(synch_parallel);
  if (mesh.isDistributed())
    global_ids_updater = new GlobalIdsUpdater(mesh, synch_parallel);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::initMaterials() {
  AKANTU_DEBUG_IN();

  // make sure the material are instantiated
  if (!are_materials_instantiated)
    instantiateMaterials();

  /// find the first igfem material
  UInt igfem_index = 0;

  while ((dynamic_cast<MaterialIGFEM *>(materials[igfem_index]) == NULL) &&
         igfem_index <= materials.size())
    ++igfem_index;

  AKANTU_DEBUG_ASSERT(igfem_index != materials.size(),
                      "No igfem materials in the material input file");

  DefaultMaterialIGFEMSelector * igfem_mat_selector =
      dynamic_cast<DefaultMaterialIGFEMSelector *>(material_selector);
  if (igfem_mat_selector != NULL)
    igfem_mat_selector->setIGFEMFallback(igfem_index);

  SolidMechanicsModel::initMaterials();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Initialize the model, basically pre-compute the shapes, shapes derivatives
 * and jacobian
 *
 */
void SolidMechanicsModelIGFEM::initModel() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::initModel();

  registerFEEngineObject<MyFEEngineIGFEMType>("IGFEMFEEngine", mesh,
                                              spatial_dimension);
  /// insert the two feengines associated with the model in the map
  this->fe_engines_per_kind[_ek_regular] = &(this->getFEEngine());
  this->fe_engines_per_kind[_ek_igfem] = &(this->getFEEngine("IGFEMFEEngine"));

  /// add the igfem type connectivities

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType type_ghost = *gt;

    Mesh::type_iterator it = mesh.firstType(spatial_dimension, type_ghost);
    Mesh::type_iterator last = mesh.lastType(spatial_dimension, type_ghost);

    for (; it != last; ++it) {
      const Array<UInt> & connectivity = mesh.getConnectivity(*it, type_ghost);
      if (connectivity.getSize() != 0) {
        ElementType type = *it;
        Vector<ElementType> types_igfem = FEEngine::getIGFEMElementTypes(type);
        for (UInt i = 0; i < types_igfem.size(); ++i)
          mesh.addConnectivityType(types_igfem(i), type_ghost);
      }
    }
  }

  getFEEngine("IGFEMFEEngine").initShapeFunctions(_not_ghost);
  getFEEngine("IGFEMFEEngine").initShapeFunctions(_ghost);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::onElementsAdded(const Array<Element> & elements,
                                               const NewElementsEvent & event) {
  AKANTU_DEBUG_IN();

  const NewIGFEMElementsEvent * igfem_event =
      dynamic_cast<const NewIGFEMElementsEvent *>(&event);
  /// insert the new and old elements in the map
  if (igfem_event != NULL) {
    this->element_map.zero();
    const Array<Element> & old_elements = igfem_event->getOldElementsList();
    for (UInt e = 0; e < elements.getSize(); ++e) {
      this->element_map[elements(e)] = old_elements(e);
    }
  }

  /// update shape functions
  getFEEngine("IGFEMFEEngine").initShapeFunctions(_not_ghost);
  getFEEngine("IGFEMFEEngine").initShapeFunctions(_ghost);

  SolidMechanicsModel::onElementsAdded(elements, event);
  this->reassignMaterial();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::onElementsRemoved(
    const Array<Element> & element_list,
    const ElementTypeMapArray<UInt> & new_numbering,
    const RemovedElementsEvent & event) {

  this->getFEEngine("IGFEMFEEngine").initShapeFunctions(_not_ghost);
  this->getFEEngine("IGFEMFEEngine").initShapeFunctions(_ghost);
  SolidMechanicsModel::onElementsRemoved(element_list, new_numbering, event);
  if (synch_parallel)
    synch_parallel->computeAllBufferSizes(*this);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::onNodesAdded(const Array<UInt> & nodes_list,
                                            const NewNodesEvent & event) {
  AKANTU_DEBUG_IN();

  const NewIGFEMNodesEvent * igfem_event =
      dynamic_cast<const NewIGFEMNodesEvent *>(&event);
  // update the node type
  if (igfem_event != NULL) {
    intersector_sphere.updateNodeType(
        nodes_list, igfem_event->getNewNodePerElem(),
        igfem_event->getElementType(), igfem_event->getGhostType());
  }

  UInt nb_nodes = mesh.getNbNodes();

  if (real_displacement)
    real_displacement->resize(nb_nodes);
  if (real_force)
    real_force->resize(nb_nodes);
  if (real_residual)
    real_residual->resize(nb_nodes);

  if (mesh.isDistributed())
    mesh.getGlobalNodesIds().resize(mesh.getNbNodes());

  if (displacement)
    displacement->resize(nb_nodes);
  if (mass)
    mass->resize(nb_nodes);
  if (velocity)
    velocity->resize(nb_nodes);
  if (acceleration)
    acceleration->resize(nb_nodes);
  if (force)
    force->resize(nb_nodes);
  if (residual)
    residual->resize(nb_nodes);
  if (blocked_dofs)
    blocked_dofs->resize(nb_nodes);

  if (previous_displacement)
    previous_displacement->resize(nb_nodes);
  if (increment_acceleration)
    increment_acceleration->resize(nb_nodes);
  if (increment)
    increment->resize(nb_nodes);

  if (current_position)
    current_position->resize(nb_nodes);

  std::vector<Material *>::iterator mat_it;
  for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->onNodesAdded(nodes_list, event);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::onNodesRemoved(const Array<UInt> & nodes_list,
                                              const Array<UInt> & new_numbering,
                                              const RemovedNodesEvent & event) {
  if (real_displacement)
    mesh.removeNodesFromArray(*real_displacement, new_numbering);
  if (real_force)
    mesh.removeNodesFromArray(*real_force, new_numbering);
  if (real_residual)
    mesh.removeNodesFromArray(*real_residual, new_numbering);

  // communicate global connectivity for slave nodes
  if (global_ids_updater)
    global_ids_updater->updateGlobalIDs(
        mesh.getNbNodes() - intersector_sphere.getNbStandardNodes());

  SolidMechanicsModel::onNodesRemoved(nodes_list, new_numbering, event);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::addDumpGroupFieldToDumper(
    const std::string & dumper_name, const std::string & field_id,
    const std::string & group_name, ElementKind element_kind,
    bool padding_flag) {
  AKANTU_DEBUG_IN();

  ElementKind _element_kind = element_kind;
  if (dumper_name == "igfem elements") {
    _element_kind = _ek_igfem;
  }

  SolidMechanicsModel::addDumpGroupFieldToDumper(
      dumper_name, field_id, group_name, _element_kind, padding_flag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::onDump() {
  this->computeValuesOnEnrichedNodes();
  this->flattenAllRegisteredInternals(_ek_igfem);
  SolidMechanicsModel::onDump();
}

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER

dumpers::Field * SolidMechanicsModelIGFEM::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag, const UInt & spatial_dimension,
    ElementKind kind) {

  dumpers::Field * field = NULL;

  if (kind != _ek_igfem)
    field = SolidMechanicsModel::createElementalField(
        field_name, group_name, padding_flag, spatial_dimension, kind);

  else {

    if (field_name == "partitions")
      field =
          mesh.createElementalField<UInt, dumpers::IGFEMElementPartitionField>(
              mesh.getConnectivities(), group_name, spatial_dimension, kind);
    else if (field_name == "material_index")
      field =
          mesh.createElementalField<UInt, Vector, dumpers::IGFEMElementalField>(
              material_index, group_name, spatial_dimension, kind);
    else {
      // this copy of field_name is used to compute derivated data such as
      // strain and von mises stress that are based on grad_u and stress
      std::string field_name_copy(field_name);

      if (field_name == "strain" || field_name == "Green strain" ||
          field_name == "principal strain" ||
          field_name == "principal Green strain")
        field_name_copy = "grad_u";
      else if (field_name == "Von Mises stress")
        field_name_copy = "stress";

      bool is_internal = this->isInternal(field_name_copy, kind);

      if (is_internal) {
        ElementTypeMap<UInt> nb_data_per_elem =
            this->getInternalDataPerElem(field_name_copy, kind);
        ElementTypeMapArray<Real> & internal_flat =
            this->flattenInternal(field_name_copy, kind);
        field =
            mesh.createElementalField<Real, dumpers::IGFEMInternalMaterialField>(
                internal_flat, group_name, spatial_dimension, kind,
                nb_data_per_elem);
        if (field_name == "strain") {
          dumpers::ComputeStrain<false> * foo =
              new dumpers::ComputeStrain<false>(*this);
          field = dumpers::FieldComputeProxy::createFieldCompute(field, *foo);
        } else if (field_name == "Von Mises stress") {
          dumpers::ComputeVonMisesStress * foo =
              new dumpers::ComputeVonMisesStress(*this);
          field = dumpers::FieldComputeProxy::createFieldCompute(field, *foo);
        } else if (field_name == "Green strain") {
          dumpers::ComputeStrain<true> * foo =
              new dumpers::ComputeStrain<true>(*this);
          field = dumpers::FieldComputeProxy::createFieldCompute(field, *foo);
        } else if (field_name == "principal strain") {
          dumpers::ComputePrincipalStrain<false> * foo =
              new dumpers::ComputePrincipalStrain<false>(*this);
          field = dumpers::FieldComputeProxy::createFieldCompute(field, *foo);
        } else if (field_name == "principal Green strain") {
          dumpers::ComputePrincipalStrain<true> * foo =
              new dumpers::ComputePrincipalStrain<true>(*this);
          field = dumpers::FieldComputeProxy::createFieldCompute(field, *foo);
        }

        /// treat the paddings
        if (padding_flag) {
          if (field_name == "stress") {
            if (spatial_dimension == 2) {
              dumpers::StressPadder<2> * foo =
                  new dumpers::StressPadder<2>(*this);
              field =
                  dumpers::FieldComputeProxy::createFieldCompute(field, *foo);
            }
          } else if (field_name == "strain" || field_name == "Green strain") {
            if (spatial_dimension == 2) {
              dumpers::StrainPadder<2> * foo =
                  new dumpers::StrainPadder<2>(*this);
              field =
                  dumpers::FieldComputeProxy::createFieldCompute(field, *foo);
            }
          }
        }
        // homogenize the field
        dumpers::ComputeFunctorInterface * foo =
            dumpers::HomogenizerProxy::createHomogenizer(*field);

        field = dumpers::FieldComputeProxy::createFieldCompute(field, *foo);
      }
    }
  }
  //  }
  return field;
}

/* -------------------------------------------------------------------------- */

dumpers::Field *
SolidMechanicsModelIGFEM::createNodalFieldReal(const std::string & field_name,
                                               const std::string & group_name,
                                               bool padding_flag) {

  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["real_displacement"] = real_displacement;

  dumpers::Field * field = NULL;
  if (padding_flag)
    field = mesh.createNodalField(real_nodal_fields[field_name], group_name, 3);
  else
    field = mesh.createNodalField(real_nodal_fields[field_name], group_name);

  if (field == NULL)
    return SolidMechanicsModel::createNodalFieldReal(field_name, group_name,
                                                     padding_flag);

  return field;
}

#else
/* -------------------------------------------------------------------------- */
dumpers::Field * SolidMechanicsModelIGFEM::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag, const UInt & spatial_dimension,
    ElementKind kind) {
  return NULL;
}

/* -------------------------------------------------------------------------- */

dumpers::Field *
SolidMechanicsModelIGFEM::createNodalFieldReal(const std::string & field_name,
                                               const std::string & group_name,
                                               bool padding_flag) {
  return NULL;
}

#endif
/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::computeValuesOnEnrichedNodes() {

  for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
    for (UInt s = 0; s < spatial_dimension; ++s)
      (*real_displacement)(n, s) = (*displacement)(n, s);
  }

  Element element;
  Vector<Real> real_coords(spatial_dimension);
  Vector<Real> interpolated(spatial_dimension);
  Array<Real>::const_vector_iterator r_displ_it =
      this->real_displacement->begin(spatial_dimension);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    element.ghost_type = *gt;
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, *gt, _ek_igfem);
    Mesh::type_iterator last = mesh.lastType(spatial_dimension, *gt, _ek_igfem);
    for (; it != last; ++it) {
      element.type = *it;
      UInt nb_element = mesh.getNbElement(*it, *gt);
      if (!nb_element)
        continue;
      UInt * elem_val = mesh.getConnectivity(*it, *gt).storage();
      UInt nb_nodes_per_element = mesh.getNbNodesPerElement(*it);
      Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);
      Matrix<Real> displ_val(spatial_dimension, nb_nodes_per_element);

      UInt nb_enriched_nodes = IGFEMHelper::getNbEnrichedNodes(*it);
      UInt nb_parent_nodes = IGFEMHelper::getNbParentNodes(*it);
      for (UInt el = 0; el < nb_element; ++el) {
        element.element = el;
        /// get the node coordinates of the element
        mesh.extractNodalValuesFromElement(
            mesh.getNodes(), nodes_coord.storage(),
            elem_val + el * nb_nodes_per_element, nb_nodes_per_element,
            spatial_dimension);

        /// get the displacement values at the nodes of the element
        mesh.extractNodalValuesFromElement(
            *(this->displacement), displ_val.storage(),
            elem_val + el * nb_nodes_per_element, nb_nodes_per_element,
            spatial_dimension);

        for (UInt i = 0; i < nb_enriched_nodes; ++i) {
          /// coordinates of enriched node
          real_coords = nodes_coord(nb_parent_nodes + i);
          /// global index of the enriched node
          UInt idx = elem_val[el * nb_nodes_per_element + nb_parent_nodes + i];
          /// compute the real displacement value
          this->getFEEngine("IGFEMFEEngine")
              .interpolate(real_coords, displ_val, interpolated, element);
          r_displ_it[idx] = interpolated;
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::transferInternalValues(
    const ID & internal, std::vector<Element> & new_elements,
    Array<Real> & added_quads, Array<Real> & internal_values) {

  /// @todo sort the new elements by their corresponding old element type and
  /// old material!!!

  /// get the number of elements for which iternals need to be transfered
  UInt nb_new_elements = new_elements.size();
  UInt nb_new_quads = added_quads.getSize() / nb_new_elements;

  Array<Real>::const_matrix_iterator quad_coords =
      added_quads.begin_reinterpret(this->spatial_dimension, nb_new_quads,
                                    nb_new_elements);
  UInt nb_internal_component = internal_values.getNbComponent();
  Array<Real>::matrix_iterator internal_val = internal_values.begin_reinterpret(
      nb_internal_component, nb_new_quads, nb_new_elements);
  Vector<Real> default_values(nb_internal_component, 0.);

  for (UInt e = 0; e < nb_new_elements; ++e, ++quad_coords, ++internal_val) {
    Element new_element = new_elements[e];
    Element old_element = this->element_map[new_element];
    UInt mat_idx = (this->material_index(
        old_element.type, old_element.ghost_type))(old_element.element);
    Material & old_material = *(this->materials[mat_idx]);
    old_material.extrapolateInternal(internal, old_element, *quad_coords,
                                     *internal_val);
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::applyEigenGradU(
    const Matrix<Real> & prescribed_eigen_grad_u, const ID & material_name,
    const GhostType ghost_type) {

  AKANTU_DEBUG_ASSERT(prescribed_eigen_grad_u.size() ==
                          spatial_dimension * spatial_dimension,
                      "The prescribed grad_u is not of the good size");
  std::vector<Material *>::iterator mat_it;
  for (mat_it = this->materials.begin(); mat_it != this->materials.end();
       ++mat_it) {
    MaterialIGFEM * mat_igfem = dynamic_cast<MaterialIGFEM *>(*mat_it);
    if (mat_igfem != NULL)
      mat_igfem->applyEigenGradU(prescribed_eigen_grad_u, material_name,
                                 ghost_type);
    else if ((*mat_it)->getName() == material_name)
      (*mat_it)->applyEigenGradU(prescribed_eigen_grad_u, ghost_type);
  }
}

} // namespace akantu
