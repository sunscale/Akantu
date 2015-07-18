/**
 * @file   solid_mechanics_model_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  solid mechanics model for IGFEM analysis
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_igfem.hh"
#include "dumpable_inline_impl.hh"
#include "material_igfem.hh"
#include "group_manager_inline_impl.cc"
#include "igfem_helper.hh"
#ifdef AKANTU_USE_IOHELPER
#  include "dumper_paraview.hh"
#  include "dumper_igfem_material_internal_field.hh"
#  include "dumper_igfem_element_partition.hh"
#  include "dumper_igfem_elemental_field.hh"
#  include "dumper_material_padders.hh"
#endif

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

const SolidMechanicsModelIGFEMOptions default_solid_mechanics_model_igfem_options(_static,
										  false,
										  false);

SolidMechanicsModelIGFEM::SolidMechanicsModelIGFEM(Mesh & mesh,
						   UInt dim,
						   const ID & id,
						   const MemoryID & memory_id) :
  SolidMechanicsModel(mesh, dim, id, memory_id),
  IGFEMEnrichment(mesh) {
  AKANTU_DEBUG_IN();

  delete material_selector;
  material_selector = new DefaultMaterialIGFEMSelector(*this);

  this->registerEventHandler(*this);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("igfem elements", id);
  this->mesh.addDumpMeshToDumper("igfem elements",
				 mesh, spatial_dimension, _not_ghost, _ek_igfem);
#endif


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModelIGFEM::~SolidMechanicsModelIGFEM() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::initFull(const ModelOptions & options) {
  AKANTU_DEBUG_IN();


  /// intialize the IGFEM enrichment
  this->initialize();

  const SolidMechanicsModelIGFEMOptions & smmc_options =
    dynamic_cast<const SolidMechanicsModelIGFEMOptions &>(options);

  this->moving_interface = smmc_options.moving_interface;

  SolidMechanicsModel::initFull(options);

  // set the initial condition to 0
  real_force->clear();
  real_displacement->clear();
  igfem_nodes->copy(this->mesh.getNodes());

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
  std::stringstream sstr_rdisp; sstr_rdisp << id << ":real_displacement";
  std::stringstream sstr_rforc; sstr_rforc << id << ":real_force";
  std::stringstream sstr_rresi; sstr_rresi << id << ":real_residual";
  std::stringstream sstr_inodes; sstr_inodes << id << ":igfem_nodes";

  real_displacement = &(alloc<Real>(sstr_rdisp.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  real_force        = &(alloc<Real>(sstr_rforc.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  real_residual     = &(alloc<Real>(sstr_rresi.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  igfem_nodes     = &(alloc<Real>(sstr_inodes.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));

  SolidMechanicsModel::initArrays(); 

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::initMaterials() {
  AKANTU_DEBUG_IN();

  // make sure the material are instantiated
  if(!are_materials_instantiated) instantiateMaterials();

  /// find the first igfem material
  UInt igfem_index = 0;

  while ((dynamic_cast<MaterialIGFEM *>(materials[igfem_index]) == NULL)
	 && igfem_index <= materials.size())
    ++igfem_index;

  AKANTU_DEBUG_ASSERT(igfem_index != materials.size(),
		      "No igfem materials in the material input file");

  (dynamic_cast<DefaultMaterialIGFEMSelector *>(material_selector))->setIGFEMFallback(igfem_index);

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

  registerFEEngineObject<MyFEEngineIGFEMType>("IGFEMFEEngine", mesh, spatial_dimension);

  /// add the igfem type connectivities

  ElementType type = _not_defined;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType type_ghost = *gt;

    Mesh::type_iterator it   = mesh.firstType(spatial_dimension, type_ghost);
    Mesh::type_iterator last = mesh.lastType(spatial_dimension, type_ghost);

    for (; it != last; ++it) {
      const Array<UInt> & connectivity = mesh.getConnectivity(*it, type_ghost);
      if (connectivity.getSize() != 0) {
	type = *it;
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
void SolidMechanicsModelIGFEM::onElementsAdded(const Array<Element> & doubled_elements,
					       const NewElementsEvent & event) {
  AKANTU_DEBUG_IN();

  UInt nb_new_elements = doubled_elements.getSize();
  Array<Element> element_list(nb_new_elements);

  /// update shape functions
  getFEEngine("IGFEMFEEngine").initShapeFunctions(_not_ghost);
  getFEEngine("IGFEMFEEngine").initShapeFunctions(_ghost);

  for (UInt e = 0; e < nb_new_elements; ++e)
    element_list(e) = doubled_elements(e, 0);

  SolidMechanicsModel::onElementsAdded(element_list, event);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::onElementsRemoved(const Array<Element> & element_list,
						 const ElementTypeMapArray<UInt> & new_numbering,
						 const RemovedElementsEvent & event) {

  this->getFEEngine("IGFEMFEEngine").initShapeFunctions(_not_ghost);
  this->getFEEngine("IGFEMFEEngine").initShapeFunctions(_ghost);
  SolidMechanicsModel::onElementsRemoved(element_list, new_numbering, event);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::onNodesAdded(const Array<UInt> & nodes_list,
					    const NewNodesEvent & event) {
  AKANTU_DEBUG_IN();

  UInt nb_new_nodes = nodes_list.getSize();
  UInt nb_nodes = mesh.getNbNodes();

  if(real_displacement) real_displacement->resize(nb_nodes);
  if(real_force) real_force->resize(nb_nodes);
  if(real_residual) real_residual->resize(nb_nodes);
  if(igfem_nodes) {igfem_nodes->resize(nb_nodes);
    for (UInt n = 0; n < nb_new_nodes; ++n) {
      UInt new_node = nodes_list(n);
      for (UInt dim = 0; dim < this->spatial_dimension; ++dim)
	(*igfem_nodes)(new_node, dim) = 0.;
    }
  }
  SolidMechanicsModel::onNodesAdded(nodes_list, event);  

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::onNodesRemoved(const Array<UInt> & element_list,
					      const Array<UInt> & new_numbering,
					      const RemovedNodesEvent & event) {
  if(real_displacement) mesh.removeNodesFromArray(*real_displacement, new_numbering);
  if(real_force        ) mesh.removeNodesFromArray(*real_force        , new_numbering);
  if(real_residual    ) mesh.removeNodesFromArray(*real_residual    , new_numbering);
  if(igfem_nodes    ) mesh.removeNodesFromArray(*igfem_nodes    , new_numbering);

  SolidMechanicsModel::onNodesRemoved(element_list, new_numbering, event);

}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::addDumpGroupFieldToDumper(const std::string & dumper_name,
							 const std::string & field_id,
							 const std::string & group_name,
							 const ElementKind & element_kind,
							 bool padding_flag) {
  AKANTU_DEBUG_IN();

  ElementKind _element_kind = element_kind;
  if (dumper_name == "igfem elements") {
    _element_kind = _ek_igfem;
  }

  SolidMechanicsModel::addDumpGroupFieldToDumper(dumper_name, 
						 field_id,
						 group_name,
						 _element_kind,
						 padding_flag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::onDump(){
  this->computeValuesOnEnrichedNodes();
  this->flattenAllRegisteredInternals(_ek_igfem);
  SolidMechanicsModel::onDump();
}

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER

dumper::Field * SolidMechanicsModelIGFEM
::createElementalField(const std::string & field_name,
		       const std::string & group_name,
		       bool padding_flag,
		       const UInt & spatial_dimension,
		       const ElementKind & kind) {
 
  dumper::Field * field = NULL;

  if (kind != _ek_igfem)
    field = SolidMechanicsModel::createElementalField(field_name, group_name, padding_flag, spatial_dimension, kind);

  else {
 

    if(field_name == "partitions")
      field = mesh.createElementalField<UInt, dumper::IGFEMElementPartitionField>(mesh.getConnectivities(),group_name,spatial_dimension, kind);
    else if(field_name == "material_index")
      field = mesh.createElementalField<UInt, Vector, dumper::IGFEMElementalField >(material_index,group_name,spatial_dimension,kind);
    else {
      // this copy of field_name is used to compute derivated data such as
      // strain and von mises stress that are based on grad_u and stress
      std::string field_name_copy(field_name);

      if (field_name == "strain"
	  || field_name == "Green strain"
	  || field_name == "principal strain"
	  || field_name == "principal Green strain")
	field_name_copy = "grad_u";
      else if (field_name == "Von Mises stress")
	field_name_copy = "stress";

      bool is_internal = this->isInternal(field_name_copy,kind);

      if (is_internal) {
	ElementTypeMap<UInt> nb_data_per_elem = this->getInternalDataPerElem(field_name_copy,kind);
	ElementTypeMapArray<Real> & internal_flat = this->flattenInternal(field_name_copy,kind);
	field = mesh.createElementalField<Real, dumper::IGFEMInternalMaterialField>(internal_flat,
										    group_name,
										    spatial_dimension,kind,nb_data_per_elem);
	if (field_name == "strain"){
		dumper::ComputeStrain<false> * foo = new dumper::ComputeStrain<false>(*this);
		field = dumper::FieldComputeProxy::createFieldCompute(field,*foo);
	} else if (field_name == "Von Mises stress") {
		dumper::ComputeVonMisesStress * foo = new dumper::ComputeVonMisesStress(*this);
		field = dumper::FieldComputeProxy::createFieldCompute(field,*foo);
	} else if (field_name == "Green strain") {
		dumper::ComputeStrain<true> * foo = new dumper::ComputeStrain<true>(*this);
		field = dumper::FieldComputeProxy::createFieldCompute(field,*foo);
	} else if (field_name == "principal strain") {
		dumper::ComputePrincipalStrain<false> * foo = new dumper::ComputePrincipalStrain<false>(*this);
		field = dumper::FieldComputeProxy::createFieldCompute(field,*foo);
	} else if (field_name == "principal Green strain") {
		dumper::ComputePrincipalStrain<true> * foo = new dumper::ComputePrincipalStrain<true>(*this);
		field = dumper::FieldComputeProxy::createFieldCompute(field,*foo);
	}

	/// treat the paddings
	if (padding_flag){
		if (field_name == "stress"){
		  if (spatial_dimension == 2) {
		    dumper::StressPadder<2> * foo = new dumper::StressPadder<2>(*this);
		    field = dumper::FieldComputeProxy::createFieldCompute(field,*foo);
		  }
		} else if (field_name == "strain" || field_name == "Green strain"){
		  if (spatial_dimension == 2) {
		    dumper::StrainPadder<2> * foo = new dumper::StrainPadder<2>(*this);
		    field = dumper::FieldComputeProxy::createFieldCompute(field,*foo);
		  }
		}
	}
	// homogenize the field
	dumper::ComputeFunctorInterface * foo =
	  dumper::HomogenizerProxy::createHomogenizer(*field);

	field = dumper::FieldComputeProxy::createFieldCompute(field,*foo);
      }
    }
  }
    //  }
  return field;
}

/* -------------------------------------------------------------------------- */

dumper::Field * SolidMechanicsModelIGFEM::createNodalFieldReal(const std::string & field_name,
							  const std::string & group_name,
							  bool padding_flag) {

  std::map<std::string,Array<Real>* > real_nodal_fields;
  real_nodal_fields["real_displacement"             ] = real_displacement;

  dumper::Field * field = NULL;
  if (padding_flag)
    field = mesh.createNodalField(real_nodal_fields[field_name],group_name, 3);
  else
    field = mesh.createNodalField(real_nodal_fields[field_name],group_name);
  
  if (field == NULL)
    return SolidMechanicsModel::createNodalFieldReal(field_name, group_name, padding_flag);
  
  return field;
}

#endif

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelIGFEM::computeValuesOnEnrichedNodes() {

  for (UInt n = 0; n < mesh.getNbNodes(); ++ n) {
    for (UInt s = 0; s < spatial_dimension; ++s)
      (*real_displacement)(n,s) = (*displacement)(n,s);
  }


  Element element;
  Vector<Real> real_coords(spatial_dimension);
  Vector<Real> interpolated(spatial_dimension);
  Array<Real>::const_vector_iterator r_displ_it = this->real_displacement->begin(spatial_dimension);
  
  for (ghost_type_t::iterator gt = ghost_type_t::begin(); gt != ghost_type_t::end(); ++gt) {
    element.ghost_type = *gt;
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, *gt, _ek_igfem);
    Mesh::type_iterator last  = mesh.lastType(spatial_dimension, *gt, _ek_igfem);
    for(;it != last; ++it) {
      element.type = *it;
      UInt nb_element = mesh.getNbElement(*it, *gt);
      if (!nb_element) continue;    
      UInt * elem_val = mesh.getConnectivity(*it, *gt).storage();
      UInt nb_nodes_per_element = mesh.getNbNodesPerElement(*it);
      Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);
      Matrix<Real> displ_val(spatial_dimension, nb_nodes_per_element);
 
      UInt nb_enriched_nodes = IGFEMHelper::getNbEnrichedNodes(*it);
      UInt nb_parent_nodes = IGFEMHelper::getNbParentNodes(*it);
      for (UInt el = 0; el < nb_element; ++el) {
	element.element = el;
	/// get the node coordinates of the element
	mesh.extractNodalValuesFromElement(mesh.getNodes(),
					   nodes_coord.storage(),
					   elem_val + el * nb_nodes_per_element,
					   nb_nodes_per_element,
					   spatial_dimension);

	/// get the displacement values at the nodes of the element
	mesh.extractNodalValuesFromElement(*(this->displacement),
					   displ_val.storage(),
					   elem_val + el * nb_nodes_per_element,
					   nb_nodes_per_element,
					   spatial_dimension);

	for (UInt i = 0; i < nb_enriched_nodes; ++i) {
	  /// coordinates of enriched node
	  real_coords = nodes_coord(nb_parent_nodes + i);
	  /// global index of the enriched node
	  UInt idx = elem_val[el * nb_nodes_per_element + nb_parent_nodes + i];
	  /// compute the real displacement value
	  this->getFEEngine("IGFEMFEEngine").interpolate(real_coords, displ_val,
							 interpolated, element);
	  r_displ_it[idx] = interpolated;	 
	}
      }
    }
  }
}


__END_AKANTU__
