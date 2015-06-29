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
#ifdef AKANTU_USE_IOHELPER
#  include "dumper_paraview.hh"
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
  SolidMechanicsModel(mesh, dim, id, memory_id) {
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

  const SolidMechanicsModelIGFEMOptions & smmc_options =
    dynamic_cast<const SolidMechanicsModelIGFEMOptions &>(options);

  this->moving_interface = smmc_options.moving_interface;

  SolidMechanicsModel::initFull(options);

  // set the initial condition to 0
  real_force->clear();
  real_displacement->clear();

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

  real_displacement = &(alloc<Real>(sstr_rdisp.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  real_force        = &(alloc<Real>(sstr_rforc.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  real_residual     = &(alloc<Real>(sstr_rresi.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));

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

  getFEEngine("IGFEMFEEngine").initShapeFunctions(_not_ghost);
  getFEEngine("IGFEMFEEngine").initShapeFunctions(_ghost);


  AKANTU_DEBUG_OUT();
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
  /// this->flattenAllRegisteredInternals(_ek_cohesive);
  SolidMechanicsModel::onDump();
}

/* -------------------------------------------------------------------------- */


__END_AKANTU__
