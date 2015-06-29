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
#ifndef __AKANTU_SOLID_MECHANICS_MODEL_IGFEM_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_IGFEM_HH__
#include "solid_mechanics_model.hh"
#include "solid_mechanics_model_event_handler.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
struct SolidMechanicsModelIGFEMOptions : public SolidMechanicsModelOptions {
  SolidMechanicsModelIGFEMOptions(AnalysisMethod analysis_method = _static,
				  bool no_init_materials = false,
				  bool moving_interface = false) :
    SolidMechanicsModelOptions(analysis_method, no_init_materials),
    moving_interface(moving_interface) { }
  bool moving_interface;
};

extern const SolidMechanicsModelIGFEMOptions default_solid_mechanics_model_igfem_options;

/* -------------------------------------------------------------------------- */
/* Solid Mechanics Model for IGFEM analysis                                   */
/* -------------------------------------------------------------------------- */
class SolidMechanicsModelIGFEM : public SolidMechanicsModel,
				 public SolidMechanicsModelEventHandler{
public:
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  typedef FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_igfem> MyFEEngineIGFEMType;

  SolidMechanicsModelIGFEM(Mesh & mesh,
			   UInt spatial_dimension = _all_dimensions,
			   const ID & id = "solid_mechanics_model_igfem",
			   const MemoryID & memory_id = 0);

  virtual ~SolidMechanicsModelIGFEM();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize the cohesive model
  virtual void initFull(const ModelOptions & options = default_solid_mechanics_model_igfem_options);

  /// initialize the model
  virtual void initModel();

  /// initialize igfem material
  virtual void initMaterials();

  ///allocate all vectors
  virtual void initArrays();

private:

  /// compute the real values of displacement, force, etc. on the enriched nodes
  void computeRealNodalFields();

  /* ------------------------------------------------------------------------ */
  /* Mesh Event Handler inherited members                                     */
  /* ------------------------------------------------------------------------ */

protected:

  virtual void onNodesAdded  (const Array<UInt> & nodes_list,
			      const NewNodesEvent & event) {};
  virtual void onElementsAdded  (const Array<Element> & nodes_list,
				 const NewElementsEvent & event) {};

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:

  virtual void onDump();

  virtual void addDumpGroupFieldToDumper(const std::string & dumper_name,
					 const std::string & field_id,
					 const std::string & group_name,
					 const ElementKind & element_kind,
					 bool padding_flag);
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// real displacements array
  Array <Real> *real_displacement;
  /// real forces array
  Array <Real> *real_force;
  /// real residuals array
  Array <Real> *real_residual;
  /// 
  bool moving_interface; 
};

/* -------------------------------------------------------------------------- */
/* IGFEMMaterialSelector                                                      */
/* -------------------------------------------------------------------------- */

class DefaultMaterialIGFEMSelector : public DefaultMaterialSelector {
public:
  DefaultMaterialIGFEMSelector(const SolidMechanicsModelIGFEM & model) :
    DefaultMaterialSelector(model.getMaterialByElement()),
    fallback_value_igfem(0) { }

  virtual UInt operator()(const Element & element) {
    if(Mesh::getKind(element.type) == _ek_igfem) 
      return fallback_value_igfem;
    else
      return DefaultMaterialSelector::operator()(element);
  }

  void setIGFEMFallback(UInt f) { this->fallback_value_igfem = f; }

private:
  UInt fallback_value_igfem;  
};
__END_AKANTU__

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_IGFEM_HH__ */


