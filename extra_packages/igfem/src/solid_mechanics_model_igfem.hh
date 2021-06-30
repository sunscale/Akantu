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
#ifndef AKANTU_SOLID_MECHANICS_MODEL_IGFEM_HH_
#define AKANTU_SOLID_MECHANICS_MODEL_IGFEM_HH_
#include "global_ids_updater.hh"
#include "igfem_enrichment.hh"
#include "solid_mechanics_model.hh"
#include "solid_mechanics_model_event_handler.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
struct SolidMechanicsModelIGFEMOptions : public SolidMechanicsModelOptions {
  SolidMechanicsModelIGFEMOptions(AnalysisMethod analysis_method = _static,
                                  bool no_init_materials = false)
      : SolidMechanicsModelOptions(analysis_method, no_init_materials) {}
};

extern const SolidMechanicsModelIGFEMOptions
    default_solid_mechanics_model_igfem_options;

/* -------------------------------------------------------------------------- */
/* Solid Mechanics Model for IGFEM analysis                                   */
/* -------------------------------------------------------------------------- */
class SolidMechanicsModelIGFEM : public SolidMechanicsModel,
                                 public SolidMechanicsModelEventHandler,
                                 public IGFEMEnrichment {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_igfem>
      MyFEEngineIGFEMType;
  typedef std::map<Element, Element> ElementMap;
  typedef std::map<ElementKind, FEEngine *> FEEnginesPerKindMap;

  SolidMechanicsModelIGFEM(Mesh & mesh,
                           UInt spatial_dimension = _all_dimensions,
                           const ID & id = "solid_mechanics_model_igfem");

  virtual ~SolidMechanicsModelIGFEM();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the cohesive model
  virtual void initFull(const ModelOptions & options =
                            default_solid_mechanics_model_igfem_options);

  /// initialize the model
  virtual void initModel();

  /// initialize igfem material
  virtual void initMaterials();

  /// register the tags associated with the parallel synchronizer
  virtual void initParallel(MeshPartition * partition,
                            DataAccessor * data_accessor = NULL);

  /// allocate all vectors
  virtual void initArrays();

  /// transfer internals from old to new elements
  void transferInternalValues(const ID & internal,
                              std::vector<Element> & new_elements,
                              Array<Real> & added_quads,
                              Array<Real> & internal_values);

  /// compute the barycenter for a sub-element
  inline void getSubElementBarycenter(UInt element, UInt sub_element,
                                      ElementType type,
                                      Vector<Real> & barycenter,
                                      GhostType ghost_type) const;

  /// apply a constant eigen_grad_u on all quadrature points of a given material
  virtual void applyEigenGradU(const Matrix<Real> & prescribed_eigen_grad_u,
                               const ID & material_name,
                               const GhostType ghost_type = _not_ghost);

private:
  /// compute the real values of displacement, force, etc. on the enriched nodes
  void computeValuesOnEnrichedNodes();

  /* ------------------------------------------------------------------------ */
  /* Mesh Event Handler inherited members                                     */
  /* ------------------------------------------------------------------------ */

protected:
  virtual void onNodesAdded(const Array<UInt> & nodes_list,
                            const NewNodesEvent & event);
  virtual void onNodesRemoved(const Array<UInt> & element_list,
                              const Array<UInt> & new_numbering,
                              const RemovedNodesEvent & event);
  virtual void onElementsAdded(const Array<Element> & nodes_list,
                               const NewElementsEvent & event);
  virtual void
  onElementsRemoved(const Array<Element> & element_list,
                    const ElementTypeMapArray<UInt> & new_numbering,
                    const RemovedElementsEvent & event);

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  virtual void onDump();

  virtual void addDumpGroupFieldToDumper(const std::string & dumper_name,
                                         const std::string & field_id,
                                         const std::string & group_name,
                                         ElementKind element_kind,
                                         bool padding_flag);

  virtual dumpers::Field * createElementalField(const std::string & field_name,
                                               const std::string & group_name,
                                               bool padding_flag,
                                               const UInt & spatial_dimension,
                                               ElementKind kind);

  virtual dumpers::Field * createNodalFieldReal(const std::string & field_name,
                                               const std::string & group_name,
                                               bool padding_flag);

  /* --------------------------------------------------------------------------
   */
  /* Accessors */
  /* --------------------------------------------------------------------------
   */
public:
  /// get the fe-engines per kind
  AKANTU_GET_MACRO(FEEnginesPerKind, fe_engines_per_kind,
                   const FEEnginesPerKindMap &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// real displacements array
  Array<Real> * real_displacement;
  /// real forces array
  Array<Real> * real_force;
  /// real residuals array
  Array<Real> * real_residual;
  /// map between and new elements (needed when the interface is moving)
  ElementMap element_map;
  /// global connectivity ids updater
  GlobalIdsUpdater * global_ids_updater;
  /// map between element kind and corresponding FEEngine object
  FEEnginesPerKindMap fe_engines_per_kind;
};

/* -------------------------------------------------------------------------- */
/* IGFEMMaterialSelector                                                      */
/* -------------------------------------------------------------------------- */

class DefaultMaterialIGFEMSelector : public DefaultMaterialSelector {
public:
  DefaultMaterialIGFEMSelector(const SolidMechanicsModelIGFEM & model)
      : DefaultMaterialSelector(model.getMaterialByElement()),
        fallback_value_igfem(0) {}

  virtual UInt operator()(const Element & element) {
    if (Mesh::getKind(element.type) == _ek_igfem)
      return fallback_value_igfem;
    else
      return DefaultMaterialSelector::operator()(element);
  }

  void setIGFEMFallback(UInt f) { this->fallback_value_igfem = f; }

protected:
  UInt fallback_value_igfem;
};

} // namespace akantu

#if defined(AKANTU_INCLUDE_INLINE_IMPL)
#include "solid_mechanics_model_igfem_inline_impl.hh"
#endif

#endif /* AKANTU_SOLID_MECHANICS_MODEL_IGFEM_HH_ */
