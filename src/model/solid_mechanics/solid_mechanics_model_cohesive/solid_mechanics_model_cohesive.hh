/**
 * @file   solid_mechanics_model_cohesive.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Mon Feb 05 2018
 *
 * @brief  Solid mechanics model for cohesive elements
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "cohesive_element_inserter.hh"
#include "material_selector_cohesive.hh"
#include "random_internal_field.hh" // included to have the specialization of
                                    // ParameterTyped::operator Real()
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH_
#define AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH_

/* -------------------------------------------------------------------------- */
namespace akantu {
class FacetSynchronizer;
class FacetStressSynchronizer;
class ElementSynchronizer;
} // namespace akantu

namespace akantu {

/* -------------------------------------------------------------------------- */
struct FacetsCohesiveIntegrationOrderFunctor {
  template <ElementType type, ElementType cohesive_type =
                                  CohesiveFacetProperty<type>::cohesive_type>
  struct _helper {
    static constexpr int get() {
      return ElementClassProperty<cohesive_type>::polynomial_degree;
    }
  };

  template <ElementType type> struct _helper<type, _not_defined> {
    static constexpr int get() {
      return ElementClassProperty<type>::polynomial_degree;
    }
  };

  template <ElementType type> static inline constexpr int getOrder() {
    return _helper<type>::get();
  }
};

/* -------------------------------------------------------------------------- */
/* Solid Mechanics Model for Cohesive elements                                */
/* -------------------------------------------------------------------------- */
class SolidMechanicsModelCohesive : public SolidMechanicsModel,
                                    public SolidMechanicsModelEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  class NewCohesiveNodesEvent : public NewNodesEvent {
  public:
    AKANTU_GET_MACRO_NOT_CONST(OldNodesList, old_nodes, Array<UInt> &);
    AKANTU_GET_MACRO(OldNodesList, old_nodes, const Array<UInt> &);

  protected:
    Array<UInt> old_nodes;
  };

  using MyFEEngineCohesiveType =
      FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_cohesive>;
  using MyFEEngineFacetType =
      FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular,
                       FacetsCohesiveIntegrationOrderFunctor>;

  SolidMechanicsModelCohesive(Mesh & mesh, UInt dim = _all_dimensions,
                              const ID & id = "solid_mechanics_model_cohesive",
                              const MemoryID & memory_id = 0);

  ~SolidMechanicsModelCohesive() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// initialize the cohesive model
  void initFullImpl(const ModelOptions & options) override;

public:
  /// set the value of the time step
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

  /// assemble the residual for the explicit scheme
  void assembleInternalForces() override;

  /// function to perform a stress check on each facet and insert
  /// cohesive elements if needed (returns the number of new cohesive
  /// elements)
  UInt checkCohesiveStress();

  /// interpolate stress on facets
  void interpolateStress();

  /// update automatic insertion after a change in the element inserter
  void updateAutomaticInsertion();

  /// insert intrinsic cohesive elements
  void insertIntrinsicElements();

  // template <SolveConvergenceMethod cmethod, SolveConvergenceCriteria
  // criteria> bool solveStepCohesive(Real tolerance, Real & error, UInt
  // max_iteration = 100,
  //                        bool load_reduction = false,
  //                        Real tol_increase_factor = 1.0,
  //                        bool do_not_factorize = false);

protected:
  /// initialize stress interpolation
  void initStressInterpolation();

  /// initialize the model
  void initModel() override;

  /// initialize cohesive material
  void initMaterials() override;

  /// init facet filters for cohesive materials
  void initFacetFilter();

  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

private:
  /// insert cohesive elements along a given physical surface of the mesh
  void insertElementsFromMeshData(const std::string & physical_name);

  /// initialize completely the model for extrinsic elements
  void initAutomaticInsertion();

  /// compute facets' normals
  void computeNormals();

  /// resize facet stress
  void resizeFacetStress();

  /// init facets_check array
  void initFacetsCheck();

  /* ------------------------------------------------------------------------ */
  /* Mesh Event Handler inherited members                                     */
  /* ------------------------------------------------------------------------ */

protected:
  void onNodesAdded(const Array<UInt> & new_nodes,
                    const NewNodesEvent & event) override;
  void onElementsAdded(const Array<Element> & element_list,
                       const NewElementsEvent & event) override;

  /* ------------------------------------------------------------------------ */
  /* SolidMechanicsModelEventHandler inherited members                        */
  /* ------------------------------------------------------------------------ */
public:
  void afterSolveStep(bool converged = true) override;

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  void onDump() override;

  void addDumpGroupFieldToDumper(const std::string & dumper_name,
                                 const std::string & field_id,
                                 const std::string & group_name,
                                 ElementKind element_kind,
                                 bool padding_flag) override;

public:
  /// register the tags associated with the parallel synchronizer for
  /// cohesive elements
  // void initParallel(MeshPartition * partition,
  //                DataAccessor * data_accessor = NULL,
  //                bool extrinsic = false);

protected:
  //void synchronizeGhostFacetsConnectivity();

  void updateCohesiveSynchronizers();
  void updateFacetStressSynchronizer();

  friend class CohesiveElementInserter;

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  UInt getNbData(const Array<Element> & elements,
                 const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
                  const SynchronizationTag & tag) override;

protected:
  UInt getNbQuadsForFacetCheck(const Array<Element> & elements) const;

  template <typename T>
  void packFacetStressDataHelper(const ElementTypeMapArray<T> & data_to_pack,
                                 CommunicationBuffer & buffer,
                                 const Array<Element> & elements) const;

  template <typename T>
  void unpackFacetStressDataHelper(ElementTypeMapArray<T> & data_to_unpack,
                                   CommunicationBuffer & buffer,
                                   const Array<Element> & elements) const;

  template <typename T, bool pack_helper>
  void packUnpackFacetStressDataHelper(ElementTypeMapArray<T> & data_to_pack,
                                       CommunicationBuffer & buffer,
                                       const Array<Element> & element) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get facet mesh
  AKANTU_GET_MACRO(MeshFacets, mesh.getMeshFacets(), const Mesh &);

  /// get stress on facets vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(StressOnFacets, facet_stress, Real);

  /// get facet material
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(FacetMaterial, facet_material, UInt);

  /// get facet material
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(FacetMaterial, facet_material, UInt);

  /// get facet material
  AKANTU_GET_MACRO(FacetMaterial, facet_material,
                   const ElementTypeMapArray<UInt> &);

  /// @todo THIS HAS TO BE CHANGED
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Tangents, tangents, Real);

  /// get element inserter
  AKANTU_GET_MACRO_NOT_CONST(ElementInserter, *inserter,
                             CohesiveElementInserter &);

  /// get is_extrinsic boolean
  AKANTU_GET_MACRO(IsExtrinsic, is_extrinsic, bool);

  /// get cohesive elements synchronizer
  AKANTU_GET_MACRO_NOT_CONST(CohesiveSynchronizer, *cohesive_synchronizer,
                   ElementSynchronizer &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  friend class CohesiveMeshGlobalDataUpdater;

  /// @todo store tangents when normals are computed:
  ElementTypeMapArray<Real> tangents;

  /// stress on facets on the two sides by quadrature point
  ElementTypeMapArray<Real> facet_stress;

  /// material to use if a cohesive element is created on a facet
  ElementTypeMapArray<UInt> facet_material;

  bool is_extrinsic{false};

  /// cohesive element inserter
  std::unique_ptr<CohesiveElementInserter> inserter;

  /// facet stress synchronizer
  std::unique_ptr<ElementSynchronizer> facet_stress_synchronizer;

  /// cohesive elements synchronizer
  std::unique_ptr<ElementSynchronizer> cohesive_synchronizer;
};

} // namespace akantu

#include "solid_mechanics_model_cohesive_inline_impl.hh"

#endif /* AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH_ */
