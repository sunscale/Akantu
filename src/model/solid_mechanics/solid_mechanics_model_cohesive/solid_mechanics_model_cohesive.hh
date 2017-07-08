/**
 * @file   solid_mechanics_model_cohesive.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Mon Dec 14 2015
 *
 * @brief  Solid mechanics model for cohesive elements
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#ifndef __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH__

#include "cohesive_element_inserter.hh"
#include "material_selector_cohesive.hh"
#include "solid_mechanics_model.hh"
#include "solid_mechanics_model_event_handler.hh"

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
#include "facet_stress_synchronizer.hh"
#include "facet_synchronizer.hh"
#endif
/* -------------------------------------------------------------------------- */

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
struct SolidMechanicsModelCohesiveOptions : public SolidMechanicsModelOptions {
  SolidMechanicsModelCohesiveOptions(
      AnalysisMethod analysis_method = _explicit_lumped_mass,
      bool extrinsic = false)
      : SolidMechanicsModelOptions(analysis_method), extrinsic(extrinsic) {}

  template <typename... pack>
  SolidMechanicsModelCohesiveOptions(use_named_args_t, pack &&... _pack)
      : SolidMechanicsModelOptions(
            OPTIONAL_NAMED_ARG(analysis_method, _explicit_lumped_mass),
            OPTIONAL_NAMED_ARG(is_extrinsic, false)) {}

  bool extrinsic{false};
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

  typedef FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_cohesive>
      MyFEEngineCohesiveType;
  typedef FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular,
                           FacetsCohesiveIntegrationOrderFunctor>
      MyFEEngineFacetType;

  SolidMechanicsModelCohesive(Mesh & mesh,
                              UInt spatial_dimension = _all_dimensions,
                              const ID & id = "solid_mechanics_model_cohesive",
                              const MemoryID & memory_id = 0);

  virtual ~SolidMechanicsModelCohesive();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// set the value of the time step
  void setTimeStep(Real time_step);

  /// assemble the residual for the explicit scheme
  virtual void assembleInternalForces();

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// function to perform a stress check on each facet and insert
  /// cohesive elements if needed (returns the number of new cohesive
  /// elements)
  UInt checkCohesiveStress();

  /// interpolate stress on facets
  void interpolateStress();

  /// initialize the cohesive model
  void
  initFull(const ModelOptions & options = SolidMechanicsModelCohesiveOptions());

  template <typename P, typename T, typename... pack>
  void initFull(named_argument::param_t<P, T &&> && first, pack &&... _pack) {
    this->initFull(
        SolidMechanicsModelCohesiveOptions{use_named_args, first, _pack...});
  }

  /// initialize the model
  void initModel();

  /// initialize cohesive material
  void initMaterials();

  /// init facet filters for cohesive materials
  void initFacetFilter();

  /// limit the cohesive element insertion to a given area
  void limitInsertion(BC::Axis axis, Real first_limit, Real second_limit);

  /// update automatic insertion after a change in the element inserter
  void updateAutomaticInsertion();

  /// insert intrinsic cohesive elements
  void insertIntrinsicElements();

  // synchronize facets physical data before insertion along physical surfaces
  void synchronizeInsertionData();

  // template <SolveConvergenceMethod cmethod, SolveConvergenceCriteria
  // criteria> bool solveStepCohesive(Real tolerance, Real & error, UInt
  // max_iteration = 100,
  //                        bool load_reduction = false,
  //                        Real tol_increase_factor = 1.0,
  //                        bool do_not_factorize = false);

  /// initialize stress interpolation
  void initStressInterpolation();

private:
  /// initialize cohesive material with intrinsic insertion (by default)
  void initIntrinsicCohesiveMaterials(UInt cohesive_index);

  ///  initialize cohesive material with intrinsic insertion (if physical
  ///  surfaces are precised)
  void initIntrinsicCohesiveMaterials(std::string cohesive_surfaces);

  /// insert cohesive elements along a given physical surface of the mesh
  void insertElementsFromMeshData(std::string physical_name);

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
  virtual void onNodesAdded(const Array<UInt> & nodes_list,
                            const NewNodesEvent & event);
  virtual void onElementsAdded(const Array<Element> & nodes_list,
                               const NewElementsEvent & event);

  /* ------------------------------------------------------------------------ */
  /* SolidMechanicsModelEventHandler inherited members                        */
  /* ------------------------------------------------------------------------ */
public:
  virtual void onEndSolveStep(const AnalysisMethod & method);

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

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// @todo store tangents when normals are computed:
  ElementTypeMapArray<Real> tangents;

  /// stress on facets on the two sides by quadrature point
  ElementTypeMapArray<Real> facet_stress;

  /// material to use if a cohesive element is created on a facet
  ElementTypeMapArray<UInt> facet_material;

  bool is_extrinsic;

  /// cohesive element inserter
  CohesiveElementInserter * inserter;

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
#include "solid_mechanics_model_cohesive_parallel.hh"
#endif
};

/* -------------------------------------------------------------------------- */
/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const SolidMechanicsModelCohesive & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "solid_mechanics_model_cohesive_inline_impl.cc"

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH__ */
