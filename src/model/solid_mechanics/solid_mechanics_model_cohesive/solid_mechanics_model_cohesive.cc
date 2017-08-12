/**
 * @file   solid_mechanics_model_cohesive.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Wed Jan 13 2016
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
#include "solid_mechanics_model_cohesive.hh"
#include "aka_iterators.hh"
#include "dumpable_inline_impl.hh"
#include "material_cohesive.hh"
#include "shape_cohesive.hh"
#ifdef AKANTU_USE_IOHELPER
#include "dumper_paraview.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */

namespace akantu {

const SolidMechanicsModelCohesiveOptions
    default_solid_mechanics_model_cohesive_options(_explicit_lumped_mass,
                                                   false);

/* -------------------------------------------------------------------------- */

SolidMechanicsModelCohesive::SolidMechanicsModelCohesive(
    Mesh & mesh, UInt dim, const ID & id, const MemoryID & memory_id)
    : SolidMechanicsModel(mesh, dim, id, memory_id), tangents("tangents", id),
      facet_stress("facet_stress", id), facet_material("facet_material", id) {
  AKANTU_DEBUG_IN();

  inserter = NULL;

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  facet_synchronizer = NULL;
  facet_stress_synchronizer = NULL;
  cohesive_element_synchronizer = NULL;
  global_connectivity = NULL;
#endif

  delete material_selector;
  material_selector = new DefaultMaterialCohesiveSelector(*this);

  this->registerEventHandler(*this, _ehp_solid_mechanics_model_cohesive);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("cohesive elements", id);
  this->mesh.addDumpMeshToDumper("cohesive elements", mesh,
                                 Model::spatial_dimension, _not_ghost,
                                 _ek_cohesive);
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModelCohesive::~SolidMechanicsModelCohesive() {
  AKANTU_DEBUG_IN();

  delete inserter;

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  delete cohesive_element_synchronizer;
  delete facet_synchronizer;
  delete facet_stress_synchronizer;
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::setTimeStep(Real time_step) {
  SolidMechanicsModel::setTimeStep(time_step);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.getDumper("cohesive elements").setTimeStep(time_step);
#endif
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initFull(const ModelOptions & options) {
  AKANTU_DEBUG_IN();

  const SolidMechanicsModelCohesiveOptions & smmc_options =
      dynamic_cast<const SolidMechanicsModelCohesiveOptions &>(options);

  this->is_extrinsic = smmc_options.extrinsic;

  if (!inserter)
    inserter = new CohesiveElementInserter(mesh, is_extrinsic,
                                           id + ":cohesive_element_inserter");

  SolidMechanicsModel::initFull(options);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initMaterials() {
  AKANTU_DEBUG_IN();

  // make sure the material are instantiated
  if (!are_materials_instantiated)
    instantiateMaterials();

  /// find the first cohesive material
  UInt cohesive_index = 0;

  while ((dynamic_cast<MaterialCohesive *>(materials[cohesive_index].get()) ==
          nullptr) &&
         cohesive_index <= materials.size())
    ++cohesive_index;

  AKANTU_DEBUG_ASSERT(cohesive_index != materials.size(),
                      "No cohesive materials in the material input file");

  material_selector->setFallback(cohesive_index);

  // set the facet information in the material in case of dynamic insertion
  if (is_extrinsic) {
    const Mesh & mesh_facets = inserter->getMeshFacets();
    facet_material.initialize(mesh_facets, _spatial_dimension =
                                               Model::spatial_dimension - 1);
    // mesh_facets.initElementTypeMapArray(facet_material, 1,
    //                                     spatial_dimension - 1);

    Element element;
    for (auto ghost_type : ghost_types) {
      element.ghost_type = ghost_type;
      for (auto & type :
           mesh_facets.elementTypes(Model::spatial_dimension - 1, ghost_type)) {
        element.type = type;
        Array<UInt> & f_material = facet_material(type, ghost_type);
        UInt nb_element = mesh_facets.getNbElement(type, ghost_type);
        f_material.resize(nb_element);
        f_material.set(cohesive_index);
        for (UInt el = 0; el < nb_element; ++el) {
          element.element = el;
          UInt mat_index = (*material_selector)(element);
          f_material(el) = mat_index;
          MaterialCohesive & mat =
              dynamic_cast<MaterialCohesive &>(*materials[mat_index]);
          mat.addFacet(element);
        }
      }
    }
    SolidMechanicsModel::initMaterials();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
    if (facet_synchronizer != NULL)
      inserter->initParallel(facet_synchronizer, cohesive_element_synchronizer);
//      inserter->initParallel(facet_synchronizer, synch_parallel);
#endif
    initAutomaticInsertion();
  } else {
    // TODO think of something a bit mor consistant than just coding the first
    // thing that comes in Fabian's head....
    typedef ParserSection::const_section_iterator const_section_iterator;
    std::pair<const_section_iterator, const_section_iterator> sub_sections =
        this->parser->getSubSections(_st_mesh);

    if (sub_sections.first != sub_sections.second) {
      std::string cohesive_surfaces =
          sub_sections.first->getParameter("cohesive_surfaces");
      this->initIntrinsicCohesiveMaterials(cohesive_surfaces);
    } else {
      this->initIntrinsicCohesiveMaterials(cohesive_index);
    }
  }

  AKANTU_DEBUG_OUT();
} // namespace akantu

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initIntrinsicCohesiveMaterials(
    std::string cohesive_surfaces) {

  AKANTU_DEBUG_IN();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  if (facet_synchronizer != NULL)
    inserter->initParallel(facet_synchronizer, cohesive_element_synchronizer);
//    inserter->initParallel(facet_synchronizer, synch_parallel);
#endif
  std::istringstream split(cohesive_surfaces);
  std::string physname;
  while (std::getline(split, physname, ',')) {
    AKANTU_DEBUG_INFO(
        "Pre-inserting cohesive elements along facets from physical surface: "
        << physname);
    insertElementsFromMeshData(physname);
  }

  synchronizeInsertionData();

  SolidMechanicsModel::initMaterials();

  if (is_default_material_selector)
    delete material_selector;
  material_selector = new MeshDataMaterialCohesiveSelector(*this);

  // UInt nb_new_elements =
  inserter->insertElements();
  // if (nb_new_elements > 0) {
  //   this->reinitializeSolver();
  // }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::synchronizeInsertionData() {

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  if (facet_synchronizer != NULL) {
    facet_synchronizer->asynchronousSynchronize(*inserter, _gst_ce_groups);
    facet_synchronizer->waitEndSynchronize(*inserter, _gst_ce_groups);
  }
#endif
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initIntrinsicCohesiveMaterials(
    UInt cohesive_index) {

  AKANTU_DEBUG_IN();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    Mesh::type_iterator first =
        mesh.firstType(Model::spatial_dimension, *gt, _ek_cohesive);
    Mesh::type_iterator last =
        mesh.lastType(Model::spatial_dimension, *gt, _ek_cohesive);

    for (; first != last; ++first) {
      Array<UInt> & mat_indexes = this->material_index(*first, *gt);
      Array<UInt> & mat_loc_num = this->material_local_numbering(*first, *gt);
      mat_indexes.set(cohesive_index);
      mat_loc_num.clear();
    }
  }
#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  if (facet_synchronizer != NULL)
    inserter->initParallel(facet_synchronizer, cohesive_element_synchronizer);
//    inserter->initParallel(facet_synchronizer, synch_parallel);
#endif

  SolidMechanicsModel::initMaterials();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Initialize the model,basically it  pre-compute the shapes, shapes derivatives
 * and jacobian
 *
 */
void SolidMechanicsModelCohesive::initModel() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::initModel();

  registerFEEngineObject<MyFEEngineCohesiveType>("CohesiveFEEngine", mesh,
                                                 Model::spatial_dimension);

  /// add cohesive type connectivity
  ElementType type = _not_defined;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType type_ghost = *gt;

    Mesh::type_iterator it =
        mesh.firstType(Model::spatial_dimension, type_ghost);
    Mesh::type_iterator last =
        mesh.lastType(Model::spatial_dimension, type_ghost);

    for (; it != last; ++it) {
      const Array<UInt> & connectivity = mesh.getConnectivity(*it, type_ghost);
      if (connectivity.getSize() != 0) {
        type = *it;
        ElementType type_facet = Mesh::getFacetType(type);
        ElementType type_cohesive =
            FEEngine::getCohesiveElementType(type_facet);
        mesh.addConnectivityType(type_cohesive, type_ghost);
      }
    }
  }

  AKANTU_DEBUG_ASSERT(type != _not_defined, "No elements in the mesh");

  getFEEngine("CohesiveFEEngine").initShapeFunctions(_not_ghost);
  getFEEngine("CohesiveFEEngine").initShapeFunctions(_ghost);

  registerFEEngineObject<MyFEEngineFacetType>(
      "FacetsFEEngine", mesh.getMeshFacets(), Model::spatial_dimension - 1);

  if (is_extrinsic) {
    getFEEngine("FacetsFEEngine").initShapeFunctions(_not_ghost);
    getFEEngine("FacetsFEEngine").initShapeFunctions(_ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::limitInsertion(BC::Axis axis,
                                                 Real first_limit,
                                                 Real second_limit) {
  AKANTU_DEBUG_IN();
  inserter->setLimit(axis, first_limit, second_limit);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::insertIntrinsicElements() {
  AKANTU_DEBUG_IN();
  // UInt nb_new_elements =
  inserter->insertIntrinsicElements();
  // if (nb_new_elements > 0) {
  //   this->reinitializeSolver();
  // }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::insertElementsFromMeshData(
    std::string physname) {
  AKANTU_DEBUG_IN();

  UInt material_index = SolidMechanicsModel::getMaterialIndex(physname);
  inserter->insertIntrinsicElements(physname, material_index);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initAutomaticInsertion() {
  AKANTU_DEBUG_IN();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  if (facet_stress_synchronizer != NULL) {
    DataAccessor * data_accessor = this;
    const ElementTypeMapArray<UInt> & rank_to_element =
        synch_parallel->getPrankToElement();

    facet_stress_synchronizer->updateFacetStressSynchronizer(
        *inserter, rank_to_element, *data_accessor);
  }
#endif

  facet_stress.initialize(inserter->getMeshFacets(),
                          _nb_component = 2 * Model::spatial_dimension *
                                          Model::spatial_dimension,
                          _spatial_dimension = Model::spatial_dimension - 1);

  // inserter->getMeshFacets().initElementTypeMapArray(
  //     facet_stress, 2 * spatial_dimension * spatial_dimension,
  //     spatial_dimension - 1);

  resizeFacetStress();

  /// compute normals on facets
  computeNormals();

  initStressInterpolation();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::updateAutomaticInsertion() {
  AKANTU_DEBUG_IN();

  inserter->limitCheckFacets();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  if (facet_stress_synchronizer != NULL) {
    DataAccessor * data_accessor = this;
    const ElementTypeMapArray<UInt> & rank_to_element =
        synch_parallel->getPrankToElement();

    facet_stress_synchronizer->updateFacetStressSynchronizer(
        *inserter, rank_to_element, *data_accessor);
  }
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initStressInterpolation() {
  Mesh & mesh_facets = inserter->getMeshFacets();

  /// compute quadrature points coordinates on facets
  Array<Real> & position = mesh.getNodes();

  ElementTypeMapArray<Real> quad_facets("quad_facets", id);
  quad_facets.initialize(mesh_facets, _nb_component = Model::spatial_dimension,
                         _spatial_dimension = Model::spatial_dimension - 1);
  // mesh_facets.initElementTypeMapArray(quad_facets, Model::spatial_dimension,
  //                                     Model::spatial_dimension - 1);

  getFEEngine("FacetsFEEngine")
      .interpolateOnIntegrationPoints(position, quad_facets);

  /// compute elements quadrature point positions and build
  /// element-facet quadrature points data structure
  ElementTypeMapArray<Real> elements_quad_facets("elements_quad_facets", id);

  elements_quad_facets.initialize(
      mesh, _nb_component = Model::spatial_dimension,
      _spatial_dimension = Model::spatial_dimension);
  // mesh.initElementTypeMapArray(elements_quad_facets,
  // Model::spatial_dimension,
  //                              Model::spatial_dimension);

  for (auto elem_gt : ghost_types) {
    for (auto & type : mesh.elementTypes(Model::spatial_dimension, elem_gt)) {
      UInt nb_element = mesh.getNbElement(type, elem_gt);
      if (nb_element == 0)
        continue;

      /// compute elements' quadrature points and list of facet
      /// quadrature points positions by element
      Array<Element> & facet_to_element =
          mesh_facets.getSubelementToElement(type, elem_gt);
      UInt nb_facet_per_elem = facet_to_element.getNbComponent();

      Array<Real> & el_q_facet = elements_quad_facets(type, elem_gt);

      ElementType facet_type = Mesh::getFacetType(type);

      UInt nb_quad_per_facet =
          getFEEngine("FacetsFEEngine").getNbIntegrationPoints(facet_type);

      el_q_facet.resize(nb_element * nb_facet_per_elem * nb_quad_per_facet);

      for (UInt el = 0; el < nb_element; ++el) {
        for (UInt f = 0; f < nb_facet_per_elem; ++f) {
          Element global_facet_elem = facet_to_element(el, f);
          UInt global_facet = global_facet_elem.element;
          GhostType facet_gt = global_facet_elem.ghost_type;
          const Array<Real> & quad_f = quad_facets(facet_type, facet_gt);

          for (UInt q = 0; q < nb_quad_per_facet; ++q) {
            for (UInt s = 0; s < Model::spatial_dimension; ++s) {
              el_q_facet(el * nb_facet_per_elem * nb_quad_per_facet +
                             f * nb_quad_per_facet + q,
                         s) = quad_f(global_facet * nb_quad_per_facet + q, s);
            }
          }
        }
      }
    }
  }

  /// loop over non cohesive materials
  for (UInt m = 0; m < materials.size(); ++m) {
    try {
      MaterialCohesive & mat __attribute__((unused)) =
          dynamic_cast<MaterialCohesive &>(*materials[m]);
    } catch (std::bad_cast &) {
      /// initialize the interpolation function
      materials[m]->initElementalFieldInterpolation(elements_quad_facets);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  // f_int += f_int_cohe
  for (auto & material : this->materials) {
    try {
      MaterialCohesive & mat = dynamic_cast<MaterialCohesive &>(*material);
      mat.computeTraction(_not_ghost);
    } catch (std::bad_cast & bce) {
    }
  }

  SolidMechanicsModel::assembleInternalForces();

  if (isDefaultSolverExplicit()) {
    for (auto & material : materials) {
      try {
        MaterialCohesive & mat = dynamic_cast<MaterialCohesive &>(*material);
        mat.computeEnergies();
      } catch (std::bad_cast & bce) {
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::computeNormals() {
  AKANTU_DEBUG_IN();

  Mesh & mesh_facets = this->inserter->getMeshFacets();
  this->getFEEngine("FacetsFEEngine")
      .computeNormalsOnIntegrationPoints(_not_ghost);

  /**
   *  @todo store tangents while computing normals instead of
   *  recomputing them as follows:
   */
  /* ------------------------------------------------------------------------ */
  UInt tangent_components =
      Model::spatial_dimension * (Model::spatial_dimension - 1);

  tangents.initialize(mesh_facets, _nb_component = tangent_components,
                      _spatial_dimension = Model::spatial_dimension - 1);
  // mesh_facets.initElementTypeMapArray(tangents, tangent_components,
  //                                     Model::spatial_dimension - 1);

  for (auto facet_type :
       mesh_facets.elementTypes(Model::spatial_dimension - 1)) {
    const Array<Real> & normals =
        this->getFEEngine("FacetsFEEngine")
            .getNormalsOnIntegrationPoints(facet_type);

    Array<Real> & tangents = this->tangents(facet_type);

    Math::compute_tangents(normals, tangents);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::interpolateStress() {

  ElementTypeMapArray<Real> by_elem_result("temporary_stress_by_facets", id);

  for (auto & material : materials) {
    try {
      MaterialCohesive & mat __attribute__((unused)) =
          dynamic_cast<MaterialCohesive &>(*material);
    } catch (std::bad_cast &) {
      /// interpolate stress on facet quadrature points positions
      material->interpolateStressOnFacets(facet_stress, by_elem_result);
    }
  }

#if defined(AKANTU_DEBUG_TOOLS)
  debug::element_manager.printData(
      debug::_dm_model_cohesive, "Interpolated stresses before", facet_stress);
#endif

  this->synchronize(_gst_smmc_facets_stress);

#if defined(AKANTU_DEBUG_TOOLS)
  debug::element_manager.printData(debug::_dm_model_cohesive,
                                   "Interpolated stresses", facet_stress);
#endif
}

/* -------------------------------------------------------------------------- */
UInt SolidMechanicsModelCohesive::checkCohesiveStress() {
  interpolateStress();

  for (auto & mat : materials) {
    try {
      MaterialCohesive & mat_cohesive = dynamic_cast<MaterialCohesive &>(*mat);
      /// check which not ghost cohesive elements are to be created
      mat_cohesive.checkInsertion();
    } catch (std::bad_cast &) {
    }
  }

  /// communicate data among processors
  this->synchronize(_gst_smmc_facets);

  /// insert cohesive elements
  UInt nb_new_elements = inserter->insertElements();

  // if (nb_new_elements > 0) {
  //   this->reinitializeSolver();
  // }

  return nb_new_elements;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::onElementsAdded(
    const Array<Element> & element_list, const NewElementsEvent & event) {
  AKANTU_DEBUG_IN();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  updateCohesiveSynchronizers();
#endif

  SolidMechanicsModel::onElementsAdded(element_list, event);

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  if (cohesive_element_synchronizer != NULL)
    cohesive_element_synchronizer->computeAllBufferSizes(*this);
#endif

  if (is_extrinsic)
    resizeFacetStress();

  ///  if (method != _explicit_lumped_mass) {
  ///    this->initSolver();
  ///  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::onNodesAdded(const Array<UInt> & new_nodes,
                                               const NewNodesEvent & event) {
  AKANTU_DEBUG_IN();

  // Array<UInt> nodes_list(nb_new_nodes);

  // for (UInt n = 0; n < nb_new_nodes; ++n)
  //   nodes_list(n) = doubled_nodes(n, 1);
  SolidMechanicsModel::onNodesAdded(new_nodes, event);

  UInt new_node, old_node;

  try {
    const auto & cohesive_event =
        dynamic_cast<const CohesiveNewNodesEvent &>(event);
    const auto & old_nodes = cohesive_event.getOldNodesList();

    auto copy = [this, &new_node, &old_node](auto & arr) {
      for (UInt s = 0; s < spatial_dimension; ++s) {
        arr(new_node, s) = arr(old_node, s);
      }
    };

    for (auto && pair : zip(new_nodes, old_nodes)) {
      std::tie(new_node, old_node) = pair;

      copy(*displacement);
      copy(*velocity);
      copy(*acceleration);
      copy(*blocked_dofs);

      if (current_position)
        copy(*current_position);

      if (previous_displacement)
        copy(*previous_displacement);
    }

    // if (this->getDOFManager().hasMatrix("M")) {
    //   this->assembleMass(old_nodes);
    // }

    // if (this->getDOFManager().hasLumpedMatrix("M")) {
    //   this->assembleMassLumped(old_nodes);
    // }

  } catch (std::bad_cast &) {
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::onEndSolveStep(const AnalysisMethod &) {

  AKANTU_DEBUG_IN();

  /*
   * This is required because the Cauchy stress is the stress measure that
   * is used to check the insertion of cohesive elements
   */
  for (auto & mat : materials) {
    if (mat->isFiniteDeformation())
      mat->computeAllCauchyStresses(_not_ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::printself(std::ostream & stream,
                                            int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "SolidMechanicsModelCohesive [" << std::endl;

  SolidMechanicsModel::printself(stream, indent + 1);

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::resizeFacetStress() {
  AKANTU_DEBUG_IN();

  Mesh & mesh_facets = inserter->getMeshFacets();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;

    Mesh::type_iterator it =
        mesh_facets.firstType(Model::spatial_dimension - 1, ghost_type);
    Mesh::type_iterator end =
        mesh_facets.lastType(Model::spatial_dimension - 1, ghost_type);
    for (; it != end; ++it) {
      ElementType type = *it;

      UInt nb_facet = mesh_facets.getNbElement(type, ghost_type);

      UInt nb_quadrature_points = getFEEngine("FacetsFEEngine")
                                      .getNbIntegrationPoints(type, ghost_type);

      UInt new_size = nb_facet * nb_quadrature_points;

      facet_stress(type, ghost_type).resize(new_size);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::addDumpGroupFieldToDumper(
    const std::string & dumper_name, const std::string & field_id,
    const std::string & group_name, const ElementKind & element_kind,
    bool padding_flag) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = Model::spatial_dimension;
  ElementKind _element_kind = element_kind;
  if (dumper_name == "cohesive elements") {
    _element_kind = _ek_cohesive;
  } else if (dumper_name == "facets") {
    spatial_dimension = Model::spatial_dimension - 1;
  }
  SolidMechanicsModel::addDumpGroupFieldToDumper(dumper_name, field_id,
                                                 group_name, spatial_dimension,
                                                 _element_kind, padding_flag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::onDump() {
  this->flattenAllRegisteredInternals(_ek_cohesive);
  SolidMechanicsModel::onDump();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
