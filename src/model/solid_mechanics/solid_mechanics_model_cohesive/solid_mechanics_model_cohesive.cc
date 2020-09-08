/**
 * @file   solid_mechanics_model_cohesive.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Wed Feb 21 2018
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
#include "solid_mechanics_model_cohesive.hh"
#include "aka_iterators.hh"
#include "cohesive_element_inserter.hh"
#include "element_synchronizer.hh"
#include "facet_synchronizer.hh"
#include "fe_engine_template.hh"
#include "global_ids_updater.hh"
#include "integrator_gauss.hh"
#include "material_cohesive.hh"
#include "mesh_accessor.hh"
#include "mesh_global_data_updater.hh"
#include "parser.hh"
#include "shape_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include "dumpable_inline_impl.hh"
#ifdef AKANTU_USE_IOHELPER
#include "dumper_iohelper_paraview.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */

namespace akantu {

class CohesiveMeshGlobalDataUpdater : public MeshGlobalDataUpdater {
public:
  CohesiveMeshGlobalDataUpdater(SolidMechanicsModelCohesive & model)
      : model(model), mesh(model.getMesh()),
        global_ids_updater(model.getMesh(), *model.cohesive_synchronizer) {}

  /* ------------------------------------------------------------------------ */
  std::tuple<UInt, UInt>
  updateData(NewNodesEvent & nodes_event,
             NewElementsEvent & elements_event) override {
    auto *cohesive_nodes_event =
        dynamic_cast<CohesiveNewNodesEvent *>(&nodes_event);
    if (cohesive_nodes_event == nullptr) {
      return std::make_tuple(nodes_event.getList().size(),
                             elements_event.getList().size());
    }

    /// update nodes type
    auto & new_nodes = cohesive_nodes_event->getList();
    auto & old_nodes = cohesive_nodes_event->getOldNodesList();

    auto local_nb_new_nodes = new_nodes.size();
    auto nb_new_nodes = local_nb_new_nodes;

    if (mesh.isDistributed()) {
      MeshAccessor mesh_accessor(mesh);
      auto & nodes_flags = mesh_accessor.getNodesFlags();
      auto nb_old_nodes = nodes_flags.size();
      nodes_flags.resize(nb_old_nodes + local_nb_new_nodes);

      for (auto && data : zip(old_nodes, new_nodes)) {
        UInt old_node;
        UInt new_node;
        std::tie(old_node, new_node) = data;
        nodes_flags(new_node) = nodes_flags(old_node);
      }

      model.updateCohesiveSynchronizers(elements_event);
      nb_new_nodes = global_ids_updater.updateGlobalIDs(new_nodes.size());
    }

    Vector<UInt> nb_new_stuff = {nb_new_nodes, elements_event.getList().size()};
    const auto & comm = mesh.getCommunicator();
    comm.allReduce(nb_new_stuff, SynchronizerOperation::_sum);

    if (nb_new_stuff(1) > 0) {
      mesh.sendEvent(elements_event);
    }

    if (nb_new_stuff(0) > 0) {
      mesh.sendEvent(nodes_event);
      // mesh.sendEvent(global_ids_updater.getChangedNodeEvent());
    }

    return std::make_tuple(nb_new_stuff(0), nb_new_stuff(1));
  }

private:
  SolidMechanicsModelCohesive & model;
  Mesh & mesh;
  GlobalIdsUpdater global_ids_updater;
};

/* -------------------------------------------------------------------------- */
SolidMechanicsModelCohesive::SolidMechanicsModelCohesive(
    Mesh & mesh, UInt dim, const ID & id, const MemoryID & memory_id)
    : SolidMechanicsModel(mesh, dim, id, memory_id,
                          ModelType::_solid_mechanics_model_cohesive),
      tangents("tangents", id), facet_stress("facet_stress", id),
      facet_material("facet_material", id) {
  AKANTU_DEBUG_IN();

  registerFEEngineObject<MyFEEngineCohesiveType>("CohesiveFEEngine", mesh,
                                                 Model::spatial_dimension);

  auto && tmp_material_selector =
      std::make_shared<DefaultMaterialCohesiveSelector>(*this);

  tmp_material_selector->setFallback(this->material_selector);
  this->material_selector = tmp_material_selector;

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("cohesive elements", id);
  this->mesh.addDumpMeshToDumper("cohesive elements", mesh,
                                 Model::spatial_dimension, _not_ghost,
                                 _ek_cohesive);
#endif

  if (this->mesh.isDistributed()) {
    /// create the distributed synchronizer for cohesive elements
    this->cohesive_synchronizer = std::make_unique<ElementSynchronizer>(
        mesh, "cohesive_distributed_synchronizer");

    auto & synchronizer = mesh.getElementSynchronizer();
    this->cohesive_synchronizer->split(synchronizer, [](auto && el) {
      return Mesh::getKind(el.type) == _ek_cohesive;
    });

    this->registerSynchronizer(*cohesive_synchronizer,
                               SynchronizationTag::_material_id);
    this->registerSynchronizer(*cohesive_synchronizer,
                               SynchronizationTag::_smm_stress);
    this->registerSynchronizer(*cohesive_synchronizer,
                               SynchronizationTag::_smm_boundary);
  }

  this->inserter = std::make_unique<CohesiveElementInserter>(
      this->mesh, id + ":cohesive_element_inserter");

  registerFEEngineObject<MyFEEngineFacetType>(
      "FacetsFEEngine", mesh.getMeshFacets(), Model::spatial_dimension - 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModelCohesive::~SolidMechanicsModelCohesive() = default;

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::setTimeStep(Real time_step,
                                              const ID & solver_id) {
  SolidMechanicsModel::setTimeStep(time_step, solver_id);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.getDumper("cohesive elements").setTimeStep(time_step);
#endif
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initFullImpl(const ModelOptions & options) {
  AKANTU_DEBUG_IN();

  const auto & smmc_options =
      aka::as_type<SolidMechanicsModelCohesiveOptions>(options);

  this->is_extrinsic = smmc_options.is_extrinsic;

  inserter->setIsExtrinsic(is_extrinsic);

  if (mesh.isDistributed()) {
    auto & mesh_facets = inserter->getMeshFacets();
    auto & synchronizer =
        aka::as_type<FacetSynchronizer>(mesh_facets.getElementSynchronizer());

    // synchronizeGhostFacetsConnectivity();

    /// create the facet synchronizer for extrinsic simulations
    if (is_extrinsic) {
      facet_stress_synchronizer = std::make_unique<ElementSynchronizer>(
          synchronizer, id + ":facet_stress_synchronizer");
      facet_stress_synchronizer->swapSendRecv();
      this->registerSynchronizer(*facet_stress_synchronizer,
                                 SynchronizationTag::_smmc_facets_stress);
    }
  }

  MeshAccessor mesh_accessor(mesh);
  mesh_accessor.registerGlobalDataUpdater(
      std::make_unique<CohesiveMeshGlobalDataUpdater>(*this));

  ParserSection section;
  bool is_empty;
  std::tie(section, is_empty) = this->getParserSection();

  if (not is_empty) {
    auto inserter_section =
        section.getSubSections(ParserType::_cohesive_inserter);
    if (inserter_section.begin() != inserter_section.end()) {
      inserter->parseSection(*inserter_section.begin());
    }
  }

  SolidMechanicsModel::initFullImpl(options);

  AKANTU_DEBUG_OUT();
} // namespace akantu

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initMaterials() {
  AKANTU_DEBUG_IN();

  // make sure the material are instantiated
  if (not are_materials_instantiated) {
    instantiateMaterials();
  }

  /// find the first cohesive material
  UInt cohesive_index = UInt(-1);

  for (auto && material : enumerate(materials)) {
    if (dynamic_cast<MaterialCohesive *>(std::get<1>(material).get()) !=
        nullptr) {
      cohesive_index = std::get<0>(material);
      break;
    }
  }

  if (cohesive_index == UInt(-1)) {
    AKANTU_EXCEPTION("No cohesive materials in the material input file");
  }

  material_selector->setFallback(cohesive_index);

  // set the facet information in the material in case of dynamic insertion
  // to know what material to call for stress checks

  const Mesh & mesh_facets = inserter->getMeshFacets();
  facet_material.initialize(
      mesh_facets, _spatial_dimension = spatial_dimension - 1,
      _with_nb_element = true,
      _default_value = material_selector->getFallbackValue());

  for_each_element(
      mesh_facets,
      [&](auto && element) {
        auto mat_index = (*material_selector)(element);
        auto & mat = aka::as_type<MaterialCohesive>(*materials[mat_index]);
        facet_material(element) = mat_index;
        if (is_extrinsic) {
          mat.addFacet(element);
        }
      },
      _spatial_dimension = spatial_dimension - 1, _ghost_type = _not_ghost);

  SolidMechanicsModel::initMaterials();

  if (is_extrinsic) {
    this->initAutomaticInsertion();
  } else {
    this->insertIntrinsicElements();
  }

  AKANTU_DEBUG_OUT();
} // namespace akantu

/* -------------------------------------------------------------------------- */
/**
 * Initialize the model,basically it  pre-compute the shapes, shapes derivatives
 * and jacobian
 */
void SolidMechanicsModelCohesive::initModel() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::initModel();

  /// add cohesive type connectivity
  ElementType type = _not_defined;
  for (auto && type_ghost : ghost_types) {
    for (const auto & tmp_type :
         mesh.elementTypes(spatial_dimension, type_ghost)) {
      const auto & connectivity = mesh.getConnectivity(tmp_type, type_ghost);
      if (connectivity.empty()) {
        continue;
      }

      type = tmp_type;
      auto type_facet = Mesh::getFacetType(type);
      auto type_cohesive = FEEngine::getCohesiveElementType(type_facet);
      mesh.addConnectivityType(type_cohesive, type_ghost);
    }
  }
  AKANTU_DEBUG_ASSERT(type != _not_defined, "No elements in the mesh");

  getFEEngine("CohesiveFEEngine").initShapeFunctions(_not_ghost);
  getFEEngine("CohesiveFEEngine").initShapeFunctions(_ghost);

  if (is_extrinsic) {
    getFEEngine("FacetsFEEngine").initShapeFunctions(_not_ghost);
    getFEEngine("FacetsFEEngine").initShapeFunctions(_ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::insertIntrinsicElements() {
  AKANTU_DEBUG_IN();
  inserter->insertIntrinsicElements();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initAutomaticInsertion() {
  AKANTU_DEBUG_IN();

  this->inserter->limitCheckFacets();
  this->updateFacetStressSynchronizer();
  this->resizeFacetStress();

  /// compute normals on facets
  this->computeNormals();

  this->initStressInterpolation();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::updateAutomaticInsertion() {
  AKANTU_DEBUG_IN();

  this->inserter->limitCheckFacets();
  this->updateFacetStressSynchronizer();

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
    for (const auto & type : mesh.elementTypes(Model::spatial_dimension, elem_gt)) {
      UInt nb_element = mesh.getNbElement(type, elem_gt);
      if (nb_element == 0) {
        continue;
      }

      /// compute elements' quadrature points and list of facet
      /// quadrature points positions by element
      const auto & facet_to_element =
          mesh_facets.getSubelementToElement(type, elem_gt);
      auto & el_q_facet = elements_quad_facets(type, elem_gt);

      auto facet_type = Mesh::getFacetType(type);
      auto nb_quad_per_facet =
          getFEEngine("FacetsFEEngine").getNbIntegrationPoints(facet_type);
      auto nb_facet_per_elem = facet_to_element.getNbComponent();

      // small hack in the loop to skip boundary elements, they are silently
      // initialized to NaN to see if this causes problems
      el_q_facet.resize(nb_element * nb_facet_per_elem * nb_quad_per_facet,
                        std::numeric_limits<Real>::quiet_NaN());

      for (auto && data :
           zip(make_view(facet_to_element),
               make_view(el_q_facet, spatial_dimension, nb_quad_per_facet))) {
        const auto & global_facet = std::get<0>(data);
        auto & el_q = std::get<1>(data);

        if (global_facet == ElementNull) {
          continue;
        }

        Matrix<Real> quad_f =
            make_view(quad_facets(global_facet.type, global_facet.ghost_type),
                      spatial_dimension, nb_quad_per_facet)
                .begin()[global_facet.element];

        el_q = quad_f;

        // for (UInt q = 0; q < nb_quad_per_facet; ++q) {
        //   for (UInt s = 0; s < Model::spatial_dimension; ++s) {
        //     el_q_facet(el * nb_facet_per_elem * nb_quad_per_facet +
        //                    f * nb_quad_per_facet + q,
        //                s) = quad_f(global_facet * nb_quad_per_facet + q,
        //                s);
        //   }
        // }
        //}
      }
    }
  }

  /// loop over non cohesive materials
  for (auto && material : materials) {
    if (aka::is_of_type<MaterialCohesive>(material)) {
      continue;
    }
    /// initialize the interpolation function
    material->initElementalFieldInterpolation(elements_quad_facets);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  // f_int += f_int_cohe
  for (auto & material : this->materials) {
    try {
      auto & mat = aka::as_type<MaterialCohesive>(*material);
      mat.computeTraction(_not_ghost);
    } catch (std::bad_cast & bce) {
    }
  }

  SolidMechanicsModel::assembleInternalForces();

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
    if (not aka::is_of_type<MaterialCohesive>(material)) {
      /// interpolate stress on facet quadrature points positions
      material->interpolateStressOnFacets(facet_stress, by_elem_result);
    }
  }

  this->synchronize(SynchronizationTag::_smmc_facets_stress);
}

/* -------------------------------------------------------------------------- */
UInt SolidMechanicsModelCohesive::checkCohesiveStress() {
  AKANTU_DEBUG_IN();

  if (not is_extrinsic) {
    AKANTU_EXCEPTION(
        "This function can only be used for extrinsic cohesive elements");
  }

  interpolateStress();

  for (auto & mat : materials) {
    if (aka::is_of_type<MaterialCohesive>(mat)) {
      /// check which not ghost cohesive elements are to be created
      auto * mat_cohesive = aka::as_type<MaterialCohesive>(mat.get());
      mat_cohesive->checkInsertion();
    }
  }

  /// communicate data among processors
  // this->synchronize(SynchronizationTag::_smmc_facets);

  /// insert cohesive elements
  UInt nb_new_elements = inserter->insertElements();

  // if (nb_new_elements > 0) {
  //   this->reinitializeSolver();
  // }

  AKANTU_DEBUG_OUT();

  return nb_new_elements;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::onElementsAdded(
    const Array<Element> & element_list, const NewElementsEvent & event) {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::onElementsAdded(element_list, event);

  if (is_extrinsic) {
    resizeFacetStress();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::onNodesAdded(const Array<UInt> & new_nodes,
                                               const NewNodesEvent & event) {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::onNodesAdded(new_nodes, event);

  const CohesiveNewNodesEvent * cohesive_event;
  if ((cohesive_event = dynamic_cast<const CohesiveNewNodesEvent *>(&event)) ==
      nullptr) {
    return;
  }

  const auto & old_nodes = cohesive_event->getOldNodesList();

  auto copy = [this, &new_nodes, &old_nodes](auto & arr) {
    UInt new_node;
    UInt old_node;

    auto view = make_view(arr, spatial_dimension);
    auto begin = view.begin();

    for (auto && pair : zip(new_nodes, old_nodes)) {
      std::tie(new_node, old_node) = pair;

      auto old_ = begin + old_node;
      auto new_ = begin + new_node;

      *new_ = *old_;
    }
  };

  copy(*displacement);
  copy(*blocked_dofs);

  if (velocity) {
    copy(*velocity);
  }

  if (acceleration) {
    copy(*acceleration);
  }

  if (current_position) {
    copy(*current_position);
  }

  if (previous_displacement) {
    copy(*previous_displacement);
  }

  // if (external_force)
  //   copy(*external_force);
  // if (internal_force)
  //   copy(*internal_force);

  if (displacement_increment) {
    copy(*displacement_increment);
  }

  copy(getDOFManager().getSolution("displacement"));
  // this->assembleMassLumped();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::afterSolveStep(bool converged) {
  AKANTU_DEBUG_IN();

  /*
   * This is required because the Cauchy stress is the stress measure that
   * is used to check the insertion of cohesive elements
   */
  if (converged) {
    for (auto & mat : materials) {
      if (mat->isFiniteDeformation()) {
        mat->computeAllCauchyStresses(_not_ghost);
      }
    }
  }

  SolidMechanicsModel::afterSolveStep(converged);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::printself(std::ostream & stream,
                                            int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "SolidMechanicsModelCohesive ["
         << "\n";
  SolidMechanicsModel::printself(stream, indent + 2);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::resizeFacetStress() {
  AKANTU_DEBUG_IN();

  this->facet_stress.initialize(getFEEngine("FacetsFEEngine"),
                                _nb_component =
                                    2 * spatial_dimension * spatial_dimension,
                                _spatial_dimension = spatial_dimension - 1);

  // for (auto && ghost_type : ghost_types) {
  //   for (const const auto & type :
  //        mesh_facets.elementTypes(spatial_dimension - 1, ghost_type)) {
  //     UInt nb_facet = mesh_facets.getNbElement(type, ghost_type);

  //     UInt nb_quadrature_points = getFEEngine("FacetsFEEngine")
  //                                     .getNbIntegrationPoints(type,
  //                                     ghost_type);

  //     UInt new_size = nb_facet * nb_quadrature_points;

  //     facet_stress(type, ghost_type).resize(new_size);
  //   }
  // }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::addDumpGroupFieldToDumper(
    const std::string & dumper_name, const std::string & field_id,
    const std::string & group_name, ElementKind element_kind,
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
