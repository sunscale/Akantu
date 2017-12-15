/**
 * @file   solid_mechanics_model_cohesive_parallel.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 *
 * @brief  Functions for parallel cohesive elements
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */
/* -------------------------------------------------------------------------- */
#include "element_synchronizer.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "solid_mechanics_model_tmpl.hh"
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#include <type_traits>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::synchronizeGhostFacetsConnectivity() {
  AKANTU_DEBUG_IN();

  const Communicator & comm = mesh.getCommunicator();
  Int psize = comm.getNbProc();

  if (psize > 1) {
    /// get global connectivity for not ghost facets
    global_connectivity =
        new ElementTypeMapArray<UInt>("global_connectivity", id);

    auto & mesh_facets = inserter->getMeshFacets();

    global_connectivity->initialize(
        mesh_facets, _spatial_dimension = spatial_dimension - 1,
        _with_nb_element = true, _with_nb_nodes_per_element = true,
        _element_kind = _ek_regular,
        _ghost_type = _not_ghost);

    mesh_facets.getGlobalConnectivity(*global_connectivity);

    /// communicate
    synchronize(_gst_smmc_facets_conn);

    /// flip facets
    MeshUtils::flipFacets(mesh_facets, *global_connectivity, _ghost);

    delete global_connectivity;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::updateCohesiveSynchronizers() {
  /// update synchronizers if needed

  if (not mesh.isDistributed())
    return;

  auto & mesh_facets = inserter->getMeshFacets();
  auto & facet_synchronizer = mesh_facets.getElementSynchronizer();
  const auto & cfacet_synchronizer = facet_synchronizer;

  // update the cohesive element synchronizer
  cohesive_synchronizer->updateSchemes(
      [&](auto && scheme, auto && proc, auto && direction) {
        auto & facet_scheme =
            cfacet_synchronizer.getCommunications().getScheme(proc, direction);

        for (auto && facet : facet_scheme) {
          const auto & connected_elements = mesh_facets.getElementToSubelement(
              facet.type,
              facet.ghost_type)(facet.element); // slow access here
          const auto & cohesive_element = connected_elements[1];

          auto && cohesive_type = FEEngine::getCohesiveElementType(facet.type);
          auto old_nb_cohesive_elements =
              mesh.getNbElement(cohesive_type, facet.ghost_type);
          old_nb_cohesive_elements -=
            mesh.getData<UInt>("facet_to_double", facet.type,
                               facet.ghost_type).size();

          if (cohesive_element.kind() == _ek_cohesive and
              cohesive_element.element >= old_nb_cohesive_elements) {
            scheme.push_back(cohesive_element);
          }
        }
      });

  const auto & comm = mesh.getCommunicator();
  auto && my_rank = comm.whoAmI();

  // update the facet stress synchronizer
  facet_stress_synchronizer->updateSchemes([&](auto && scheme, auto && proc,
                                               auto && /*direction*/) {
    auto it_element = scheme.begin();
    for (auto && element : scheme) {
      auto && facet_check = inserter->getCheckFacets(
          element.type, element.ghost_type)(element.element); // slow access
                                                              // here

      if (facet_check) {
        auto && connected_elements = mesh_facets.getElementToSubelement(
            element.type, element.ghost_type)(element.element); // slow access
                                                                // here
        auto && rank_left = facet_synchronizer.getRank(connected_elements[0]);
        auto && rank_right = facet_synchronizer.getRank(connected_elements[1]);

        // keep element if the element is still a boundary element between two
        // processors
        if ((rank_left == Int(proc) and rank_right == my_rank) or
            (rank_left == my_rank and rank_right == Int(proc))) {
          *it_element = element;
          ++it_element;
        }
      }
    }
    scheme.resize(it_element - scheme.begin());
  });
}

/* -------------------------------------------------------------------------- */
template <typename T>
void SolidMechanicsModelCohesive::packFacetStressDataHelper(
    const ElementTypeMapArray<T> & data_to_pack, CommunicationBuffer & buffer,
    const Array<Element> & elements) const {
  packUnpackFacetStressDataHelper<T, true>(
      const_cast<ElementTypeMapArray<T> &>(data_to_pack), buffer, elements);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void SolidMechanicsModelCohesive::unpackFacetStressDataHelper(
    ElementTypeMapArray<T> & data_to_unpack, CommunicationBuffer & buffer,
    const Array<Element> & elements) const {
  packUnpackFacetStressDataHelper<T, false>(data_to_unpack, buffer, elements);
}

/* -------------------------------------------------------------------------- */
template <typename T, bool pack_helper>
void SolidMechanicsModelCohesive::packUnpackFacetStressDataHelper(
    ElementTypeMapArray<T> & data_to_pack, CommunicationBuffer & buffer,
    const Array<Element> & elements) const {
  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  UInt nb_quad_per_elem = 0;
  UInt sp2 = spatial_dimension * spatial_dimension;
  UInt nb_component = sp2 * 2;
  bool element_rank = false;
  Mesh & mesh_facets = inserter->getMeshFacets();

  Array<T> * vect = nullptr;
  Array<std::vector<Element>> * element_to_facet = nullptr;

  auto & fe_engine = this->getFEEngine("FacetsFEEngine");
  for (auto && el : elements) {
    if (el.type == _not_defined)
      AKANTU_EXCEPTION("packUnpackFacetStressDataHelper called with wrong inputs");

    if (el.type != current_element_type ||
        el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type = el.ghost_type;
      vect = &data_to_pack(el.type, el.ghost_type);

      element_to_facet =
          &(mesh_facets.getElementToSubelement(el.type, el.ghost_type));

      nb_quad_per_elem =
          fe_engine.getNbIntegrationPoints(el.type, el.ghost_type);
    }

    if (pack_helper)
      element_rank =
          (*element_to_facet)(el.element)[0].ghost_type != _not_ghost;
    else
      element_rank =
          (*element_to_facet)(el.element)[0].ghost_type == _not_ghost;

    for (UInt q = 0; q < nb_quad_per_elem; ++q) {
      Vector<T> data(vect->storage() +
                         (el.element * nb_quad_per_elem + q) * nb_component +
                         element_rank * sp2,
                     sp2);

      if (pack_helper)
        buffer << data;
      else
        buffer >> data;
    }
  }
}

/* -------------------------------------------------------------------------- */
UInt SolidMechanicsModelCohesive::getNbQuadsForFacetCheck(
    const Array<Element> & elements) const {
  UInt nb_quads = 0;
  UInt nb_quad_per_facet = 0;

  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  auto & fe_engine = this->getFEEngine("FacetsFEEngine");
  for(auto & el : elements) {
    if (el.type != current_element_type ||
        el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type = el.ghost_type;

      nb_quad_per_facet = fe_engine
                              .getNbIntegrationPoints(el.type, el.ghost_type);
    }

    nb_quads += nb_quad_per_facet;
  }

  return nb_quads;
}

/* -------------------------------------------------------------------------- */
UInt SolidMechanicsModelCohesive::getNbData(
    const Array<Element> & elements, const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  if (elements.size() == 0)
    return 0;

  /// regular element case
  if (elements(0).kind() == _ek_regular) {

    switch (tag) {
    case _gst_smmc_facets: {
      size += elements.size() * sizeof(bool);
      break;
    }
    case _gst_smmc_facets_conn: {
      UInt nb_nodes = Mesh::getNbNodesPerElementList(elements);
      size += nb_nodes * sizeof(UInt);
      break;
    }
    case _gst_smmc_facets_stress: {
      UInt nb_quads = getNbQuadsForFacetCheck(elements);
      size += nb_quads * spatial_dimension * spatial_dimension * sizeof(Real);
      break;
    }
    default: { size += SolidMechanicsModel::getNbData(elements, tag); }
    }
  }
  /// cohesive element case
  else if (elements(0).kind() == _ek_cohesive) {

    switch (tag) {
    case _gst_material_id: {
      size += elements.size() * sizeof(UInt);
      break;
    }
    case _gst_smm_boundary: {
      UInt nb_nodes_per_element = 0;

      for (auto && el : elements) {
        nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
      }

      // force, displacement, boundary
      size += nb_nodes_per_element * spatial_dimension *
              (2 * sizeof(Real) + sizeof(bool));
      break;
    }
    default:
      break;
    }

    if (tag != _gst_material_id && tag != _gst_smmc_facets) {
      splitByMaterial(elements, [&](auto && mat, auto && elements) {
        size += mat.getNbData(elements, tag);
      });
    }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::packData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  if (elements.size() == 0)
    return;

  if (elements(0).kind() == _ek_regular) {
    switch (tag) {
    case _gst_smmc_facets: {
      packElementalDataHelper(inserter->getInsertionFacetsByElement(), buffer,
                              elements, false, getFEEngine());
      break;
    }
    case _gst_smmc_facets_conn: {
      packElementalDataHelper(*global_connectivity, buffer, elements, false,
                              getFEEngine());
      break;
    }
    case _gst_smmc_facets_stress: {
      packFacetStressDataHelper(facet_stress, buffer, elements);
      break;
    }
    default: { SolidMechanicsModel::packData(buffer, elements, tag); }
    }
  } else if (elements(0).kind() == _ek_cohesive) {
    switch (tag) {
    case _gst_material_id: {
      packElementalDataHelper(material_index, buffer, elements, false,
                              getFEEngine("CohesiveFEEngine"));
      break;
    }
    case _gst_smm_boundary: {
      packNodalDataHelper(*internal_force, buffer, elements, mesh);
      packNodalDataHelper(*velocity, buffer, elements, mesh);
      packNodalDataHelper(*blocked_dofs, buffer, elements, mesh);
      break;
    }
    default: {}
    }

    if (tag != _gst_material_id && tag != _gst_smmc_facets) {
      splitByMaterial(elements, [&](auto && mat, auto && elements) {
        mat.packData(buffer, elements, tag);
      });
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::unpackData(CommunicationBuffer & buffer,
                                             const Array<Element> & elements,
                                             const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  if (elements.size() == 0)
    return;

  if (elements(0).kind() == _ek_regular) {
    switch (tag) {
    case _gst_smmc_facets: {
      unpackElementalDataHelper(inserter->getInsertionFacetsByElement(), buffer,
                                elements, false, getFEEngine());
      break;
    }
    case _gst_smmc_facets_conn: {
      unpackElementalDataHelper(*global_connectivity, buffer, elements, false,
                                getFEEngine());
      break;
    }
    case _gst_smmc_facets_stress: {
      unpackFacetStressDataHelper(facet_stress, buffer, elements);
      break;
    }
    default: { SolidMechanicsModel::unpackData(buffer, elements, tag); }
    }
  } else if (elements(0).kind() == _ek_cohesive) {
    switch (tag) {
    case _gst_material_id: {
      unpackElementalDataHelper(material_index, buffer, elements, false,
                                getFEEngine("CohesiveFEEngine"));
      break;
    }
    case _gst_smm_boundary: {
      unpackNodalDataHelper(*internal_force, buffer, elements, mesh);
      unpackNodalDataHelper(*velocity, buffer, elements, mesh);
      unpackNodalDataHelper(*blocked_dofs, buffer, elements, mesh);
      break;
    }
    default: {}
    }

    if (tag != _gst_material_id && tag != _gst_smmc_facets) {
      splitByMaterial(elements, [&](auto && mat, auto && elements) {
        mat.unpackData(buffer, elements, tag);
      });
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
