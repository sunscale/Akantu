/**
 * @file   solid_mechanics_model_cohesive.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Fri Sep 05 2014
 *
 * @brief  Solid mechanics model for cohesive elements
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include <algorithm>
#include "shape_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "dumpable_inline_impl.hh"
#include "material_cohesive.hh"

#ifdef AKANTU_USE_IOHELPER
#  include "dumper_paraview.hh"
#endif

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

const SolidMechanicsModelCohesiveOptions default_solid_mechanics_model_cohesive_options(_explicit_lumped_mass,
											false,
											false,
											true);

/* -------------------------------------------------------------------------- */

SolidMechanicsModelCohesive::SolidMechanicsModelCohesive(Mesh & mesh,
							 UInt dim,
							 const ID & id,
							 const MemoryID & memory_id) :
  SolidMechanicsModel(mesh, dim, id, memory_id),
  tangents("tangents", id),
  stress_on_facet("stress_on_facet", id),
  facet_stress("facet_stress", id),
  facet_material("facet_material", id) {
  AKANTU_DEBUG_IN();

  inserter = NULL;

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  facet_synchronizer = NULL;
  facet_stress_synchronizer = NULL;
  cohesive_distributed_synchronizer = NULL;
  global_connectivity = NULL;
#endif

  delete material_selector;
  material_selector = new DefaultMaterialCohesiveSelector(*this);

  this->registerEventHandler(*this);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("cohesive elements", id);
  this->mesh.addDumpMeshToDumper("cohesive elements",
			    mesh, spatial_dimension, _not_ghost, _ek_cohesive);
#endif

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
SolidMechanicsModelCohesive::~SolidMechanicsModelCohesive() {
  AKANTU_DEBUG_IN();

  delete inserter;

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  delete cohesive_distributed_synchronizer;
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

  this->stress_interpolation = smmc_options.stress_interpolation;
  this->is_extrinsic = smmc_options.extrinsic;

  if (!inserter)
    inserter = new CohesiveElementInserter(mesh, is_extrinsic, synch_parallel,
					   id+":cohesive_element_inserter");

  SolidMechanicsModel::initFull(options);

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  if (facet_synchronizer != NULL)
    inserter->initParallel(facet_synchronizer);
#endif

  if (is_extrinsic)
    initAutomaticInsertion();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initMaterials() {
  AKANTU_DEBUG_IN();

  // make sure the material are instantiated
  if(!are_materials_instantiated) instantiateMaterials();

  /// find the first cohesive material
  UInt cohesive_index = 0;

  while ((dynamic_cast<MaterialCohesive *>(materials[cohesive_index]) == NULL)
	 && cohesive_index <= materials.size())
    ++cohesive_index;

  AKANTU_DEBUG_ASSERT(cohesive_index != materials.size(),
		      "No cohesive materials in the material input file");

  material_selector->setFallback(cohesive_index);

  // set the facet information in the material in case of dynamic insertion
  if (is_extrinsic) {
    const Mesh & mesh_facets = inserter->getMeshFacets();
    mesh_facets.initElementTypeMapArray(facet_material, 1, spatial_dimension - 1);

    Element element;
    for (ghost_type_t::iterator gt = ghost_type_t::begin(); gt != ghost_type_t::end(); ++gt) {
      element.ghost_type = *gt;
      Mesh::type_iterator first = mesh_facets.firstType(spatial_dimension - 1, *gt);
      Mesh::type_iterator last  = mesh_facets.lastType(spatial_dimension - 1, *gt);
      for(;first != last; ++first) {
	element.type = *first;
	Array<UInt> & f_material = facet_material(*first, *gt);
	UInt nb_element = mesh_facets.getNbElement(*first, *gt);
	f_material.resize(nb_element);
	f_material.set(cohesive_index);
	for (UInt el = 0; el < nb_element; ++el) {
	  element.element = el;
	  UInt mat_index = (*material_selector)(element);
	  f_material(el) = mat_index;
	  MaterialCohesive & mat = dynamic_cast<MaterialCohesive &>(*materials[mat_index]);
	  mat.addFacet(element);
	}
      }
    }
  } else {
    for (ghost_type_t::iterator gt = ghost_type_t::begin(); gt != ghost_type_t::end(); ++gt) {
      Mesh::type_iterator first = mesh.firstType(spatial_dimension, *gt, _ek_cohesive);
      Mesh::type_iterator last  = mesh.lastType(spatial_dimension, *gt, _ek_cohesive);

      for(;first != last; ++first) {
	Array<UInt> & el_id_by_mat = element_index_by_material(*first, *gt);
	Vector<UInt> el_mat(2); el_mat(0) = cohesive_index; el_mat(1) = 0;
	el_id_by_mat.set(el_mat);
      }
    }
  }

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

  registerFEEngineObject<MyFEEngineCohesiveType>("CohesiveFEEngine", mesh, spatial_dimension);

  /// add cohesive type connectivity
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
	ElementType type_facet = Mesh::getFacetType(type);
	ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);
	mesh.addConnectivityType(type_cohesive, type_ghost);
      }
    }
  }

  AKANTU_DEBUG_ASSERT(type != _not_defined, "No elements in the mesh");

  getFEEngine("CohesiveFEEngine").initShapeFunctions(_not_ghost);
  getFEEngine("CohesiveFEEngine").initShapeFunctions(_ghost);

  registerFEEngineObject<MyFEEngineType>("FacetsFEEngine",
					 mesh.getMeshFacets(),
					 spatial_dimension - 1);

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
  inserter->insertIntrinsicElements();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initAutomaticInsertion() {
  AKANTU_DEBUG_IN();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  if (facet_stress_synchronizer != NULL) {
    DataAccessor * data_accessor = this;
    const ElementTypeMapArray<UInt> & rank_to_element = synch_parallel->getPrankToElement();

    facet_stress_synchronizer->updateFacetStressSynchronizer(*inserter,
							     rank_to_element,
							     *data_accessor);
  }
#endif

  inserter->getMeshFacets().initElementTypeMapArray(facet_stress,
						    2 * spatial_dimension * spatial_dimension,
						    spatial_dimension - 1);

  resizeFacetStress();

  /// compute normals on facets
  computeNormals();

  /// allocate stress_on_facet to store element stress on facets
  mesh.initElementTypeMapArray(stress_on_facet,
			      spatial_dimension * spatial_dimension,
			      spatial_dimension);

  Mesh::type_iterator it   = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last = mesh.lastType(spatial_dimension);

  for (; it != last; ++it) {
    ElementType type = *it;
    UInt nb_element = mesh.getNbElement(type);

    UInt nb_facet_per_elem = Mesh::getNbFacetsPerElement(type);
    ElementType type_facet = Mesh::getFacetType(type);
    UInt nb_quad_per_facet = getFEEngine("FacetsFEEngine").getNbQuadraturePoints(type_facet);

    stress_on_facet(type).resize(nb_quad_per_facet * nb_facet_per_elem * nb_element);
  }

  if (stress_interpolation)
    initStressInterpolation();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::updateAutomaticInsertion() {
  AKANTU_DEBUG_IN();

  inserter->limitCheckFacets();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
    if (facet_stress_synchronizer != NULL) {
      DataAccessor * data_accessor = this;
      const ElementTypeMapArray<UInt> & rank_to_element = synch_parallel->getPrankToElement();

      facet_stress_synchronizer->updateFacetStressSynchronizer(*inserter,
							       rank_to_element,
							       *data_accessor);
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
  mesh_facets.initElementTypeMapArray(quad_facets,
				     spatial_dimension, spatial_dimension - 1);

  getFEEngine("FacetsFEEngine").interpolateOnQuadraturePoints(position, quad_facets);

  /// compute elements quadrature point positions and build
  /// element-facet quadrature points data structure
  ElementTypeMapArray<Real> elements_quad_facets("elements_quad_facets", id);

  mesh.initElementTypeMapArray(elements_quad_facets,
			      spatial_dimension,
			      spatial_dimension);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType elem_gt = *gt;

    Mesh::type_iterator it   = mesh.firstType(spatial_dimension, elem_gt);
    Mesh::type_iterator last = mesh.lastType(spatial_dimension, elem_gt);

    for (; it != last; ++it) {
      ElementType type = *it;
      UInt nb_element = mesh.getNbElement(type, elem_gt);
      if (nb_element == 0) continue;

      /// compute elements' quadrature points and list of facet
      /// quadrature points positions by element
      Array<Element> & facet_to_element = mesh_facets.getSubelementToElement(type,
									     elem_gt);
      UInt nb_facet_per_elem = facet_to_element.getNbComponent();

      Array<Real> & el_q_facet = elements_quad_facets(type, elem_gt);

      ElementType facet_type = Mesh::getFacetType(type);

      UInt nb_quad_per_facet =
	getFEEngine("FacetsFEEngine").getNbQuadraturePoints(facet_type);

      el_q_facet.resize(nb_element * nb_facet_per_elem * nb_quad_per_facet);

      for (UInt el = 0; el < nb_element; ++el) {
	for (UInt f = 0; f < nb_facet_per_elem; ++f) {
	  Element global_facet_elem = facet_to_element(el, f);
	  UInt global_facet = global_facet_elem.element;
	  GhostType facet_gt = global_facet_elem.ghost_type;
	  const Array<Real> & quad_f = quad_facets(facet_type, facet_gt);

	  for (UInt q = 0; q < nb_quad_per_facet; ++q) {
	    for (UInt s = 0; s < spatial_dimension; ++s) {
	      el_q_facet(el * nb_facet_per_elem * nb_quad_per_facet
			 + f * nb_quad_per_facet + q, s)
		= quad_f(global_facet * nb_quad_per_facet + q, s);
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
    } catch(std::bad_cast&) {
      /// initialize the interpolation function
      materials[m]->initElementalFieldInterpolation(elements_quad_facets);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::updateResidual(bool need_initialize) {
  AKANTU_DEBUG_IN();

  if (need_initialize) initializeUpdateResidualData();

  // f -= fint
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    try {
      MaterialCohesive & mat = dynamic_cast<MaterialCohesive &>(**mat_it);
      mat.computeTraction(_not_ghost);
    } catch (std::bad_cast & bce) { }
  }

  SolidMechanicsModel::updateResidual(false);

  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    try {
      MaterialCohesive & mat = dynamic_cast<MaterialCohesive &>(**mat_it);
      mat.computeEnergies();
    } catch (std::bad_cast & bce) { }
  }


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::computeNormals() {
  AKANTU_DEBUG_IN();

  Mesh & mesh_facets = inserter->getMeshFacets();
  getFEEngine("FacetsFEEngine").computeNormalsOnControlPoints(_not_ghost);

  /**
   *  @todo store tangents while computing normals instead of
   *  recomputing them as follows:
   */
  /* ------------------------------------------------------------------------ */
  UInt tangent_components = spatial_dimension * (spatial_dimension - 1);

  mesh_facets.initElementTypeMapArray(tangents,
				     tangent_components,
				     spatial_dimension - 1);

  Mesh::type_iterator it   = mesh_facets.firstType(spatial_dimension - 1);
  Mesh::type_iterator last = mesh_facets.lastType(spatial_dimension - 1);

  for (; it != last; ++it) {
    ElementType facet_type = *it;

    const Array<Real> & normals =
      getFEEngine("FacetsFEEngine").getNormalsOnQuadPoints(facet_type);

    UInt nb_quad = normals.getSize();

    Array<Real> & tang = tangents(facet_type);
    tang.resize(nb_quad);

    Real * normal_it = normals.storage();
    Real * tangent_it = tang.storage();

    /// compute first tangent
    for (UInt q = 0; q < nb_quad; ++q) {

      /// if normal is orthogonal to xy plane, arbitrarly define tangent
      if ( Math::are_float_equal(Math::norm2(normal_it), 0) )
	tangent_it[0] = 1;
      else
	Math::normal2(normal_it, tangent_it);

      normal_it += spatial_dimension;
      tangent_it += tangent_components;
    }

    /// compute second tangent (3D case)
    if (spatial_dimension == 3) {
      normal_it = normals.storage();
      tangent_it = tang.storage();

      for (UInt q = 0; q < nb_quad; ++q) {
	Math::normal3(normal_it, tangent_it, tangent_it + spatial_dimension);
	normal_it += spatial_dimension;
	tangent_it += tangent_components;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::averageStressOnFacets(UInt material_index) {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it   = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last = mesh.lastType(spatial_dimension);

  /// loop over element type
  for (; it != last; ++it) {
    ElementType type = *it;
    UInt nb_element = mesh.getNbElement(type);
    UInt nb_quad_per_element = getFEEngine().getNbQuadraturePoints(type);
    const Array<Real> & stress = materials[material_index]->getStress(type);
    Array<Real> & s_on_facet = stress_on_facet(type);

    UInt nb_facet_per_elem = Mesh::getNbFacetsPerElement(type);
    ElementType type_facet = Mesh::getFacetType(type);
    UInt nb_quad_per_facet = getFEEngine("FacetsFEEngine").getNbQuadraturePoints(type_facet);
    UInt nb_facet_quad_per_elem = nb_quad_per_facet * nb_facet_per_elem;

    Array<Real>::const_iterator<Matrix<Real> > stress_it
      = stress.begin(spatial_dimension, spatial_dimension);
    Array<Real>::iterator<Matrix<Real> > s_on_facet_it
      = s_on_facet.begin(spatial_dimension, spatial_dimension);

    Matrix<Real> average_stress(spatial_dimension, spatial_dimension);

    for (UInt el = 0; el < nb_element; ++el) {

      /// compute average stress
      average_stress.clear();

      for (UInt q = 0; q < nb_quad_per_element; ++q, ++stress_it)
	average_stress += *stress_it;

      average_stress /= nb_quad_per_element;

      /// store average stress
      for (UInt q = 0; q < nb_facet_quad_per_elem; ++q, ++s_on_facet_it)
	*s_on_facet_it = average_stress;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::fillStressOnFacet() {
  AKANTU_DEBUG_IN();

  Mesh & mesh_facets = inserter->getMeshFacets();
  UInt sp2 = spatial_dimension * spatial_dimension;
  UInt sp4 = sp2 * 2;

  /// loop over materials
  for (UInt m = 0; m < materials.size(); ++m) {
    try {
      MaterialCohesive & mat __attribute__((unused)) =
	dynamic_cast<MaterialCohesive &>(*materials[m]);
    } catch(std::bad_cast&) {

      if (stress_interpolation)
	/// interpolate stress on facet quadrature points positions
	materials[m]->interpolateStress(stress_on_facet);
      else
	averageStressOnFacets(m);

      GhostType ghost_type = _not_ghost;

      Mesh::type_iterator it   = mesh.firstType(spatial_dimension, ghost_type);
      Mesh::type_iterator last = mesh.lastType(spatial_dimension, ghost_type);

      /// loop over element type
      for (; it != last; ++it) {
	ElementType type = *it;
	UInt nb_element = mesh.getNbElement(type, ghost_type);
	if (nb_element == 0) continue;

	Array<Real> & stress_on_f = stress_on_facet(type, ghost_type);

	/// store the interpolated stresses on the facet_stress vector
	Array<Element> & facet_to_element =
	  mesh_facets.getSubelementToElement(type, ghost_type);

	UInt nb_facet_per_elem = facet_to_element.getNbComponent();

	Array<Element>::iterator<Vector<Element> > facet_to_el_it =
	  facet_to_element.begin(nb_facet_per_elem);

	Array<Real>::iterator< Matrix<Real> > stress_on_f_it =
	  stress_on_f.begin(spatial_dimension, spatial_dimension);

	ElementType facet_type = _not_defined;
	GhostType facet_gt = _casper;
	Array<std::vector<Element> > * element_to_facet = NULL;
	Array<Real> * f_stress = NULL;
	Array<bool> * f_check = NULL;
	UInt nb_quad_per_facet = 0;
	UInt element_rank = 0;

	Element current_el(type, 0, ghost_type);

	for (UInt el = 0; el < nb_element; ++el, ++facet_to_el_it) {
	  current_el.element = el;
	  for (UInt f = 0; f < nb_facet_per_elem; ++f) {
	    Element global_facet_elem = (*facet_to_el_it)(f);
	    UInt global_facet = global_facet_elem.element;

	    if (facet_type != global_facet_elem.type ||
		facet_gt != global_facet_elem.ghost_type) {
	      facet_type = global_facet_elem.type;
	      facet_gt = global_facet_elem.ghost_type;

	      element_to_facet =
		&( mesh_facets.getElementToSubelement(facet_type, facet_gt) );
	      f_stress = &( facet_stress(facet_type, facet_gt) );

	      nb_quad_per_facet =
		getFEEngine("FacetsFEEngine").getNbQuadraturePoints(facet_type, facet_gt);

	      f_check = &( inserter->getCheckFacets(facet_type, facet_gt) );
	    }

	    if (!(*f_check)(global_facet)) {
	      stress_on_f_it += nb_quad_per_facet;
	      continue;
	    }

	    for (UInt q = 0; q < nb_quad_per_facet; ++q, ++stress_on_f_it) {

	      element_rank = (*element_to_facet)(global_facet)[0] != current_el;

	      Matrix<Real> facet_stress_local(f_stress->storage()
					      + (global_facet * nb_quad_per_facet + q) * sp4
					      + element_rank * sp2,
					      spatial_dimension,
					      spatial_dimension);

	      facet_stress_local = *stress_on_f_it;
	    }
	  }
	}
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::reassignMaterial() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::reassignMaterial();

  std::vector< Array<Element> > element_to_add   (materials.size());
  std::vector< Array<Element> > element_to_remove(materials.size());

  Element element;
  for (ghost_type_t::iterator gt = ghost_type_t::begin(); gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;
    element.ghost_type = ghost_type;

    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, ghost_type, _ek_cohesive);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, ghost_type, _ek_cohesive);
    for(; it != end; ++it) {

      ElementType type = *it;
      element.type = type;
      element.kind = Mesh::getKind(type);

      UInt nb_element = mesh.getNbElement(type, ghost_type);

      Array<UInt> & el_index_by_mat = element_index_by_material(type, ghost_type);

      for (UInt el = 0; el < nb_element; ++el) {
	element.element = el;

	UInt old_material = el_index_by_mat(el, 0);
	UInt new_material = (*material_selector)(element);

	if(old_material != new_material) {
	  element_to_add   [new_material].push_back(element);
	  element_to_remove[old_material].push_back(element);
	}
      }
    }
  }

  std::vector<Material *>::iterator mat_it;
  UInt mat_index = 0;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it, ++mat_index) {
    (*mat_it)->removeElements(element_to_remove[mat_index]);
    (*mat_it)->addElements   (element_to_add[mat_index]);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::checkCohesiveStress() {
  AKANTU_DEBUG_IN();

  fillStressOnFacet();

#if defined(AKANTU_DEBUG_TOOLS)
  debug::element_manager.printData(debug::_dm_model_cohesive,
				   "Interpolated stresses before",
				   facet_stress);
#endif

  synch_registry->synchronize(_gst_smmc_facets_stress);

#if defined(AKANTU_DEBUG_TOOLS)
  debug::element_manager.printData(debug::_dm_model_cohesive,
				   "Interpolated stresses",
				   facet_stress);
#endif

  for (UInt m = 0; m < materials.size(); ++m) {
    try {
      MaterialCohesive & mat_cohesive = dynamic_cast<MaterialCohesive &>(*materials[m]);
      /// check which not ghost cohesive elements are to be created
      mat_cohesive.checkInsertion();
    } catch(std::bad_cast&) { }
  }

  /// communicate data among processors
  synch_registry->synchronize(_gst_smmc_facets);

  /// insert cohesive elements
  inserter->insertElements();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::onElementsAdded(const Array<Element> & element_list,
						  const NewElementsEvent & event) {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::onElementsAdded(element_list, event);

  /// update shape functions
  getFEEngine("CohesiveFEEngine").initShapeFunctions(_not_ghost);
  getFEEngine("CohesiveFEEngine").initShapeFunctions(_ghost);

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  updateCohesiveSynchronizers();
#endif

  if (is_extrinsic) resizeFacetStress();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::onNodesAdded(const Array<UInt> & doubled_nodes,
					       __attribute__((unused)) const NewNodesEvent & event) {
  AKANTU_DEBUG_IN();

  UInt nb_new_nodes = doubled_nodes.getSize();
  Array<UInt> nodes_list(nb_new_nodes);

  for (UInt n = 0; n < nb_new_nodes; ++n)
    nodes_list(n) = doubled_nodes(n, 1);

  SolidMechanicsModel::onNodesAdded(nodes_list, event);

  for (UInt n = 0; n < nb_new_nodes; ++n) {

    UInt old_node = doubled_nodes(n, 0);
    UInt new_node = doubled_nodes(n, 1);

    for (UInt dim = 0; dim < spatial_dimension; ++dim) {
      (*displacement)(new_node, dim) = (*displacement)(old_node, dim);
      (*velocity)    (new_node, dim) = (*velocity)    (old_node, dim);
      (*acceleration)(new_node, dim) = (*acceleration)(old_node, dim);
      (*blocked_dofs)(new_node, dim) = (*blocked_dofs)(old_node, dim);

      if (current_position)
	(*current_position)(new_node, dim) = (*current_position)(old_node, dim);

      if (increment_acceleration)
	(*increment_acceleration)(new_node, dim)
	  = (*increment_acceleration)(old_node, dim);

      if (increment)
	(*increment)(new_node, dim) = (*increment)(old_node, dim);

      if (previous_displacement)
	(*previous_displacement)(new_node, dim)
	  = (*previous_displacement)(old_node, dim);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::onEndSolveStep(const AnalysisMethod & method) {

  AKANTU_DEBUG_IN();

  /******************************************************************************
    This is required because the Cauchy stress is the stress measure that is used
     to check the insertion of cohesive elements
  ******************************************************************************/

  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    if(mat.isFiniteDeformation())
      mat.computeAllCauchyStresses(_not_ghost);
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "SolidMechanicsModelCohesive [" << std::endl;

  SolidMechanicsModel::printself(stream, indent + 1);

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::resizeFacetStress() {
  AKANTU_DEBUG_IN();

  Mesh & mesh_facets = inserter->getMeshFacets();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {
    GhostType ghost_type = *gt;

    Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, ghost_type);
    Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, ghost_type);
    for(; it != end; ++it) {
      ElementType type = *it;

      UInt nb_facet = mesh_facets.getNbElement(type, ghost_type);

      UInt nb_quadrature_points =
	getFEEngine("FacetsFEEngine").getNbQuadraturePoints(type, ghost_type);

      UInt new_size = nb_facet * nb_quadrature_points;

      facet_stress(type, ghost_type).resize(new_size);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::addDumpGroupFieldToDumper(const std::string & dumper_name,
							    const std::string & field_id,
							    const std::string & group_name,
							    const ElementKind & element_kind,
							    bool padding_flag) {
  AKANTU_DEBUG_IN();

  ElementKind _element_kind = element_kind;
  if (dumper_name == "cohesive elements") {
    _element_kind = _ek_cohesive;
  }
  SolidMechanicsModel::addDumpGroupFieldToDumper(dumper_name, 
						 field_id,
						 group_name,
						 _element_kind,
						 padding_flag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::onDump(){
  this->flattenAllRegisteredInternals(_ek_cohesive);
  SolidMechanicsModel::onDump();
}

/* -------------------------------------------------------------------------- */



__END_AKANTU__
