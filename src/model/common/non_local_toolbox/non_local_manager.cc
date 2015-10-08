/**
 * @file   non_local_manager.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Sep 21 15:32:10 2015
 *
 * @brief  Implementation of non-local manager
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "non_local_manager.hh"
#include "non_local_neighborhood.hh"
#include "material_non_local.hh"
#include "base_weight_function.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
NonLocalManager::NonLocalManager(SolidMechanicsModel & model, 
				const ID & id,
				const MemoryID & memory_id) : 
  Memory(id, memory_id),
  Parsable(_st_neighborhoods, id),
  model(model),
  quad_positions("quad_positions", id),
  volumes("volumes", id),
  spatial_dimension(this->model.getSpatialDimension()),
  compute_stress_calls(0){ 
  Mesh & mesh = this->model.getMesh();

  /// initialize the element type map arrays
  mesh.initElementTypeMapArray(quad_positions, spatial_dimension, spatial_dimension, false, _ek_regular, true);

  /// parse the neighborhood information from the input file
  const Parser & parser = getStaticParser();

  /// iterate over all the non-local sections and store them in a map
  std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
    weight_sect = parser.getSubSections(_st_non_local);
  Parser::const_section_iterator it = weight_sect.first;
  for (; it != weight_sect.second; ++it) {
    const ParserSection & section = *it;
    ID name = section.getName();
    this->weight_function_types[name] = section;    
  }
}

/* -------------------------------------------------------------------------- */
NonLocalManager::~NonLocalManager() {

  /// delete neighborhoods
  NeighborhoodMap::iterator it;
  for (it = neighborhoods.begin(); it != neighborhoods.end(); ++it) {
    if(it->second) delete it->second;
  }

  /// delete non-local variables
  std::map<ID, NonLocalVariable *>::iterator it_variables;
  for (it_variables = non_local_variables.begin(); it_variables != non_local_variables.end(); ++it_variables) {
    if(it_variables->second) delete it_variables->second;
  }

  std::map<ID, ElementTypeMapReal *>::iterator it_internals;
  for (it_internals = weight_function_internals.begin(); it_internals != weight_function_internals.end(); ++it_internals) {
    if(it_internals->second) delete it_internals->second;
  }

}

/* -------------------------------------------------------------------------- */
void NonLocalManager::setJacobians(const FEEngine & fe_engine, const ElementKind & kind) {
  Mesh & mesh = this->model.getMesh();
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, gt, kind);
    Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, gt, kind);
    for(; it != last_type; ++it) {
      jacobians(*it, gt) = &fe_engine.getIntegratorInterface().getJacobians(*it, gt);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::createNeighborhood(const ID & weight_func, const ID & neighborhood_id) {

  AKANTU_DEBUG_IN();

  
  const ParserSection & section = this->weight_function_types[weight_func];
  const ID weight_func_type = section.getOption();
  /// create new neighborhood for given ID
  std::stringstream sstr; sstr << id << ":neighborhood:" << neighborhood_id;

  if (weight_func_type == "base_wf")
    neighborhoods[neighborhood_id] = new NonLocalNeighborhood<BaseWeightFunction>(*this, this->quad_positions, sstr.str());
  else if (weight_func_type == "remove_wf")
    neighborhoods[neighborhood_id] = new NonLocalNeighborhood<RemoveDamagedWeightFunction>(*this, this->quad_positions, sstr.str());                                              
  else if (weight_func_type == "stress_wf")
    neighborhoods[neighborhood_id] = new NonLocalNeighborhood<StressBasedWeightFunction>(*this, this->quad_positions, sstr.str());
  else if (weight_func_type == "damage_wf")
    neighborhoods[neighborhood_id] = new NonLocalNeighborhood<DamagedWeightFunction>(*this, this->quad_positions, sstr.str());
  else
    AKANTU_EXCEPTION("error in weight function type provided in material file");

  neighborhoods[neighborhood_id]->parseSection(section);
  neighborhoods[neighborhood_id]->initNeighborhood();

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
void NonLocalManager::createNeighborhoodSynchronizers() {

  NeighborhoodMap::iterator it;
  for (it = neighborhoods.begin(); it != neighborhoods.end(); ++it) {
    it->second->createGridSynchronizer();
  }

}

/* -------------------------------------------------------------------------- */
void NonLocalManager::flattenInternal(ElementTypeMapReal & internal_flat,
				      const GhostType & ghost_type,
				      const ElementKind & kind) {
  
  const ID field_name = internal_flat.getName();
  for (UInt m = 0; m < this->non_local_materials.size(); ++m) {
    Material & material = *(this->non_local_materials[m]);
    if (material.isInternal(field_name, kind)) 
      material.flattenInternal(field_name, internal_flat, ghost_type, kind);
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::averageInternals(const GhostType & ghost_type) {
  /// update the weights of the weight function
  if (ghost_type == _not_ghost) 
    this->computeWeights();

  /// loop over all neighborhoods and compute the non-local variables
  NeighborhoodMap::iterator neighborhood_it = neighborhoods.begin();
  NeighborhoodMap::iterator neighborhood_end = neighborhoods.end(); 
  for (; neighborhood_it != neighborhood_end; ++neighborhood_it) {
    /// loop over all the non-local variables of the given neighborhood
    std::map<ID, NonLocalVariable *>::iterator non_local_variable_it = non_local_variables.begin();
    std::map<ID, NonLocalVariable *>::iterator non_local_variable_end = non_local_variables.end();
    for(; non_local_variable_it != non_local_variable_end; ++non_local_variable_it) {
      NonLocalVariable * non_local_var = non_local_variable_it->second;
      neighborhood_it->second->weightedAverageOnNeighbours(non_local_var->local, non_local_var->non_local, non_local_var->nb_component, ghost_type);

    }
  }    
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::init(){

  /// exchange the missing ghosts for the non-local neighborhoods
  this->createNeighborhoodSynchronizers();

  /// insert the ghost quadrature points of the non-local materials into the non-local neighborhoods
  for(UInt m = 0; m < this->non_local_materials.size(); ++m) {
    switch (spatial_dimension) {
    case 1:
      dynamic_cast<MaterialNonLocal<1> &>(*(this->non_local_materials[m])).insertQuadsInNeighborhoods(_ghost); 
    break;      
    case 2:
      dynamic_cast<MaterialNonLocal<2> &>(*(this->non_local_materials[m])).insertQuadsInNeighborhoods(_ghost); 
    break;      
    case 3:
      dynamic_cast<MaterialNonLocal<3> &>(*(this->non_local_materials[m])).insertQuadsInNeighborhoods(_ghost); 
    break;      
    }

  }

  this->setJacobians(this->model.getFEEngine(), _ek_regular);
  this->initNonLocalVariables();
  this->updatePairLists();
  this->initElementTypeMap(1, volumes);
  this->computeWeights();
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::initNonLocalVariables(){
  /// loop over all the non-local variables
  std::map<ID, NonLocalVariable *>::iterator non_local_variable_it = non_local_variables.begin();
  std::map<ID, NonLocalVariable *>::iterator non_local_variable_end = non_local_variables.end();

  for(; non_local_variable_it != non_local_variable_end; ++non_local_variable_it) {
    NonLocalVariable & variable = *(non_local_variable_it->second);
    this->initElementTypeMap(variable.nb_component, variable.non_local);      
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::initElementTypeMap(UInt nb_component, ElementTypeMapReal & element_map) {
  Mesh & mesh = this->model.getMesh();
  /// need to resize the arrays
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt, _ek_regular);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, _ek_regular);
    for(; it != end; ++it) {
      ElementType el_type = *it;
      UInt nb_element = mesh.getNbElement(*it, gt);
      UInt nb_quads = this->model.getFEEngine().getNbQuadraturePoints(*it, gt);
      if (!element_map.exists(el_type, gt)) {
	element_map.alloc(nb_element * nb_quads,
			  nb_component,
			  el_type,
			  gt);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::distributeInternals(ElementKind kind) {

  /// loop over all the non-local variables and copy back their values into the materials
  std::map<ID, NonLocalVariable *>::iterator non_local_variable_it = non_local_variables.begin();
  std::map<ID, NonLocalVariable *>::iterator non_local_variable_end = non_local_variables.end();
  for(; non_local_variable_it != non_local_variable_end; ++non_local_variable_it) {
    NonLocalVariable * non_local_var = non_local_variable_it->second;
    const ID field_name = non_local_var->non_local.getName();
    /// loop over all the materials
    for (UInt m = 0; m < this->non_local_materials.size(); ++m) {
      if (this->non_local_materials[m]->isInternal(field_name, kind))

	switch (spatial_dimension) {
	case 1:
	  dynamic_cast<MaterialNonLocal<1> &>(*(this->non_local_materials[m])).updateNonLocalInternals(non_local_var->non_local, field_name, non_local_var->nb_component); 
	  break;      
	case 2:
	  dynamic_cast<MaterialNonLocal<2> &>(*(this->non_local_materials[m])).updateNonLocalInternals(non_local_var->non_local, field_name, non_local_var->nb_component);
	  break;      
	case 3:
	  dynamic_cast<MaterialNonLocal<3> &>(*(this->non_local_materials[m])).updateNonLocalInternals(non_local_var->non_local, field_name, non_local_var->nb_component);
	  break;      
	}
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::computeAllNonLocalStresses() {

  /// update the flattened version of the internals
  std::map<ID, NonLocalVariable *>::iterator non_local_variable_it = non_local_variables.begin();
  std::map<ID, NonLocalVariable *>::iterator non_local_variable_end = non_local_variables.end();

  for(; non_local_variable_it != non_local_variable_end; ++non_local_variable_it) {
    non_local_variable_it->second->local.clear();
    non_local_variable_it->second->non_local.clear();
    for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
      GhostType ghost_type = (GhostType) gt;
      this->flattenInternal(non_local_variable_it->second->local, ghost_type, _ek_regular);
    }
  }

  this->volumes.clear();  

  /// loop over all neighborhoods and compute the non-local variables
  NeighborhoodMap::iterator neighborhood_it = neighborhoods.begin();
  NeighborhoodMap::iterator neighborhood_end = neighborhoods.end(); 
  for (; neighborhood_it != neighborhood_end; ++neighborhood_it) {
    neighborhood_it->second->getSynchronizerRegistry().asynchronousSynchronize(_gst_mnl_for_average);
  }

  this->averageInternals(_not_ghost);

  AKANTU_DEBUG_INFO("Wait distant non local stresses");

  /// loop over all neighborhoods and compute the non-local variables
  neighborhood_it = neighborhoods.begin();
  for (; neighborhood_it != neighborhood_end; ++neighborhood_it) {
    neighborhood_it->second->getSynchronizerRegistry().waitEndSynchronize(_gst_mnl_for_average);
  }

  this->averageInternals(_ghost);

  /// copy the results in the materials
  this->distributeInternals(_ek_regular);
    /// loop over all the materials and update the weights
    for (UInt m = 0; m < this->non_local_materials.size(); ++m) {
	switch (spatial_dimension) {
	case 1:
	  dynamic_cast<MaterialNonLocal<1> &>(*(this->non_local_materials[m])).computeNonLocalStresses(_not_ghost); break;      
	case 2:
	  dynamic_cast<MaterialNonLocal<2> &>(*(this->non_local_materials[m])).computeNonLocalStresses(_not_ghost); break;      
	case 3:
	  dynamic_cast<MaterialNonLocal<3> &>(*(this->non_local_materials[m])).computeNonLocalStresses(_not_ghost); break;      
	}
    }
    
    ++this->compute_stress_calls;
}


__END_AKANTU__
