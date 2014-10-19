/**
 * @file   model_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 25 2010
 * @date last modification: Tue Jul 29 2014
 *
 * @brief  inline implementation of the model class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
inline SynchronizerRegistry & Model::getSynchronizerRegistry(){
  AKANTU_DEBUG_ASSERT(synch_registry,"synchronizer registry not initialized:"
		      << " did you call createSynchronizerRegistry ?");
  return  *synch_registry;
}

/* -------------------------------------------------------------------------- */
template <typename FEEngineClass>
inline FEEngineClass & Model::getFEEngineClassBoundary(std::string name) {
  AKANTU_DEBUG_IN();

  if (name == "") name = default_fem;

  FEEngineMap::const_iterator it_boun = fems_boundary.find(name);

  FEEngineClass * tmp_fem_boundary;

  if (it_boun == fems_boundary.end()){
    AKANTU_DEBUG_INFO("Creating FEEngine boundary " << name);

    FEEngineMap::const_iterator it = fems.find(name);
    AKANTU_DEBUG_ASSERT(it != fems.end(), "The FEEngine " << name << " is not registered");

    UInt spatial_dimension = it->second->getElementDimension();
    std::stringstream sstr; sstr << id << ":fem_boundary:" << name;

    tmp_fem_boundary = new FEEngineClass(it->second->getMesh(),
				    spatial_dimension-1,
				    sstr.str(),
				    memory_id);
    fems_boundary[name] = tmp_fem_boundary;
  } else {
    tmp_fem_boundary = dynamic_cast<FEEngineClass *>(it_boun->second);
  }

  AKANTU_DEBUG_OUT();
  return *tmp_fem_boundary;
}


/* -------------------------------------------------------------------------- */
template <typename FEEngineClass>
inline FEEngineClass & Model::getFEEngineClass(std::string name) const{
  AKANTU_DEBUG_IN();

  if (name == "") name = default_fem;

  FEEngineMap::const_iterator it = fems.find(name);
  AKANTU_DEBUG_ASSERT(it != fems.end(), "The FEEngine " << name << " is not registered");

  AKANTU_DEBUG_OUT();
  return dynamic_cast<FEEngineClass &>(*(it->second));
}

/* -------------------------------------------------------------------------- */

inline void Model::unRegisterFEEngineObject(const std::string & name){

  FEEngineMap::iterator it = fems.find(name);
  AKANTU_DEBUG_ASSERT(it != fems.end(), "FEEngine object with name "
		      << name << " was not found");

  delete((*it).second);
  fems.erase(it);
  if (!fems.empty())
    default_fem = (*fems.begin()).first;
}

/* -------------------------------------------------------------------------- */

template <typename FEEngineClass>
inline void Model::registerFEEngineObject(const std::string & name,
				     Mesh & mesh,
				     UInt spatial_dimension){
  if (fems.size() == 0) default_fem = name;

#ifndef AKANTU_NDEBUG
  FEEngineMap::iterator it = fems.find(name);
  AKANTU_DEBUG_ASSERT(it == fems.end(), "FEEngine object with name "
		      << name << " was already created");
#endif

  std::stringstream sstr; sstr << id << ":fem:" << name;
  fems[name] = new FEEngineClass(mesh, spatial_dimension, sstr.str(), memory_id);

  // MeshUtils::buildFacets(fems[name]->getMesh());

  // std::stringstream sstr2; sstr2 << id << ":fem_boundary:" << name;
  // fems_boundary[name] = new FEEngineClass(mesh, spatial_dimension-1, sstr2.str(), memory_id);
}

/* -------------------------------------------------------------------------- */
inline FEEngine & Model::getFEEngine(const ID & name) const{
  AKANTU_DEBUG_IN();
  ID tmp_name = name;
  if (name == "") tmp_name = default_fem;

  FEEngineMap::const_iterator it = fems.find(tmp_name);

  AKANTU_DEBUG_ASSERT(it != fems.end(),
		      "The FEEngine " << tmp_name << " is not registered");

  AKANTU_DEBUG_OUT();
  return *(it->second);
}


/* -------------------------------------------------------------------------- */
inline FEEngine & Model::getFEEngineBoundary(const ID & name){
  AKANTU_DEBUG_IN();

  ID tmp_name = name;
  if (name == "") tmp_name = default_fem;

  FEEngineMap::const_iterator it = fems_boundary.find(tmp_name);
  AKANTU_DEBUG_ASSERT(it != fems_boundary.end(),
		      "The FEEngine boundary  " << tmp_name << " is not registered");
  AKANTU_DEBUG_ASSERT(it->second != NULL,
		      "The FEEngine boundary " << tmp_name << " was not created");

  AKANTU_DEBUG_OUT();
  return *(it->second);
}
/* -------------------------------------------------------------------------- */
/// @todo : should merge with a single function which handles local and global
inline void Model::changeLocalEquationNumberForPBC(std::map<UInt,UInt> & pbc_pair,
					    UInt dimension){
  for (std::map<UInt,UInt>::iterator it = pbc_pair.begin();
       it != pbc_pair.end();++it) {
    Int node_master = (*it).second;
    Int node_slave = (*it).first;
    for (UInt i = 0; i < dimension; ++i) {
      (*dof_synchronizer->getLocalDOFEquationNumbersPointer())
	(node_slave*dimension+i) = dimension*node_master+i;
      (*dof_synchronizer->getGlobalDOFEquationNumbersPointer())
	(node_slave*dimension+i) = dimension*node_master+i;
    }
  }
}
/* -------------------------------------------------------------------------- */
inline bool Model::isPBCSlaveNode(const UInt node) const {
  // if no pbc is defined, is_pbc_slave_node is of size zero
  if (is_pbc_slave_node.getSize() == 0)
    return false;
  else
    return is_pbc_slave_node(node);
}

/* -------------------------------------------------------------------------- */
inline UInt Model::getNbQuadraturePoints(const Array<Element> & elements,
					 const ID & fem_id) const {
  UInt nb_quad = 0;
  Array<Element>::const_iterator<Element> it  = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    nb_quad += getFEEngine(fem_id).getNbQuadraturePoints(el.type,
						    el.ghost_type);
  }
  return nb_quad;
}

/* -------------------------------------------------------------------------- */

