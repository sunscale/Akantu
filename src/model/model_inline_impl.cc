/**
 * @file   model_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 25 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  inline implementation of the model class
 *
 * @section LICENSE
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
#include "model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MODEL_INLINE_IMPL_CC__
#define __AKANTU_MODEL_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename FEEngineClass>
inline FEEngineClass & Model::getFEEngineClassBoundary(std::string name) {
  if (name == "")
    name = default_fem;

  auto it_boun = fems_boundary.find(name);

  if (it_boun == fems_boundary.end()) {
    AKANTU_DEBUG_INFO("Creating FEEngine boundary " << name);

    auto it = fems.find(name);
    if (it == fems.end()) {
      AKANTU_EXCEPTION("The FEEngine " << name << " is not registered");
    }

    auto spatial_dimension = it->second->getElementDimension();
    fems_boundary[name] = std::make_unique<FEEngineClass>(
        it->second->getMesh(), spatial_dimension - 1,
        id + ":fem_boundary:" + name, memory_id);
  }

  return aka::as_type<FEEngineClass>(*fems_boundary[name]);
}

/* -------------------------------------------------------------------------- */
template <typename FEEngineClass>
inline FEEngineClass & Model::getFEEngineClass(std::string name) const {
  if (name == "")
    name = default_fem;

  auto it = fems.find(name);
  if (it == fems.end()) {
    AKANTU_EXCEPTION("The FEEngine " << name << " is not registered");
  }

  return aka::as_type<FEEngineClass>(*(it->second));
}

/* -------------------------------------------------------------------------- */
inline void Model::unRegisterFEEngineObject(const std::string & name) {
  auto it = fems.find(name);
  if (it == fems.end()) {
    AKANTU_EXCEPTION("FEEngine object with name " << name << " was not found");
  }

  fems.erase(it);
  if (not fems.empty() and default_fem == name)
    default_fem = (*fems.begin()).first;
}

/* -------------------------------------------------------------------------- */
template <typename FEEngineClass>
inline void Model::registerFEEngineObject(const std::string & name, Mesh & mesh,
                                          UInt spatial_dimension) {
  if (fems.size() == 0)
    default_fem = name;

  auto it = fems.find(name);
  if (it != fems.end()) {
    AKANTU_EXCEPTION("FEEngine object with name " << name
                                                  << " was already created");
  }

  fems[name] = std::make_unique<FEEngineClass>(
      mesh, spatial_dimension, id + ":fem:" + name + std::to_string(memory_id),
      memory_id);
}

/* -------------------------------------------------------------------------- */
inline FEEngine & Model::getFEEngine(const ID & name) const {
  ID tmp_name = (name == "") ? default_fem : name;

  auto it = fems.find(tmp_name);

  if (it == fems.end()) {
    AKANTU_EXCEPTION("The FEEngine " << tmp_name << " is not registered");
  }
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
inline FEEngine & Model::getFEEngineBoundary(const ID & name) {
  ID tmp_name = (name == "") ? default_fem : name;

  auto it = fems_boundary.find(tmp_name);
  if (it == fems_boundary.end()) {
    AKANTU_EXCEPTION("The FEEngine boundary  " << tmp_name
                                               << " is not registered");
  }
  AKANTU_DEBUG_ASSERT(it->second != nullptr, "The FEEngine boundary "
                                                 << tmp_name
                                                 << " was not created");
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Model::allocNodalField(Array<T> *& array, UInt nb_component,
                            const ID & name) {
  if (array)
    return;

  UInt nb_nodes = mesh.getNbNodes();
  array = &(alloc<T>(id + ":" + name, nb_nodes, nb_component, T()));
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Model::allocNodalField(std::unique_ptr<Array<T>> & array,
                            UInt nb_component, const ID & name) const {
  if (array)
    return;

  UInt nb_nodes = mesh.getNbNodes();
  array =
      std::make_unique<Array<T>>(nb_nodes, nb_component, T(), id + ":" + name);
}

/* -------------------------------------------------------------------------- */
inline UInt Model::getNbIntegrationPoints(const Array<Element> & elements,
                                          const ID & fem_id) const {
  UInt nb_quad = 0;
  for (auto && el : elements) {
    nb_quad +=
        getFEEngine(fem_id).getNbIntegrationPoints(el.type, el.ghost_type);
  }
  return nb_quad;
}

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* __AKANTU_MODEL_INLINE_IMPL_CC__ */
