/**
 * @file   embedded_interface_model.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Mon Mar 9 2015
 * @date last modification: Mon Mar 9 2015
 *
 * @brief  Model of Solid Mechanics with embedded interfaces
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

#include "aka_common.hh"

#include "embedded_interface_model.hh"
#include "material_reinforcement.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

EmbeddedInterfaceModel::EmbeddedInterfaceModel(Mesh & mesh,
                                               UInt spatial_dimension,
                                               const ID & id,
                                               const MemoryID & memory_id) :
  SolidMechanicsModel(mesh, spatial_dimension, id, memory_id),
  interface_mesh(NULL),
  interface_container(mesh)
{}

EmbeddedInterfaceModel::~EmbeddedInterfaceModel()
{}

void EmbeddedInterfaceModel::initModel() {
  SolidMechanicsModel::initModel();

  registerFEEngineObject<MyFEEngineType>("EmbeddedInterfaceFEEngine", *interface_mesh, spatial_dimension);

}

void EmbeddedInterfaceModel::initMaterials() {
  SolidMechanicsModel::initMaterials();
}

void EmbeddedInterfaceModel::initInterface(const std::list<Interface> & interface_list) {

}

void EmbeddedInterfaceModel::assembleStiffnessMatrix() {

}


__END_AKANTU__

