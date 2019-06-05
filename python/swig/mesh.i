/**
 * @file   mesh.i
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Dec 12 2014
 * @date last modification: Wed Jan 13 2016
 *
 * @brief  mesh wrapper
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
%{
#include "mesh.hh"
#include "node_group.hh"
#include "solid_mechanics_model.hh"
#include "python_functor.hh"
#include "mesh_utils.hh"
#include "aka_bbox.hh"
#include "mesh_accessor.hh"
#include "communicator.hh"

using akantu::IntegrationPoint;
using akantu::Vector;
using akantu::ElementTypeMapArray;
using akantu::MatrixProxy;
using akantu::Matrix;
using akantu::UInt;
using akantu::Real;
using akantu::Array;
using akantu::BBox;
using akantu::Communicator;
using akantu::SolidMechanicsModel;
%}


namespace akantu {
  %ignore NewNodesEvent;
  %ignore RemovedNodesEvent;
  %ignore NewElementsEvent;
  %ignore RemovedElementsEvent;
  %ignore MeshEventHandler;
  %ignore MeshEvent< UInt >;
  %ignore MeshEvent< Element >;
  %ignore Mesh::extractNodalCoordinatesFromPBCElement;
  %ignore Mesh::getGroupDumer;
  %ignore Mesh::getFacetLocalConnectivity;
  %ignore Mesh::getAllFacetTypes;
  %ignore Mesh::getCommunicator;
  %ignore Mesh::getConnectivities;
  %ignode Mesh::getBBox;
  %ignore GroupManager::getElementGroups;
  %ignore Dumpable::addDumpFieldExternalReal;
}

print_self(Mesh)

// Swig considers enums to be ints, and it creates a conflict with two versions of getNbElement()
%rename(getNbElementByDimension)
akantu::Mesh::getNbElement(const UInt spatial_dimension = _all_dimensions,
                           const GhostType& ghost_type = _not_ghost,
                           const ElementKind& kind = _ek_not_defined) const;

%extend akantu::Mesh {
  PyObject * getElementGroups(){
    return akantu::PythonFunctor::convertToPython($self->getElementGroups());
  }

  PyObject * getAllConnectivities(){
    return akantu::PythonFunctor::convertToPython($self->getConnectivities());
  }

  void resizeMesh(UInt nb_nodes, UInt nb_element, const ElementType & type) {
    Array<Real> & nodes = const_cast<Array<Real> &>($self->getNodes());
    nodes.resize(nb_nodes);

    $self->addConnectivityType(type);
    Array<UInt> & connectivity = const_cast<Array<UInt> &>($self->getConnectivity(type));
    connectivity.resize(nb_element);
  }

  Array<Real> & getNodalDataReal(const ID & name, UInt nb_components = 1) {
    auto && data = $self->getNodalData<Real>(name, nb_components);
    data.resize($self->getNbNodes());
    return data;
  }

  bool hasDataReal(const ID & name,
                   const ElementType & type) {
    return $self->hasData<Real>(name, type);
  }

  Array<Real> & getElementalDataReal(const ID & name,
                                     const ElementType & type,
                                     UInt nb_components = 1) {
    auto && data = $self->getElementalDataArrayAlloc<Real>(name, type,
                                                           akantu::_not_ghost,
                                                           nb_components);
    data.resize($self->getNbElement(type, akantu::_not_ghost));
    return data;
  }

  Array<UInt> & getElementalDataUInt(const ID & name,
                                     const ElementType & type,
                                     UInt nb_components = 1) {
    auto && data = $self->getElementalDataArrayAlloc<akantu::UInt>(name, type,
                                                           akantu::_not_ghost,
                                                           nb_components);
    data.resize($self->getNbElement(type, akantu::_not_ghost));
    return data;
  }

  Array<Real> & computeBarycenters(const ElementType & type) {
    auto dim = $self->getSpatialDimension();
    auto && data = $self->getElementalDataArrayAlloc<akantu::Real>("barycenters", type,
                                                           akantu::_not_ghost, dim);
    auto nb_el = data.size();
    auto total_nb_el = $self->getNbElement(type, akantu::_not_ghost);

    data.resize(total_nb_el);

    auto bary_it = make_view(data, dim).begin() + nb_el;
    for (auto el = nb_el; el < total_nb_el; ++el) {
      $self->getBarycenter(akantu::Element{type, el, akantu::_not_ghost},
                           *bary_it);
      ++bary_it;
    }
    return data;
  }

  void ready() {
    akantu::MeshAccessor ma(* $self);
    ma.makeReady();
  }
}

%extend akantu::GroupManager {
  void createGroupsFromStringMeshData(const std::string & dataset_name) {
    if (dataset_name == "physical_names"){
      AKANTU_EXCEPTION("Deprecated behavior: no need to call 'createGroupsFromStringMeshData' for physical names");
    }
    $self->createGroupsFromMeshData<std::string>(dataset_name);
  }

  void createGroupsFromUIntMeshData(const std::string & dataset_name) {
    $self->createGroupsFromMeshData<akantu::UInt>(dataset_name);
  }
}

%extend akantu::NodeGroup {
    akantu::Array<akantu::Real> & getGroupedNodes(akantu::Array<akantu::Real, true> & surface_array, Mesh & mesh) {
    auto && group_node = $self->getNodes();
    auto && full_array = mesh.getNodes();
    surface_array.resize(group_node.size());

    for (UInt i = 0; i < group_node.size(); ++i) {
      for (UInt cmp = 0; cmp < full_array.getNbComponent(); ++cmp) {
        surface_array(i, cmp) = full_array(group_node(i), cmp);
      }
    }

    akantu::Array<akantu::Real> & res(surface_array);
    return res;
  }

  akantu::Array<akantu::Real> & getGroupedArray(akantu::Array<akantu::Real, true> & surface_array,
                                                akantu::SolidMechanicsModel & model, int type) {
    akantu::Array<akantu::Real> * full_array;

    switch (type) {

    case 0 : full_array = new akantu::Array<akantu::Real>(model.getDisplacement());
      break;
    case 1 : full_array = new akantu::Array<akantu::Real>(model.getVelocity());
      break;
    case 2 : full_array = new akantu::Array<akantu::Real>(model.getExternalForce());
      break;
    }
    akantu::Array<akantu::UInt> group_node = $self->getNodes();
    surface_array.resize(group_node.size());

    for (UInt i = 0; i < group_node.size(); ++i) {
      for (UInt cmp = 0; cmp < full_array->getNbComponent(); ++cmp) {

        surface_array(i,cmp) = (*full_array)(group_node(i),cmp);
      }
    }

    akantu::Array<akantu::Real> & res(surface_array);
    return res;
  }
}

%include "group_manager.hh"
%include "node_group.hh"
%include "dumper_iohelper.hh"
%include "dumpable_iohelper.hh"
%include "element_group.hh"
%include "mesh.hh"
%include "mesh_utils.hh"
%include "aka_bbox.hh"

namespace akantu{
%extend Dumpable {
    void addDumpFieldExternalReal(const std::string & field_id,
                                  const Array<Real> & field){
      $self->addDumpFieldExternal<Real>(field_id,field);
    }
  }
 }
