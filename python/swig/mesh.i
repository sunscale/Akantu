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
//#include "dumpable_inline_impl.hh"

using akantu::IntegrationPoint;
using akantu::Vector;
using akantu::ElementTypeMapArray;
using akantu::MatrixProxy;
using akantu::Matrix;
using akantu::UInt;
using akantu::Real;
using akantu::Array;
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
  %ignore GroupManager::getElementGroups;
}

print_self(Mesh)

// Swig considers enums to be ints, and it creates a conflict with two versions of getNbElement()
%rename(getNbElementByDimension)
akantu::Mesh::getNbElement(const UInt spatial_dimension = _all_dimensions,
                           const GhostType& ghost_type = _not_ghost,
                           const ElementKind& kind = _ek_not_defined) const;

%extend akantu::Mesh {

  void resizeMesh(UInt nb_nodes, UInt nb_element, const ElementType & type) {
    Array<Real> & nodes = const_cast<Array<Real> &>($self->getNodes());
    nodes.resize(nb_nodes);

    $self->addConnectivityType(type);
    Array<UInt> & connectivity = const_cast<Array<UInt> &>($self->getConnectivity(type));
    connectivity.resize(nb_element);
  }

#if defined(AKANTU_COHESIVE_ELEMENT)
  Array<Real> & getCohesiveBarycenter(SpacialDirection dir) {
    UInt spatial_dimension = $self->getSpatialDimension();
    ElementTypeMapArray<Real> & barycenter =
        $self->registerData<Real>("barycenter");
    $self->initElementTypeMapArray(barycenter, 1, spatial_dimension, false,
                                   akantu::_ek_cohesive, true);
    akantu::ElementType type = *($self->firstType(
        spatial_dimension, akantu::_not_ghost, akantu::_ek_cohesive));

    Vector<Real> bary(spatial_dimension);
    Array<Real> & bary_coh = barycenter(type);
    for (UInt i = 0; i < $self->getNbElement(type); ++i) {
      bary.clear();
      $self->getBarycenter(i, type, bary.storage());
      bary_coh(i) = bary(dir);
    }
    return bary_coh;
  }
#endif
}

%extend akantu::GroupManager {
  void createGroupsFromStringMeshData(const std::string & dataset_name) {
    $self->createGroupsFromMeshData<std::string>(dataset_name);
  }

  void createGroupsFromUIntMeshData(const std::string & dataset_name) {
    $self->createGroupsFromMeshData<akantu::UInt>(dataset_name);
  }
}

%extend akantu::NodeGroup {
  akantu::Array<akantu::Real> & getGroupedNodes(akantu::Array<akantu::Real, true> & surface_array, Mesh & mesh) {
    akantu::Array<akantu::UInt> group_node = $self->getNodes();
    akantu::Array<akantu::Real> & full_array = mesh.getNodes();
    surface_array.resize(group_node.size());

    for (UInt i = 0; i < group_node.size(); ++i) {
      for (UInt cmp = 0; cmp < full_array.getNbComponent(); ++cmp) {

        surface_array(i,cmp) = full_array(group_node(i),cmp);
      }
    }

    akantu::Array<akantu::Real> & res(surface_array);
    return res;
  }

  akantu::Array<akantu::Real> & getGroupedArray(akantu::Array<akantu::Real, true> & surface_array, akantu::SolidMechanicsModel & model, int type) {
    akantu::Array<akantu::Real> * full_array;

    switch (type) {

    case 0 : full_array = new akantu::Array<akantu::Real>(model.getDisplacement());
      break;
    case 1 : full_array = new akantu::Array<akantu::Real>(model.getVelocity());
      break;
    case 2 : full_array = new akantu::Array<akantu::Real>(model.getForce());
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
%include "element_group.hh"
%include "node_group.hh"
%include "dumper_iohelper.hh"
%include "dumpable_iohelper.hh"
%include "mesh.hh"

namespace akantu{
%extend Dumpable {
    void addDumpFieldExternalReal(const std::string & field_id,
                                  const Array<Real> & field){
      $self->addDumpFieldExternal<Real>(field_id,field);
    }
  }
 }
