%{
#include "mesh.hh"
#include "node_group.hh"
#include "solid_mechanics_model.hh"

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
}

print_self(Mesh)


%extend akantu::Mesh {

  void resizeMesh(UInt nb_nodes, UInt nb_element, const ElementType & type) {

    Array<Real> & nodes = const_cast<Array<Real> &>($self->getNodes());
    nodes.resize(nb_nodes);

    $self->addConnectivityType(type);
    Array<UInt> & connectivity = const_cast<Array<UInt> &>($self->getConnectivity(type));
    connectivity.resize(nb_element);
    


  }

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
    surface_array.resize(group_node.getSize());

    for (UInt i = 0; i < group_node.getSize(); ++i) {
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
    surface_array.resize(group_node.getSize());
    
    for (UInt i = 0; i < group_node.getSize(); ++i) {
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

%include "mesh.hh"


