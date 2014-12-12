%{
#include "mesh.hh"

using akantu::Vector;
using akantu::ElementTypeMapArray;
using akantu::MatrixProxy;
using akantu::Matrix;
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
}

print_self(Mesh)

%extend akantu::GroupManager {
  void createGroupsFromStringMeshData(const std::string & dataset_name) {
    $self->createGroupsFromMeshData<std::string>(dataset_name);
  }

  void createGroupsFromUIntMeshData(const std::string & dataset_name) {
    $self->createGroupsFromMeshData<akantu::UInt>(dataset_name);
  }
}

%include "group_manager.hh"


%include "mesh.hh"


namespace akantu {
  %template(MeshEventUInt) MeshEvent<UInt>;
  %template(MeshEventElement) MeshEvent<Element>;
}
