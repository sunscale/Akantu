namespace akantu {
  %ignore MeshPartition::getPartitions;
  %ignore MeshPartition::getPartition;
  %ignore MeshPartition::getGhostPartitionCSR;
}

%include "mesh_partition.hh"
%include "mesh_utils.hh"

