/**
 * @file   test_solid_mechanics_model_pbc_parallel.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Fri Apr 13 16:31:38 2012
 *
 * @brief  test if pbc works in parallel if partition is strips
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
#include <limits>
#include <fstream>
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"
#endif //AKANTU_USE_IOHELPER

#ifdef AKANTU_USE_IOHELPER
akantu::ElementType type = akantu::_quadrangle_4;
iohelper::ElemType paraview_type = iohelper::QUAD1;
#endif //AKANTU_USE_IOHELPER

#ifdef AKANTU_USE_IOHELPER
static void paraviewInit(iohelper::Dumper & dumper, const akantu::SolidMechanicsModel & model);
static void paraviewDump(iohelper::Dumper & dumper);
#endif

akantu::Array<akantu::Real> proc_rank(0,1);

int main(int argc, char *argv[])
{
  akantu::debug::setDebugLevel(akantu::dblWarning);
  akantu::initialize(argc, argv);

  akantu::StaticCommunicator & comm = 
    akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm->getNbProc();
  akantu::Int prank = comm->whoAmI();

  akantu::UInt spatial_dimension = 2;
  //  akantu::UInt max_steps = 1;
  akantu::Real time_factor = 0.2;

  akantu::Mesh mesh(spatial_dimension);

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    akantu::MeshIOMSH mesh_io;
    mesh_io.read("square_structured.msh", mesh);
    partition = new akantu::MeshPartitionScotch(mesh, spatial_dimension);
    partition->setNbPartition(psize);

    // create the partition
    akantu::Array<akantu::Int> part_tab(0,1);
    mesh.computeBoundingBox();
    akantu::Real rank_border = 0.5 * (mesh.getYMax() - mesh.getYMin());
    rank_border += 1e-10;

    akantu::Mesh::type_iterator it  = mesh.firstType(spatial_dimension);
    akantu::Mesh::type_iterator end = mesh.lastType(spatial_dimension);
    for(; it != end; ++it) {
      akantu::ElementType c_type = *it;
      akantu::UInt nb_element = mesh.getNbElement(*it);
      
      for (akantu::UInt el=0; el<nb_element; ++el) {
	akantu::Real barycenter[spatial_dimension];
	mesh.getBarycenter(el,c_type,barycenter);
	if (barycenter[1] > rank_border)
	  part_tab.push_back(0);
	else
	  part_tab.push_back(1);
      }
    }
    
    partition->fillPartitionInformation(mesh,part_tab.storage());
  } else {
    mesh.computeBoundingBox();
  }

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);
  model->initParallel(partition);

  /// model initialization
  model->initArrays();
  
  /// set vectors to 0
  model->getForce().clear();
  model->getVelocity().clear();
  model->getAcceleration().clear();
  model->getDisplacement().clear();

  model->initExplicit();
  model->initModel();
  model->readMaterials("material.dat");
  model->initMaterials();

  if(prank == 0)
    std::cout << model->getMaterial(0) << std::endl;

  model->setPBC(1,0,0);
  model->initPBC();
  model->assembleMassLumped();

  akantu::UInt nb_element = mesh.getNbElement(type);
  akantu::UInt nb_quads = model->getFEEngine().getNbQuadraturePoints(type);
  for(akantu::UInt i=0; i<nb_element * nb_quads; ++i) 
    proc_rank.push_back(prank);

  /// boundary conditions
  akantu::UInt nb_nodes = model->getFEEngine().getMesh().getNbNodes();
  akantu::Real eps = 1e-16;
  akantu::Array<akantu::Real> & coords = const_cast<akantu::Array<akantu::Real> & >(model->getFEEngine().getMesh().getNodes());
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    // block top and bottom nodes
    if(std::abs(coords(i,1)-mesh.getYMax()) <= eps 
       || std::abs(coords(i,1)-mesh.getYMin()) <= eps) {
      model->getBlockedDOFs().storage()[spatial_dimension*i + 1] = true;
    }
    // correct coordinates (gmsh's unprecision)
    for (akantu::UInt d=0; d<spatial_dimension; ++d) {
      akantu::Real cor = std::floor(10 * coords(i,d) + 0.5) / 10.;
      coords(i,d) = cor;
      //std::cout << cor << " ";
    }
    //std::cout << std::endl;
  }
  //std::cout << std::endl;
  model->synchronizeBoundaries();

  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  if(prank == 0)  
    std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);

#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
  model->updateResidual();
  iohelper::DumperParaview dumper;
  paraviewInit(dumper, *model);
#endif //AKANTU_USE_IOHELPER

  // modify displacements
  akantu::Array<akantu::Real> & displacement = model->getDisplacement();
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    displacement(i,1) = std::abs(coords(i,1) - mesh.getYMin()) * 0.0001;
  }
  model->synchronizeBoundaries();

#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
  model->updateResidual();
  model->synchronizeResidual();
  paraviewDump(dumper);
#endif //AKANTU_USE_IOHELPER

  // test (traction at top and bottom boundary should be 2.826923077e7)
  // therefore the nodal residual should be 2.826923077e6
  akantu::Real solution = 2.826923077e6;
  akantu::Real adm_error = 1e-3;
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    akantu::Real trac = std::abs(model->getResidual().storage()[spatial_dimension*i + 1]);
    if((std::abs(coords(i,1)-mesh.getYMax()) <= eps 
       || std::abs(coords(i,1)-mesh.getYMin()) <= eps) &&
       std::abs(trac - solution) > adm_error) {
      std::cerr << "Boundary residual in y direction is " << trac << 
	" but should be " << solution << "!!!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  akantu::finalize();
  if(prank == 0)
    std::cout << "Test successful!" << std::endl;
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
/* iohelper::Dumper vars                                                                */
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER
void paraviewInit(iohelper::Dumper & dumper, const akantu::SolidMechanicsModel & model) {
  akantu::StaticCommunicator * comm = 
    akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm->getNbProc();
  akantu::Int prank = comm->whoAmI();

  akantu::UInt spatial_dimension = model.getSpatialDimension();
  akantu::UInt nb_nodes = model.getFEEngine().getMesh().getNbNodes();
  akantu::UInt nb_element = model.getFEEngine().getMesh().getNbElement(type);
  dumper.SetParallelContext(prank, psize);
  dumper.SetPoints(model.getFEEngine().getMesh().getNodes().storage(),
		   spatial_dimension, nb_nodes, "pbc_parallel");
  dumper.SetConnectivity((int *)model.getFEEngine().getMesh().getConnectivity(type).storage(),
			 paraview_type, nb_element, iohelper::C_MODE);
  dumper.AddNodeDataField(model.getDisplacement().storage(),
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model.getVelocity().storage(),
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model.getAcceleration().storage(),
			  spatial_dimension, "acceleration");
  dumper.AddNodeDataField(model.getResidual().storage(),
			  spatial_dimension, "force");
  dumper.AddNodeDataField(model.getMass().storage(),
			  spatial_dimension, "mass");
  dumper.AddNodeDataField(model.getForce().storage(),
			  spatial_dimension, "applied_force");
  dumper.AddElemDataField(model.getMaterial(0).getStrain(type).storage(),
   			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model.getMaterial(0).getStress(type).storage(),
   			  spatial_dimension*spatial_dimension, "stress");
  dumper.AddElemDataField(proc_rank.storage(), 1, "rank");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */
void paraviewDump(iohelper::Dumper & dumper) {
  dumper.Dump();
}
#endif
