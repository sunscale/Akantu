/**
 * @file   test_solid_mechanics_model.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Sep 03 15:56:56 2010
 *
 * @brief  test of the class SolidMechanicsModel
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
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "fe_engine.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  UInt max_steps = 1;
  Real epot, ekin;

  Mesh mesh(2);
  MeshIOMSH mesh_io;
  mesh_io.read("triangle.msh", mesh);
  mesh.getBlockedDOFs().createBoundariesFromGeometry();

  SolidMechanicsModel * model = new SolidMechanicsModel(mesh);

  /// model initialization
  model->initArrays();
  UInt nb_nodes = model->getFEEngine().getMesh().getNbNodes();
  memset(model->getForce().storage(),        0, 2*nb_nodes*sizeof(Real));
  memset(model->getVelocity().storage(),     0, 2*nb_nodes*sizeof(Real));
  memset(model->getAcceleration().storage(), 0, 2*nb_nodes*sizeof(Real));
  memset(model->getDisplacement().storage(), 0, 2*nb_nodes*sizeof(Real));

  model->readMaterials("material.dat");
  model->initMaterials();
  model->initModel();

  Real time_step = model->getStableTimeStep();
  model->setTimeStep(time_step/10.);

  model->assembleMass();

  std::cout << *model << std::endl;

  /// boundary conditions
  // Real eps = 1e-16;
  // for (UInt i = 0; i < nb_nodes; ++i) {
  //   model->getDisplacement().storage()[2*i] = model->getFEEngine().getMesh().getNodes().storage()[2*i] / 100.;

  //   if(model->getFEEngine().getMesh().getNodes().storage()[2*i] <= eps) {
  //     model->getBlockedDOFs().storage()[2*i    ] = true;
  //     if(model->getFEEngine().getMesh().getNodes().storage()[2*i + 1] <= eps)
  // 	model->getBlockedDOFs().storage()[2*i + 1] = true;
  //   }
  //   if(model->getFEEngine().getMesh().getNodes().storage()[2*i + 1] <= eps) {
  //     model->getBlockedDOFs().storage()[2*i + 1] = true;
  //   }

  // }

  // Boundary condition (Neumann)
  Matrix<Real> stress(2,2);
  stress.eye(1e3);
  model.applyBC(BC::Neumann::FromHigherDim(stress), "0");

  // const Mesh::ConnectivityTypeList & type_list = fem_boundary.getMesh().getConnectivityTypeList();
  // Mesh::ConnectivityTypeList::const_iterator it;
  // for(it = type_list.begin(); it != type_list.end(); ++it) {
  //   if(Mesh::getSpatialDimension(*it) != fem_boundary.getElementDimension()) continue;

  //   //    ElementType facet_type = Mesh::getFacetElementType(*it);
  //   UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
  //   UInt nb_quad              = FEEngine::getNbQuadraturePoints(*it);


  //   UInt nb_element;
  //   const Array<Real> * shapes;
  //   Array<Real> quad_coords(0,2,"quad_coords");
  //   const Array<Real> * normals_on_quad;

  //   nb_element   = fem_boundary.getMesh().getNbElement(*it);
  //   fem_boundary.interpolateOnQuadraturePoints(mesh.getNodes(), quad_coords, 2, _segment_2);
  //   normals_on_quad = &(fem_boundary.getNormalsOnQuadPoints(*it));

  //   shapes       = &(fem_boundary.getShapes(*it));

  //   Array<Real> * sigma_funct = new Array<Real>(nb_element, 4*nb_quad, "myfunction");
  //   Array<Real> * funct = new Array<Real>(nb_element, 2*nb_quad, "myfunction");

  //   Real * sigma_funct_val = sigma_funct->storage();
  //   Real * shapes_val = shapes->storage();

  //   /// compute t * \phi_i for each nodes of each element
  //   for (UInt el = 0; el < nb_element; ++el) {
  //     for (UInt q = 0; q < nb_quad; ++q) {
  // 	trac(quad_coords.storage()+el*nb_quad*2+q,sigma_funct_val);
  // 	sigma_funct_val += 4;
  //     }
  //   }

  //   Math::matrix_vector(2,2,*sigma_funct,*normals_on_quad,*funct);
  //   funct->extendComponentsInterlaced(nb_nodes_per_element,2);

  //   Real * funct_val = funct->storage();
  //   for (UInt el = 0; el < nb_element; ++el) {
  //     for (UInt q = 0; q < nb_quad; ++q) {
  // 	for (UInt n = 0; n < nb_nodes_per_element; ++n) {
  // 	  *funct_val++ *= *shapes_val;
  // 	  *funct_val++ *= *shapes_val++;
  // 	}
  //     }
  //   }


  //   Array<Real> * int_funct = new Array<Real>(nb_element, 2*nb_nodes_per_element,
  // 						    "inte_funct");
  //   fem_boundary.integrate(*funct, *int_funct, 2*nb_nodes_per_element, *it);
  //   delete funct;

  //   fem_boundary.assembleArray(*int_funct,model->getForce(), 2, *it);
  //   delete int_funct;
  // }


  //  model->getDisplacement().storage()[1] = 0.1;


#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  dumper.SetMode(iohelper::TEXT);

  dumper.SetPoints(model->getFEEngine().getMesh().getNodes().storage(), 2, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model->getFEEngine().getMesh().getConnectivity(_triangle_3).storage(),
			 iohelper::TRIANGLE1, model->getFEEngine().getMesh().getNbElement(_triangle_3), iohelper::C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().storage(), 2, "displacements");
  dumper.AddNodeDataField(model->getVelocity().storage(), 2, "velocity");
  dumper.AddNodeDataField(model->getForce().storage(), 2, "force");
  dumper.AddNodeDataField(model->getResidual().storage(), 2, "residual");
  dumper.AddElemDataField(model->getMaterial(0).getStrain(_triangle_3).storage(), 4, "strain");
  dumper.AddElemDataField(model->getMaterial(0).getStress(_triangle_3).storage(), 4, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("force", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
#endif //AKANTU_USE_IOHELPER

  model->setPotentialEnergyFlagOn();
  for(UInt s = 0; s < max_steps; ++s) {
    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    epot = model->getPotentialEnergy();
    ekin = model->getKineticEnergy();

    std::cout << s << " " << epot << " " << ekin << " " << epot + ekin
	      << std::endl;

#ifdef AKANTU_USE_IOHELPER
    if(s % 10 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  return EXIT_SUCCESS;
}



