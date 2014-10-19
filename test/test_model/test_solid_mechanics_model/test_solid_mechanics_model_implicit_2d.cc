/**
 * @file   test_solid_mechanics_model_implicit_2d.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Mar 03 2011
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  test of traction in implicit
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

/* -------------------------------------------------------------------------- */
#include <limits>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"

#define bar_length 1
#define bar_height 1

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  debug::setDebugLevel(dblWarning);
  initialize("material.dat", argc, argv);

  UInt spatial_dimension = 2;

  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  MeshPartition * partition = NULL;
  if(prank == 0) {
    mesh.read("square_implicit2.msh");
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    //   partition->reorder();
    partition->partitionate(psize);
  }
  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initParallel(partition);
  delete partition;

  model.initFull(SolidMechanicsModelOptions(_static));

  if (prank == 0) std::cout << model.getMaterial("steel") << std::endl;

  /// boundary conditions
  const  Array<Real> & position = mesh.getNodes();
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & displacment = model.getDisplacement();

  UInt nb_nodes = model.getFEEngine().getMesh().getNbNodes();
  for (UInt n = 0; n < nb_nodes; ++n) {
    if(position(n,0) < Math::getTolerance()) boundary(n,0) = true;
    if(position(n,1) < Math::getTolerance()) boundary(n,1) = true;

    if(std::abs(position(n,0) - bar_length) < Math::getTolerance()) {
      boundary(n,0) = true;
      displacment(n,0) = 0.1;
    }
  }

  model.setBaseName("implicit_2d");
  model.addDumpField("displacement");
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("force"       );
  model.addDumpField("residual"    );
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );

  model.updateResidual();
  model.dump();

  model.solveStep<_scm_newton_raphson_tangent_modified, _scc_increment>(1e-3, 100);
  model.dump();

  finalize();

  return EXIT_SUCCESS;
}
