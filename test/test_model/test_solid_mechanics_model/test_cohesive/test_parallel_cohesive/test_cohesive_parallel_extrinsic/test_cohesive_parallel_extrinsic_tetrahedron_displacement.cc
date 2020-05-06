/**
 * @file   test_cohesive_parallel_extrinsic_tetrahedron_displacement.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Nov 05 2014
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Displacement test for 3D cohesive elements
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#include "dumper_paraview.hh"
#include "material_cohesive_linear.hh"
#include "solid_mechanics_model_cohesive.hh"

#ifdef AKANTU_USE_IOHELPER
#include "dumper_paraview.hh"
#endif

/* -------------------------------------------------------------------------- */
using namespace akantu;

bool checkDisplacement(SolidMechanicsModelCohesive & model, ElementType type,
                       std::ofstream & error_output, UInt step,
                       bool barycenters);

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  debug::setDebugLevel(dblWarning);

  const UInt max_steps = 500;
  Math::setTolerance(1.e-12);

  UInt spatial_dimension = 3;
  ElementType type = _tetrahedron_10;

  Mesh mesh(spatial_dimension);

  const auto & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  akantu::MeshPartition * partition = NULL;
  if (prank == 0) {
    // Read the mesh
    mesh.read("tetrahedron.msh");

    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    //    debug::setDebugLevel(dblDump);
    partition->partitionate(psize);
    //    debug::setDebugLevel(dblWarning);
  }

  SolidMechanicsModelCohesive model(mesh);

  model.initParallel(partition, NULL, true);

  // debug::setDebugLevel(dblDump);
  // std::cout << mesh << std::endl;
  // debug::setDebugLevel(dblWarning);

  model.initFull(
      SolidMechanicsModelCohesiveOptions(_explicit_lumped_mass, true));

  /* ------------------------------------------------------------------------ */
  /* Facet part                                                               */
  /* ------------------------------------------------------------------------ */

  // Array<Real> limits(spatial_dimension, 2);
  // limits(0, 0) = -0.01;
  // limits(0, 1) = 0.01;
  // limits(1, 0) = -100;
  // limits(1, 1) = 100;
  // limits(2, 0) = -100;
  // limits(2, 1) = 100;

  // model.enableFacetsCheckOnArea(limits);

  /* ------------------------------------------------------------------------ */
  /* End of facet part                                                        */
  /* ------------------------------------------------------------------------ */

  // debug::setDebugLevel(dblDump);
  // std::cout << mesh_facets << std::endl;
  // debug::setDebugLevel(dblWarning);

  Real time_step = model.getStableTimeStep() * 0.1;
  model.setTimeStep(time_step);
  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();

  Array<Real> & position = mesh.getNodes();
  Array<Real> & velocity = model.getVelocity();
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & displacement = model.getDisplacement();
  //  const Array<Real> & residual = model.getResidual();

  UInt nb_nodes = mesh.getNbNodes();

  /// boundary conditions
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (position(n, 0) > 0.99 || position(n, 0) < -0.99) {
      for (UInt dim = 0; dim < spatial_dimension; ++dim) {
        boundary(n, dim) = true;
      }
    }

    if (position(n, 0) > 0.99 || position(n, 0) < -0.99) {
      for (UInt dim = 0; dim < spatial_dimension; ++dim) {
        boundary(n, dim) = true;
      }
    }
  }

  // #if defined (AKANTU_DEBUG_TOOLS)
  //   Vector<Real> facet_center(spatial_dimension);
  //   facet_center(0) =  0;
  //   facet_center(1) = -0.16666667;
  //   facet_center(2) =  0.5;

  //   debug::element_manager.setMesh(mesh);
  //   debug::element_manager.addModule(debug::_dm_material_cohesive);
  //   debug::element_manager.addModule(debug::_dm_debug_tools);
  //   //debug::element_manager.addModule(debug::_dm_integrator);
  // #endif

  /// initial conditions
  Real loading_rate = 1;
  Real disp_update = loading_rate * time_step;
  for (UInt n = 0; n < nb_nodes; ++n) {
    velocity(n, 0) = loading_rate * position(n, 0);
    velocity(n, 1) = loading_rate * position(n, 0);
  }

  model.synchronizeBoundaries();
  model.updateResidual();

  std::stringstream paraview_output;
  paraview_output << "extrinsic_parallel_tetrahedron_" << psize;

  model.setBaseName(paraview_output.str());
  model.addDumpFieldVector("displacement");
  model.addDumpFieldVector("velocity");
  model.addDumpFieldVector("acceleration");
  model.addDumpFieldVector("residual");
  model.addDumpFieldTensor("stress");
  model.addDumpFieldTensor("grad_u");
  model.addDumpField("partitions");
  //  model.getDumper().getDumper().setMode(iohelper::BASE64);
  model.dump();

  model.setBaseNameToDumper("cohesive elements",
                            paraview_output.str() + "_cohesive_elements");
  model.addDumpFieldVectorToDumper("cohesive elements", "displacement");
  model.addDumpFieldToDumper("cohesive elements", "damage");
  model.dump("cohesive elements");

  std::stringstream error_stream;
  error_stream << "error"
               << ".csv";
  std::ofstream error_output;
  error_output.open(error_stream.str().c_str());
  error_output << "# Step, Average, Max, Min" << std::endl;

  if (checkDisplacement(model, type, error_output, 0, true)) {
  }

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {

    /// update displacement on extreme nodes
    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      if (position(n, 0) > 0.99 || position(n, 0) < -0.99) {
        displacement(n, 0) += disp_update * position(n, 0);
        displacement(n, 1) += disp_update * position(n, 0);
      }
    }

    model.checkCohesiveStress();

    model.solveStep();

    if (s % 100 == 0) {
      if (prank == 0)
        std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }
  }

  model.dump();
  model.dump("cohesive elements");

  if (!checkDisplacement(model, type, error_output, max_steps, false)) {
    finalize();
    return EXIT_FAILURE;
  }

  finalize();
  return EXIT_SUCCESS;
}

bool checkDisplacement(SolidMechanicsModelCohesive & model, ElementType type,
                       std::ofstream & error_output, UInt step,
                       bool barycenters) {

  Mesh & mesh = model.getMesh();
  UInt spatial_dimension = mesh.getSpatialDimension();
  const Array<UInt> & connectivity = mesh.getConnectivity(type);
  const Array<Real> & displacement = model.getDisplacement();
  UInt nb_element = mesh.getNbElement(type);
  UInt nb_nodes_per_elem = Mesh::getNbNodesPerElement(type);

  const auto & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  if (psize == 1) {
    std::stringstream displacement_file;
    displacement_file << "displacement/displacement_" << std::setfill('0')
                      << std::setw(6) << step;
    std::ofstream displacement_output;
    displacement_output.open(displacement_file.str().c_str());

    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt n = 0; n < nb_nodes_per_elem; ++n) {
        UInt node = connectivity(el, n);

        for (UInt dim = 0; dim < spatial_dimension; ++dim) {
          displacement_output << std::setprecision(15)
                              << displacement(node, dim) << " ";
        }
        displacement_output << std::endl;
      }
    }

    displacement_output.close();

    if (barycenters) {
      std::stringstream barycenter_file;
      barycenter_file << "displacement/barycenters";
      std::ofstream barycenter_output;
      barycenter_output.open(barycenter_file.str().c_str());

      Element element(type, 0);
      Vector<Real> bary(spatial_dimension);

      for (UInt el = 0; el < nb_element; ++el) {
        element.element = el;
        mesh.getBarycenter(element, bary);

        for (UInt dim = 0; dim < spatial_dimension; ++dim) {
          barycenter_output << std::setprecision(15) << bary(dim) << " ";
        }
        barycenter_output << std::endl;
      }

      barycenter_output.close();
    }
  } else {

    if (barycenters)
      return true;

    /// read data
    std::stringstream displacement_file;
    displacement_file << "displacement/displacement_" << std::setfill('0')
                      << std::setw(6) << step;
    std::ifstream displacement_input;
    displacement_input.open(displacement_file.str().c_str());

    Array<Real> displacement_serial(0, spatial_dimension);
    Vector<Real> disp_tmp(spatial_dimension);

    while (displacement_input.good()) {
      for (UInt i = 0; i < spatial_dimension; ++i)
        displacement_input >> disp_tmp(i);

      displacement_serial.push_back(disp_tmp);
    }

    std::stringstream barycenter_file;
    barycenter_file << "displacement/barycenters";
    std::ifstream barycenter_input;
    barycenter_input.open(barycenter_file.str().c_str());

    Array<Real> barycenter_serial(0, spatial_dimension);

    while (barycenter_input.good()) {
      for (UInt dim = 0; dim < spatial_dimension; ++dim)
        barycenter_input >> disp_tmp(dim);

      barycenter_serial.push_back(disp_tmp);
    }

    Element element(type, 0);
    Vector<Real> bary(spatial_dimension);

    Array<Real>::iterator<Vector<Real>> it;
    Array<Real>::iterator<Vector<Real>> begin =
        barycenter_serial.begin(spatial_dimension);
    Array<Real>::iterator<Vector<Real>> end =
        barycenter_serial.end(spatial_dimension);

    Array<Real>::const_iterator<Vector<Real>> disp_it;
    Array<Real>::iterator<Vector<Real>> disp_serial_it;

    Vector<Real> difference(spatial_dimension);
    Array<Real> error;

    /// compute error
    for (UInt el = 0; el < nb_element; ++el) {
      element.element = el;
      mesh.getBarycenter(element, bary);

      /// find element
      for (it = begin; it != end; ++it) {
        UInt matched_dim = 0;

        while (matched_dim < spatial_dimension &&
               Math::are_float_equal(bary(matched_dim), (*it)(matched_dim)))
          ++matched_dim;

        if (matched_dim == spatial_dimension)
          break;
      }

      if (it == end) {
        std::cout << "Element barycenter not found!" << std::endl;
        return false;
      }

      UInt matched_el = it - begin;

      disp_serial_it = displacement_serial.begin(spatial_dimension) +
                       matched_el * nb_nodes_per_elem;

      for (UInt n = 0; n < nb_nodes_per_elem; ++n, ++disp_serial_it) {
        UInt node = connectivity(el, n);
        if (!mesh.isLocalOrMasterNode(node))
          continue;

        disp_it = displacement.begin(spatial_dimension) + node;

        difference = *disp_it;
        difference -= *disp_serial_it;

        error.push_back(difference.norm());
      }
    }

    /// compute average error
    Real average_error = std::accumulate(error.begin(), error.end(), 0.);
    comm.allReduce(&average_error, 1, _so_sum);

    UInt error_size = error.getSize();
    comm.allReduce(&error_size, 1, _so_sum);

    average_error /= error_size;

    /// compute maximum and minimum
    Real max_error = *std::max_element(error.begin(), error.end());
    comm.allReduce(&max_error, 1, _so_max);

    Real min_error = *std::min_element(error.begin(), error.end());
    comm.allReduce(&min_error, 1, _so_min);

    /// output data
    if (prank == 0) {
      error_output << step << ", " << average_error << ", " << max_error << ", "
                   << min_error << std::endl;
    }

    if (max_error > 1.e-9) {
      std::cout << "Displacement error is too big!" << std::endl;
      return false;
    }
  }

  return true;
}
