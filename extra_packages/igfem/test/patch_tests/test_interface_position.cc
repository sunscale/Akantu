/**
 * @file   test_interface_position.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  patch test for interface close to standard nodes
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "dumpable_inline_impl.hh"
#include "solid_mechanics_model_igfem.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

Real computeL2Error(SolidMechanicsModelIGFEM & model,
                    ElementTypeMapReal & error_per_element);

int main(int argc, char * argv[]) {

  initialize("material_test_interface_position.dat", argc, argv);
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  /// create a mesh and read the regular elements from the mesh file
  /// mesh creation
  const UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  akantu::MeshPartition * partition = NULL;
  if (prank == 0) {
    mesh.read("test_interface_position.msh");
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  /// model creation
  SolidMechanicsModelIGFEM model(mesh);
  model.initParallel(partition);
  delete partition;
  model.initFull();

  /// add fields that should be dumped
  model.setBaseName("regular_elements");
  model.setBaseNameToDumper("igfem elements", "igfem elements");
  model.addDumpField("material_index");
  model.addDumpFieldVector("displacement");
  model.addDumpField("blocked_dofs");
  model.addDumpField("stress");
  model.addDumpField("partitions");
  model.addDumpFieldToDumper("igfem elements", "lambda");
  model.addDumpFieldVectorToDumper("igfem elements", "displacement");
  model.addDumpFieldVectorToDumper("igfem elements", "real_displacement");
  model.addDumpFieldToDumper("igfem elements", "blocked_dofs");
  model.addDumpFieldToDumper("igfem elements", "material_index");
  model.addDumpFieldToDumper("igfem elements", "stress");
  model.addDumpFieldToDumper("igfem elements", "partitions");
  /// dump mesh before the IGFEM interface is created
  model.dump();
  model.dump("igfem elements");

  /// create the interace:
  UInt nb_standard_nodes = mesh.getNbNodes();
  std::list<SK::Sphere_3> sphere_list;
  SK::Sphere_3 sphere_1(SK::Point_3(0., 0., 0.), 0.25 * 0.25);
  sphere_list.push_back(sphere_1);
  model.registerGeometryObject(sphere_list, "inside");
  model.update();

  /// dump mesh after the IGFEM interface is created
  model.dump();
  model.dump("igfem elements");

  /// apply the boundary conditions: left and bottom side on rollers
  /// imposed displacement along right side
  mesh.computeBoundingBox();
  const Vector<Real> & lower_bounds = mesh.getLowerBounds();
  const Vector<Real> & upper_bounds = mesh.getUpperBounds();
  Real bottom = lower_bounds(1);
  Real left = lower_bounds(0);
  Real right = upper_bounds(0);
  Real eps = std::abs((right - left) * 1e-6);
  const Array<Real> & pos = mesh.getNodes();
  Array<Real> & disp = model.getDisplacement();
  Array<bool> & boun = model.getBlockedDOFs();

  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    if (std::abs(pos(i, 1) - bottom) < eps) {
      boun(i, 1) = true;
      disp(i, 1) = 0.0;
    }
    if (std::abs(pos(i, 0) - left) < eps) {
      boun(i, 0) = true;
      disp(i, 0) = 0.0;
    }
    if (std::abs(pos(i, 0) - right) < eps) {
      boun(i, 0) = true;
      disp(i, 0) = 1.0;
    }
  }

  /// compute the volume of the mesh
  Real int_volume = 0.;
  std::map<ElementKind, FEEngine *> fe_engines = model.getFEEnginesPerKind();
  std::map<ElementKind, FEEngine *>::const_iterator fe_it = fe_engines.begin();
  for (; fe_it != fe_engines.end(); ++fe_it) {
    ElementKind kind = fe_it->first;
    FEEngine & fe_engine = *(fe_it->second);
    Mesh::type_iterator it =
        mesh.firstType(spatial_dimension, _not_ghost, kind);
    Mesh::type_iterator last_type =
        mesh.lastType(spatial_dimension, _not_ghost, kind);
    for (; it != last_type; ++it) {
      ElementType type = *it;
      Array<Real> Volume(mesh.getNbElement(type) *
                             fe_engine.getNbIntegrationPoints(type),
                         1, 1.);
      int_volume += fe_engine.integrate(Volume, type);
    }
  }

  comm.allReduce(&int_volume, 1, _so_sum);
  if (prank == 0)
    if (!Math::are_float_equal(int_volume, 4)) {
      finalize();
      std::cout << "Error in area computation of the 2D mesh" << std::endl;
      return EXIT_FAILURE;
    }

  /// solve the system
  model.assembleStiffnessMatrix();
  Real error = 0;
  bool converged = false;
  bool factorize = false;

  converged = model.solveStep<_scm_newton_raphson_tangent,
                              SolveConvergenceCriteria::_increment>(
      1e-12, error, 2, factorize);
  if (!converged) {
    std::cout << "The solver did not converge!!! The error is: " << error
              << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  /// store the error on each element for visualization
  ElementTypeMapReal error_per_element("error_per_element");
  mesh.addDumpFieldExternal("error_per_element", error_per_element,
                            spatial_dimension, _not_ghost, _ek_regular);
  mesh.addDumpFieldExternalToDumper("igfem elements", "error_per_element",
                                    error_per_element, spatial_dimension,
                                    _not_ghost, _ek_igfem);
  mesh.initElementTypeMapArray(error_per_element, 1, spatial_dimension, false,
                               _ek_regular, true);
  mesh.initElementTypeMapArray(error_per_element, 1, spatial_dimension, false,
                               _ek_igfem, true);

  Real L2_error = computeL2Error(model, error_per_element);
  comm.allReduce(&L2_error, 1, _so_sum);

  if (prank == 0) {
    std::cout << "Error: " << L2_error << std::endl;
    if (L2_error > 1e-13) {
      finalize();
      std::cout << "The patch test did not pass!!!!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  /// dump the deformed mesh
  model.dump();
  model.dump("igfem elements");

  /* --------------------------------------------------------------------------
   */
  /// move the interface very close the standard nodes, but far enough
  /// to not cut trough the standard nodes
  model.moveInterface(0.5 * (1 - 1e-9));
  model.dump();
  model.dump("igfem elements");
  UInt nb_igfem_triangle_4 = mesh.getNbElement(_igfem_triangle_4, _not_ghost);
  UInt nb_igfem_triangle_5 = mesh.getNbElement(_igfem_triangle_5, _not_ghost);
  comm.allReduce(&nb_igfem_triangle_4, 1, _so_sum);
  comm.allReduce(&nb_igfem_triangle_5, 1, _so_sum);

  if (prank == 0) {
    if ((nb_igfem_triangle_4 != 0) || (nb_igfem_triangle_5 != 8)) {
      std::cout << "something went wrong in the interface creation"
                << std::endl;
      finalize();
      return EXIT_FAILURE;
    }
  }

  if ((psize == 0) && (mesh.getNbNodes() - nb_standard_nodes != 8)) {
    std::cout << "something went wrong in the interface node creation"
              << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  converged = model.solveStep<_scm_newton_raphson_tangent,
                              SolveConvergenceCriteria::_increment>(
      1e-12, error, 2, factorize);
  if (!converged) {
    std::cout << "The solver did not converge!!! The error is: " << error
              << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  L2_error = computeL2Error(model, error_per_element);
  comm.allReduce(&L2_error, 1, _so_sum);
  if (prank == 0) {
    std::cout << "Error: " << L2_error << std::endl;
    if (L2_error > 1e-13) {
      finalize();
      std::cout << "The patch test did not pass!!!!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  /// dump the new interface
  model.dump();
  model.dump("igfem elements");

  /* --------------------------------------------------------------------------
   */
  /// move the interface so that it cuts through the standard nodes
  model.moveInterface((0.5 * (1 - 1e-10)));
  nb_igfem_triangle_4 = mesh.getNbElement(_igfem_triangle_4, _not_ghost);
  nb_igfem_triangle_5 = mesh.getNbElement(_igfem_triangle_5, _not_ghost);
  comm.allReduce(&nb_igfem_triangle_4, 1, _so_sum);
  comm.allReduce(&nb_igfem_triangle_5, 1, _so_sum);
  if (prank == 0) {
    if ((nb_igfem_triangle_4 != 8) || (nb_igfem_triangle_5 != 0)) {
      std::cout << "something went wrong in the interface creation"
                << std::endl;
      finalize();
      return EXIT_FAILURE;
    }
  }

  if ((psize == 0) && (mesh.getNbNodes() - nb_standard_nodes != 4)) {
    std::cout << "something went wrong in the interface node creation"
              << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  converged = model.solveStep<_scm_newton_raphson_tangent,
                              SolveConvergenceCriteria::_increment>(
      1e-12, error, 2, factorize);
  if (!converged) {
    std::cout << "The solver did not converge!!! The error is: " << error
              << std::endl;
    finalize();
    return EXIT_FAILURE;
  }
  L2_error = computeL2Error(model, error_per_element);
  comm.allReduce(&L2_error, 1, _so_sum);
  if (prank == 0) {
    std::cout << "Error: " << L2_error << std::endl;
    if (L2_error > 1e-13) {
      finalize();
      std::cout << "The patch test did not pass!!!!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  /// dump the new interface
  model.dump();
  model.dump("igfem elements");

  finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
Real computeL2Error(SolidMechanicsModelIGFEM & model,
                    ElementTypeMapReal & error_per_element) {
  Real error = 0;
  Real normalization = 0;

  Mesh & mesh = model.getMesh();
  UInt spatial_dimension = mesh.getSpatialDimension();
  ElementTypeMapReal quad_coords("quad_coords");
  GhostType ghost_type = _not_ghost;
  const std::map<ElementKind, FEEngine *> & fe_engines =
      model.getFEEnginesPerKind();
  std::map<ElementKind, FEEngine *>::const_iterator fe_it = fe_engines.begin();
  for (; fe_it != fe_engines.end(); ++fe_it) {
    ElementKind kind = fe_it->first;
    FEEngine & fe_engine = *(fe_it->second);
    mesh.initElementTypeMapArray(quad_coords, spatial_dimension,
                                 spatial_dimension, false, kind, true);
    fe_engine.computeIntegrationPointsCoordinates(quad_coords);
    Mesh::type_iterator it =
        mesh.firstType(spatial_dimension, ghost_type, kind);
    Mesh::type_iterator last_type =
        mesh.lastType(spatial_dimension, ghost_type, kind);
    for (; it != last_type; ++it) {
      ElementType type = *it;
      UInt nb_elements = mesh.getNbElement(type, ghost_type);
      UInt nb_quads = fe_engine.getNbIntegrationPoints(type);
      /// interpolate the displacement at the quadrature points
      Array<Real> displ_on_quads(nb_quads * nb_elements, spatial_dimension,
                                 "displ_on_quads");
      Array<Real> quad_coords(nb_quads * nb_elements, spatial_dimension,
                              "quad_coords");
      fe_engine.interpolateOnIntegrationPoints(
          model.getDisplacement(), displ_on_quads, spatial_dimension, type);
      fe_engine.computeIntegrationPointsCoordinates(quad_coords, type,
                                                    ghost_type);
      Array<Real> & el_error = error_per_element(type, ghost_type);
      el_error.resize(nb_elements);
      Array<Real>::const_vector_iterator displ_it =
          displ_on_quads.begin(spatial_dimension);
      Array<Real>::const_vector_iterator coord_it =
          quad_coords.begin(spatial_dimension);
      Vector<Real> error_vec(spatial_dimension);
      for (UInt e = 0; e < nb_elements; ++e) {
        Vector<Real> error_per_quad(nb_quads);
        Vector<Real> normalization_per_quad(nb_quads);
        for (UInt q = 0; q < nb_quads; ++q, ++displ_it, ++coord_it) {
          Real exact = 0.5 * (*coord_it)(0) + 0.5;
          error_vec = *displ_it;
          error_vec(0) -= exact;
          error_per_quad(q) = error_vec.dot(error_vec);
          normalization_per_quad(q) = std::abs(exact) * std::abs(exact);
          ///  std::cout << error_vec  << std::endl;
        }
        /// integrate the error in the element and the corresponding
        /// normalization
        Real int_error =
            fe_engine.integrate(error_per_quad, type, e, ghost_type);
        error += int_error;
        el_error(e) = std::sqrt(int_error);
        normalization +=
            fe_engine.integrate(normalization_per_quad, type, e, ghost_type);
      }
    }
  }
  return (std::sqrt(error) / std::sqrt(normalization));
}
