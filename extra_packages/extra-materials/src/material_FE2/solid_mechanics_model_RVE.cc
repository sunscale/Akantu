/**
 * @file   solid_mechanics_model_RVE.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Jan 13 15:32:35 2016
 *
 * @brief  Implementation of SolidMechanicsModelRVE
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
#include "solid_mechanics_model_RVE.hh"
#include "material_damage_iterative.hh"
#ifdef AKANTU_USE_MUMPS
#include "solver_mumps.hh"
#endif
#ifdef AKANTU_USE_PETSC
#include "solver_petsc.hh"
#include "petsc_matrix.hh"
#endif
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
SolidMechanicsModelRVE::SolidMechanicsModelRVE(Mesh & mesh, bool use_RVE_mat_selector,
					       UInt nb_gel_pockets, UInt dim, 
					       const ID & id, const MemoryID & memory_id) :
  SolidMechanicsModel(mesh, dim, id, memory_id),
  volume(0.),
  use_RVE_mat_selector(use_RVE_mat_selector),
  static_communicator_dummy(StaticCommunicator::getStaticCommunicatorDummy()),
  nb_gel_pockets(nb_gel_pockets),
  nb_dumps(0) {
  AKANTU_DEBUG_IN();
  /// create node groups for PBCs
  mesh.createGroupsFromMeshData<std::string>("physical_names");
  /// find the four corner nodes of the RVE
  findCornerNodes();

  /// remove the corner nodes from the surface node groups:
  /// This most be done because corner nodes a not periodic
  mesh.getElementGroup("top").removeNode(corner_nodes(2));
  mesh.getElementGroup("top").removeNode(corner_nodes(3));
  mesh.getElementGroup("left").removeNode(corner_nodes(3));
  mesh.getElementGroup("left").removeNode(corner_nodes(0));
  mesh.getElementGroup("bottom").removeNode(corner_nodes(1));
  mesh.getElementGroup("bottom").removeNode(corner_nodes(0));
  mesh.getElementGroup("right").removeNode(corner_nodes(2));
  mesh.getElementGroup("right").removeNode(corner_nodes(1));

  const ElementGroup & bottom = mesh.getElementGroup("bottom");
  bottom_nodes.insert( bottom.node_begin(), bottom.node_end() );

  const ElementGroup & left = mesh.getElementGroup("left");
  left_nodes.insert( left.node_begin(), left.node_end() );

  /// enforce periodicity on the displacement fluctuations
  SurfacePair surface_pair_1 = std::make_pair("top","bottom");
  SurfacePair surface_pair_2 = std::make_pair("right","left");
  SurfacePairList surface_pairs_list;
  surface_pairs_list.push_back(surface_pair_1);
  surface_pairs_list.push_back(surface_pair_2);
  this->setPBC(surface_pairs_list);
 AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModelRVE::~SolidMechanicsModelRVE() {
  delete static_communicator_dummy;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::initFull(const ModelOptions & options) {
  AKANTU_DEBUG_IN();
  SolidMechanicsModel::initFull(options);

  this->initMaterials();

  /// compute the volume of the RVE
  FEEngine * fem = this->fems["SolidMechanicsFEEngine"];
  GhostType gt = _not_ghost;
  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt, _ek_not_defined);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, _ek_not_defined);
  for(; it != end; ++it) {
    const ElementType element_type = *it;
    Array<Real> Volume(this->mesh.getNbElement(element_type) * this->fems["SolidMechanicsFEEngine"]->getNbIntegrationPoints(element_type), 1, 1.);
    this->volume = fem->integrate(Volume, element_type);
  }

  std::cout << "The volume of the RVE is " << this->volume << std::endl; 

  /// dumping
  std::stringstream base_name;
  base_name << this->id; // << this->memory_id - 1;
  this->setBaseName       (base_name.str());
  this->addDumpFieldVector("displacement");
  this->addDumpField      ("stress"      );
  this->addDumpField      ("grad_u"      );
  this->addDumpField      ("eigen_grad_u"      );
  this->addDumpField      ("blocked_dofs"      );
  this->addDumpField      ("material_index"      );
  this->addDumpField      ("damage"      );
  this->addDumpField      ("Sc");
  this->addDumpField      ("force");
  this->addDumpField      ("equivalent_stress");
  this->addDumpField      ("residual");

  this->dump();
  this->nb_dumps +=1;
 AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::applyBoundaryConditions(const Matrix<Real> & displacement_gradient) {
  AKANTU_DEBUG_IN();
  /// get the position of the nodes
  const Array<Real> & pos = mesh.getNodes();
  /// storage for the coordinates of a given node and the displacement that will be applied
  Vector<Real> x(spatial_dimension);
  Vector<Real> appl_disp(spatial_dimension);

  /// fix top right node
  UInt node = this->corner_nodes(2);
  x(0) = pos(node,0); x(1) = pos(node,1);
  appl_disp.mul<false>(displacement_gradient,x);
  (*this->blocked_dofs)(node,0) = true; (*this->displacement)(node,0) = appl_disp(0);
  (*this->blocked_dofs)(node,1) = true; (*this->displacement)(node,1) = appl_disp(1);
  // (*this->blocked_dofs)(node,0) = true; (*this->displacement)(node,0) = 0.;
  // (*this->blocked_dofs)(node,1) = true; (*this->displacement)(node,1) = 0.;

  /// apply Hx at all the other corner nodes; H: displ. gradient
  node = this->corner_nodes(0);
  x(0) = pos(node,0); x(1) = pos(node,1);
  appl_disp.mul<false>(displacement_gradient,x);
  (*this->blocked_dofs)(node,0) = true; (*this->displacement)(node,0) = appl_disp(0);
  (*this->blocked_dofs)(node,1) = true; (*this->displacement)(node,1) = appl_disp(1);

  node = this->corner_nodes(1);
  x(0) = pos(node,0); x(1) = pos(node,1);
  appl_disp.mul<false>(displacement_gradient,x);
  (*this->blocked_dofs)(node,0) = true; (*this->displacement)(node,0) = appl_disp(0);
  (*this->blocked_dofs)(node,1) = true; (*this->displacement)(node,1) = appl_disp(1);

  node = this->corner_nodes(3);
  x(0) = pos(node,0); x(1) = pos(node,1);
  appl_disp.mul<false>(displacement_gradient,x);
  (*this->blocked_dofs)(node,0) = true; (*this->displacement)(node,0) = appl_disp(0);
  (*this->blocked_dofs)(node,1) = true; (*this->displacement)(node,1) = appl_disp(1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::findCornerNodes() {
  AKANTU_DEBUG_IN();
  mesh.computeBoundingBox();

  // find corner nodes
  const Array<Real> & position = mesh.getNodes();
  const Vector<Real> & lower_bounds = mesh.getLowerBounds();
  const Vector<Real> & upper_bounds = mesh.getUpperBounds();

  AKANTU_DEBUG_ASSERT(spatial_dimension == 2, "This is 2D only!");
  corner_nodes.resize(4);
  corner_nodes.set(UInt(-1));

  for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
    // node 1
    if (Math::are_float_equal(position(n, 0), lower_bounds(0)) &&
	Math::are_float_equal(position(n, 1), lower_bounds(1))) {
      corner_nodes(0) = n;
    }
    // node 2
    else if (Math::are_float_equal(position(n, 0), upper_bounds(0)) &&
	     Math::are_float_equal(position(n, 1), lower_bounds(1))) {
      corner_nodes(1) = n;
    }
    // node 3
    else if (Math::are_float_equal(position(n, 0), upper_bounds(0)) &&
	     Math::are_float_equal(position(n, 1), upper_bounds(1))) {
      corner_nodes(2) = n;
    }
    // node 4
    else if (Math::are_float_equal(position(n, 0), lower_bounds(0)) &&
	     Math::are_float_equal(position(n, 1), upper_bounds(1))) {
      corner_nodes(3) = n;
    }
  }

  for (UInt i = 0; i < corner_nodes.getSize(); ++i) {
    if (corner_nodes(i) == UInt(-1))
      AKANTU_DEBUG_ERROR("The corner node " << i + 1 << " wasn't found");
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::advanceASR(const Matrix<Real> & prestrain) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(spatial_dimension == 2, "This is 2D only!");

  /// apply the new eigenstrain
  GhostType gt = _not_ghost;

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt, _ek_not_defined);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, _ek_not_defined);
  for(; it != end; ++it) {
    const ElementType element_type = *it;
    Array<Real> & prestrain_vect = const_cast<Array<Real> &>(this->getMaterial("gel").getInternal<Real>("eigen_grad_u")(element_type, gt));
    Array<Real>::iterator< Matrix<Real> > prestrain_it = prestrain_vect.begin(spatial_dimension, spatial_dimension);
    Array<Real>::iterator< Matrix<Real> > prestrain_end = prestrain_vect.end(spatial_dimension, spatial_dimension);

    for (; prestrain_it != prestrain_end; ++prestrain_it) 
      (*prestrain_it) = prestrain;
  }


  /// advance the damage
  MaterialDamageIterative<2> & mat_paste = dynamic_cast<MaterialDamageIterative<2> & >(*this->materials[1]);
  MaterialDamageIterative<2> & mat_aggregate = dynamic_cast<MaterialDamageIterative<2> & >(*this->materials[0]); 
  UInt nb_damaged_elements = 0;
  Real max_eq_stress_aggregate = 0;
  Real max_eq_stress_paste = 0;
  Real error = 0;
  bool converged = false;

  do {   	
    
    converged = this->solveStep<_scm_newton_raphson_tangent, _scc_increment>(1e-6, error, 2, false, *static_communicator_dummy);
    AKANTU_DEBUG_ASSERT(converged, "Did not converge");
    std::cout << "the error is " << error << std::endl;
    /// compute damage 
    max_eq_stress_aggregate = mat_aggregate.getNormMaxEquivalentStress();
    max_eq_stress_paste = mat_paste.getNormMaxEquivalentStress();
    
    nb_damaged_elements = 0;
    if (max_eq_stress_aggregate > max_eq_stress_paste)
      nb_damaged_elements = mat_aggregate.updateDamage();
    else if (max_eq_stress_aggregate < max_eq_stress_paste)
      nb_damaged_elements = mat_paste.updateDamage();
    else
      nb_damaged_elements = (mat_paste.updateDamage() + mat_aggregate.updateDamage()); 

    std::cout << "the number of damaged elements is " << nb_damaged_elements << std::endl;
  } while (nb_damaged_elements);

  if (this->nb_dumps % 10 == 0) {
    this->dump();
  }
  this->nb_dumps += 1;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
Real SolidMechanicsModelRVE::averageTensorField(UInt row_index, UInt col_index, const ID & field_type) {
  AKANTU_DEBUG_IN();
  FEEngine * fem = this->fems["SolidMechanicsFEEngine"];
  Real average = 0;
  GhostType gt = _not_ghost;

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt, _ek_not_defined);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, _ek_not_defined);
  for(; it != end; ++it) {
    const ElementType element_type = *it;
    if (field_type == "stress") {
      for(UInt m = 0; m < this->materials.size(); ++m) {
	const Array<Real> & stress_vec = this->materials[m]->getStress(element_type);
	const Array<UInt> & elem_filter = this->materials[m]->getElementFilter(element_type);
	Array<Real> int_stress_vec(elem_filter.getSize(), spatial_dimension*spatial_dimension, "int_of_stress");

	fem->integrate(stress_vec, int_stress_vec, spatial_dimension*spatial_dimension, element_type, _not_ghost, elem_filter);

	for (UInt k = 0; k < elem_filter.getSize(); ++k)
	  average += int_stress_vec(k, row_index * spatial_dimension + col_index); //3 is the value for the yy (in 3D, the value is 4)
      }
    }

    else if (field_type == "strain") {
      for(UInt m = 0; m < this->materials.size(); ++m) {
	const Array<Real> & gradu_vec = this->materials[m]->getGradU(element_type);
	const Array<UInt> & elem_filter = this->materials[m]->getElementFilter(element_type);
	Array<Real> int_gradu_vec(elem_filter.getSize(), spatial_dimension*spatial_dimension, "int_of_gradu");

	fem->integrate(gradu_vec, int_gradu_vec, spatial_dimension*spatial_dimension, element_type, _not_ghost, elem_filter);

	for (UInt k = 0; k < elem_filter.getSize(); ++k)
	  /// averaging is done only for normal components, so stress and strain are equal
	  average += 0.5 * (int_gradu_vec(k, row_index * spatial_dimension + col_index) + int_gradu_vec(k, col_index * spatial_dimension + row_index));
      }
    }

    else if (field_type == "eigen_grad_u") {
      for(UInt m = 0; m < this->materials.size(); ++m) {
	const Array<Real> & eigen_gradu_vec = this->materials[m]->getInternal<Real>("eigen_grad_u")(element_type);
	const Array<UInt> & elem_filter = this->materials[m]->getElementFilter(element_type);
	Array<Real> int_eigen_gradu_vec(elem_filter.getSize(), spatial_dimension*spatial_dimension, "int_of_gradu");

	fem->integrate(eigen_gradu_vec, int_eigen_gradu_vec, spatial_dimension*spatial_dimension, element_type, _not_ghost, elem_filter);

	for (UInt k = 0; k < elem_filter.getSize(); ++k)
	  /// averaging is done only for normal components, so stress and strain are equal
	  average += int_eigen_gradu_vec(k, row_index * spatial_dimension + col_index);
      }
    }

    else
      AKANTU_DEBUG_ERROR("Averaging not implemented for this field!!!");
  }
  return average/this->volume;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::homogenizeStiffness(Matrix<Real> & C_macro) {
  AKANTU_DEBUG_IN();
  const UInt dim = 2;
  AKANTU_DEBUG_ASSERT(this->spatial_dimension == dim, "Is only implemented for 2D!!!");

  /// apply three independent loading states to determine C
  /// 1. eps_el = (1;0;0) 2. eps_el = (0,1,0) 3. eps_el = (0,0,0.5)

  /// clear the eigenstrain
  GhostType gt = _not_ghost;
  Mesh::type_iterator it  = mesh.firstType(dim, gt, _ek_not_defined);
  Mesh::type_iterator end = mesh.lastType(dim, gt, _ek_not_defined);
  Matrix<Real> zero_eigengradu(dim, dim, 0.);
  for(; it != end; ++it) {
    const ElementType element_type = *it;
    Array<Real> & prestrain_vect = const_cast<Array<Real> &>(this->getMaterial("gel").getInternal<Real>("eigen_grad_u")(element_type, gt));
    Array<Real>::iterator< Matrix<Real> > prestrain_it = prestrain_vect.begin(spatial_dimension, spatial_dimension);
    Array<Real>::iterator< Matrix<Real> > prestrain_end = prestrain_vect.end(spatial_dimension, spatial_dimension);

    for (; prestrain_it != prestrain_end; ++prestrain_it) 
      (*prestrain_it) = zero_eigengradu;
  }

  /// storage for results of 3 different loading states
  UInt voigt_size = VoigtHelper<dim>::size;
  Matrix<Real> stresses(voigt_size, voigt_size, 0.);
  Matrix<Real> strains(voigt_size, voigt_size, 0.);
  Matrix<Real> H(dim, dim, 0.);

  /// save the damage state before fillin up cracks
  ElementTypeMapReal saved_damage("saved_damage");
  mesh.initElementTypeMapArray(saved_damage, 1, spatial_dimension, false, _ek_regular, true);
  saved_damage.clear();
  /// fill the cracks for the tension test
  //  this->fillCracks(saved_damage);
  
  /// virtual test 1:
  H(0,0) = 0.01;
  this->performVirtualTesting(H, stresses, strains, 0);

  /// virtual test 2:
  H.clear();
  H(1,1) = 0.01;
  this->performVirtualTesting(H, stresses, strains, 1);

  /// virtual test 3:
  H.clear();
  H(0,1) = 0.01;
  this->performVirtualTesting(H, stresses, strains, 2);

  /// drain cracks
  //this->drainCracks(saved_damage);
  /// compute effective stiffness
  Matrix<Real> eps_inverse(voigt_size, voigt_size);
  eps_inverse.inverse(strains);
  C_macro.mul<false, false>(stresses, eps_inverse);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::performVirtualTesting(const Matrix<Real> & H, Matrix<Real> & eff_stresses, Matrix<Real> & eff_strains, const UInt test_no) {
  AKANTU_DEBUG_IN();
  this->applyBoundaryConditions(H);

  /// solve system
  this->assembleStiffnessMatrix();
  Real error = 0;
  bool converged= this->solveStep<_scm_newton_raphson_tangent_not_computed, _scc_increment>(1e-6, error, 2, false, *static_communicator_dummy);
  std::cout << "error in tension test " << error << std::endl;
  AKANTU_DEBUG_ASSERT(converged, "Did not converge");

  /// get average stress and strain
  eff_stresses(0, test_no) = this->averageTensorField(0,0, "stress");
  eff_strains(0, test_no) = this->averageTensorField(0,0, "strain");
  eff_stresses(1, test_no) = this->averageTensorField(1,1, "stress");
  eff_strains(1, test_no) = this->averageTensorField(1,1, "strain");
  eff_stresses(2, test_no) = this->averageTensorField(1,0, "stress");
  eff_strains(2, test_no) = 2. * this->averageTensorField(1,0, "strain");
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::homogenizeEigenGradU(Matrix<Real> & eigen_gradu_macro) {
  AKANTU_DEBUG_IN();
  eigen_gradu_macro(0,0) = this->averageTensorField(0,0, "eigen_grad_u");
  eigen_gradu_macro(1,1) = this->averageTensorField(1,1, "eigen_grad_u");
  eigen_gradu_macro(0,1) = this->averageTensorField(0,1, "eigen_grad_u");
  eigen_gradu_macro(1,0) = this->averageTensorField(1,0, "eigen_grad_u");
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::initMaterials() {
  AKANTU_DEBUG_IN();

  // make sure the material are instantiated
  if(!are_materials_instantiated) instantiateMaterials();

  if (use_RVE_mat_selector) {
    const Vector<Real> & lowerBounds = mesh.getLowerBounds();
    const Vector<Real> & upperBounds = mesh.getUpperBounds();
    Real bottom  = lowerBounds(1);
    Real top = upperBounds(1);
    Real box_size = std::abs(top-bottom);
    Real eps = box_size * 1e-6;


    if(is_default_material_selector) delete material_selector;
    material_selector = new GelMaterialSelector(*this, box_size, "gel", this->nb_gel_pockets, eps);
    is_default_material_selector = false;
  }

  SolidMechanicsModel::initMaterials();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::initSolver(__attribute__((unused)) SolverOptions & options) {
 AKANTU_DEBUG_IN();
#if !defined(AKANTU_USE_MUMPS) && !defined(AKANTU_USE_PETSC)// or other solver in the future \todo add AKANTU_HAS_SOLVER in CMake
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#else
  UInt nb_global_nodes = mesh.getNbGlobalNodes();

  delete jacobian_matrix;
  std::stringstream sstr; sstr << id << ":jacobian_matrix";

#ifdef AKANTU_USE_PETSC
  jacobian_matrix = new PETScMatrix(nb_global_nodes * spatial_dimension, _symmetric, sstr.str(), memory_id);
#else
jacobian_matrix = new SparseMatrix(nb_global_nodes * spatial_dimension, _unsymmetric, sstr.str(), memory_id, 1);
#endif //AKANTU_USE PETSC
  jacobian_matrix->buildProfile(mesh, *dof_synchronizer, spatial_dimension);

  if (!isExplicit()) {
    delete stiffness_matrix;
    std::stringstream sstr_sti; sstr_sti << id << ":stiffness_matrix";
#ifdef AKANTU_USE_PETSC
    stiffness_matrix = new SparseMatrix(nb_global_nodes * spatial_dimension, _symmetric, sstr.str(), memory_id, 1);
    stiffness_matrix->buildProfile(mesh, *dof_synchronizer, spatial_dimension);
#else
    stiffness_matrix = new SparseMatrix(*jacobian_matrix, sstr_sti.str(), memory_id);
#endif //AKANTU_USE_PETSC
  }

  delete solver;
  std::stringstream sstr_solv; sstr_solv << id << ":solver";
#ifdef AKANTU_USE_PETSC
  solver = new SolverPETSc(*jacobian_matrix, sstr_solv.str(), memory_id);
#elif defined(AKANTU_USE_MUMPS)
  solver = new SolverMumps(*jacobian_matrix, sstr_solv.str(), memory_id);
  dof_synchronizer->initScatterGatherCommunicationScheme();
#else
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#endif //AKANTU_USE_MUMPS

  SolverMumpsOptions opt(SolverMumpsOptions::_serial_split);

  if(solver)
    solver->initialize(opt);
#endif //AKANTU_HAS_SOLVER
 AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::initArrays() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = mesh.getNbNodes();
  std::stringstream sstr_disp; sstr_disp << id << ":displacement";
  //  std::stringstream sstr_mass; sstr_mass << id << ":mass";
  std::stringstream sstr_velo; sstr_velo << id << ":velocity";
  std::stringstream sstr_acce; sstr_acce << id << ":acceleration";
  std::stringstream sstr_forc; sstr_forc << id << ":force";
  std::stringstream sstr_resi; sstr_resi << id << ":residual";
  std::stringstream sstr_boun; sstr_boun << id << ":blocked_dofs";

  displacement = &(alloc<Real>(sstr_disp.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  //  mass         = &(alloc<Real>(sstr_mass.str(), nb_nodes, spatial_dimension, 0));
  velocity     = &(alloc<Real>(sstr_velo.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  acceleration = &(alloc<Real>(sstr_acce.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  force        = &(alloc<Real>(sstr_forc.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  residual     = &(alloc<Real>(sstr_resi.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  blocked_dofs = &(alloc<bool>(sstr_boun.str(), nb_nodes, spatial_dimension, false));

  std::stringstream sstr_curp; sstr_curp << id << ":current_position";
  current_position = &(alloc<Real>(sstr_curp.str(), 0, spatial_dimension, REAL_INIT_VALUE));

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt, _ek_not_defined);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, _ek_not_defined);
    for(; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it, gt);
      material_index.alloc(nb_element, 1, *it, gt);
      material_local_numbering.alloc(nb_element, 1, *it, gt);
    }
  }
 
  dof_synchronizer = new DOFSynchronizer(mesh, spatial_dimension,
					 *static_communicator_dummy);
  dof_synchronizer->initLocalDOFEquationNumbers();
  dof_synchronizer->initGlobalDOFEquationNumbers();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::fillCracks(ElementTypeMapReal & saved_damage) {
  const Material & mat_gel = this->getMaterial("gel");
  const Real E_gel = mat_gel.getParam<Real>("E");
  Real E_homogenized = 0.;
  GhostType gt = _not_ghost;
  for (UInt m = 0; m < this->getNbMaterials(); ++m) {
    Material & mat = this->getMaterial(m); 
    if (mat.getName() == "gel" || mat.getName() == "FE2_mat")
      continue;
    const Real E = mat.getParam<Real>("E");
    InternalField<Real> & damage = mat.getInternal<Real>("damage");
    Mesh::type_iterator it =
      mesh.firstType(spatial_dimension, gt, _ek_regular);
    Mesh::type_iterator end =
      mesh.lastType(spatial_dimension, gt, _ek_regular);
    for (; it != end; ++it) {
      const ElementType element_type = *it;
      const Array<UInt> & elem_filter = mat.getElementFilter(element_type, gt);
      if (!elem_filter.getSize())
	continue;
      Array<Real> & damage_vec = damage(element_type, gt);
      Array<Real> & saved_damage_vec = saved_damage(element_type, gt);
      for (UInt i = 0; i < damage_vec.getSize(); ++i) {
	saved_damage_vec(elem_filter(i)) = damage_vec(i);
	E_homogenized = (E_gel - E) * damage_vec(i) + E;
	damage_vec(i) = 1. -(E_homogenized/E);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::drainCracks(const ElementTypeMapReal & saved_damage) {
  GhostType gt = _not_ghost; 
  for (UInt m = 0; m < this->getNbMaterials(); ++m) {
    Material & mat = this->getMaterial(m); 
    if (mat.getName() == "gel" || mat.getName() == "FE2_mat")
      continue;
    else {
      InternalField<Real> & damage = mat.getInternal<Real>("damage");
      Mesh::type_iterator it =
	mesh.firstType(spatial_dimension, gt, _ek_regular);
      Mesh::type_iterator end =
	mesh.lastType(spatial_dimension, gt, _ek_regular);
      for (; it != end; ++it) {
	const ElementType element_type = *it;
	const Array<UInt> & elem_filter = mat.getElementFilter(element_type, gt);
	if (!elem_filter.getSize())
	  continue;
	Array<Real> & damage_vec = damage(element_type, gt);
	const Array<Real> & saved_damage_vec = saved_damage(element_type, gt);
	for (UInt i = 0; i < damage_vec.getSize(); ++i) {
	  damage_vec(i) = saved_damage_vec(elem_filter(i));
	}
      }
    }
  }
}

} // namespace akantu
