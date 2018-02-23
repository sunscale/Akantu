/**
 * @file   dof_manager_petsc.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Oct 07 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  DOFManaterPETSc is the PETSc implementation of the DOFManager
 *
 * @section LICENSE
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
#include "dof_manager_petsc.hh"

#include "cppargparse.hh"

#if defined(AKANTU_USE_MPI)
#include "mpi_type_wrapper.hh"
#include "static_communicator.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <petscsys.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

#if not defined(PETSC_CLANGUAGE_CXX)
/// small hack to use the c binding of petsc when the cxx binding does notation
/// exists
int aka_PETScError(int ierr) {
  CHKERRQ(ierr);
  return 0;
}
#endif

UInt DOFManagerPETSc::petsc_dof_manager_instances = 0;

/// Error handler to make PETSc errors caught by Akantu
#if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5
static PetscErrorCode PETScErrorHandler(MPI_Comm, int line, const char * dir,
                                        const char * file,
                                        PetscErrorCode number,
                                        PetscErrorType type,
                                        const char * message, void *) {
  AKANTU_ERROR("An error occured in PETSc in file \""
               << file << ":" << line << "\" - PetscErrorCode " << number
               << " - \"" << message << "\"");
}
#else
static PetscErrorCode PETScErrorHandler(MPI_Comm, int line, const char * func,
                                        const char * dir, const char * file,
                                        PetscErrorCode number,
                                        PetscErrorType type,
                                        const char * message, void *) {
  AKANTU_ERROR("An error occured in PETSc in file \""
               << file << ":" << line << "\" - PetscErrorCode " << number
               << " - \"" << message << "\"");
}
#endif

/* -------------------------------------------------------------------------- */
DOFManagerPETSc::DOFManagerPETSc(const ID & id, const MemoryID & memory_id)
    : DOFManager(id, memory_id) {

// check if the akantu types and PETSc one are consistant
#if __cplusplus > 199711L
  static_assert(sizeof(Int) == sizeof(PetscInt),
                "The integer type of Akantu does not match the one from PETSc");
  static_assert(sizeof(Real) == sizeof(PetscReal),
                "The integer type of Akantu does not match the one from PETSc");
#else
  AKANTU_DEBUG_ASSERT(
      sizeof(Int) == sizeof(PetscInt),
      "The integer type of Akantu does not match the one from PETSc");
  AKANTU_DEBUG_ASSERT(
      sizeof(Real) == sizeof(PetscReal),
      "The integer type of Akantu does not match the one from PETSc");
#endif

  if (this->petsc_dof_manager_instances == 0) {
#if defined(AKANTU_USE_MPI)
    StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
    const StaticCommunicatorMPI & mpi_st_comm =
        dynamic_cast<const StaticCommunicatorMPI &>(
            comm.getRealStaticCommunicator());

    this->mpi_communicator =
        mpi_st_comm.getMPITypeWrapper().getMPICommunicator();
#else
    this->mpi_communicator = PETSC_COMM_SELF;
#endif

    cppargparse::ArgumentParser & argparser = getStaticArgumentParser();
    int & argc = argparser.getArgC();
    char **& argv = argparser.getArgV();

    PetscErrorCode petsc_error = PetscInitialize(&argc, &argv, NULL, NULL);

    if (petsc_error != 0) {
      AKANTU_ERROR("An error occured while initializing Petsc (PetscErrorCode "
                   << petsc_error << ")");
    }

    PetscPushErrorHandler(PETScErrorHandler, NULL);
    this->petsc_dof_manager_instances++;
  }

  VecCreate(PETSC_COMM_WORLD, &this->residual);
  VecCreate(PETSC_COMM_WORLD, &this->solution);
}

/* -------------------------------------------------------------------------- */
DOFManagerPETSc::~DOFManagerPETSc() {
  PetscErrorCode ierr;
  ierr = VecDestroy(&(this->residual));
  CHKERRXX(ierr);

  ierr = VecDestroy(&(this->solution));
  CHKERRXX(ierr);

  this->petsc_dof_manager_instances--;
  if (this->petsc_dof_manager_instances == 0) {
    PetscFinalize();
  }
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                                   DOFSupportType & support_type) {
  DOFManager::registerDOFs(dof_id, dofs_array, support_type);

  PetscErrorCode ierr;

  PetscInt current_size;
  ierr = VecGetSize(this->residual, &current_size);
  CHKERRXX(ierr);

  if (current_size == 0) { // first time vector is set
    PetscInt local_size = this->pure_local_system_size;
    ierr = VecSetSizes(this->residual, local_size, PETSC_DECIDE);
    CHKERRXX(ierr);

    ierr = VecSetFromOptions(this->residual);
    CHKERRXX(ierr);

#ifndef AKANTU_NDEBUG
    PetscInt global_size;
    ierr = VecGetSize(this->residual, &global_size);
    CHKERRXX(ierr);
    AKANTU_DEBUG_ASSERT(this->system_size == UInt(global_size),
                        "The local value of the system size does not match the "
                        "one determined by PETSc");
#endif
    PetscInt start_dof, end_dof;
    VecGetOwnershipRange(this->residual, &start_dof, &end_dof);

    PetscInt * global_indices = new PetscInt[local_size];
    global_indices[0] = start_dof;

    for (PetscInt d = 0; d < local_size; d++)
      global_indices[d + 1] = global_indices[d] + 1;

// To be change if we switch to a block definition
#if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5
    ISLocalToGlobalMappingCreate(this->communicator, 1, local_size,
                                 global_indices, PETSC_COPY_VALUES,
                                 &this->is_ltog);

#else
    ISLocalToGlobalMappingCreate(this->communicator, local_size, global_indices,
                                 PETSC_COPY_VALUES, &this->is_ltog);
#endif

    VecSetLocalToGlobalMapping(this->residual, this->is_ltog);
    delete[] global_indices;

    ierr = VecDuplicate(this->residual, &this->solution);
    CHKERRXX(ierr);

  } else { // this is an update of the object already created
    AKANTU_TO_IMPLEMENT();
  }

  /// set the solution to zero
  // ierr = VecZeroEntries(this->solution);
  // CHKERRXX(ierr);
}

/* -------------------------------------------------------------------------- */
/**
 * This function creates the non-zero pattern of the PETSc matrix. In
 * PETSc the parallel matrix is partitioned across processors such
 * that the first m0 rows belong to process 0, the next m1 rows belong
 * to process 1, the next m2 rows belong to process 2 etc.. where
 * m0,m1,m2,.. are the input parameter 'm'. i.e each processor stores
 * values corresponding to [m x N] submatrix
 * (http://www.mcs.anl.gov/petsc/).
 * @param mesh mesh discretizing the domain we want to analyze
 * @param dof_synchronizer dof synchronizer that maps the local
 * dofs to the global dofs and the equation numbers, i.e., the
 * position at which the dof is assembled in the matrix
 */

// void SparseMatrixPETSc::buildProfile(const Mesh & mesh,
//                                      const DOFSynchronizer &
//                                      dof_synchronizer,
//                                      UInt nb_degree_of_freedom) {
//   AKANTU_DEBUG_IN();

//   // clearProfile();
//   this->dof_synchronizer = &const_cast<DOFSynchronizer &>(dof_synchronizer);
//   this->setSize();
//   PetscErrorCode ierr;

//   /// resize arrays to store the number of nonzeros in each row
//   this->d_nnz.resize(this->local_size);
//   this->o_nnz.resize(this->local_size);

//   /// set arrays to zero everywhere
//   this->d_nnz.set(0);
//   this->o_nnz.set(0);

//   // if(irn_jcn_to_k) delete irn_jcn_to_k;
//   // irn_jcn_to_k = new std::map<std::pair<UInt, UInt>, UInt>;

//   coordinate_list_map::iterator irn_jcn_k_it;

//   Int * eq_nb_val = dof_synchronizer.getGlobalDOFEquationNumbers().storage();
//   UInt nb_global_dofs = dof_synchronizer.getNbGlobalDOFs();
//   Array<Int> index_pair(2);

//   /// Loop over all the ghost types
//   for (ghost_type_t::iterator gt = ghost_type_t::begin();
//        gt != ghost_type_t::end(); ++gt) {
//     const GhostType & ghost_type = *gt;
//     Mesh::type_iterator it =
//         mesh.firstType(mesh.getSpatialDimension(), ghost_type,
//         _ek_not_defined);
//     Mesh::type_iterator end =
//         mesh.lastType(mesh.getSpatialDimension(), ghost_type,
//         _ek_not_defined);
//     for (; it != end; ++it) {
//       UInt nb_element = mesh.getNbElement(*it, ghost_type);
//       UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
//       UInt size_mat = nb_nodes_per_element * nb_degree_of_freedom;

//       UInt * conn_val = mesh.getConnectivity(*it, ghost_type).storage();
//       Int * local_eq_nb_val =
//           new Int[nb_degree_of_freedom * nb_nodes_per_element];

//       for (UInt e = 0; e < nb_element; ++e) {
//         Int * tmp_local_eq_nb_val = local_eq_nb_val;
//         for (UInt i = 0; i < nb_nodes_per_element; ++i) {
//           UInt n = conn_val[i];
//           for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
//             /**
//              * !!!!!!Careful!!!!!! This is a ugly fix. @todo this is a
//              * very ugly fix, because the offset for the global
//              * equation number, where the dof will be assembled, is
//              * hardcoded. In the future a class dof manager has to be
//              * added to Akantu to handle the mapping between the dofs
//              * and the equation numbers
//              *
//              */

//             *tmp_local_eq_nb_val++ =
//                 eq_nb_val[n * nb_degree_of_freedom + d] -
//                 (dof_synchronizer.isPureGhostDOF(n * nb_degree_of_freedom +
//                 d)
//                      ? nb_global_dofs
//                      : 0);
//           }
//         }

//         for (UInt i = 0; i < size_mat; ++i) {
//           Int c_irn = local_eq_nb_val[i];
//           UInt j_start = 0;
//           for (UInt j = j_start; j < size_mat; ++j) {
//             Int c_jcn = local_eq_nb_val[j];
//             index_pair(0) = c_irn;
//             index_pair(1) = c_jcn;
//             AOApplicationToPetsc(this->petsc_matrix_wrapper->ao, 2,
//                                  index_pair.storage());
//             if (index_pair(0) >= first_global_index &&
//                 index_pair(0) < first_global_index + this->local_size) {
//               KeyCOO irn_jcn = keyPETSc(c_irn, c_jcn);
//               irn_jcn_k_it = irn_jcn_k.find(irn_jcn);

//               if (irn_jcn_k_it == irn_jcn_k.end()) {
//                 irn_jcn_k[irn_jcn] = nb_non_zero;

//                 // check if node is slave node
//                 if (index_pair(1) >= first_global_index &&
//                     index_pair(1) < first_global_index + this->local_size)
//                   this->d_nnz(index_pair(0) - first_global_index) += 1;
//                 else
//                   this->o_nnz(index_pair(0) - first_global_index) += 1;
//                 nb_non_zero++;
//               }
//             }
//           }
//         }
//         conn_val += nb_nodes_per_element;
//       }

//       delete[] local_eq_nb_val;
//     }
//   }

//   // /// for pbc @todo correct it for parallel
//   // if(StaticCommunicator::getStaticCommunicator().getNbProc() == 1) {
//   //   for (UInt i = 0; i < size; ++i) {
//   //    KeyCOO irn_jcn = key(i, i);
//   //    irn_jcn_k_it = irn_jcn_k.find(irn_jcn);
//   //    if(irn_jcn_k_it == irn_jcn_k.end()) {
//   //      irn_jcn_k[irn_jcn] = nb_non_zero;
//   //      irn.push_back(i + 1);
//   //      jcn.push_back(i + 1);
//   //      nb_non_zero++;
//   //    }
//   //   }
//   // }

//   // std::string mat_type;
//   // mat_type.resize(20, 'a');
//   // std::cout << "MatType: " << mat_type << std::endl;
//   // const char * mat_type_ptr = mat_type.c_str();
//   MatType type;
//   MatGetType(this->petsc_matrix_wrapper->mat, &type);
//   /// std::cout << "the matrix type is: " << type << std::endl;
//   /**
//    * PETSc will use only one of the following preallocation commands
//    * depending on the matrix type and ignore the rest. Note that for
//    * the block matrix format a block size of 1 is used. This might
//    * result in bad performance. @todo For better results implement
//    * buildProfile() with larger block size.
//    *
//    */
//   /// build profile:
//   if (strcmp(type, MATSEQAIJ) == 0) {
//     ierr = MatSeqAIJSetPreallocation(this->petsc_matrix_wrapper->mat, 0,
//                                      d_nnz.storage());
//     CHKERRXX(ierr);
//   } else if ((strcmp(type, MATMPIAIJ) == 0)) {
//     ierr = MatMPIAIJSetPreallocation(this->petsc_matrix_wrapper->mat, 0,
//                                      d_nnz.storage(), 0, o_nnz.storage());
//     CHKERRXX(ierr);
//   } else {
//     AKANTU_ERROR("The type " << type
//                                    << " of PETSc matrix is not handled by"
//                                    << " akantu in the preallocation step");
//   }

//   // ierr =  MatSeqSBAIJSetPreallocation(this->petsc_matrix_wrapper->mat, 1,
//   //                                  0, d_nnz.storage()); CHKERRXX(ierr);

//   if (this->sparse_matrix_type == _symmetric) {
//     /// set flag for symmetry to enable ICC/Cholesky preconditioner
//     ierr = MatSetOption(this->petsc_matrix_wrapper->mat, MAT_SYMMETRIC,
//                         PETSC_TRUE);
//     CHKERRXX(ierr);
//     /// set flag for symmetric positive definite
//     ierr = MatSetOption(this->petsc_matrix_wrapper->mat, MAT_SPD,
//     PETSC_TRUE);
//     CHKERRXX(ierr);
//   }
//   /// once the profile has been build ignore any new nonzero locations
//   ierr = MatSetOption(this->petsc_matrix_wrapper->mat,
//                       MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
//   CHKERRXX(ierr);

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */

} // akantu
