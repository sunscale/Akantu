/**
 * @file   solver_mumps.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  implem of SolverMumps class
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
 * @section DESCRIPTION
 *
 * @subsection Ctrl_param Control parameters
 *
 * ICNTL(1),
 * ICNTL(2),
 * ICNTL(3) : output streams for error, diagnostics, and global messages
 *
 * ICNTL(4) : verbose level : 0 no message - 4 all messages
 *
 * ICNTL(5) : type of matrix, 0 assembled, 1 elementary
 *
 * ICNTL(6) : control  the permutation and scaling(default 7)  see mumps doc for
 * more information
 *
 * ICNTL(7) : determine  the pivot  order (default  7) see  mumps doc  for more
 * information
 *
 * ICNTL(8) : describe the scaling method used
 *
 * ICNTL(9) : 1 solve A x = b, 0 solve At x = b
 *
 * ICNTL(10) : number of iterative refinement when NRHS = 1
 *
 * ICNTL(11) : > 0 return statistics
 *
 * ICNTL(12) : only used for SYM = 2, ordering strategy
 *
 * ICNTL(13) :
 *
 * ICNTL(14) : percentage of increase of the estimated working space
 *
 * ICNTL(15-17) : not used
 *
 * ICNTL(18) : only  used if ICNTL(5) = 0, 0 matrix  centralized, 1 structure on
 * host and mumps  give the mapping, 2 structure on  host and distributed matrix
 * for facto, 3 distributed matrix
 *
 * ICNTL(19) : > 0, Shur complement returned
 *
 * ICNTL(20) : 0 rhs dense, 1 rhs sparse
 *
 * ICNTL(21) : 0 solution in rhs, 1 solution distributed in ISOL_loc and SOL_loc
 * allocated by user
 *
 * ICNTL(22) : 0 in-core, 1 out-of-core
 *
 * ICNTL(23) : maximum memory allocatable by mumps pre proc
 *
 * ICNTL(24) : controls the detection of "null pivot rows"
 *
 * ICNTL(25) :
 *
 * ICNTL(26) :
 *
 * ICNTL(27) :
 *
 * ICNTL(28) : 0 automatic choice, 1 sequential analysis, 2 parallel analysis
 *
 * ICNTL(29) : 0 automatic choice, 1 PT-Scotch, 2 ParMetis
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"

#if defined(AKANTU_USE_MPI)
#  include "static_communicator_mpi.hh"
#  include "mpi_type_wrapper.hh"
#endif

#include "solver_mumps.hh"
#include "dof_synchronizer.hh"


/* -------------------------------------------------------------------------- */
// static std::ostream & operator <<(std::ostream & stream, const DMUMPS_STRUC_C & _this) {
//   stream << "DMUMPS Data [" << std::endl;
//   stream << " + job          : " << _this.job          << std::endl;
//   stream << " + par          : " << _this.par          << std::endl;
//   stream << " + sym          : " << _this.sym          << std::endl;
//   stream << " + comm_fortran : " << _this.comm_fortran << std::endl;
//   stream << " + nz           : " << _this.nz           << std::endl;
//   stream << " + irn          : " << _this.irn          << std::endl;
//   stream << " + jcn          : " << _this.jcn          << std::endl;
//   stream << " + nz_loc       : " << _this.nz_loc       << std::endl;
//   stream << " + irn_loc      : " << _this.irn_loc      << std::endl;
//   stream << " + jcn_loc      : " << _this.jcn_loc      << std::endl;
//   stream << "]";
//   return stream;
// }

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
SolverMumps::SolverMumps(SparseMatrix & matrix,
                         const ID & id,
                         const MemoryID & memory_id) :
  Solver(matrix, id, memory_id), is_mumps_data_initialized(false), rhs_is_local(true) {
  AKANTU_DEBUG_IN();

#ifdef AKANTU_USE_MPI
  parallel_method = SolverMumpsOptions::_fully_distributed;
#else //AKANTU_USE_MPI
  parallel_method = SolverMumpsOptions::_not_parallel;
#endif //AKANTU_USE_MPI

  CommunicatorEventHandler & comm_event_handler = *this;

  communicator.registerEventHandler(comm_event_handler);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolverMumps::~SolverMumps() {
  AKANTU_DEBUG_IN();
  this->destroyMumpsData();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::destroyMumpsData() {
  AKANTU_DEBUG_IN();
  
  if(this->is_mumps_data_initialized) {
    this->mumps_data.job = _smj_destroy; // destroy
    dmumps_c(&this->mumps_data);
    this->is_mumps_data_initialized = false;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::onCommunicatorFinalize(const StaticCommunicator & comm) {
  AKANTU_DEBUG_IN();

  try{
    destroyMumpsData();
  } catch(...) {}

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::initMumpsData() {
  // Default Scaling
  icntl(8) = 77;

  // Assembled matrix
  icntl(5) = 0;

  /// Default centralized dense second member
  icntl(20) = 0;
  icntl(21) = 0;

  icntl(28) = 0; //automatic choice for analysis analysis

  switch(this->parallel_method) {
  case SolverMumpsOptions::_fully_distributed:
      icntl(18) = 3; //fully distributed

      this->mumps_data.nz_loc  = matrix->getNbNonZero();
      this->mumps_data.irn_loc = matrix->getIRN().storage();
      this->mumps_data.jcn_loc = matrix->getJCN().storage();
      break;
  case SolverMumpsOptions::_not_parallel:
  case SolverMumpsOptions::_master_slave_distributed:
    icntl(18) = 0; //centralized

    if(prank == 0) {
      this->mumps_data.nz  = matrix->getNbNonZero();
      this->mumps_data.irn = matrix->getIRN().storage();
      this->mumps_data.jcn = matrix->getJCN().storage();
    } else {
      this->mumps_data.nz  = 0;
      this->mumps_data.irn = NULL;
      this->mumps_data.jcn = NULL;
    }
    break;
  default:
    AKANTU_DEBUG_ERROR("This case should not happen!!");
  }
}

/* -------------------------------------------------------------------------- */
void SolverMumps::initialize(SolverOptions & options) {
  AKANTU_DEBUG_IN();

  if(SolverMumpsOptions * opt = dynamic_cast<SolverMumpsOptions *>(&options)) {
    this->parallel_method = opt->parallel_method;
  }

  this->mumps_data.par = 1; // The host is part of computations

  switch(this->parallel_method) {
  case SolverMumpsOptions::_not_parallel: break;
  case SolverMumpsOptions::_master_slave_distributed:
    this->mumps_data.par = 0; // The host is not part of the computations
  case SolverMumpsOptions::_fully_distributed:
#ifdef AKANTU_USE_MPI
    const StaticCommunicatorMPI & mpi_st_comm = dynamic_cast<const StaticCommunicatorMPI &>(communicator.getRealStaticCommunicator());
    this->mumps_data.comm_fortran = MPI_Comm_c2f(mpi_st_comm.getMPITypeWrapper().getMPICommunicator());
#endif
    break;
  }


  this->mumps_data.sym = 2 * (matrix->getSparseMatrixType() == _symmetric);
  this->prank = communicator.whoAmI();

  this->mumps_data.job = _smj_initialize; //initialize
  dmumps_c(&this->mumps_data);

  this->is_mumps_data_initialized = true;

  /* ------------------------------------------------------------------------ */
  UInt size = matrix->getSize();

  if(prank == 0) {
    std::stringstream sstr_rhs; sstr_rhs << id << ":rhs";
    this->rhs = &(alloc<Real>(sstr_rhs.str(), size, 1, 0.));
  } else {
    this->rhs = NULL;
  }

  this->mumps_data.nz_alloc = 0;
  this->mumps_data.n        = size;
  /* ------------------------------------------------------------------------ */
  // Output setup
  if(AKANTU_DEBUG_TEST(dblTrace)) {
    icntl(1) = 6;
    icntl(2) = 2;
    icntl(3) = 2;
    icntl(4) = 4;
  } else {
    /// No outputs
    icntl(1) = 6; // error output
    icntl(2) = 0; // dignostics output
    icntl(3) = 0; // informations
    icntl(4) = 0; // no outputs
  }

  if(AKANTU_DEBUG_TEST(dblDump)) {
    strcpy(this->mumps_data.write_problem, "mumps_matrix.mtx");
  }

  this->analysis();

//  icntl(14) = 80;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::setRHS(const Array<Real> & rhs) {
  if(prank == 0) {

    std::copy(rhs.storage(), rhs.storage() + this->rhs->getSize(), this->rhs->storage());

    //    DebugLevel dbl = debug::getDebugLevel();
//    debug::setDebugLevel(dblError);
//    matrix->getDOFSynchronizer().gather(rhs, 0, this->rhs);
//    debug::setDebugLevel(dbl);

  } else {
    this->matrix->getDOFSynchronizer().gather(rhs, 0);
  }
}

/* -------------------------------------------------------------------------- */
void SolverMumps::analysis() {
  AKANTU_DEBUG_IN();

  initMumpsData();

  this->mumps_data.job = _smj_analyze; //analyze
  dmumps_c(&this->mumps_data);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::factorize() {
  AKANTU_DEBUG_IN();

  if(parallel_method == SolverMumpsOptions::_fully_distributed)
    this->mumps_data.a_loc  = this->matrix->getA().storage();
  else {
    if(prank == 0)
      this->mumps_data.a  = this->matrix->getA().storage();
  }

  if(prank == 0) {
    this->mumps_data.rhs = this->rhs->storage();
  }


  this->mumps_data.job = _smj_factorize; // factorize
  dmumps_c(&this->mumps_data);

  this->printError();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::solve() {
  AKANTU_DEBUG_IN();

  if(prank == 0) {
    this->mumps_data.rhs = this->rhs->storage();
  }

  this->mumps_data.job = _smj_solve; // solve
  dmumps_c(&this->mumps_data);

  this->printError();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::solve(Array<Real> & solution) {
  AKANTU_DEBUG_IN();

  this->solve();

  if(prank == 0) {
//    matrix->getDOFSynchronizer().scatter(solution, 0, this->rhs);
    std::copy(this->rhs->storage(), this->rhs->storage() + this->rhs->getSize(), solution.storage());

  } else {
    this->matrix->getDOFSynchronizer().scatter(solution, 0);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::printError() {
  if(info(1) != 0) {
    communicator.allReduce(&info(1), 1, _so_min);
    switch(info(1)) {     
    case -10: AKANTU_DEBUG_ERROR("The matrix is singular"); break;
    case  -9: {
      icntl(14) += 10;
      if(icntl(14) != 90) {
	//std::cout << "Dynamic memory increase of 10%" << std::endl;
      	this->analysis();
      	this->factorize();
      	this->solve();
      } else {
	AKANTU_DEBUG_ERROR("The MUMPS workarray is too small INFO(2)=" << info(2) << "No further increase possible"); break;
      }
    }
    default:
      AKANTU_DEBUG_ERROR("Error in mumps during solve process, check mumps user guide INFO(1) ="
                         << info(1));
    }
  }
}


__END_AKANTU__
