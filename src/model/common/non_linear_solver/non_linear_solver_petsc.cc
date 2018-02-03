#include "non_linear_solver_petsc.hh"
#include "mpi_communicator_data.hh"

namespace akantu {

NonLinearSolverPETSc::NonLinearSolverPETSc(
    DOFManagerPETSc & dof_manager,
    const NonLinearSolverType & non_linear_solver_type, const ID & id,
    UInt memory_id)
    : NonLinearSolver(dof_manager, non_linear_solver_type, id, memory_id),
      dof_manager(dof_manager) {
  std::unordered_map<NonLinearSolverType, std::string>
      petsc_non_linear_solver_types{{_nls_newton_raphson, "newtonls"},
                                    {_nls_linear, "ksponly"},
                                    {_nls_gmres, "ngmres"},
                                    {_nls_bfgs, "qn"},
                                    {_nls_cg, "ncg"}};

  for(const auto & pair : petsc_non_linear_solver_types) {
    supported_type.insert(pair.first);
  }

  this->checkIfTypeIsSupported();

  const auto & mpi_data = dynamic_cast<const MPICommunicatorData &>(
      communicator.getCommunicatorData());
  MPI_Comm mpi_comm = mpi_data.getMPICommunicator();

  SNESCreate(mpi_comm, &snes);

  auto it = petsc_non_linear_solver_types.find(non_linear_solver_type);
  if(it != petsc_non_linear_solver_types.end()) {
    SNESSetType(snes, it->second);
  }
}

void NonLinearSolverPETSc::parseSection(const ParserSection & section) {
  auto parameters = section.getParameters();
  for(auto && param : range(parameters.first, parameters.second)) {
    PetscOptionSetValue(NULL, parm.getName().c_str(), parm.getValue().c_str());
  }

  SNESSetFromOptions(snes);

  PetscOptionClear(NULL);
}

} // akantu
