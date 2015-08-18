/**
 * @file   dof_manager.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Jul 22 11:43:43 2015
 *
 * @brief  Class handling the different types of dofs
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
#include "aka_memory.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DOF_MANAGER_HH__
#define __AKANTU_DOF_MANAGER_HH__

__BEGIN_AKANTU__

class DOFManager : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DOFManager(const Mesh & mesh, const ID & id = "dof_manager",
             const MemoryID & memory_id = 0);
  virtual ~DOFManager();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void registerDOFs(const ID & dof_id, Array<Real> & dofs_array);

  /// Get the part of the solution corresponding to the dof_id
  virtual void getSolution(const ID & dof_id, Array<Real> & solution_array) = 0;

  /// Assemble an array to the global residual array
  virtual void assembleToResidual(const ID & dof_id,
                                  const Array<Real> & array_to_assemble,
                                  Real scale_factor = 1.) = 0;

  /**
   * Assemble elementary values to a local array of the size nb_nodes *
   * nb_dof_per_node. The dof number is implicitly considered as
   * conn(el, n) * nb_nodes_per_element + d.
   * With 0 < n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalArrayLocalArray(
      const Array<Real> & elementary_vect, Array<Real> & array_assembeled,
      const ElementType & type, const GhostType & ghost_type, Real scale_factor,
      const Array<UInt> & filter_elements);

  /**
   * Assemble elementary values to the global residual array. The dof number is
   * implicitly considered as conn(el, n) * nb_nodes_per_element + d.
   * With 0 < n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalArrayResidual(
      const ID & dof_id, const Array<Real> & elementary_vect,
      const ElementType & type, const GhostType & ghost_type, Real scale_factor,
      const Array<UInt> & filter_elements);

  /**
   * Assemble elementary values to the global residual array. The dof number is
   * implicitly considered as conn(el, n) * nb_nodes_per_element + d.  With 0 <
   * n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void
  assembleElementalMatricesToMatrix(const ID & matrix_id, const ID & dof_id,
                                    const Array<Real> & elemental_mat) = 0;

  /// notation fully defined yet...
  virtual void assemblePreassembledMatrix(const ID & matrix_id,
                                          const ID & dof_id_m,
                                          const ID & dof_id_n,
                                          const Matrix<Real> & matrix) = 0;

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// fill a Vector with the equation numbers corresponding to the given
  /// connectivity
  inline void extractElementEquationNumber(
      const Array<UInt> & equation_numbers, const Vector<UInt> & connectivity,
      UInt nb_degree_of_freedom, Vector<UInt> & local_equation_number);

  /// register a matrix
  void registerSparseMatrix(const ID & matrix_id, SparseMatrix & matrix);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the equation numbers corresponding to a dof ID
  Array<Real> & getEquationNumbers(const ID & dof_id);

  const Array<Real> & getDOFs(const ID & id) const;

  /// Get an instance of a new SparseMatrix
  virtual SparseMatrix & getNewMatrix(const ID & matrix_id,
                                      const MatrixType & matrix_type);

  /// Get an instance of a new SparseMatrix as a copy of the SparseMatrix
  /// matrix_to_copy_id
  virtual SparseMatrix & getNewMatrix(const ID & matrix_id,
                                      const ID & matrix_to_copy_id);

  /// Get the reference of an existing matrix
  SparseMatrix & getMatrix(const ID & matrix_id);

  AKANTU_GET_MACRO(SystemSize, this->system_size, UInt);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// store a reference to the dof arrays
  std::map<ID, Array<Real> *> dofs;
  /// equation numbers corresponding to the dofglobalids arrays
  std::map<ID, Array<UInt> *> equation_numbers;

  /// list of sparse matrices that where created
  std::map<ID, SparseMatrix *> matrices;

  /// reference to the underlying mesh
  const Mesh & mesh;

  /// Total number of degrees of freedom
  UInt local_system_size;

  /// Total number of degrees of freedom
  UInt system_size;
};

__END_AKANTU__

#endif /* __AKANTU_DOF_MANAGER_HH__ */
