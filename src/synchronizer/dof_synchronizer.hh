/**
 * @file   dof_synchronizer.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 17 2011
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Synchronize Array of DOFs
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_array.hh"
#include "aka_common.hh"
#include "synchronizer_impl.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
class Mesh;
class DOFManagerDefault;
} // namespace akantu

#ifndef __AKANTU_DOF_SYNCHRONIZER_HH__
#define __AKANTU_DOF_SYNCHRONIZER_HH__

namespace akantu {

class DOFSynchronizer : public SynchronizerImpl<UInt> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DOFSynchronizer(DOFManagerDefault & dof_manager,
                  const ID & id = "dof_synchronizer", MemoryID memory_id = 0);
  ~DOFSynchronizer() override;

  virtual void registerDOFs(const ID & dof_id);
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  template <typename T>
  /// Gather the DOF value on the root proccessor
  void gather(const Array<T> & to_gather, Array<T> & gathered);

  /// Gather the DOF value on the root proccessor
  template <typename T> void gather(const Array<T> & to_gather);

  /// Scatter a DOF Array form root to all processors
  template <typename T>
  void scatter(Array<T> & scattered, const Array<T> & to_scatter);

  /// Scatter a DOF Array form root to all processors
  template <typename T> void scatter(Array<T> & scattered);

  /// Uses the synchronizer to perform a reduction on the vector
  template <template <class> class Op, typename T>
  void reduceSynchronize(Array<T> & array) const;

  template <typename T> void synchronize(Array<T> & vector) const;

  void onNodesAdded(const Array<UInt> & nodes);

protected:
  /// check if dof changed set on at least one processor
  bool hasChanged();

  /// init the scheme for scatter and gather operation, need extra memory
  void initScatterGatherCommunicationScheme();

  Int getRank(const UInt & /*node*/) const final { AKANTU_TO_IMPLEMENT(); }

private:
  /// Root processor for scatter/gather operations
  Int root;

  /// information on the dofs
  DOFManagerDefault & dof_manager;

  /// dofs from root
  Array<UInt> root_dofs;

  /// Dofs received from slaves proc (only on master)
  std::map<UInt, Array<UInt>> master_receive_dofs;

  bool dof_changed;
};
} // namespace akantu

#include "dof_synchronizer_inline_impl.cc"

#endif /* __AKANTU_DOF_SYNCHRONIZER_HH__ */
