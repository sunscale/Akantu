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
  void onNodesAdded(const Array<UInt> & nodes);

protected:
  Int getRank(const UInt & /*node*/) const final { AKANTU_TO_IMPLEMENT(); }

  /// list the entities to send to root process
  void fillEntityToSend(Array<UInt> & dofs_to_send) override;

  inline UInt canScatterSize() override;
  inline UInt gatheredSize() override;

  inline UInt localToGlobalEntity(const UInt & local) override;

private:
  /// information on the dofs
  DOFManagerDefault & dof_manager;
};

} // namespace akantu

#include "dof_synchronizer_inline_impl.cc"

#endif /* __AKANTU_DOF_SYNCHRONIZER_HH__ */
