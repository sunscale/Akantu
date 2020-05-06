/**
 * @file   test_dof_data_accessor.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Dec 09 2014
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  data accessor class for testing the
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

#include "aka_common.hh"
#include "data_accessor.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

class TestDOFAccessor : public DataAccessor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  inline TestDOFAccessor(const Array<Int> & global_dof_equation_numbers);

  /* ------------------------------------------------------------------------ */
  /* Ghost Synchronizer inherited members                                     */
  /* ------------------------------------------------------------------------ */
protected:
  inline UInt getNbDataForDOFs(const Array<UInt> & dofs,
                               SynchronizationTag tag) const;
  inline void packDOFData(CommunicationBuffer & buffer,
                          const Array<UInt> & dofs,
                          SynchronizationTag tag) const;
  inline void unpackDOFData(CommunicationBuffer & buffer,
                            const Array<UInt> & dofs, SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  const Array<Int> & global_dof_equation_numbers;
};

/* -------------------------------------------------------------------------- */
/* TestDOFSynchronizer implementation */
/* -------------------------------------------------------------------------- */
inline TestDOFAccessor::TestDOFAccessor(
    const Array<Int> & global_dof_equation_numbers)
    : global_dof_equation_numbers(global_dof_equation_numbers) {}

inline UInt TestDOFAccessor::getNbDataForDOFs(const Array<UInt> & dofs,
                                              __attribute__((unused))
                                              SynchronizationTag tag) const {
  if (dofs.size())
    // return Mesh::getSpatialDimension(elements(0).type) * sizeof(Real) *
    // elements.size();
    return sizeof(Int) * dofs.size();
  else
    return 0;
}

inline void TestDOFAccessor::packDOFData(CommunicationBuffer & buffer,
                                         const Array<UInt> & dofs,
                                         __attribute__((unused))
                                         SynchronizationTag tag) const {
  Array<UInt>::const_scalar_iterator bit = dofs.begin();
  Array<UInt>::const_scalar_iterator bend = dofs.end();
  for (; bit != bend; ++bit) {
    buffer << this->global_dof_equation_numbers[*bit];
  }
}

inline void TestDOFAccessor::unpackDOFData(CommunicationBuffer & buffer,
                                           const Array<UInt> & dofs,
                                           __attribute__((unused))
                                           SynchronizationTag tag) {
  Array<UInt>::const_scalar_iterator bit = dofs.begin();
  Array<UInt>::const_scalar_iterator bend = dofs.end();
  for (; bit != bend; ++bit) {
    Int global_dof_eq_nb_local = global_dof_equation_numbers[*bit];
    Int global_dof_eq_nb = 0;
    buffer >> global_dof_eq_nb;
    std::cout << *bit << global_dof_eq_nb_local << std::endl;
    Real tolerance = Math::getTolerance();
    if (!(std::abs(global_dof_eq_nb - global_dof_eq_nb_local) <= tolerance))
      AKANTU_ERROR(
          "Unpacking an unknown value for the dof: "
          << *bit << "(global_dof_equation_number = " << global_dof_eq_nb_local
          << " and buffer = " << global_dof_eq_nb << ") - tag: " << tag);
  }
}

} // namespace akantu
