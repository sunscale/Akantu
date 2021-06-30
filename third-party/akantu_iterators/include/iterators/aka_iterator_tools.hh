/**
 * @file   aka_iterator_tools.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  jeu déc 12 2019
 *
 * @brief A Documented file.
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * akantu-iterators is free  software: you can redistribute it and/or  modify it
 * under the terms  of the  GNU Lesser  General Public  License as  published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * akantu-iterators is  distributed in the  hope that it  will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public
 * License  for more details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with akantu-iterators. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */

#ifndef AKA_ITERATOR_TOOLS_H
#define AKA_ITERATOR_TOOLS_H

#ifndef AKANTU_ITERATORS_NAMESPACE
#define AKANTU_ITERATORS_NAMESPACE akantu
#endif

namespace AKANTU_ITERATORS_NAMESPACE {

/* -------------------------------------------------------------------------- */
namespace iterators {

  namespace details {
    template <bool enable> struct CopyAssignmentEnabler {};
    template <> struct CopyAssignmentEnabler<false> {
      CopyAssignmentEnabler() = default;
      CopyAssignmentEnabler(const CopyAssignmentEnabler &) = default;
      CopyAssignmentEnabler(CopyAssignmentEnabler &&) = default;
      auto operator=(const CopyAssignmentEnabler &)
          -> CopyAssignmentEnabler & = delete;
      auto operator=(CopyAssignmentEnabler &&)
          -> CopyAssignmentEnabler & = default;
    };

    template <bool enable> struct MoveAssignmentEnabler {};
    template <> struct MoveAssignmentEnabler<false> {
      MoveAssignmentEnabler() = default;
      MoveAssignmentEnabler(const MoveAssignmentEnabler &) = default;
      MoveAssignmentEnabler(MoveAssignmentEnabler &&) = default;
      auto operator=(const MoveAssignmentEnabler &)
          -> MoveAssignmentEnabler & = delete;
      auto operator=(MoveAssignmentEnabler &&)
          -> MoveAssignmentEnabler & = default;
    };

  } // namespace details
} // namespace iterators
} // namespace AKANTU_ITERATORS_NAMESPACE

#endif // AKA_ITERATOR_TOOLS_H
