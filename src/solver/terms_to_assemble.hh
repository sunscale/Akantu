/**
 * @file   terms_to_assemble.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  List of terms to assemble to a matrix
 *
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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_TERMS_TO_ASSEMBLE_HH_
#define AKANTU_TERMS_TO_ASSEMBLE_HH_

namespace akantu {

class TermsToAssemble {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TermsToAssemble() = default;
  virtual ~TermsToAssemble() = default;

  class TermToAssemble {
  public:
    TermToAssemble(UInt i, UInt j) : _i(i), _j(j), val(0.) {}
    inline TermToAssemble & operator=(Real val) {
      this->val = val;
      return *this;
    }
    inline TermToAssemble operator+=(Real val) {
      this->val += val;
      return *this;
    }
    inline operator Real() const { return val; }
    inline UInt i() const { return _i; }
    inline UInt j() const { return _j; }

  private:
    UInt _i, _j;
    Real val;
  };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline TermToAssemble & operator()(UInt i, UInt j) {
    terms.emplace_back(i, j);
    return terms.back();
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
private:
  using TermsContainer = std::vector<TermToAssemble>;

public:
  using const_terms_iterator = TermsContainer::const_iterator;

  const_terms_iterator begin() const { return terms.begin(); }
  const_terms_iterator end() const { return terms.end(); }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  TermsContainer terms;
};

} // namespace akantu

#endif /* AKANTU_TERMS_TO_ASSEMBLE_HH_ */
