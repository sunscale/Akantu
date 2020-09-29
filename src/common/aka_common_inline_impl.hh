/**
 * @file   aka_common_inline_impl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  inline implementations of common akantu type descriptions
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
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <iostream>
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
/// standard output stream operator for GhostType
inline std::ostream & operator<<(std::ostream & stream, GhostType type) {
  switch (type) {
  case _not_ghost:
    stream << "not_ghost";
    break;
  case _ghost:
    stream << "ghost";
    break;
  case _casper:
    stream << "Casper the friendly ghost";
    break;
  }
  return stream;
}

/* -------------------------------------------------------------------------- */
inline std::string to_lower(const std::string & str) {
  std::string lstr = str;
  std::transform(lstr.begin(), lstr.end(), lstr.begin(), (int (*)(int))tolower);
  return lstr;
}

namespace {
  template <typename pred>
  inline std::string trim_p(const std::string & to_trim, pred && p) {
    std::string trimed = to_trim;
    auto && not_ = [&](auto && a) { return not p(a); };

    // left trim
    trimed.erase(trimed.begin(),
                 std::find_if(trimed.begin(), trimed.end(), not_));
    // right trim
    trimed.erase(std::find_if(trimed.rbegin(), trimed.rend(), not_).base(),
                 trimed.end());
    return trimed;
  }

} // namespace

/* -------------------------------------------------------------------------- */
inline std::string trim(const std::string & to_trim) {
  return trim_p(to_trim, [&](auto && a) { return std::isspace(a); });
}

inline std::string trim(const std::string & to_trim, char c) {
  return trim_p(to_trim, [&c](auto && a) { return (a == c); });
}

/* -------------------------------------------------------------------------- */
template <typename T> std::string printMemorySize(UInt size) {
  Real real_size = size * sizeof(T);

  UInt mult = 0;
  if (real_size != 0) {
    mult = (std::log(real_size) / std::log(2)) / 10;
  }

  std::stringstream sstr;

  real_size /= Real(1 << (10 * mult));
  sstr << std::setprecision(2) << std::fixed << real_size;

  std::string size_prefix;
  switch (mult) {
  case 0:
    sstr << "";
    break;
  case 1:
    sstr << "Ki";
    break;
  case 2:
    sstr << "Mi";
    break;
  case 3:
    sstr << "Gi";
    break; // I started on this type of machines
           // (32bit computers) (Nicolas)
  case 4:
    sstr << "Ti";
    break;
  case 5:
    sstr << "Pi";
    break;
  case 6:
    sstr << "Ei";
    break; // theoritical limit of RAM of the current
           // computers in 2014 (64bit computers) (Nicolas)
  case 7:
    sstr << "Zi";
    break;
  case 8:
    sstr << "Yi";
    break;
  default:
    AKANTU_ERROR(
        "The programmer in 2014 didn't thought so far (even wikipedia does not "
        "go further)."
        << " You have at least 1024 times more than a yobibit of RAM!!!"
        << " Just add the prefix corresponding in this switch case.");
  }

  sstr << "Byte";

  return sstr.str();
}

} // namespace akantu
