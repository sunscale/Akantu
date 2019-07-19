/**
 * @file   aka_array_printer.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  mer jun 19 2019
 *
 * @brief A Documented file.
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
#include "aka_array.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_ARRAY_PRINTER_HH__
#define __AKANTU_AKA_ARRAY_PRINTER_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class container, bool no_explicit = true> class ArrayPrinter {
public:
  ArrayPrinter(const container & cont) : cont(cont) {}

  void printself(std::ostream & stream, int indent = 0) const {
    std::string space(indent, AKANTU_INDENT);

    stream << space << "{";
    for (UInt i = 0; i < this->cont.size(); ++i) {
      stream << this->cont[i];
      if (i != this->cont.size() - 1)
        stream << ", ";
    }
    stream << "}" << std::endl;
  }

private:
  const container & cont;
};

/* -------------------------------------------------------------------------- */
template <class T> class ArrayPrinter<Array<T>> {
public:
  ArrayPrinter(const Array<T> & cont) : cont(cont) {}

  void printself(std::ostream & stream, int indent = 0) const {
    std::string space(indent, AKANTU_INDENT);

    stream << space << "{";
    for (UInt i = 0; i < this->cont.size(); ++i) {
      stream << "{";
      for (UInt j = 0; j < this->cont.getNbComponent(); ++j) {
        stream << this->cont(i, j);
        if (j != this->cont.getNbComponent() - 1)
          stream << ", ";
      }
      stream << "}";
      if (i != this->cont.size() - 1)
        stream << ", ";
    }
    stream << "}" << std::endl;
  }

private:
  const Array<T> & cont;
};

template <class container>
decltype(auto) make_printer(const container & array) {
  return ArrayPrinter<container>(array);
}

/* -------------------------------------------------------------------------- */
template <class T>
inline std::ostream & operator<<(std::ostream & stream,
                                 const ArrayPrinter<T> & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* __AKANTU_AKA_ARRAY_PRINTER_HH__ */
