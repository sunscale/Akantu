/**
 * @file   dumper_padding_helper.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Mon Jun 19 2017
 *
 * @brief  Padding helper for field iterators
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_DUMPER_PADDING_HELPER_HH__
#define __AKANTU_DUMPER_PADDING_HELPER_HH__
/* -------------------------------------------------------------------------- */
#include "dumper_compute.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */

class PadderInterface {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  PadderInterface() {
    padding_m = 0;
    padding_n = 0;
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:
  void setPadding(UInt m, UInt n = 0) {
    padding_m = m;
    padding_n = n;
  }

  virtual UInt getPaddedDim(UInt nb_data) { return nb_data; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

public:
  /// padding informations
  UInt padding_n, padding_m;
};

/* -------------------------------------------------------------------------- */

template <class input_type, class output_type>
class PadderGeneric : public ComputeFunctor<input_type, output_type>,
                      public PadderInterface {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  PadderGeneric() : PadderInterface() {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:
  inline output_type pad(const input_type & in,
                         __attribute__((unused)) UInt nb_data) {
    return in; // trick due to the fact that IOHelper padds the vectors (avoid a
               // copy of data)
  }
};
/* -------------------------------------------------------------------------- */

template <class T>
class PadderGeneric<Vector<T>, Matrix<T>>
    : public ComputeFunctor<Vector<T>, Matrix<T>>, public PadderInterface {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  inline Matrix<T> pad(const Vector<T> & _in, UInt nrows, UInt ncols,
                       UInt nb_data) {
    Matrix<T> in(_in.storage(), nrows, ncols);

    if (padding_m <= nrows && padding_n * nb_data <= ncols)
      return in;
    else {
      Matrix<T> ret(padding_m, padding_n * nb_data);
      UInt nb_cols_per_data = in.cols() / nb_data;
      for (UInt d = 0; d < nb_data; ++d)
        for (UInt i = 0; i < in.rows(); ++i)
          for (UInt j = 0; j < nb_cols_per_data; ++j)
            ret(i, j + d * padding_n) = in(i, j + d * nb_cols_per_data);
      return ret;
    }
  }
};

/* -------------------------------------------------------------------------- */

} // namespace akantu
__END_AKANTU_DUMPER__

#endif /* __AKANTU_DUMPER_PADDING_HELPER_HH__ */
