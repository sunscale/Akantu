/**
 * @file   aka_fwd.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Apr 13 2012
 * @date last modification: Wed Oct 25 2017
 *
 * @brief  File containing forward declarations in akantu.
 * This file helps if circular #include would be needed because two classes
 * refer both to each other. This file usually does not need any modification.
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
#ifndef AKANTU_FWD_HH_
#define AKANTU_FWD_HH_

namespace cppargparse {
class ArgumentParser;
}

namespace akantu {
// forward declaration
template <int dim, class model_type> struct ContactData;

template <typename T> class Matrix;
template <typename T> class Vector;
template <typename T> class Tensor3;

template <typename T, bool is_scal = aka::is_scalar<T>::value> class Array;
template <typename T, typename SupportType = ElementType>
class ElementTypeMapArray;

template <class T> class SpatialGrid;

// Model element
template <class ModelPolicy> class ModelElement;

extern const Array<UInt> empty_filter;

class Parser;
class ParserSection;

extern Parser static_parser; // NOLINT

extern cppargparse::ArgumentParser static_argparser; // NOLINT

class Mesh;
class SparseMatrix;
} // namespace akantu

#endif /* AKANTU_FWD_HH_ */
