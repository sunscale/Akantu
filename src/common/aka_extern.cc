/**
 * @file   aka_extern.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Thu Nov 19 2015
 *
 * @brief  initialisation of all global variables
 * to insure the order of creation
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "aka_common.hh"
#include "aka_math.hh"
#include "aka_random_generator.hh"
#include "cppargparse.hh"
#include "parser.hh"
#include "communication_tag.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <limits>
/* -------------------------------------------------------------------------- */
#if defined(AKANTU_DEBUG_TOOLS)
#include "aka_debug_tools.hh"
#endif

namespace akantu {

/** \todo write function to get this
 *   values from the environment or a config file
 */

/* -------------------------------------------------------------------------- */
/* error.hpp variables                                                        */
/* -------------------------------------------------------------------------- */
namespace debug {
/// standard output for debug messages
std::ostream * _akantu_debug_cout = &std::cerr;

/// standard output for normal messages
std::ostream & _akantu_cout = std::cout;

/// parallel context used in debug messages
std::string _parallel_context = "";

Debugger debugger;

#if defined(AKANTU_DEBUG_TOOLS)
DebugElementManager element_manager;
#endif
}

/* -------------------------------------------------------------------------- */
/// list of ghost iterable types
ghost_type_t ghost_types(_casper);

/* -------------------------------------------------------------------------- */
/// Paser for commandline arguments
::cppargparse::ArgumentParser static_argparser;

/// Parser containing the information parsed by the input file given to initFull
Parser static_parser;

bool Parser::permissive_parser = false;
/* -------------------------------------------------------------------------- */
Real Math::tolerance = std::numeric_limits<Real>::epsilon();
/* -------------------------------------------------------------------------- */
const UInt _all_dimensions = UInt(-1);
/* -------------------------------------------------------------------------- */
const Array<UInt> empty_filter(0, 1, "empty_filter");
/* -------------------------------------------------------------------------- */
template <> long int RandGenerator<Real>::_seed = 0;
template <>
long int RandGenerator<bool>::_seed =
    0; // useless just defined due to a template instantiation
template <> long int RandGenerator<UInt>::_seed = 0;
template <> long int RandGenerator<Int>::_seed = 0;
#if not defined(_WIN32)
template <> long int Rand48Generator<Real>::_seed = 0;
#endif
/* -------------------------------------------------------------------------- */
int Tag::max_tag = 0;
/* -------------------------------------------------------------------------- */

}  // akantu
