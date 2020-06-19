/**
 * @file   aka_extern.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  initialisation of all global variables
 * to insure the order of creation
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
#include "aka_array.hh"
#include "aka_common.hh"
#include "aka_math.hh"
#include "aka_named_argument.hh"
#include "aka_random_generator.hh"
#include "communication_tag.hh"
#include "cppargparse.hh"
#include "parser.hh"
#include "solid_mechanics_model.hh"
#if defined(AKANTU_COHESIVE_ELEMENT)
#include "solid_mechanics_model_cohesive.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <limits>
/* -------------------------------------------------------------------------- */
#if defined(AKANTU_DEBUG_TOOLS)
#include "aka_debug_tools.hh"
#endif

namespace akantu {

/* -------------------------------------------------------------------------- */
/* error.hpp variables                                                        */
/* -------------------------------------------------------------------------- */
namespace debug {
  /** \todo write function to get this
   *   values from the environment or a config file
   */
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
} // namespace debug

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
Real Math::tolerance = 1e2 * std::numeric_limits<Real>::epsilon();

/* -------------------------------------------------------------------------- */
const UInt _all_dimensions [[gnu::unused]] = UInt(-1);

/* -------------------------------------------------------------------------- */
const Array<UInt> empty_filter(0, 1, "empty_filter");

/* -------------------------------------------------------------------------- */
template <> long int RandomGenerator<UInt>::_seed = 5489u;
template <> std::default_random_engine RandomGenerator<UInt>::generator(5489u);
/* -------------------------------------------------------------------------- */
int Tag::max_tag = 0;

/* -------------------------------------------------------------------------- */

} // namespace akantu
