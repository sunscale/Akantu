/**
 * @file   aka_common.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Mon Feb 05 2018
 *
 * @brief  Initialization of global variables
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
#include "aka_random_generator.hh"
#include "aka_static_memory.hh"
#include "communicator.hh"

#include "cppargparse.hh"
#include "parser.hh"

#include "communication_tag.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <ctime>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void initialize(int & argc, char **& argv) {
  AKANTU_DEBUG_IN();

  initialize("", argc, argv);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void initialize(const std::string & input_file, int & argc, char **& argv) {
  AKANTU_DEBUG_IN();
  StaticMemory::getStaticMemory();
  Communicator & comm = Communicator::getStaticCommunicator(argc, argv);

  Tag::setMaxTag(comm.getMaxTag());

  debug::debugger.setParallelContext(comm.whoAmI(), comm.getNbProc());
  debug::setDebugLevel(dblError);

  static_argparser.setParallelContext(comm.whoAmI(), comm.getNbProc());
  static_argparser.setExternalExitFunction(debug::exit);
  static_argparser.addArgument("--aka_input_file", "Akantu's input file", 1,
                               cppargparse::_string, std::string());
  static_argparser.addArgument(
      "--aka_debug_level",
      std::string("Akantu's overall debug level") +
          std::string(" (0: error, 1: exceptions, 4: warnings, 5: info, ..., "
                      "100: dump") +
          std::string(" more info on levels can be foind in aka_error.hh)"),
      1, cppargparse::_integer, (long int)(dblWarning));

  static_argparser.addArgument(
      "--aka_print_backtrace",
      "Should Akantu print a backtrace in case of error", 0,
      cppargparse::_boolean, false, true);

  static_argparser.addArgument("--aka_seed", "The seed to use on prank 0", 1,
                               cppargparse::_integer);

  static_argparser.parse(argc, argv, cppargparse::_remove_parsed);

  std::string infile = static_argparser["aka_input_file"];
  if (infile == "")
    infile = input_file;
  debug::debugger.printBacktrace(static_argparser["aka_print_backtrace"]);

  if ("" != infile) {
    readInputFile(infile);
  }

  long int seed;
  if (static_argparser.has("aka_seed")) {
    seed = static_argparser["aka_seed"];
  } else {
    seed =
        static_parser.getParameter("seed", time(nullptr), _ppsc_current_scope);
  }

  seed *= (comm.whoAmI() + 1);
  RandomGenerator<UInt>::seed(seed);

  long int dbl_level = static_argparser["aka_debug_level"];
  debug::setDebugLevel(DebugLevel(dbl_level));

  AKANTU_DEBUG_INFO("Random seed set to " << seed);

  std::atexit(finalize);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void finalize() {
  AKANTU_DEBUG_IN();

  // if (StaticCommunicator::isInstantiated()) {
  //   StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  //   delete &comm;
  // }

  if (StaticMemory::isInstantiated()) {
    delete &(StaticMemory::getStaticMemory());
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void readInputFile(const std::string & input_file) {
  static_parser.parse(input_file);
}

/* -------------------------------------------------------------------------- */
cppargparse::ArgumentParser & getStaticArgumentParser() {
  return static_argparser;
}

/* -------------------------------------------------------------------------- */
Parser & getStaticParser() { return static_parser; }

/* -------------------------------------------------------------------------- */
const ParserSection & getUserParser() {
  return *(static_parser.getSubSections(ParserType::_user).first);
}

std::unique_ptr<Communicator> Communicator::static_communicator;

std::ostream & operator<<(std::ostream & stream, NodeFlag flag) {
  using under = std::underlying_type_t<NodeFlag>;
  int digits = std::log(std::numeric_limits<under>::max() + 1) / std::log(16);
  std::ios_base::fmtflags ff;
  ff = stream.flags();
  auto value = static_cast<std::common_type_t<under, unsigned int>>(flag);
  stream << "0x" << std::hex << std::setw(digits) << std::setfill('0') << value;
  stream.flags(ff);
  return stream;
}

} // namespace akantu
