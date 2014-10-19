/**
 * @file   aka_common.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  Initialization of global variables
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "aka_static_memory.hh"
#include "static_communicator.hh"
#include "static_solver.hh"
#include "aka_random_generator.hh"

#include "parser.hh"
#include "cppargparse.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void initialize(int & argc, char ** & argv) {
  AKANTU_DEBUG_IN();

  initialize("", argc, argv);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void initialize(const std::string & input_file, int & argc, char ** & argv) {
  AKANTU_DEBUG_IN();
  StaticMemory::getStaticMemory();
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator(argc, argv);
  debug::debugger.setParallelContext(comm.whoAmI(), comm.getNbProc());
  debug::initSignalHandler();

  static_argparser.setParallelContext(comm.whoAmI(), comm.getNbProc());
  static_argparser.setExternalExitFunction(debug::exit);
  static_argparser.addArgument("--aka_input_file", "Akantu's input file",
			       1, cppargparse::_string, std::string());
  static_argparser.addArgument("--aka_debug_level", std::string("Akantu's overall debug level") +
			       std::string(" (0: error, 1: exceptions, 4: warnings, 5: info, ..., 100: dump,") +
			       std::string(" more info on levels can be foind in aka_error.hh)"),
			       1,
			       cppargparse::_integer, int(dblWarning));

  static_argparser.addArgument("--aka_print_backtrace",
			       "Should Akantu print a backtrace in case of error",
			       0,
			       cppargparse::_boolean, false, true);

  static_argparser.parse(argc, argv, cppargparse::_remove_parsed);

  
  std::string infile = static_argparser["aka_input_file"];
  if(infile == "") infile = input_file;

  debug::setDebugLevel(dblError);

  if ("" != infile) {
    static_parser.parse(infile);
  }

  long int seed;
  try {
    seed = static_parser.getParameter("seed", _ppsc_current_scope);
  } catch (debug::Exception &) {
    seed = time(NULL);
  }

  int dbl_level = static_argparser["aka_debug_level"];
  debug::setDebugLevel(DebugLevel(dbl_level));
  debug::debugger.printBacktrace(static_argparser["aka_print_backtrace"]);

  seed *= (comm.whoAmI() + 1);
  Rand48Generator<Real>::seed(seed);
  RandGenerator<Real>::seed(seed);
  AKANTU_DEBUG_INFO("Random seed set to " << seed);

  /// initialize external solvers 
  StaticSolver::getStaticSolver().initialize(argc, argv);

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
void finalize() {
  AKANTU_DEBUG_IN();

  /// finalize external solvers 
  StaticSolver::getStaticSolver().finalize();

  if(StaticMemory::isInstantiated()) delete &(StaticMemory::getStaticMemory());
  if(StaticCommunicator::isInstantiated()) {
    StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
    comm.barrier();
    delete &comm;
  }


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
cppargparse::ArgumentParser & getStaticArgumentParser() {
  return static_argparser;
}

/* -------------------------------------------------------------------------- */
Parser & getStaticParser() {
  return static_parser;
}

/* -------------------------------------------------------------------------- */
const ParserSection & getUserParser() {
  return *(static_parser.getSubSections(_st_user).first);
}

__END_AKANTU__
