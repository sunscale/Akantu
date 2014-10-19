/**
 * @file   cppargparse.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Apr 03 2014
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  Get the commandline options and store them as short, long and others
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <map>
#include <string>
#include <iostream>
#include <vector>


#ifndef __CPPARGPARSE_HH__
#define __CPPARGPARSE_HH__

namespace cppargparse {

enum ArgumentType {
  _string,
  _integer,
  _float,
  _boolean
};

enum ArgumentNargs {
  _one_if_possible = -1,
  _at_least_one = -2,
  _any = -3
};

enum ParseFlags {
  _no_flags           = 0x0,
  _stop_on_not_parsed = 0x1,
  _remove_parsed      = 0x2
};

inline ParseFlags operator|(const ParseFlags & a, const ParseFlags & b) {
  ParseFlags tmp = ParseFlags(int(a) | int(b));
  return tmp;
}

class ArgumentParser {
public:
  class Argument {
  public:
    Argument() : name(std::string()) {}
    virtual ~Argument() {}
    virtual void printself(std::ostream & stream) const = 0;
    template<class T> operator T() const;
    std::string name;
  };

  /// constructor
  ArgumentParser();

  /// destroy everything
  ~ArgumentParser();

  /// add an argument with a description
  void addArgument(const std::string & name_or_flag,
		     const std::string & help,
		     int nargs = 1,
		     ArgumentType type = _string);

  /// add an argument with an help and a default value
  template<class T>
  void addArgument(const std::string & name_or_flag,
		     const std::string & help,
		     int nargs,
		     ArgumentType type,
		     T def);

  /// add an argument with an help and a default + const value
  template<class T>
  void addArgument(const std::string & name_or_flag,
		     const std::string & help,
		     int nargs,
		     ArgumentType type,
		     T def,
		     T cons);

  /// parse argc, argv
  void parse(int & argc, char ** & argv,
	     int flags = _stop_on_not_parsed, bool parse_help = true);

  /// print the content in the stream
  void printself(std::ostream & stream) const;

  /// print the help text
  void print_help(std::ostream & stream = std::cout) const;

  /// print the usage text
  void print_usage(std::ostream & stream = std::cout) const;

  /// set an external function to replace the exit function from the stdlib
  void setExternalExitFunction(void (*external_exit)(int)) {
    this->external_exit = external_exit;
  }

  /// accessor for a registered argument that was parsed, throw an exception if
  /// the argument does not exist or was not set (parsed or default value)
  const Argument & operator[](const std::string & name) const;

  bool has(const std::string &) const;
  
  /// set the parallel context to avoid multiple help messages in multiproc/thread cases
  void setParallelContext(int prank, int psize);

public:
  class _Argument;
  template<class T> class ArgumentStorage;

private:
  _Argument & _addArgument(const std::string & name_or_flag,
			     const std::string & description,
			     int nargs,
			     ArgumentType type);

  void _exit(const std::string & msg = "", int status = 0);
  bool checkType(ArgumentType type, const std::string & value) const;

  void print_usage_nargs(std::ostream & stream, const _Argument & argument) const;
  void print_help_argument(std::ostream & stream, const _Argument & argument) const;

private:
  typedef std::map<std::string, _Argument *> _Arguments;
  typedef std::map<std::string, Argument *> Arguments;
  typedef std::map<std::string, _Argument *> ArgumentKeyMap;
  typedef std::vector<_Argument *> PositionalArgument;

  _Arguments arguments;

  ArgumentKeyMap     key_args;
  PositionalArgument pos_args;
  Arguments success_parsed;

  std::string program_name;
  void (*external_exit)(int);

  int prank, psize;
};

}

inline std::ostream & operator<<(std::ostream & stream, const cppargparse::ArgumentParser & argparse) {
  argparse.printself(stream);
  return stream;
}

#endif /* __CPPARGPARSE_HH__ */

#include "cppargparse_tmpl.hh"
