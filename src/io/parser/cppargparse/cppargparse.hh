/**
 * @file   cppargparse.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Apr 03 2014
 * @date last modification: Sun Dec 03 2017
 *
 * @brief  Get the commandline options and store them as short, long and others
 *
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

/* -------------------------------------------------------------------------- */
#include <iostream>
#include <map>
#include <string>
#include <vector>

#ifndef __CPPARGPARSE_HH__
#define __CPPARGPARSE_HH__

/* -------------------------------------------------------------------------- */
namespace cppargparse {

/// define the types of the arguments
enum ArgumentType { _string, _integer, _float, _boolean };

/// Defines how many arguments to expect
enum ArgumentNargs { _one_if_possible = -1, _at_least_one = -2, _any = -3 };

/// Flags for the parse function of ArgumentParser
enum ParseFlags {
  _no_flags = 0x0,           ///< Default behavior
  _stop_on_not_parsed = 0x1, ///< Stop on unknown arguments
  _remove_parsed = 0x2       ///< Remove parsed arguments from argc argv
};

/// Helps to combine parse flags
inline ParseFlags operator|(const ParseFlags & a, const ParseFlags & b) {
  auto tmp = ParseFlags(int(a) | int(b));
  return tmp;
}

/* -------------------------------------------------------------------------- */

/**
 * ArgumentParser is a class that mimics the Python argparse module
 */
class ArgumentParser {
public:
  /// public definition of an argument
  class Argument {
  public:
    Argument() : name(std::string()) {}
    virtual ~Argument() = default;
    virtual void printself(std::ostream & stream) const = 0;
    template <class T> operator T() const;
    std::string name;
  };

  /// constructor
  ArgumentParser();

  /// destroy everything
  ~ArgumentParser();

  /// add an argument with a description
  void addArgument(const std::string & name_or_flag, const std::string & help,
                   int nargs = 1, ArgumentType type = _string);

  /// add an argument with an help and a default value
  template <class T>
  void addArgument(const std::string & name_or_flag, const std::string & help,
                   int nargs, ArgumentType type, T def);

  /// add an argument with an help and a default + const value
  template <class T>
  void addArgument(const std::string & name_or_flag, const std::string & help,
                   int nargs, ArgumentType type, T def, T cons);

  /// parse argc, argv
  void parse(int & argc, char **& argv, int flags = _stop_on_not_parsed,
             bool parse_help = true);

  /// get the last argc parsed
  int & getArgC() { return *(this->argc); }

  /// get the last argv parsed
  char **& getArgV() { return *(this->argv); }

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

  /// is the argument present
  bool has(const std::string &) const;

  /// set the parallel context to avoid multiple help messages in
  /// multiproc/thread cases
  void setParallelContext(int prank, int psize);

public:
  /// Internal class describing the arguments
  struct _Argument;
  /// Stores that value of an argument
  template <class T> class ArgumentStorage;

private:
  /// Internal function to be used by the public addArgument
  _Argument & _addArgument(const std::string & name_or_flag,
                           const std::string & description, int nargs,
                           ArgumentType type);

  void _exit(const std::string & msg = "", int status = 0);
  bool checkType(ArgumentType type, const std::string & value) const;

  /// function to help to print help
  void print_usage_nargs(std::ostream & stream,
                         const _Argument & argument) const;
  /// function to help to print help
  void print_help_argument(std::ostream & stream,
                           const _Argument & argument) const;

private:
  /// public arguments storage
  using Arguments = std::map<std::string, Argument *>;
  /// internal arguments storage
  using _Arguments = std::map<std::string, _Argument *>;
  /// association key argument
  using ArgumentKeyMap = std::map<std::string, _Argument *>;
  /// position arguments
  using PositionalArgument = std::vector<_Argument *>;

  /// internal storage of arguments declared by the user
  _Arguments arguments;
  /// list of arguments successfully parsed
  Arguments success_parsed;
  /// keys associated to arguments
  ArgumentKeyMap key_args;
  /// positional arguments
  PositionalArgument pos_args;

  /// program name
  std::string program_name;

  /// exit function to use
  void (*external_exit)(int){nullptr};

  /// Parallel context, rank and size of communicator
  int prank{0}, psize{1};

  /// The last argc parsed (those are the modified version after parse)
  int * argc;

  /// The last argv parsed (those are the modified version after parse)
  char *** argv;
};

inline std::ostream & operator<<(std::ostream & stream,
                                 const ArgumentParser & argparse) {
  argparse.printself(stream);
  return stream;
}

} // namespace cppargparse

#endif /* __CPPARGPARSE_HH__ */

#include "cppargparse_tmpl.hh"
