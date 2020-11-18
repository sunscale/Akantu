/**
 * @file   dumper_text.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue May 14 2013
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  header for dumper text
 *
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * IOHelper is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef IOHELPER_DUMPER_TEXT_H_
#define IOHELPER_DUMPER_TEXT_H_
/* -------------------------------------------------------------------------- */
#include <map>
#include <string>
#include "dumper.hh"
#include "file_manager.hh"
/* -------------------------------------------------------------------------- */

namespace iohelper {

/** Class DumperText
 * Implementation of a dumper to text file
 */

class DumperText : public Dumper, public Visitor {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  
public:
  
  DumperText(TextDumpMode md = _tdm_space, 
	     const std::string & prefix = "./",
	     bool file_per_time_step = false);
  ~DumperText() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  void dump(const std::string & name, UInt count) override;
  void setEmbeddedValue(const std::string & /*name*/,
                        int /*value*/) {};

  void dumpDescription(char descr_sep = ' ') override;
  virtual void dumpFieldDescription(char descr_sep = ' ');
  virtual void dumpTimeDescription(char descr_sep = ' ');
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */

  void setDataSubDirectory(const std::string & name);
  std::string getDataSubDirectory();
  void setDataFileExtensions(const std::string & ext);
  void setDumpMode(const TextDumpMode & mode);

  void setPrecision(UInt prec) { this->precision = prec; };
 
  //! visitor system
  template <typename T> void visitField(T & visited);
  template <typename T> void visitVariable(T & visited);
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  using FileMap = std::map<std::string, File *>;

  /**
   * another base name so that you will never understand how iohelper
   * does what it does (now useless so it is commented but not removed
   * for historical reasons)
   */
  //  std::string yet_another_base_name;

  /// name of simulation, time and field description options
  std::string sim_name;
  std::string time_name;
  std::string field_name;
  char description_sep;

  /// this is a separator !!
  TextDumpMode dump_mode;
  char separator;
  char comment;
  UInt precision;
  bool file_per_time_step;

  bool is_first_dump;

  FileMap file_map;
};

/* -------------------------------------------------------------------------- */
template <typename T>
void DumperText::visitField(T & visited) {
  File file;
  
  if (this->file_per_time_step || this->is_first_dump) {
    file.open(this->getAbsoluteFilePath(this->getBaseName() + "_" + visited.getName(),
					"data_fields"),
	      std::fstream::out);
  }
  else {
    file.open(this->getAbsoluteFilePath(this->getBaseName() + "_" + visited.getName(),
					"data_fields"),
	      std::fstream::out | std::fstream::app);
  }
  file << std::scientific << std::setprecision(this->precision);

  typename T::iterator it = visited.begin();
  typename T::iterator end = visited.end();
  
  UInt dim = visited.getDim();
  
  for (; it != end; ++it) {
    for (UInt i=0; i<dim; ++i) {
      if (i != 0) {
        file << this->separator;
      }
      file << (*it)[i];
    }
    file << std::endl;
  }

  file << std::endl;  
  file.close();
}

/* -------------------------------------------------------------------------- */
template <typename T>
void DumperText::visitVariable(T & visited) {

  // only root rank dumps variables
  if (this->my_rank != this->root_rank) {
    return;
  }

  const std::string & name = visited.getName();

  auto it = this->file_map.find(name);
  auto end = this->file_map.end();

  File * file;
  if (it == end) {
    auto * new_file = new File;
    new_file->open(this->getAbsoluteFilePath(this->getBaseName() + "_" + name,
					     "data_variables"), 
		   std::fstream::out);
    this->file_map[name] = new_file;
    file = new_file;
    (*file) << std::scientific << std::setprecision(this->precision);
  }
  else {
    file = it->second;
  }

  UInt dim = visited.getDim();
  
  //  File file = it->second;
  for (UInt i=0; i<dim; ++i) {
    if (i != 0) {
      (*file) << this->separator;
    }
    (*file) << (*visited)[i];
  }
  (*file) << std::endl;
  
}

/* -------------------------------------------------------------------------- */
}

#endif /* IOHELPER_DUMPER_TEXT_H_ */
