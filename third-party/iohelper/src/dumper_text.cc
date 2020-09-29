/**
 * @file   dumper_text.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue May 14 2013
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  implementation for text dumper
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
#include <sys/stat.h>
#include <iomanip>

#include "dumper_text.hh"
/* -------------------------------------------------------------------------- */
#if defined(__INTEL_COMPILER)
/// remark #981: operands are evaluated in unspecified order
#pragma warning ( disable : 981 )
#endif //defined(__INTEL_COMPILER)

namespace iohelper {

/* -------------------------------------------------------------------------- */
DumperText::DumperText(TextDumpMode md, const std::string & prefix, 
		       bool file_per_time_step) : 
  Dumper(prefix), sim_name("sim_description"), 
  time_name("time_description"), field_name("field_description"), 
  description_sep(' '), dump_mode(md), separator(' '), comment('#'), 
  precision(6), file_per_time_step(file_per_time_step), is_first_dump(true) {

  // fields and variables dump options
  if (file_per_time_step) {
    this->registerDumpOptions("data_fields", "", ".out", _df_counter | _df_proc_id);
  } else {
    this->registerDumpOptions("data_fields", "", ".out", _df_proc_id);
  }

  this->registerDumpOptions("data_variables", "", ".out");

  // description dump options
  this->registerDumpOptions(time_name, "", ".time");
  this->registerDumpOptions(sim_name, "", ".info");
  this->registerDumpOptions(field_name, "", ".fields");

  this->setDumpMode(md);
}


/* -------------------------------------------------------------------------- */
DumperText::~DumperText() {
  
  // only root rank has this files
  //  if (this->my_rank == this->root_rank) {
  auto it = this->file_map.begin();
  auto end = this->file_map.end();

  for (; it != end; ++it) {
    it->second->close();
    delete it->second;
  }

  file_map.clear();
  //}
  
}

/* -------------------------------------------------------------------------- */
void DumperText::setDumpMode(const TextDumpMode & mode) {

  this->dump_mode = mode;
  
  switch(mode) {
  case _tdm_space: {
    this->separator = ' ';
    this->setDataFileExtensions(".out");
    break;
  }
  case _tdm_csv: {
    this->separator = ',';
    this->setDataFileExtensions(".csv");
    break;
  }
  default: {
    IOHELPER_THROW("Unkown dump mode " << mode, _et_options_error);
  }
  }
}

/* -------------------------------------------------------------------------- */
void DumperText::dump(const std::string & name, UInt count) {
  Dumper::dump(name, count);

  //create sub directory to store files
  iohelper_mkdir(std::string(this->getAbsoluteFolderPath("data_fields")).c_str()   , 0755);
  iohelper_mkdir(std::string(this->getAbsoluteFolderPath("data_variables")).c_str(), 0755);
  
  /* dump description file */
  if (this->is_first_dump) {
    this->dumpDescription(this->description_sep);
  }

  /* nodal data */
  auto it = per_node_data.begin();
  auto end = per_node_data.end();
  for (; it != end ; ++it) {
    (*it).second->accept(*this);
  }
  
  /* element data */
  it = per_element_data.begin();
  end = per_element_data.end();
  for (; it != end ; ++it) {
    (*it).second->accept(*this);
  }

  /* global data */
  auto git = global_data.begin();
  auto gend = global_data.end();
  for (; git != gend; ++git) {
    (*git).second->accept(*this);
  }

  if (this->write_time_desc_file) {
    this->dumpTimeDescription(this->description_sep);
    this->current_time += this->time_step;
  }

  this->is_first_dump = false;
  this->incDumpStep();
}

/* -------------------------------------------------------------------------- */
void DumperText::setDataSubDirectory(const std::string & name) {
  this->getDumpOptions("data_fields").setFolder(name);
  this->getDumpOptions("data_variables").setFolder(name);
}

/* -------------------------------------------------------------------------- */
std::string DumperText::getDataSubDirectory() {
  return this->getDumpOptions("data_fields").getFolder();
}

/* -------------------------------------------------------------------------- */
void DumperText::setDataFileExtensions(const std::string & ext) {
  this->getDumpOptions("data_fields").extension    = ext;
  this->getDumpOptions("data_variables").extension = ext;
}

/* -------------------------------------------------------------------------- */
void DumperText::dumpDescription(const char descr_sep) {

  this->dumpFieldDescription(descr_sep);

  File file;
  file.open(this->getAbsoluteFilePath(this->getBaseName(), this->sim_name),
				      std::fstream::out);
  
  file << this->comment 
       << " [0]-version [1]-counter-width [2]-world-size" 
       << " [3]-proc-width [4]-file-per-time-step" << std::endl;
  file << "0-1";
  file << descr_sep << this->getDumpStepWidth();
  file << descr_sep << this->world_size;
  file << descr_sep << this->my_rank_width;
  file << descr_sep << this->file_per_time_step;
  file << std::endl << std::endl;

  file << "base_name" << descr_sep << this->getBaseName() << std::endl;

  // file << this->comment << " field description relative file path" << std::endl;
  file << this->field_name << descr_sep 
       << this->getRelativeFilePath(this->getBaseName(), this->field_name) << std::endl;

  // file << this->comment << " time description relative file path" << std::endl;
  if (this->write_time_desc_file) {
    file << this->time_name << descr_sep 
	 << this->getRelativeFilePath(this->getBaseName(), this->time_name) << std::endl;
  }

  // file << this->comment << " time description relative file path" << std::endl;
  // file << "datafolder" << descr_sep 
  //      << this->getDataSubDirectory() << std::endl;

  file.close();
}

/* -------------------------------------------------------------------------- */
void DumperText::dumpTimeDescription(const char descr_sep) {

  // only root rank dumps variables
  if (this->my_rank != this->root_rank) {
    return;
  }

  auto it = this->file_map.find(this->time_name);
  auto end = this->file_map.end();

  File * file;
  if (it == end) {
    File * new_file = new File;
    new_file->open(this->getAbsoluteFilePath(this->getBaseName(),
					     this->time_name), 
		   std::fstream::out);
    this->file_map[this->time_name] = new_file;
    file = new_file;
    (*file) << std::scientific << std::setprecision(this->precision);
  }
  else {
    file = it->second;
  }

  (*file) << this->getDumpStep();
  (*file) << descr_sep << this->current_time;
  (*file) << std::endl;
}

/* -------------------------------------------------------------------------- */
void DumperText::dumpFieldDescription(const char descr_sep) {

  File file;
  file.open(this->getAbsoluteFilePath(this->getBaseName(), 
				      this->field_name),
	    std::fstream::out);
  
  file << this->comment 
       << " [0]-field_name [1]-nodal_element_global_data [2]-data_type"
       << " [3]-dump_mode [4]-nb_components [5]-rel_file_path" << std::endl;

  /* nodal data */
  auto it = per_node_data.begin();
  auto end = per_node_data.end();
  for (; it != end ; ++it) {
    FieldInterface * fld = (*it).second;
    file << fld->getName();
    file << descr_sep << "N"; // nodal information
    file << descr_sep << fld->getDataType(); 
    file << descr_sep << this->dump_mode;
    file << descr_sep << fld->getDim();
    file << descr_sep << this->getRelativeFolderPath("data_fields");
    file << std::endl;
  }
  
  /* element data */
  it = per_element_data.begin();
  end = per_element_data.end();
  for (; it != end ; ++it) {
    FieldInterface * fld = (*it).second;
    file << fld->getName(); 
    file << descr_sep << "E"; // element information
    file << descr_sep << fld->getDataType();
    file << descr_sep << this->dump_mode;
    file << descr_sep << fld->getDim();
    file << descr_sep << this->getRelativeFolderPath("data_fields");
    file << std::endl;
  }

  /* global data */
  auto git = global_data.begin();
  auto gend = global_data.end();
  for (; git != gend; ++git) {
    VariableInterface * var = (*git).second;
    file << var->getName();
    file << descr_sep << "G"; // global information
    file << descr_sep << var->getDataType();
    file << descr_sep << this->dump_mode;
    file << descr_sep << var->getDim();
    file << descr_sep << this->getRelativeFolderPath("data_variables");
    file << std::endl;
  }
  
  file.close();

}

}

/* -------------------------------------------------------------------------- */



