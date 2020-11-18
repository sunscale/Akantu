/**
 * @file   dumper.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Mon Jun 10 2013
 *
 * @brief  implementation of main dumper
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * IOHelper is free  software: you can redistribute it and/or  modify it under
 * the terms  of the  GNU Lesser  General Public  License as  published by  the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for
 * more details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include <cmath>

#include "dumper.hh"
#include "field_inline_impl.hh"
#include "variable_inline_impl.hh"
/* -------------------------------------------------------------------------- */
#if defined(__INTEL_COMPILER)
/// remark #981: operands are evaluated in unspecified order
#pragma warning(disable : 981)
#endif // defined(__INTEL_COMPILER)

namespace iohelper {

/* -------------------------------------------------------------------------- */
Dumper::Dumper(const std::string & prefix)
    : dump_step(0), dump_step_width(4), time_step(0), current_time(0),
      write_time_desc_file(false), mode(0), world_size(-1), my_rank(-1),
      root_rank(0), my_rank_width(3), time_description_file_name("") {

  // needs to be after the definition of the directory separator
  this->setPrefix(prefix);
}

/* -------------------------------------------------------------------------- */
Dumper::Dumper(const std::string & prefix, const std::string & base_name)
    : dump_step(0), dump_step_width(4), time_step(0), current_time(0),
      write_time_desc_file(false), mode(0), world_size(-1), my_rank(-1),
      root_rank(0), my_rank_width(3), time_description_file_name("") {
  // needs to be after the definition of the directory separator
  this->setPrefix(prefix);

  setBaseName(base_name);
}

/* -------------------------------------------------------------------------- */
Dumper::~Dumper() {
  {
    auto it = per_node_data.begin();
    auto end = per_node_data.end();
    while (it != end) {
      delete (*it).second;
      ++it;
    }
  }
  {
    auto it = per_element_data.begin();
    auto end = per_element_data.end();
    while (it != end) {
      delete (*it).second;
      ++it;
    }
  }
  {
    auto it = global_data.begin();
    auto end = global_data.end();
    while (it != end) {
      delete (*it).second;
      ++it;
    }
  }
}

/* -------------------------------------------------------------------------- */
void Dumper::dump(const std::string & name, UInt count) {
  if (count != UInt(-1)) {
    dump_step = count;
  }

  if (name != "") {
    base_name = name;
  }

  if (time_description_file_name == "") {
    time_description_file_name = base_name;
  }
}

/* -------------------------------------------------------------------------- */
void Dumper::registerDumpOptions(const std::string & key,
                                 const std::string & folder,
                                 const std::string & extension,
                                 DumpFlag dump_flag) {
  DumpOptions & dos = dump_options[key];
  dos.setFolder(folder);
  dos.extension = extension;
  dos.dump_flags = dump_flag;
}

/* -------------------------------------------------------------------------- */
void Dumper::init() {
  if (world_size == -1 || my_rank == -1) {
    //    DUMP("world_size and my_rank variables are not well set: going to
    //    sequential dump");
    world_size = 1;
    my_rank = 0;
  }
}

/* -------------------------------------------------------------------------- */
void Dumper::setPoints(Real * points, int dimension, int nb,
                       const std::string & name) {
  addNodeDataField(std::string("positions"), points, dimension, nb);
  setBaseName(name);
}

/* -------------------------------------------------------------------------- */
void Dumper::setBaseName(const std::string & name) { this->base_name = name; }

/* -------------------------------------------------------------------------- */
void Dumper::setConnectivity(int * connectivity, ElemType elem_type,
                             UInt nb_elem, int mode) {
  addElemDataField(std::string("connectivities"), connectivity, elem_type,
                   nb_node_per_elem[elem_type], nb_elem);
  // ElemType * types = new ElemType[nb_elem];
  // for (UInt i = 0; i < nb_elem; ++i) {
  //   types[i] = elem_type;
  // }
  // addElemDataField(std::string("element_type"), types, 1, nb_elem);
  connectivity_mode = mode;
}

/* -------------------------------------------------------------------------- */
void Dumper::DumpOptions::setFolder(const std::string & fld) {
  this->folder = Dumper::checkDirectoryName(fld);
  ;
}

/* -------------------------------------------------------------------------- */
const std::string & Dumper::DumpOptions::getFolder() const {
  return this->folder;
}

/* -------------------------------------------------------------------------- */
std::string Dumper::checkDirectoryName(std::string fname) {
  if (fname.size() > 0 &&
      fname[fname.size() - 1] != IOHELPER_DIRECTORY_SEPARATOR) {
    fname += IOHELPER_DIRECTORY_SEPARATOR;
  }
  return fname;
}

/* -------------------------------------------------------------------------- */
Dumper::DumpOptions & Dumper::getDumpOptions(const std::string & key) {
  auto it = this->dump_options.find(key);
  if (it == this->dump_options.end())
    IOHELPER_THROW("No dump options registered under the name " << key,
                   _et_options_error);

  return it->second;
}

/* -------------------------------------------------------------------------- */
std::string Dumper::getRelativeFilePath(const std::string & name,
                                        const std::string & key, UInt proc) {
  return this->getRelativeFolderPath(key) + this->getFileName(name, key, proc);
}

/* -------------------------------------------------------------------------- */
std::string Dumper::getRelativeFilePath(const std::string & name,
                                        const std::string & key) {
  return this->getRelativeFolderPath(key) + this->getFileName(name, key);
}

/* -------------------------------------------------------------------------- */
std::string Dumper::getAbsoluteFilePath(const std::string & name,
                                        const std::string & key, UInt proc) {
  return prefix + this->getRelativeFilePath(name, key, proc);
}

/* -------------------------------------------------------------------------- */
std::string Dumper::getAbsoluteFilePath(const std::string & name,
                                        const std::string & key) {
  return prefix + this->getRelativeFilePath(name, key);
}

/* -------------------------------------------------------------------------- */
std::string Dumper::getFileName(const std::string & name,
                                const std::string & key) {
  return this->getFileName(name, key, this->my_rank);
}

/* -------------------------------------------------------------------------- */
std::string Dumper::getFileName(const std::string & name,
                                const std::string & key, UInt proc) {
  std::stringstream sstr;
  sstr << name;

  const DumpOptions & dos = this->getDumpOptions(key);

  if (dos.dump_flags & _df_counter) {
    sstr << "_";
    sstr.width(this->dump_step_width);
    sstr.fill('0');
    sstr << dump_step;
  }

  if (world_size > 1 && (dos.dump_flags & _df_proc_id)) {
    sstr << ".proc";
    sstr.width(this->my_rank_width);
    sstr.fill('0');
    sstr << proc;
  }

  sstr << dos.extension;

  return sstr.str();
}

/* -------------------------------------------------------------------------- */
std::string Dumper::getRelativeFolderPath(const std::string & key) {
  return getDumpOptions(key).getFolder();
}

/* -------------------------------------------------------------------------- */
std::string Dumper::getAbsoluteFolderPath(const std::string & key) {
  return prefix + this->getRelativeFolderPath(key);
}

/* -------------------------------------------------------------------------- */
std::string Dumper::getBaseName() { return this->base_name; }

/* -------------------------------------------------------------------------- */
void Dumper::setParallelContext(int me, int wld_size, int root) {
  this->my_rank = me;
  this->world_size = wld_size;
  this->my_rank_width = std::ceil(std::log10(this->world_size));
  this->root_rank = root;
}

/* -------------------------------------------------------------------------- */
void Dumper::printNodeDataFields() {
  auto it = per_node_data.begin();
  auto end = per_node_data.end();
  int count = 0;
  while (it != end) {
    std::cout << "Field " << ++count << " : " << it->second->getName()
              << std::endl;
    ++it;
  }
}

/* -------------------------------------------------------------------------- */
void Dumper::activateTimeDescFiles(Real delta_t, Real initial_time) {
  this->time_step = delta_t;
  this->current_time = initial_time;
  this->write_time_desc_file = true;
}

/* -------------------------------------------------------------------------- */

}
