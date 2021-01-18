/**
 * @file   dumper.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Wed Nov 13 2013
 *
 * @brief  dumper interface
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
#ifndef IOHELPER_DUMPER_H_
#define IOHELPER_DUMPER_H_
/* -------------------------------------------------------------------------- */
#include "container_array.hh"
#include "field.hh"
#include "field_interface.hh"
#include "iohelper_common.hh"
#include "variable.hh"
#include "variable_interface.hh"
#include <map>
#include <string>
/* -------------------------------------------------------------------------- */
#include "visitor.hh"
/* -------------------------------------------------------------------------- */

#if !defined(_WIN32)
#define IOHELPER_DIRECTORY_SEPARATOR '/'
#define iohelper_mkdir(path, mode) mkdir(path, mode)
#else
#define IOHELPER_DIRECTORY_SEPARATOR '\\'
#define iohelper_mkdir(path, mode) mkdir(path)
#include <io.h>
#endif

namespace iohelper {

/** Class Dumper
 * Interface of a dumper
 */

class Dumper {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  Dumper(const std::string & prefix, const std::string & base_name);
  Dumper(const std::string & prefix);
  virtual ~Dumper();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  //! dump to file
  virtual void dump(const std::string & name = std::string(),
                    UInt count = UInt(-1));
  //! initialisation of the dumper
  void init();

  //! TODO set comment
  void printNodeDataFields();

  //! dump of field information
  virtual void dumpDescription(__attribute__((unused))
                               const char descr_sep = ' '){};

public:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */

  //! give vector with coordinates
  void setPoints(Real * points, int dimension, int nb,
                 const std::string & name);

  //! set number of filtered elements
  void setNumberFilteredElements(int nb_filtered);

  //! give vector to connectivity
  virtual void setConnectivity(int * connectivity, ElemType elem_type,
                               UInt nb_elem, int mode);

  //! give vector to per node data
  template <typename T>
  void addNodeDataField(const std::string & name, T * data, UInt dimension,
                        UInt size, UInt stride = 1);

  //! give a generic container as per node data
  template <typename Cont>
  void addNodeDataField(const std::string & name, Cont & data);

  //! give vector to per element data
  template <typename T>
  void addElemDataField(const std::string & name, T * data, ElemType elem_type,
                        UInt dimension, UInt size, UInt stride = 1);

  //! give a generic container as per elem data
  template <typename Cont>
  void addElemDataField(const std::string & name, Cont & data);

  //! give a generic container as per node data
  template <typename VarType>
  void addVariable(const std::string & name, VarType & data);

  //! set mode
  virtual void setMode(int mode) { this->mode = mode; }

  //! set rank and world size params for parallel treatment
  void setParallelContext(int me, int wld_size, int root = 0);

  //! set current value for the dump step
  void setDumpStep(int s) { dump_step = s; };
  Int getDumpStep() const { return dump_step; };
  Int getDumpStepWidth() const { return dump_step_width; };

  //! increment the dumpstep
  void incDumpStep(UInt s = 1) { dump_step += s; };

  void setPrefix(const std::string & prefix) {
    this->prefix = Dumper::checkDirectoryName(prefix);
  };

  virtual void activateTimeDescFiles(Real delta_t, Real initial_time = 0.);

  void setTimeStep(Real delta_t) { this->time_step = delta_t; }
  void setCurrentTime(Real time) { this->current_time = time; }

  void setTimeDescriptionFileName(const std::string & name) {
    this->time_description_file_name = name;
  }

protected:
  static std::string checkDirectoryName(std::string fname);

  std::string getFolder(const std::string & key);

  //! get file name with relative path
  std::string getRelativeFilePath(const std::string & name,
                                  const std::string & key, UInt proc);

  std::string getRelativeFilePath(const std::string & name,
                                  const std::string & key);

  std::string getAbsoluteFilePath(const std::string & name,
                                  const std::string & key, UInt proc);

  std::string getAbsoluteFilePath(const std::string & name,
                                  const std::string & key);

  std::string getFileName(const std::string & name, const std::string & key,
                          UInt proc);

  std::string getFileName(const std::string & name, const std::string & key);

  std::string getRelativeFolderPath(const std::string & key);
  std::string getAbsoluteFolderPath(const std::string & key);

  //! get base_name
  std::string getBaseName();
  void setBaseName(const std::string & name);

public:
  enum DumpFlag { _df_no_flag = 0x0, _df_counter = 0x1, _df_proc_id = 0x2 };

protected:
  struct DumpOptions {
  private:
    std::string folder;

  public:
    std::string extension;
    DumpFlag dump_flags;

  public:
    void setFolder(const std::string & fld);
    const std::string & getFolder() const;
  };

  using DumpOptionsMap = std::map<std::string, DumpOptions>;

protected:
  void registerDumpOptions(const std::string & key, const std::string & folder,
                           const std::string & extension,
                           DumpFlag dump_flag = _df_no_flag);

  DumpOptions & getDumpOptions(const std::string & key);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  std::string base_name;
  std::string prefix;

  DumpOptionsMap dump_options;

  UInt dump_step;
  UInt dump_step_width;

protected:
  //! current delta t
  Real time_step;

  //! current time (dah!)
  Real current_time;
  //! if implemented the dumper will write a time description file
  bool write_time_desc_file;

  //! flag to produce zipped files
  UInt mode;

  using field_map = std::map<std::string, FieldInterface *>;
  using variable_map = std::map<std::string, VariableInterface *>;

  //! vector of additional per node data
  field_map per_node_data;
  //! vector of additional per element data
  field_map per_element_data;
  //! map of global variables
  variable_map global_data;
  //! for parallel runs is the total number of processors
  Int world_size;
  //! for parallel runs is rank of the process
  Int my_rank;
  //! fortran or C style connectivity indexing
  Int connectivity_mode;
  //! for parallel runs is rank of the process that should write the data to
  //! file
  Int root_rank;
  //! number of zeros to put to my_rank when creating files
  Int my_rank_width;

  //! file name for the time description files (if not set use the first base
  //! name)
  std::string time_description_file_name;
};

inline Dumper::DumpFlag operator|(const Dumper::DumpFlag & a,
                                  const Dumper::DumpFlag & b) {
  auto tmp = Dumper::DumpFlag(UInt(a) | UInt(b));
  return tmp;
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Dumper::addNodeDataField(const std::string & name, T * data,
                              UInt dimension, UInt size, UInt stride) {

  auto * cont = new ContainerArray<T>(data, dimension, size, stride);
  addNodeDataField(name, *cont);
}

/* -------------------------------------------------------------------------- */
template <typename Cont>
void Dumper::addNodeDataField(const std::string & name, Cont & data) {
  auto * test = new Field<Cont>(data, name);
  per_node_data[name] = test;
}

/* -------------------------------------------------------------------------- */
template <typename VarType>
void Dumper::addVariable(const std::string & name, VarType & data) {
  auto * vari = new Variable<VarType>(data, name);
  global_data[name] = vari;
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Dumper::addElemDataField(const std::string & name, T * data,
                              ElemType elem_type, UInt dimension, UInt size,
                              UInt stride) {

  auto * cont = new ContainerArray<T>(data, dimension, size, stride);
  cont->setElemType(elem_type);
  addElemDataField(name, *cont);
}

/* -------------------------------------------------------------------------- */
template <typename Cont>
void Dumper::addElemDataField(const std::string & name, Cont & data) {
  auto * test = new Field<Cont>(data, name);
  per_element_data[name] = test;
}

/* -------------------------------------------------------------------------- */

}

#endif /* IOHELPER_DUMPER_H_ */
