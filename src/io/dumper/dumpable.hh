/**
 * @file   dumpable.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Oct 26 2012
 * @date last modification: Fri Sep 05 2014
 *
 * @brief  Interface for object who wants to dump themselves
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
#include "element_type_map.hh"
#ifdef AKANTU_USE_IOHELPER
#include "dumper_iohelper.hh"
#endif //AKANTU_USE_IOHELPER

#ifndef __AKANTU_DUMPABLE_HH__
#define __AKANTU_DUMPABLE_HH__
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER

#include "dumper_iohelper.hh"
#include <set>

__BEGIN_AKANTU__

class Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Dumpable();
  virtual ~Dumpable();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// create a new dumper (of templated type T) and register it under dumper_name. file_name is used for construction of T. is default states if this dumper is the default dumper. 
  template<class T>
  inline void registerDumper(const std::string & dumper_name,
			     const std::string & file_name = "",
			     const bool is_default = false);

  /// register an externally created dumper
  void registerExternalDumper(DumperIOHelper & dumper,
                              const std::string & dumper_name,
                              const bool is_default = false);

  /// register a mesh to the default dumper
  void addDumpMesh(const Mesh & mesh, UInt spatial_dimension = _all_dimensions,
                   const GhostType & ghost_type = _not_ghost,
                   const ElementKind & element_kind = _ek_not_defined);

  /// register a mesh to the default identified by its name
  void addDumpMeshToDumper(const std::string & dumper_name,
                           const Mesh & mesh, UInt spatial_dimension = _all_dimensions,
                           const GhostType & ghost_type = _not_ghost,
                           const ElementKind & element_kind = _ek_not_defined);

  /// register a filtered mesh as the default dumper
  void addDumpFilteredMesh(const Mesh & mesh,
                           const ElementTypeMapArray<UInt> & elements_filter,
                           const Array<UInt> & nodes_filter,
                           UInt spatial_dimension = _all_dimensions,
                           const GhostType & ghost_type = _not_ghost,
                           const ElementKind & element_kind = _ek_not_defined);

  /// register a filtered mesh and provides a name
  void addDumpFilteredMeshToDumper(const std::string & dumper_name,
                                   const Mesh & mesh,
                                   const ElementTypeMapArray<UInt> & elements_filter,
                                   const Array<UInt> & nodes_filter,
                                   UInt spatial_dimension = _all_dimensions,
                                   const GhostType & ghost_type = _not_ghost,
                                   const ElementKind & element_kind = _ek_not_defined);

  /// to implement
  virtual void addDumpField(const std::string & field_id);
  /// to implement
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
                                    const std::string & field_id);
  
  /// add a field 
  virtual void addDumpFieldExternal(const std::string & field_id,
                                    dumper::Field * field);
  virtual void addDumpFieldExternalToDumper(const std::string & dumper_name,
                                            const std::string & field_id,
                                            dumper::Field * field);

  template<typename T>
  inline void addDumpFieldExternal(const std::string & field_id,
                            const Array<T> & field);
  template<typename T>
  inline void addDumpFieldExternalToDumper(const std::string & dumper_name,
                                    const std::string & field_id,
                                    const Array<T> & field);
  template<typename T>
  inline void addDumpFieldExternal(const std::string & field_id,
                            const ElementTypeMapArray<T> & field,
                            UInt spatial_dimension = _all_dimensions,
                            const GhostType & ghost_type = _not_ghost,
                            const ElementKind & element_kind = _ek_not_defined);
  template<typename T>
  inline void addDumpFieldExternalToDumper(const std::string & dumper_name,
                                    const std::string & field_id,
                                    const ElementTypeMapArray<T> & field,
                                    UInt spatial_dimension = _all_dimensions,
                                    const GhostType & ghost_type = _not_ghost,
                                    const ElementKind & element_kind = _ek_not_defined);

  void removeDumpField(const std::string & field_id);
  void removeDumpFieldFromDumper(const std::string & dumper_name,
                                 const std::string & field_id);

  virtual void addDumpFieldVector(const std::string & field_id);
  virtual void addDumpFieldVectorToDumper(const std::string & dumper_name,
                                          const std::string & field_id);

  virtual void addDumpFieldTensor(const std::string & field_id);
  virtual void addDumpFieldTensorToDumper(const std::string & dumper_name,
                                          const std::string & field_id);

  void setDirectory(const std::string & directory);
  void setDirectoryToDumper(const std::string & dumper_name,
                            const std::string & directory);

  void setBaseName(const std::string & basename);

  void setBaseNameToDumper(const std::string & dumper_name,
                           const std::string & basename);
  void setTimeStepToDumper(Real time_step);
  void setTimeStepToDumper(const std::string & dumper_name,
                           Real time_step);

  virtual void dump();
  virtual void dump(UInt step);
  virtual void dump(Real time, UInt step);
  virtual void dump(const std::string & dumper_name);
  virtual void dump(const std::string & dumper_name, UInt step);
  virtual void dump(const std::string & dumper_name, Real time, UInt step);


public:
  void internalAddDumpFieldToDumper(const std::string & dumper_name,
                                    const std::string & field_id,
                                    dumper::Field * field);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  DumperIOHelper & getDumper();
  DumperIOHelper & getDumper(const std::string & dumper_name);

  template<class T> T & getDumper(const std::string & dumper_name);
  std::string getDefaultDumperName() const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  typedef std::map<std::string, DumperIOHelper *> DumperMap;
  typedef std::set<std::string> DumperSet;

  DumperMap dumpers;
  std::string default_dumper;

  DumperSet external_dumpers;
};

__END_AKANTU__

#else

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused"

namespace dumper {
  class Field;
}

class DumperIOHelper;
class Mesh;

/* -------------------------------------------------------------------------- */
class Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Dumpable() {};
  virtual ~Dumpable() { };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  template<class T>
  inline void registerDumper(const std::string & dumper_name,
                      const std::string & file_name = "",
                      const bool is_default = false) { }

  void registerExternalDumper(DumperIOHelper * dumper,
                              const std::string & dumper_name,
                              const bool is_default = false) { }

  void addDumpMesh(const Mesh & mesh,
                   UInt spatial_dimension = _all_dimensions,
                   const GhostType & ghost_type = _not_ghost,
                   const ElementKind & element_kind = _ek_not_defined) {
  }

  void addDumpMeshToDumper(const std::string & dumper_name,
                           const Mesh & mesh,
                           UInt spatial_dimension = _all_dimensions,
                           const GhostType & ghost_type = _not_ghost,
                           const ElementKind & element_kind = _ek_not_defined) {
  }

  void addDumpFilteredMesh(const Mesh & mesh,
                           const ElementTypeMapArray<UInt> & elements_filter,
                           const Array<UInt> & nodes_filter,
                           UInt spatial_dimension = _all_dimensions,
                           const GhostType & ghost_type = _not_ghost,
                           const ElementKind & element_kind = _ek_not_defined) {
  }

  void addDumpFilteredMeshToDumper(const std::string & dumper_name,
                                   const Mesh & mesh,
                                   const ElementTypeMapArray<UInt> & elements_filter,
                                   const Array<UInt> & nodes_filter,
                                   UInt spatial_dimension = _all_dimensions,
                                   const GhostType & ghost_type = _not_ghost,
                                   const ElementKind & element_kind = _ek_not_defined) {
  }

  virtual void addDumpField(const std::string & field_id){
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
                                    const std::string & field_id){
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  virtual void addDumpFieldExternal(const std::string & field_id,
                                    dumper::Field * field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  virtual void addDumpFieldExternalToDumper(const std::string & dumper_name,
                                    const std::string & field_id,
                                    dumper::Field * field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  template<typename T>
  void addDumpFieldExternal(const std::string & field_id,
                            const Array<T> & field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  template<typename T>
  void addDumpFieldExternalToDumper(const std::string & dumper_name,
                                    const std::string & field_id,
                                    const Array<T> & field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  template<typename T>
  void addDumpFieldExternal(const std::string & field_id,
                            const ElementTypeMapArray<T> & field,
                            UInt spatial_dimension = _all_dimensions,
                            const GhostType & ghost_type = _not_ghost,
                            const ElementKind & element_kind = _ek_not_defined) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  template<typename T>
  void addDumpFieldExternalToDumper(const std::string & dumper_name,
                                    const std::string & field_id,
                                    const ElementTypeMapArray<T> & field,
                                    UInt spatial_dimension = _all_dimensions,
                                    const GhostType & ghost_type = _not_ghost,
                                    const ElementKind & element_kind = _ek_not_defined) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void removeDumpField(const std::string & field_id) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  void removeDumpFieldFromDumper(const std::string & dumper_name,
                               const std::string & field_id) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void setDirectory(const std::string & directory) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  void setDirectoryToDumper(const std::string & dumper_name,
                            const std::string & directory) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void setBaseName(const std::string & basename) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  void setBaseNameToDumper(const std::string & dumper_name,
                           const std::string & basename) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void dump() {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  void dump(const std::string & dumper_name) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void dump(UInt step) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  void dump(const std::string & dumper_name,
            UInt step) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void dump(Real current_time,
            UInt step) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  void dump(const std::string & dumper_name,
            Real current_time,
            UInt step) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

protected:

  void internalAddDumpFieldToDumper(const std::string & dumper_name,
                                    const std::string & field_id,
                                    dumper::Field * field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

protected:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  DumperIOHelper & getDumper() {
    AKANTU_DEBUG_ERROR("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  DumperIOHelper & getDumper(const std::string & dumper_name){
    AKANTU_DEBUG_ERROR("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  template<class T>
  T & getDumper(const std::string & dumper_name) {
    AKANTU_DEBUG_ERROR("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  std::string getDefaultDumperName() {
    AKANTU_DEBUG_ERROR("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

#pragma GCC diagnostic pop

__END_AKANTU__

#endif //AKANTU_USE_IOHELPER

//#include "dumpable_inline_impl.hh"

#endif /* __AKANTU_DUMPABLE_HH__ */
