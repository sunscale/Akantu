/**
 * @file   dumpable_dummy.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Oct 26 2012
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Interface for object who wants to dump themselves
 *
 * @section LICENSE
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

#ifndef __AKANTU_DUMPABLE_DUMMY_HH__
#define __AKANTU_DUMPABLE_DUMMY_HH__
/* -------------------------------------------------------------------------- */

namespace akantu {

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
  Dumpable(){};
  virtual ~Dumpable(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  template <class T>
  inline void registerDumper(const std::string & dumper_name,
                             const std::string & file_name = "",
                             const bool is_default = false) {}

  void registerExternalDumper(std::shared_ptr<DumperIOHelper> dumper,
                              const std::string & dumper_name,
                              const bool is_default = false) {}

  void addDumpMesh(const Mesh & mesh, UInt spatial_dimension = _all_dimensions,
                   const GhostType & ghost_type = _not_ghost,
                   const ElementKind & element_kind = _ek_not_defined) {}

  void addDumpMeshToDumper(const std::string & dumper_name, const Mesh & mesh,
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

  void addDumpFilteredMeshToDumper(
      const std::string & dumper_name, const Mesh & mesh,
      const ElementTypeMapArray<UInt> & elements_filter,
      const Array<UInt> & nodes_filter,
      UInt spatial_dimension = _all_dimensions,
      const GhostType & ghost_type = _not_ghost,
      const ElementKind & element_kind = _ek_not_defined) {}

  virtual void addDumpField(const std::string & field_id) {
    AKANTU_TO_IMPLEMENT();
  }
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
                                    const std::string & field_id) {
    AKANTU_TO_IMPLEMENT();
  }

  virtual void addDumpFieldExternal(const std::string & field_id,
                                    std::shared_ptr<dumper::Field> field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }
  virtual void
  addDumpFieldExternalToDumper(const std::string & dumper_name,
                               const std::string & field_id,
                               std::shared_ptr<dumper::Field> field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }

  template <typename T>
  void addDumpFieldExternal(const std::string & field_id,
                            const Array<T> & field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }
  template <typename T>
  void addDumpFieldExternalToDumper(const std::string & dumper_name,
                                    const std::string & field_id,
                                    const Array<T> & field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }

  template <typename T>
  void
  addDumpFieldExternal(const std::string & field_id,
                       const ElementTypeMapArray<T> & field,
                       UInt spatial_dimension = _all_dimensions,
                       const GhostType & ghost_type = _not_ghost,
                       const ElementKind & element_kind = _ek_not_defined) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }
  template <typename T>
  void addDumpFieldExternalToDumper(
      const std::string & dumper_name, const std::string & field_id,
      const ElementTypeMapArray<T> & field,
      UInt spatial_dimension = _all_dimensions,
      const GhostType & ghost_type = _not_ghost,
      const ElementKind & element_kind = _ek_not_defined) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }

  void removeDumpField(const std::string & field_id) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }
  void removeDumpFieldFromDumper(const std::string & dumper_name,
                                 const std::string & field_id) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }

  void setDirecory(const std::string & directory) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }
  void setDirectoryToDumper(const std::string & dumper_name,
                            const std::string & directory) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }

  void setBaseName(const std::string & basename) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }
  void setBaseNameToDumper(const std::string & dumper_name,
                           const std::string & basename) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }

  void setTextModeToDumper(const std::string & dumper_name) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }
  void setTextModeToDumper() {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }

  void dump() {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }
  void dump(const std::string & dumper_name) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }

  void dump(UInt step) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }
  void dump(const std::string & dumper_name, UInt step) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }

  void dump(Real current_time, UInt step) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }
  void dump(const std::string & dumper_name, Real current_time, UInt step) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }

protected:
  void internalAddDumpFieldToDumper(const std::string & dumper_name,
                                    const std::string & field_id,
                                    std::shared_ptr<dumper::Field> field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on "
                         "AKANTU_USE_IOHELPER in cmake.");
  }

protected:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  DumperIOHelper & getDumper() {
    AKANTU_ERROR("No dumper activated at compilation, turn on "
                 "AKANTU_USE_IOHELPER in cmake.");
  }

  DumperIOHelper & getDumper(const std::string & dumper_name) {
    AKANTU_ERROR("No dumper activated at compilation, turn on "
                 "AKANTU_USE_IOHELPER in cmake.");
  }

  template <class T> T & getDumper(const std::string & dumper_name) {
    AKANTU_ERROR("No dumper activated at compilation, turn on "
                 "AKANTU_USE_IOHELPER in cmake.");
  }

  std::string getDefaultDumperName() {
    AKANTU_ERROR("No dumper activated at compilation, turn on "
                 "AKANTU_USE_IOHELPER in cmake.");
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

#pragma GCC diagnostic pop

} // namespace akantu

#endif /* __AKANTU_DUMPABLE_DUMMY_HH__ */
