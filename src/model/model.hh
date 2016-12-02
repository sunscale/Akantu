/**
 * @file   model.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Oct 16 2015
 *
 * @brief  Interface of a model
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "aka_memory.hh"
#include "fe_engine.hh"
#include "mesh.hh"
#include "mesh_partition.hh"
#include "mesh_utils.hh"
#include "model_solver.hh"
#include "parser.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MODEL_HH__
#define __AKANTU_MODEL_HH__

namespace akantu {
class SynchronizerRegistry;
class DataAccessorBase;
template <class Entity> class DataAccessor;
} // akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

struct ModelOptions {
  virtual ~ModelOptions() {}
};

class DumperIOHelper;

class Model : public Memory, public ModelSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  Model(Mesh & mesh, UInt spatial_dimension = _all_dimensions,
        const ID & id = "model", const MemoryID & memory_id = 0);

  virtual ~Model();

  typedef std::map<std::string, FEEngine *> FEEngineMap;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void initFull(const ModelOptions & options);

  virtual void initModel() = 0;

  // /// change local equation number so that PBC is assembled properly
  // void changeLocalEquationNumberForPBC(std::map<UInt, UInt> & pbc_pair,
  //                                      UInt dimension);
  /// function to print the containt of the class
  virtual void printself(__attribute__((unused)) std::ostream & stream,
                         __attribute__((unused)) int indent = 0) const {};

  // /// initialize the model for PBC
  // void setPBC(UInt x, UInt y, UInt z);
  // void setPBC(SurfacePairList & surface_pairs);

  virtual void initPBC();

  /// set the parser to use
  void setParser(Parser & parser);

  /* ------------------------------------------------------------------------ */
  /* Access to the dumpable interface of the boundaries                       */
  /* ------------------------------------------------------------------------ */
  /// Dump the data for a given group
  void dumpGroup(const std::string & group_name);
  void dumpGroup(const std::string & group_name,
                 const std::string & dumper_name);
  /// Dump the data for all boundaries
  void dumpGroup();
  /// Set the directory for a given group
  void setGroupDirectory(const std::string & directory,
                         const std::string & group_name);
  /// Set the directory for all boundaries
  void setGroupDirectory(const std::string & directory);
  /// Set the base name for a given group
  void setGroupBaseName(const std::string & basename,
                        const std::string & group_name);
  /// Get the internal dumper of a given group
  DumperIOHelper & getGroupDumper(const std::string & group_name);

  /* ------------------------------------------------------------------------ */
  /* Function for non local capabilities                                      */
  /* ------------------------------------------------------------------------ */
  virtual void updateDataForNonLocalCriterion(__attribute__((unused))
                                              ElementTypeMapReal & criterion) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors */
  /* ------------------------------------------------------------------------ */
public:
  /// get id of model
  AKANTU_GET_MACRO(ID, id, const ID)

  /// get the number of surfaces
  AKANTU_GET_MACRO(Mesh, mesh, Mesh &)

  /// synchronize the boundary in case of parallel run
  virtual void synchronizeBoundaries(){};

  /// return the fem object associated with a provided name
  inline FEEngine & getFEEngine(const ID & name = "") const;

  /// return the fem boundary object associated with a provided name
  virtual FEEngine & getFEEngineBoundary(const ID & name = "");

  /// register a fem object associated with name
  template <typename FEEngineClass>
  inline void registerFEEngineObject(const std::string & name, Mesh & mesh,
                                     UInt spatial_dimension);
  /// unregister a fem object associated with name
  inline void unRegisterFEEngineObject(const std::string & name);

  /// return the synchronizer registry
  SynchronizerRegistry & getSynchronizerRegistry();

  /// return the fem object associated with a provided name
  template <typename FEEngineClass>
  inline FEEngineClass & getFEEngineClass(std::string name = "") const;

  /// return the fem boundary object associated with a provided name
  template <typename FEEngineClass>
  inline FEEngineClass & getFEEngineClassBoundary(std::string name = "");

  /// get the pbc pairs
  std::map<UInt, UInt> & getPBCPairs() { return pbc_pair; };

  /// returns if node is slave in pbc
  inline bool isPBCSlaveNode(const UInt node) const {
    return false; /* TODO repair PBC*/
  }

  /// returns the array of pbc slave nodes (boolean information)
  AKANTU_GET_MACRO(IsPBCSlaveNode, is_pbc_slave_node, const Array<bool> &)

  /* ------------------------------------------------------------------------ */
  /* Pack and unpack helper functions                                         */
  /* ------------------------------------------------------------------------ */
public:
  inline UInt getNbIntegrationPoints(const Array<Element> & elements,
                                     const ID & fem_id = ID()) const;

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface (kept for convenience) and dumper relative functions  */
  /* ------------------------------------------------------------------------ */

  void setTextModeToDumper();

  virtual void addDumpGroupFieldToDumper(const std::string & field_id,
                                         dumper::Field * field,
                                         DumperIOHelper & dumper);

  virtual void addDumpField(const std::string & field_id);

  virtual void addDumpFieldVector(const std::string & field_id);

  virtual void addDumpFieldToDumper(const std::string & dumper_name,
                                    const std::string & field_id);

  virtual void addDumpFieldVectorToDumper(const std::string & dumper_name,
                                          const std::string & field_id);

  virtual void addDumpFieldTensorToDumper(const std::string & dumper_name,
                                          const std::string & field_id);

  virtual void addDumpFieldTensor(const std::string & field_id);

  virtual void setBaseName(const std::string & basename);

  virtual void setBaseNameToDumper(const std::string & dumper_name,
                                   const std::string & basename);

  virtual void addDumpGroupField(const std::string & field_id,
                                 const std::string & group_name);

  virtual void addDumpGroupFieldToDumper(const std::string & dumper_name,
                                         const std::string & field_id,
                                         const std::string & group_name,
                                         const ElementKind & element_kind,
                                         bool padding_flag);

  virtual void addDumpGroupFieldToDumper(const std::string & dumper_name,
                                         const std::string & field_id,
                                         const std::string & group_name,
                                         UInt spatial_dimension,
                                         const ElementKind & element_kind,
                                         bool padding_flag);

  virtual void removeDumpGroupField(const std::string & field_id,
                                    const std::string & group_name);
  virtual void removeDumpGroupFieldFromDumper(const std::string & dumper_name,
                                              const std::string & field_id,
                                              const std::string & group_name);

  virtual void addDumpGroupFieldVector(const std::string & field_id,
                                       const std::string & group_name);

  virtual void addDumpGroupFieldVectorToDumper(const std::string & dumper_name,
                                               const std::string & field_id,
                                               const std::string & group_name);

  virtual dumper::Field *
  createNodalFieldReal(__attribute__((unused)) const std::string & field_name,
                       __attribute__((unused)) const std::string & group_name,
                       __attribute__((unused)) bool padding_flag) {
    return NULL;
  }

  virtual dumper::Field *
  createNodalFieldUInt(__attribute__((unused)) const std::string & field_name,
                       __attribute__((unused)) const std::string & group_name,
                       __attribute__((unused)) bool padding_flag) {
    return NULL;
  }

  virtual dumper::Field *
  createNodalFieldBool(__attribute__((unused)) const std::string & field_name,
                       __attribute__((unused)) const std::string & group_name,
                       __attribute__((unused)) bool padding_flag) {
    return NULL;
  }

  virtual dumper::Field *
  createElementalField(__attribute__((unused)) const std::string & field_name,
                       __attribute__((unused)) const std::string & group_name,
                       __attribute__((unused)) bool padding_flag,
                       __attribute__((unused)) const UInt & spatial_dimension,
                       __attribute__((unused)) const ElementKind & kind) {
    return NULL;
  }

  void setDirectory(const std::string & directory);
  void setDirectoryToDumper(const std::string & dumper_name,
                            const std::string & directory);

  virtual void dump();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Mesh
  Mesh & mesh;

  /// Spatial dimension of the problem
  UInt spatial_dimension;

  /// the main fem object present in all  models
  FEEngineMap fems;

  /// the fem object present in all  models for boundaries
  FEEngineMap fems_boundary;

  /// default fem object
  std::string default_fem;

  /// pbc pairs
  std::map<UInt, UInt> pbc_pair;

  /// flag per node to know is pbc slave
  Array<bool> is_pbc_slave_node;

  /// parser to the pointer to use
  Parser * parser;
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream, const Model & _this) {
  _this.printself(stream);
  return stream;
}

} // akantu

#include "model_inline_impl.cc"

#endif /* __AKANTU_MODEL_HH__ */
