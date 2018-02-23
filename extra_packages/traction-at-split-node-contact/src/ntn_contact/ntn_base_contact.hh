/**
 * @file   ntn_base_contact.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  base contact for ntn and ntrf contact
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AST_NTN_BASE_CONTACT_HH__
#define __AST_NTN_BASE_CONTACT_HH__

/* -------------------------------------------------------------------------- */
// akantu
#include "filtered_synchronizer.hh"
#include "solid_mechanics_model.hh"

// simtools
#include "synchronized_array.hh"

namespace akantu {

class NTNBaseContact;

/* -------------------------------------------------------------------------- */
class NTNContactSynchElementFilter : public SynchElementFilter {
public:
  // constructor
  NTNContactSynchElementFilter(NTNBaseContact * contact);

  // answer to: do we need this element ?
  virtual bool operator()(const Element & e);

private:
  const NTNBaseContact * contact;
  const ElementTypeMapArray<UInt> & connectivity;
};

/* -------------------------------------------------------------------------- */
class NTNBaseContact : protected Memory, public DataAccessor, public Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNBaseContact(SolidMechanicsModel & model, const ContactID & id = "contact",
                 const MemoryID & memory_id = 0);
  virtual ~NTNBaseContact();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initializes ntn contact parallel
  virtual void initParallel();

  /// add split node
  virtual void addSplitNode(UInt node);

  /// update normals, lumped boundary, and impedance
  virtual void updateInternalData();

  /// update (compute the normals)
  virtual void updateNormals() = 0;

  /// update the lumped boundary B matrix
  virtual void updateLumpedBoundary();

  /// update the impedance matrix
  virtual void updateImpedance() = 0;

  /// compute the normal contact force
  virtual void computeContactPressure();

  /// impose the normal contact force
  virtual void applyContactPressure();

  /// register synchronizedarrays for sync
  virtual void registerSynchronizedArray(SynchronizedArrayBase & array);

  /// dump restart file
  virtual void dumpRestart(const std::string & file_name) const;

  /// read restart file
  virtual void readRestart(const std::string & file_name);

  /// compute the normal gap
  virtual void computeNormalGap(Array<Real> & gap) const = 0;

  /// compute relative normal field (only value that has to be multiplied with
  /// the normal)
  /// relative to master nodes
  virtual void
  computeRelativeNormalField(const Array<Real> & field,
                             Array<Real> & rel_normal_field) const = 0;

  /// compute relative tangential field (complet array)
  /// relative to master nodes
  virtual void
  computeRelativeTangentialField(const Array<Real> & field,
                                 Array<Real> & rel_tang_field) const = 0;

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// updateLumpedBoundary
  virtual void
  internalUpdateLumpedBoundary(const Array<UInt> & nodes,
                               const ElementTypeMapArray<UInt> & elements,
                               SynchronizedArray<Real> & boundary);

  // to find the slave_elements or master_elements
  virtual void findBoundaryElements(const Array<UInt> & interface_nodes,
                                    ElementTypeMapArray<UInt> & elements);

  /// synchronize arrays
  virtual void syncArrays(SyncChoice sync_choice);

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  inline virtual UInt getNbDataForElements(const Array<Element> & elements,
                                           SynchronizationTag tag) const;

  inline virtual void packElementData(CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      SynchronizationTag tag) const;

  inline virtual void unpackElementData(CommunicationBuffer & buffer,
                                        const Array<Element> & elements,
                                        SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
                                    const std::string & field_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Model, model, SolidMechanicsModel &)

  AKANTU_GET_MACRO(Slaves, slaves, const SynchronizedArray<UInt> &)
  AKANTU_GET_MACRO(Normals, normals, const SynchronizedArray<Real> &)
  AKANTU_GET_MACRO(ContactPressure, contact_pressure,
                   const SynchronizedArray<Real> &)
  AKANTU_GET_MACRO(LumpedBoundarySlaves, lumped_boundary_slaves,
                   const SynchronizedArray<Real> &)
  AKANTU_GET_MACRO(Impedance, impedance, const SynchronizedArray<Real> &)
  AKANTU_GET_MACRO(IsInContact, is_in_contact, const SynchronizedArray<bool> &)

  AKANTU_GET_MACRO(SlaveElements, slave_elements,
                   const ElementTypeMapArray<UInt> &)

  AKANTU_GET_MACRO(SynchronizerRegistry, synch_registry, SynchronizerRegistry *)

  /// get number of nodes that are in contact (globally, on all procs together)
  /// is_in_contact = true
  virtual UInt getNbNodesInContact() const;

  /// get index of node in either slaves or masters array
  /// if node is in neither of them, return -1
  virtual Int getNodeIndex(UInt node) const;

  /// get number of contact nodes: nodes in the system locally (on this proc)
  /// is_in_contact = true and false, because just in the system
  virtual UInt getNbContactNodes() const { return this->slaves.getSize(); };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  typedef std::set<const ElementGroup *> SurfacePtrSet;

  SolidMechanicsModel & model;

  /// array of slave nodes
  SynchronizedArray<UInt> slaves;
  /// array of normals
  SynchronizedArray<Real> normals;
  /// array indicating if nodes are in contact
  SynchronizedArray<Real> contact_pressure;
  /// array indicating if nodes are in contact
  SynchronizedArray<bool> is_in_contact;
  /// boundary matrix for slave nodes
  SynchronizedArray<Real> lumped_boundary_slaves;
  /// impedance matrix
  SynchronizedArray<Real> impedance;

  /// contact surface
  SurfacePtrSet contact_surfaces;

  /// element list for dump and lumped_boundary
  ElementTypeMapArray<UInt> slave_elements;
  CSR<Element> node_to_elements;

  /// parallelisation
  SynchronizerRegistry * synch_registry;
  Synchronizer * synchronizer;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "ntn_base_contact_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const NTNBaseContact & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* __AST_NTN_BASE_CONTACT_HH__ */
