/**
 * @file   ntn_contact.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Feb 20 15:13:23 2012
 *
 * @brief  contact for node to node discretization
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AST_NTN_CONTACT_HH__
#define __AST_NTN_CONTACT_HH__

/* -------------------------------------------------------------------------- */
// akantu
#include "solid_mechanics_model.hh"

// simtools
#include "syncronized_array.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
using namespace akantu;

/* -------------------------------------------------------------------------- */
class NTNContact : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  NTNContact(SolidMechanicsModel & model,
	     const ContactID & id = "contact",
	     const MemoryID & memory_id = 0);
  virtual ~NTNContact() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// add surface pair and pair nodes according to the surface normal
  void addSurfacePair(Surface slave, Surface master, UInt surface_normal_dir);
  
  /// add node pair
  void addNodePair(UInt slave, UInt master);

  // add node pairs from a list with pairs(*,0)=slaves and pairs(*,1)=masters
  void addNodePairs(Array<UInt> & pairs);

  /// update normals, lumped boundary, and impedance
  void updateInternalData();

  /// update (compute the normals on the master nodes)
  void updateNormals();

  /// update the lumped boundary B matrix
  void updateLumpedBoundary();

  /// update the impedance matrix
  void updateImpedance();

  /// compute the normal contact force
  void computeContactPressure();

  /// impose the normal contact force
  void applyContactPressure();

  /// register syncronizedarrays for sync
  void registerSyncronizedArray(SyncronizedArrayBase & array);

  /// dump restart file
  void dumpRestart(std::string file_name) const;

  /// read restart file
  void readRestart(std::string file_name);

  /// compute relative normal field (only value that has to be multiplied with the normal)
  /// relative to master nodes
  void computeRelativeNormalField(const Array<Real> & field,
				  Array<Real> & rel_normal_field) const;

  /// compute relative tangential field (complet array)
  /// relative to master nodes
  void computeRelativeTangentialField(const Array<Real> & field,
				      Array<Real> & rel_tang_field) const;
  
  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;
  
protected:
  /// synchronize arrays
  void syncArrays(SyncChoice sync_choice);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Model, model, SolidMechanicsModel &)

  AKANTU_GET_MACRO(Slaves,                    slaves, const SyncronizedArray<UInt> &)
  AKANTU_GET_MACRO(Masters,                  masters, const SyncronizedArray<UInt> &)
  AKANTU_GET_MACRO(Normals,                  normals, const SyncronizedArray<Real> &)
  AKANTU_GET_MACRO(ContactPressure, contact_pressure, const SyncronizedArray<Real> &)
  AKANTU_GET_MACRO(LumpedBoundary,   lumped_boundary, const SyncronizedArray<Real> &)
  AKANTU_GET_MACRO(Impedance,              impedance, const SyncronizedArray<Real> &)
  AKANTU_GET_MACRO(IsInContact,        is_in_contact, const SyncronizedArray<bool> &)

  /// get number of contact nodes
  UInt getNbContactNodes() const { return this->slaves.getSize(); };

  /// get index of node in either slaves or masters array
  /// if node is in neither of them, return -1
  Int getNodeIndex(UInt node) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  ContactID id;

  SolidMechanicsModel & model;

  /// array of slave nodes
  SyncronizedArray<UInt> slaves;
  /// array of master nodes
  SyncronizedArray<UInt> masters;
  /// array of normals on master nodes
  SyncronizedArray<Real> normals;
  /// array indicating if nodes are in contact
  SyncronizedArray<Real> contact_pressure;
  /// array indicating if nodes are in contact
  SyncronizedArray<bool> is_in_contact;
  /// boundary matrix, lumped_boundary[:,0] master nodes, lumped_boundary[:,1] slave nodes
  SyncronizedArray<Real> lumped_boundary;
  /// impedance matrix
  SyncronizedArray<Real> impedance;

  /// contact surface
  std::set<Surface> contact_surfaces;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntn_contact_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const NTNContact & _this)
{
  _this.printself(stream);
  return stream;
}

__END_SIMTOOLS__

#endif /* __AST_NTN_CONTACT_HH__ */
