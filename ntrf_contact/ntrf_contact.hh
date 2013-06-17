/**
 * @file   ntrf_contact.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Fri Nov  2 17:33:49 2012
 *
 * @brief  contact for node to rigid flat interface
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

#ifndef __AST_NTRF_CONTACT_HH__
#define __AST_NTRF_CONTACT_HH__

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "syncronized_array.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
using namespace akantu;

/* -------------------------------------------------------------------------- */
class NTRFContact : protected Memory, public Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  NTRFContact(SolidMechanicsModel & model,
	      const ContactID & id = "contact",
	      const MemoryID & memory_id = 0);
  virtual ~NTRFContact() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void setReferencePoint(Real x=0., Real y=0., Real z=0.);
  void setNormal(Real x=1., Real y=0., Real z=0.);

  /// add surface and nodes according to the surface normal
  void addSurface(const Surface & surf);

  /// add node
  void addNode(UInt node);

  // add nodes from a list
  void addNodes(Array<UInt> & nodes);

  /// update normals, lumped boundary, and impedance
  void updateInternalData();

  /// update the lumped boundary B matrix
  void updateLumpedBoundary();

  /// compute the normal contact force
  void computeContactPressure();

  /// impose the normal contact force
  void applyContactPressure();

  /// register syncronizedarrays for sync
  void registerSyncronizedArray(SyncronizedArrayBase & array);

  /// dump restart file
  void dumpRestart(const std::string & file_name) const;

  /// read restart file
  void readRestart(const std::string & file_name);

  /// compute the normal gap
  void computeNormalGap(Array<Real> & gap) const;

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
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
				    const std::string & field_id);
  //  virtual void addDumpFieldVector(const std::string & field_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Model, model, SolidMechanicsModel &)

  AKANTU_GET_MACRO(Slaves,                    slaves, const SyncronizedArray<UInt> &)
  AKANTU_GET_MACRO(ContactPressure, contact_pressure, const SyncronizedArray<Real> &)
  AKANTU_GET_MACRO(LumpedBoundary,   lumped_boundary, const SyncronizedArray<Real> &)
  AKANTU_GET_MACRO(IsInContact,        is_in_contact, const SyncronizedArray<bool> &)

  AKANTU_GET_MACRO(Elements, elements, const ByElementTypeArray<UInt> &)

  /// get number of contact nodes: nodes in the system locally (on this proc)
  UInt getNbContactNodes() const { return this->slaves.getSize(); };

  /// get number of nodes that are in contact (globally, on all procs together)
  UInt getNbNodesInContact() const;

  /// get index of node in either slaves or masters array
  /// if node is in neither of them, return -1
  Int getNodeIndex(UInt node) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  typedef std::set<const SubBoundary *> SurfacePtrSet;
  ContactID id;

  SolidMechanicsModel & model;

  /// array of slave nodes
  SyncronizedArray<UInt> slaves;
  /// array indicating if nodes are in contact
  SyncronizedArray<Real> contact_pressure;
  /// array indicating if nodes are in contact
  SyncronizedArray<bool> is_in_contact;
  /// boundary matrix, lumped_boundary[:,0] slave nodes
  SyncronizedArray<Real> lumped_boundary;

  /// reference point for rigid flat surface
  Array<Real> reference_point;
  /// outpointing normal of rigid flat surface
  Array<Real> normal;

  /// contact surface
  SurfacePtrSet contact_surfaces;

  /// element list for dump
  ByElementTypeArray<UInt> elements;
  CSR<Element> node_to_elements;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntrf_contact_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const NTRFContact & _this)
{
  _this.printself(stream);
  return stream;
}

__END_SIMTOOLS__

#endif /* __AST_NTRF_CONTACT_HH__ */
