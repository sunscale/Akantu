/**
 * @file   ntn_contact.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  contact for node to node discretization
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AST_NTN_CONTACT_HH_
#define AST_NTN_CONTACT_HH_

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_base_contact.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
class NTNContact : public NTNBaseContact {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNContact(SolidMechanicsModel & model, const ID & id = "contact");
  ~NTNContact() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// add surface pair and pair nodes according to the surface normal
  void addSurfacePair(const ID & slave, const ID & master,
                      UInt surface_normal_dir);

  /// fills the pairs vector with interface node pairs (*,0)=slaves,
  /// (*,1)=masters
  static void pairInterfaceNodes(const ElementGroup & slave_boundary,
                                 const ElementGroup & master_boundary,
                                 UInt surface_normal_dir, const Mesh & mesh,
                                 Array<UInt> & pairs);

  // add node pairs from a list with pairs(*,0)=slaves and pairs(*,1)=masters
  void addNodePairs(const Array<UInt> & pairs);

  /// add node pair
  void addSplitNode(UInt slave, UInt master) override;

  /// update (compute the normals on the master nodes)
  void updateNormals() override;

  /// update the lumped boundary B matrix
  void updateLumpedBoundary() override;

  /// update the impedance matrix
  void updateImpedance() override;

  /// impose the normal contact force
  void applyContactPressure() override;

  /// dump restart file
  void dumpRestart(const std::string & file_name) const override;

  /// read restart file
  void readRestart(const std::string & file_name) override;

  /// compute the normal gap
  void computeNormalGap(Array<Real> & gap) const override {
    this->computeRelativeNormalField(this->model.getCurrentPosition(), gap);
  };

  /// compute relative normal field (only value that has to be multiplied with
  /// the normal)
  /// relative to master nodes
  void
  computeRelativeNormalField(const Array<Real> & field,
                             Array<Real> & rel_normal_field) const override;

  /// compute relative tangential field (complet array)
  /// relative to master nodes
  void
  computeRelativeTangentialField(const Array<Real> & field,
                                 Array<Real> & rel_tang_field) const override;

  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

protected:
  /// synchronize arrays
  void syncArrays(SyncChoice sync_choice) override;

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  void addDumpFieldToDumper(const std::string & dumper_name,
                            const std::string & field_id) override;
  //  virtual void addDumpFieldVector(const std::string & field_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Masters, masters, const SynchronizedArray<UInt> &)
  AKANTU_GET_MACRO(LumpedBoundaryMasters, lumped_boundary_masters,
                   const SynchronizedArray<Real> &)

  /// get interface node pairs (*,0) are slaves, (*,1) are masters
  void getNodePairs(Array<UInt> & pairs) const;

  /// get index of node in either slaves or masters array
  /// if node is in neither of them, return -1
  Int getNodeIndex(UInt node) const override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// array of master nodes
  SynchronizedArray<UInt> masters;
  /// lumped boundary of master nodes
  SynchronizedArray<Real> lumped_boundary_masters;

  // element list for dump and lumped_boundary
  ElementTypeMapArray<UInt> master_elements;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntn_contact_inline_impl.hh"

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const NTNContact & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* AST_NTN_CONTACT_HH_ */
