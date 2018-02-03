/**
 * @file   ntrf_contact.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  contact for node to rigid flat interface
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AST_NTRF_CONTACT_HH__
#define __AST_NTRF_CONTACT_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_base_contact.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
class NTRFContact : public NTNBaseContact {
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

  // add nodes from a list
  void addNodes(Array<UInt> & nodes);

  /// update (copy the normal to all normals)
  virtual void updateNormals();

  /// update the impedance matrix
  virtual void updateImpedance();

  /// compute the normal gap
  virtual void computeNormalGap(Array<Real> & gap) const;

  /// compute relative normal field (only value that has to be multiplied with the normal)
  /// relative to master nodes
  virtual void computeRelativeNormalField(const Array<Real> & field,
					  Array<Real> & rel_normal_field) const;

  /// compute relative tangential field (complet array)
  /// relative to master nodes
  virtual void computeRelativeTangentialField(const Array<Real> & field,
					      Array<Real> & rel_tang_field) const;

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

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

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// reference point for rigid flat surface
  Vector<Real> reference_point;
  /// outpointing normal of rigid flat surface
  Vector<Real> normal;
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

} // namespace akantu

#endif /* __AST_NTRF_CONTACT_HH__ */
