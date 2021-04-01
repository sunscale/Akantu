/**
 * @file   phasefield.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Mar 2 2020
 * @date last modification: Mon Mar 2 2020
 *
 * @brief  Mother class for all phasefield laws
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
#include "aka_factory.hh"
#include "aka_memory.hh"
#include "data_accessor.hh"
#include "parsable.hh"
#include "parser.hh"
/* -------------------------------------------------------------------------- */
#include "internal_field.hh"
#include "random_internal_field.hh"
/* -------------------------------------------------------------------------- */


#ifndef __AKANTU_PHASEFIELD_HH__
#define __AKANTU_PHASEFIELD_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {
  class Model;
  class PhaseFieldModel;
  class PhaseField;
} // namespace akantu

namespace akantu {

template<typename T>
using InternalPhaseField = InternalFieldTmpl<PhaseField, T>;
  
using PhaseFieldFactory =
  Factory<PhaseField, ID, UInt, const ID &, PhaseFieldModel &, const ID &>;
  
class PhaseField : public Memory,
		   public DataAccessor<Element>,
		   public Parsable {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PhaseField(const PhaseField & phase) = delete;
  PhaseField & operator=(const PhaseField & phase) =  delete;

  /// Initialize phasefield with defaults
  PhaseField(PhaseFieldModel & model, const ID & id = "");

  /// Initialize phasefield with custom mesh & fe_engine
  PhaseField(PhaseFieldModel & model, UInt dim, const Mesh & mesh,
	     FEEngine & fe_engine, const ID & id = "");

  /// Destructor
  ~PhaseField() override;


protected:
  void initialize();

  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  template <typename T> void registerInternal(InternalPhaseField<T> & /*vect*/) {
    AKANTU_TO_IMPLEMENT();
  }

  template <typename T> void unregisterInternal(InternalPhaseField<T> & /*vect*/) {
    AKANTU_TO_IMPLEMENT();
  }

  
  /// initialize the phasefield computed parameter
  virtual void initPhaseField();

  /// 
  virtual void beforeSolveStep();

  ///
  virtual void afterSolveStep();

  /// assemble the residual for this phasefield
  virtual void assembleInternalForces(GhostType ghost_type);

  /// assemble the stiffness matrix for this phasefield
  virtual void assembleStiffnessMatrix(GhostType ghost_type);

  /// compute the driving force for this phasefield
  virtual void computeAllDrivingForces(GhostType ghost_type = _not_ghost);
  
  /// save the phi in the phi internal field if needed
  virtual void savePreviousState();
  
  /// add an element to the local mesh filter
  inline UInt addElement(const ElementType & type, UInt element,
                         const GhostType & ghost_type);
  inline UInt addElement(const Element & element);

  
  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  
protected:
  /// resize the internals arrrays
  virtual void resizeInternals();

  /// function called to updatet the internal parameters when the
  /// modifiable parameters are modified   
  virtual void updateInternalParameters();

  // constitutive law for driving force
  virtual void computeDrivingForce(const ElementType & /* el_type */,
				   GhostType /* ghost_type */ = _not_ghost) {
    AKANTU_TO_IMPLEMENT();
  }

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  template <typename T>
  inline void packElementDataHelper(const ElementTypeMapArray<T> & data_to_pack,
                                    CommunicationBuffer & buffer,
                                    const Array<Element> & elements,
                                    const ID & fem_id = ID()) const;

  template <typename T>
  inline void unpackElementDataHelper(ElementTypeMapArray<T> & data_to_unpack,
                                      CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      const ID & fem_id = ID());


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Name, name, const std::string &);

  AKANTU_GET_MACRO(Model, model, const PhaseFieldModel &)

  AKANTU_GET_MACRO(ID, Memory::getID(), const ID &);
  
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Strain, strain, Real);

  AKANTU_GET_MACRO(Strain, strain, const ElementTypeMapArray<Real> &);

  AKANTU_GET_MACRO_NOT_CONST(Strain, strain, ElementTypeMapArray<Real> &);

  
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Damage, damage, Real);

  AKANTU_GET_MACRO_NOT_CONST(Damage, damage, ElementTypeMapArray<Real> &);
  AKANTU_GET_MACRO(Damage, damage, const ElementTypeMapArray<Real> &);
  
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementFilter, element_filter, UInt);

  AKANTU_GET_MACRO(ElementFilter, element_filter,
                   const ElementTypeMapArray<UInt> &);

  template <typename T>
  const InternalPhaseField<T> & getInternal(const ID & id) const;
  template <typename T> InternalPhaseField<T> & getInternal(const ID & id);

  template <typename T>
  inline bool isInternal(const ID & id, const ElementKind & element_kind) const;

  
  template <typename T> inline void setParam(const ID & param, T value);
  inline const Parameter & getParam(const ID & param) const;

  template <typename T>
  void flattenInternal(const std::string & field_id,
                       ElementTypeMapArray<T> & internal_flat,
                       const GhostType ghost_type = _not_ghost,
                       ElementKind element_kind = _ek_not_defined) const;


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// boolean to know if the material has been initialized
  bool is_init;

  std::map<ID, InternalPhaseField<Real> *> internal_vectors_real;
  std::map<ID, InternalPhaseField<UInt> *> internal_vectors_uint;
  std::map<ID, InternalPhaseField<bool> *> internal_vectors_bool;

  
protected:
  /// Link to the fem object in the model
  FEEngine & fem;

  /// phasefield name
  std::string name;

  /// The model to whch the phasefield belong
  PhaseFieldModel & model;

  /// length scale parameter
  Real l0;

  /// critical energy release rate
  Real g_c;

  /// Young's modulus
  Real E;

  /// Poisson ratio
  Real nu;

  /// Lame's first parameter
  Real lambda;

  /// Lame's second paramter
  Real mu;
  
  /// spatial dimension
  UInt spatial_dimension;

  /// list of element handled by the phasefield
  ElementTypeMapArray<UInt> element_filter;

  /// damage arrays ordered by element types
  InternalPhaseField<Real> damage;
  
  /// phi arrays ordered by element types
  InternalPhaseField<Real> phi;
  
  /// strain arrays ordered by element types
  InternalPhaseField<Real> strain;

  /// driving force ordered by element types
  InternalPhaseField<Real> driving_force;

  /// damage energy ordered by element types
  InternalPhaseField<Real> damage_energy;

  /// damage energy density ordered by element types
  InternalPhaseField<Real> damage_energy_density;

  
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const PhaseField & _this) {
  _this.printself(stream);
  return stream;
}


} // namespace akantu


#include "phasefield_inline_impl.cc"

#include "internal_field_tmpl.hh"
#include "random_internal_field_tmpl.hh"

/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
/* Auto loop                                                                  */
/* -------------------------------------------------------------------------- */
/// This can be used to automatically write the loop on quadrature points in
/// functions such as computeStress. This macro in addition to write the loop
/// provides two tensors (matrices) sigma and grad_u

#define PHASEFIELD_PHI_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type)	\
      auto && strain_view =						\
	make_view(this->strain(el_type, ghost_type), this->spatial_dimension, \
		  this->spatial_dimension);				\
      auto phi_view =							\
	make_view(this->phi(el_type, ghost_type));			\
      auto && phi_hist_view =						\
	make_view(this->phi.previous(el_type, ghost_type));		\
									\
      for (auto && data : zip(strain_view, phi_view, phi_hist_view)) {	\
	[[gnu::unused]] Matrix<Real> & strain_q = std::get<0>(data);	\
	[[gnu::unused]] Real & phi_q = std::get<1>(data);		\
	[[gnu::unused]] Real & phi_hist_q = std::get<2>(data);

#define PHASEFIELD_PHI_QUADRATURE_POINT_LOOP_END }

#define PHASEFIELD_ENERGY_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type) \
      auto && eng_den_view =						\
	make_view(this->damage_energy_density(el_type, ghost_type));	\
      auto phi_view =							\
	make_view(this->phi(el_type, ghost_type));			\
									\
      for (auto && data : zip(eng_den_view, phi_view)) {		\
 	[[gnu::unused]] Real & dam_eng_density = std::get<0>(data);	\
	[[gnu::unused]] Real & phi_q = std::get<1>(data);

#define PHASEFIELD_ENERGY_QUADRATURE_POINT_LOOP_END }
	

#define INSTANTIATE_PHASEFIELD_ONLY(phase_name)				\
  template class phase_name<1>;						\
  template class phase_name<2>;						\
  template class phase_name<3>

#define PHASEFIELD_DEFAULT_PER_DIM_ALLOCATOR(id, phase_name)		\
  [](UInt dim, const ID &, PhaseFieldModel & model,			\
     const ID & id) -> std::unique_ptr<PhaseField> {			\
    switch (dim) {							\
    case 1:								\
    case 2:								\
    case 3:								\
      return std::make_unique<phase_name>(model, id);			\
    default:								\
      AKANTU_EXCEPTION("The dimension "					\
                       << dim << "is not a valid dimension for the phasefield "	\
                       << #id);						\
    }									\
  }

#define INSTANTIATE_PHASEFIELD(id, phase_name)				\
  static bool phasefield_is_alocated_##id [[gnu::unused]] =		\
				 PhaseFieldFactory::getInstance().registerAllocator( \
          #id, PHASEFIELD_DEFAULT_PER_DIM_ALLOCATOR(id, phase_name))


#endif
