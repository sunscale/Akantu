/**
 * @file   material_phasefield.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Aug 14 2018
 * @date last modification: Wed Aug 14 2018
 *
 * @brief  Mother class for all materials in phasefield model
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
#include "aka_voigthelper.hh"
#include "data_accessor.hh"
#include "integration_point.hh"
#include "parsable.hh"
#include "parser.hh"
/* -------------------------------------------------------------------------- */
#include "internal_field.hh"
#include "random_internal_field.hh"
/* -------------------------------------------------------------------------- */
#include "mesh_events.hh"

#ifndef __AKANTU_MATERIAL_PHASEFIELD_HH__
#define __AKANTU_MATERIAL_PHASEFIELD_HH__


namespace akantu {
class Model;
class PhaseFieldModel;
} // namespace akantu

namespace akantu {

  
/**
 * Interface of all materials
 * Prerequisites for a new material
 * - inherit from this class
 * - implement the following methods:
 * \code
 *  virtual Real getStableTimeStep(Real h, const Element & element =
 * ElementNull);
 *
 *  virtual void computeStress(ElementType el_type,
 *                             GhostType ghost_type = _not_ghost);
 *
 *  virtual void computeTangentStiffness(const ElementType & el_type,
 *                                       Array<Real> & tangent_matrix,
 *                                       GhostType ghost_type = _not_ghost);
 * \endcode
 *
 */
class MaterialPhaseField : public Memory,
			     public DataAccessor<Element>,
			     public Parsable,
			     public MeshEventHandler {
			     
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
#if __cplusplus > 199711L
  MaterialPhaseField(const MaterialPhaseField & mat) = delete;
  MaterialPhaseField & operator=(const MaterialPhaseField & mat) = delete;
#endif

  /// Initialize material with defaults
  MaterialPhaseField(PhaseFieldModel & model, const ID & id = "");

  /// Initialize material with custom mesh & fe_engine
  MaterialPhaseField(PhaseFieldModel & model, UInt dim, const Mesh & mesh,
		     FEEngine & fe_engine, const ID & id = "");

  /// Destructor
  ~MaterialPhaseField() override;

  /* ------------------------------------------------------------------------ */
  /* Function that materials can/should reimplement                           */
  /* ------------------------------------------------------------------------ */
protected:
  /// constitutive law
  virtual void computeStress(__attribute__((unused)) ElementType el_type,
                             __attribute__((unused))
                             GhostType ghost_type = _not_ghost) {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the tangent stiffness matrix
  virtual void computeTangentModuli(__attribute__((unused))
                                    const ElementType & el_type,
                                    __attribute__((unused))
                                    Array<Real> & tangent_matrix,
                                    __attribute__((unused))
                                    GhostType ghost_type = _not_ghost) {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the fracture energy
  virtual void computeFractureEnergy(ElementType el_type,
                                      GhostType ghost_type = _not_ghost);

  /// compute the fracture energy for an element
  virtual void
  computeFractureEnergyByElement(__attribute__((unused)) ElementType type,
				 __attribute__((unused)) UInt index,
				 __attribute__((unused))
				 Vector<Real> & epot_on_quad_points) {
    AKANTU_TO_IMPLEMENT();
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  template <typename T>
  void registerInternal(__attribute__((unused)) InternalField<T> & vect) {
    AKANTU_TO_IMPLEMENT();
  }

  template <typename T>
  void unregisterInternal(__attribute__((unused)) InternalField<T> & vect) {
    AKANTU_TO_IMPLEMENT();
  }

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// aseemble the residual for this material
  virtual void assembleInternalForces(GhostType ghost_type);
  
  /// assemble the damage matrix
  virtual void assembleDamageMatrix(GhostType ghost_type);
  
  /// assemble the damage gradient matrix
  virtual void assembleDamageGradMatrix(GhostType ghost_type);
  
  /// add an element to the local mesh filter
  inline UInt addElement(const ElementType & type, UInt element,
                         const GhostType & ghost_type);
  inline UInt addElement(const Element & element);

  /// add many elements at once
  void addElements(const Array<Element> & elements_to_add);

  /// remove many element at once
  void removeElements(const Array<Element> & elements_to_remove);

  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /**
   * function to initialize the elemental field interpolation
   * function by inverting the quadrature points' coordinates
   */
  void initElementalFieldInterpolation(
      const ElementTypeMapArray<Real> & interpolation_points_coordinates);

  /* ------------------------------------------------------------------------ */
  /* Common part                                                              */
  /* ------------------------------------------------------------------------ */
protected:
  /* ------------------------------------------------------------------------ */
  
  /// compute the dissipated energy by element
  void computePotentialEnergyByElements();

  /// resize the intenals arrays
  virtual void resizeInternals();

protected:
  /// assemble the residual
  template <UInt dim>
  void assembleInternalForces(GhostType ghost_type);

  /// assemble the damage matrix
  template <UInt dim>
  void assembleDamageMatrix(const ElementType & type,
			    GhostType ghost_type);

  /// assemble the damage gradient matrix
  template <UInt dim>
  void assembleDamageGradMatrix(const ElementType & type,
				GhostType ghost_type);
  /* ------------------------------------------------------------------------ */
  /* Conversion functions                                                     */
  /* ------------------------------------------------------------------------ */
public:
  template <UInt dim>
  static inline void gradUToF(const Matrix<Real> & grad_u, Matrix<Real> & F);
  
  template <UInt dim>
  static inline void gradUToEpsilon(const Matrix<Real> & grad_u,
                                    Matrix<Real> & epsilon);
  template <UInt dim>
  static inline void gradUToGreenStrain(const Matrix<Real> & grad_u,
                                        Matrix<Real> & epsilon);

protected:
  /// converts global element to local element
  inline Element convertToLocalElement(const Element & global_element) const;
  /// converts local element to global element
  inline Element convertToGlobalElement(const Element & local_element) const;

  /// converts global quadrature point to local quadrature point
  inline IntegrationPoint
  convertToLocalPoint(const IntegrationPoint & global_point) const;
  /// converts local quadrature point to global quadrature point
  inline IntegrationPoint
  convertToGlobalPoint(const IntegrationPoint & local_point) const;


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
  /* MeshEventHandler inherited members                                       */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  void onNodesAdded(const Array<UInt> &, const NewNodesEvent &) override{};
  void onNodesRemoved(const Array<UInt> &, const Array<UInt> &,
                      const RemovedNodesEvent &) override{};
  void onElementsAdded(const Array<Element> & element_list,
                       const NewElementsEvent & event) override;
  void onElementsRemoved(const Array<Element> & element_list,
                         const ElementTypeMapArray<UInt> & new_numbering,
                         const RemovedElementsEvent & event) override;
  void onElementsChanged(const Array<Element> &, const Array<Element> &,
                         const ElementTypeMapArray<UInt> &,
                         const ChangedElementsEvent &) override{};


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Name, name, const std::string &);

  AKANTU_GET_MACRO(Model, model, const PhaseFieldModel &)

  AKANTU_GET_MACRO(ID, Memory::getID(), const ID &);

  /// get first lame coefficient
  AKANTU_GET_MACRO(E, E, Real);
  /// get Second lame coefficient 
  AKANTU_GET_MACRO(Nu, nu, Real);
  /// get Griffith fracture energy
  AKANTU_GET_MACRO(Gc, gc, Real);
  /// get length scale
  AKANTU_GET_MACRO(LengthScale, l0, Real);
  // get the dimensions
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /// return the dissipated energy for the subset of elements contained by the
  /// material
  Real getDissipatedEnergy();
  /// return the dissipated energy for the provided element
  Real getDissipatedEnergy(ElementType & type, UInt index);

  /// return the energy (identified by id) for the subset of elements contained
  /// by the material
  virtual Real getEnergy(const std::string & energy_id);
  /// return the energy (identified by id) for the provided element
  virtual Real getEnergy(const std::string & energy_id, ElementType type,
                         UInt index);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementFilter, element_filter, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(GradU, gradu, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Stress, stress, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PotentialEnergy, dissipated_energy,
                                         Real);
  AKANTU_GET_MACRO(GradU, gradu, const ElementTypeMapArray<Real> &);
  AKANTU_GET_MACRO(Stress, stress, const ElementTypeMapArray<Real> &);
  AKANTU_GET_MACRO(ElementFilter, element_filter,
                   const ElementTypeMapArray<UInt> &);
  AKANTU_GET_MACRO(FEEngine, fem, FEEngine &);

  bool isNonLocal() const { return is_non_local; }

  template <typename T>
  const Array<T> & getArray(const ID & id, const ElementType & type,
                            const GhostType & ghost_type = _not_ghost) const;
  template <typename T>
  Array<T> & getArray(const ID & id, const ElementType & type,
                      const GhostType & ghost_type = _not_ghost);

  template <typename T>
  const InternalField<T> & getInternal(const ID & id) const;
  template <typename T> InternalField<T> & getInternal(const ID & id);

  template <typename T>
  inline bool isInternal(const ID & id, const ElementKind & element_kind) const;

  template <typename T>
  ElementTypeMap<UInt>
  getInternalDataPerElem(const ID & id, const ElementKind & element_kind) const;

protected:
  /// Link to the fem object in the model
  FEEngine & fem;

  /// material name
  std::string name;

  /// The model to witch the material belong
  PhaseFieldModel & model;

  /// Young's modulus
  Real E;

  /// Poisson's ratio 
  Real nu;

  /// First Lamé coefficient
  Real lambda;
  
  /// Second Lamé coefficient (shear modulus)
  Real mu;

  /// Griffith's fracture energy
  Real gc;

  /// length scale
  Real l0;
  
  /// spatial dimension
  UInt spatial_dimension;

  /// list of element handled by the material
  ElementTypeMapArray<UInt> element_filter;

  /// damage arrays ordered by element types
  InternalField<Real> damage;

  /// stresses arrays ordered by element types
  InternalField<Real> stress;

  /// eigengrad_u arrays ordered by element types
  InternalField<Real> eigengradu;

  /// grad_u arrays ordered by element types
  InternalField<Real> gradu;

  /// Green Lagrange strain (Finite deformation)
  InternalField<Real> green_strain;

  /// Second Piola-Kirchhoff stress tensor arrays ordered by element types
  /// (Finite deformation)
  InternalField<Real> piola_kirchhoff_2;

  /// dissipated energy by element
  InternalField<Real> dissipated_energy;

  /// tell if using in non local mode or not
  bool is_non_local;

  /// tell if the material need the previous stress state
  bool use_previous_stress;

  /// tell if the material need the previous strain state
  bool use_previous_gradu;

  /// elemental field interpolation coordinates
  InternalField<Real> interpolation_inverse_coordinates;

  /// elemental field interpolation points
  InternalField<Real> interpolation_points_matrices;

  /// vector that contains the names of all the internals that need to
  /// be transferred when material interfaces move
  std::vector<ID> internals_to_transfer;
  
};

} // namespace akantu
  
#endif
