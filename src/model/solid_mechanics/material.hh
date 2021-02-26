/**
 * @file   material.hh
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Mother class for all materials
 *
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
#include "solid_mechanics_model_event_handler.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_HH_
#define AKANTU_MATERIAL_HH_

/* -------------------------------------------------------------------------- */
namespace akantu {
class Model;
class SolidMechanicsModel;
class Material;
} // namespace akantu

namespace akantu {

using MaterialFactory =
    Factory<Material, ID, UInt, const ID &, SolidMechanicsModel &, const ID &>;

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
 *  virtual void computeTangentStiffness(ElementType el_type,
 *                                       Array<Real> & tangent_matrix,
 *                                       GhostType ghost_type = _not_ghost);
 * \endcode
 *
 */
class Material : public Memory,
                 public DataAccessor<Element>,
                 public Parsable,
                 public MeshEventHandler,
                 protected SolidMechanicsModelEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  Material(const Material & mat) = delete;
  Material & operator=(const Material & mat) = delete;

  /// Initialize material with defaults
  Material(SolidMechanicsModel & model, const ID & id = "");

  /// Initialize material with custom mesh & fe_engine
  Material(SolidMechanicsModel & model, UInt dim, const Mesh & mesh,
           FEEngine & fe_engine, const ID & id = "");

  /// Destructor
  ~Material() override;

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Function that materials can/should reimplement                           */
  /* ------------------------------------------------------------------------ */
protected:
  /// constitutive law
  virtual void computeStress(ElementType /* el_type */,
                             GhostType /* ghost_type */ = _not_ghost) {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the tangent stiffness matrix
  virtual void computeTangentModuli(ElementType /*el_type*/,
                                    Array<Real> & /*tangent_matrix*/,
                                    GhostType /*ghost_type*/ = _not_ghost) {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the potential energy
  virtual void computePotentialEnergy(ElementType el_type);

  /// compute the potential energy for an element
  virtual void
  computePotentialEnergyByElement(ElementType /*type*/, UInt /*index*/,
                                  Vector<Real> & /*epot_on_quad_points*/) {
    AKANTU_TO_IMPLEMENT();
  }

  virtual void updateEnergies(ElementType /*el_type*/) {}

  virtual void updateEnergiesAfterDamage(ElementType /*el_type*/) {}

  /// set the material to steady state (to be implemented for materials that
  /// need it)
  virtual void setToSteadyState(ElementType /*el_type*/,
                                GhostType /*ghost_type*/ = _not_ghost) {}

  /// function called to update the internal parameters when the modifiable
  /// parameters are modified
  virtual void updateInternalParameters() {}

public:
  /// extrapolate internal values
  virtual void extrapolateInternal(const ID & id, const Element & element,
                                   const Matrix<Real> & points,
                                   Matrix<Real> & extrapolated);

  /// compute the p-wave speed in the material
  virtual Real getPushWaveSpeed(const Element & /*element*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the s-wave speed in the material
  virtual Real getShearWaveSpeed(const Element & /*element*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// get a material celerity to compute the stable time step (default: is the
  /// push wave speed)
  virtual Real getCelerity(const Element & element) const {
    return getPushWaveSpeed(element);
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  template <typename T> void registerInternal(InternalField<T> & /*vect*/) {
    AKANTU_TO_IMPLEMENT();
  }

  template <typename T> void unregisterInternal(InternalField<T> & /*vect*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// compute the residual for this material
  //  virtual void updateResidual(GhostType ghost_type = _not_ghost);

  /// assemble the residual for this material
  virtual void assembleInternalForces(GhostType ghost_type);

  /// save the stress in the previous_stress if needed
  virtual void savePreviousState();

  /// restore the stress from previous_stress if needed
  virtual void restorePreviousState();

  /// compute the stresses for this material
  virtual void computeAllStresses(GhostType ghost_type = _not_ghost);
  // virtual void
  // computeAllStressesFromTangentModuli(GhostType ghost_type = _not_ghost);
  virtual void computeAllCauchyStresses(GhostType ghost_type = _not_ghost);

  /// set material to steady state
  void setToSteadyState(GhostType ghost_type = _not_ghost);

  /// compute the stiffness matrix
  virtual void assembleStiffnessMatrix(GhostType ghost_type);

  /// add an element to the local mesh filter
  inline UInt addElement(ElementType type, UInt element, GhostType ghost_type);
  inline UInt addElement(const Element & element);

  /// add many elements at once
  void addElements(const Array<Element> & elements_to_add);

  /// remove many element at once
  void removeElements(const Array<Element> & elements_to_remove);

  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /**
   * interpolate stress on given positions for each element by means
   * of a geometrical interpolation on quadrature points
   */
  void interpolateStress(ElementTypeMapArray<Real> & result,
                         GhostType ghost_type = _not_ghost);

  /**
   * interpolate stress on given positions for each element by means
   * of a geometrical interpolation on quadrature points and store the
   * results per facet
   */
  void interpolateStressOnFacets(ElementTypeMapArray<Real> & result,
                                 ElementTypeMapArray<Real> & by_elem_result,
                                 GhostType ghost_type = _not_ghost);

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
  static inline UInt getTangentStiffnessVoigtSize(UInt dim);

  /// compute the potential energy by element
  void computePotentialEnergyByElements();

  /// resize the intenals arrays
  virtual void resizeInternals();

  /* ------------------------------------------------------------------------ */
  /* Finite deformation functions                                             */
  /* This functions area implementing what is described in the paper of Bathe */
  /* et al, in IJNME, Finite Element Formulations For Large Deformation       */
  /* Dynamic Analysis, Vol 9, 353-386, 1975                                   */
  /* ------------------------------------------------------------------------ */
protected:
  /// assemble the residual
  template <UInt dim> void assembleInternalForces(GhostType ghost_type);

  template <UInt dim>
  void computeAllStressesFromTangentModuli(ElementType type,
                                           GhostType ghost_type);

  template <UInt dim>
  void assembleStiffnessMatrix(ElementType type, GhostType ghost_type);

  /// assembling in finite deformation
  template <UInt dim>
  void assembleStiffnessMatrixNL(ElementType type, GhostType ghost_type);

  template <UInt dim>
  void assembleStiffnessMatrixL2(ElementType type, GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Conversion functions                                                     */
  /* ------------------------------------------------------------------------ */
public:
  /// Size of the Stress matrix for the case of finite deformation see: Bathe et
  /// al, IJNME, Vol 9, 353-386, 1975
  static inline UInt getCauchyStressMatrixSize(UInt dim);

  /// Sets the stress matrix according to Bathe et al, IJNME, Vol 9, 353-386,
  /// 1975
  template <UInt dim>
  static inline void setCauchyStressMatrix(const Matrix<Real> & S_t,
                                           Matrix<Real> & sigma);

  /// write the stress tensor in the Voigt notation.
  template <UInt dim>
  static inline decltype(auto) stressToVoigt(const Matrix<Real> & stress) {
    return VoigtHelper<dim>::matrixToVoigt(stress);
  }

  /// write the strain tensor in the Voigt notation.
  template <UInt dim>
  static inline decltype(auto) strainToVoigt(const Matrix<Real> & strain) {
    return VoigtHelper<dim>::matrixToVoigtWithFactors(strain);
  }

  /// write a voigt vector to stress
  template <UInt dim>
  static inline void voigtToStress(const Vector<Real> & voigt,
                                   Matrix<Real> & stress) {
    VoigtHelper<dim>::voigtToMatrix(voigt, stress);
  }

  /// Computation of Cauchy stress tensor in the case of finite deformation from
  /// the 2nd Piola-Kirchhoff for a given element type
  template <UInt dim>
  void StoCauchy(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// Computation the Cauchy stress the 2nd Piola-Kirchhoff and the deformation
  /// gradient
  template <UInt dim>
  inline void StoCauchy(const Matrix<Real> & F, const Matrix<Real> & S,
                        Matrix<Real> & sigma, const Real & C33 = 1.0) const;

  template <UInt dim>
  static inline void gradUToF(const Matrix<Real> & grad_u, Matrix<Real> & F);

  template <UInt dim>
  static inline decltype(auto) gradUToF(const Matrix<Real> & grad_u);

  static inline void rightCauchy(const Matrix<Real> & F, Matrix<Real> & C);
  static inline void leftCauchy(const Matrix<Real> & F, Matrix<Real> & B);

  template <UInt dim>
  static inline void gradUToEpsilon(const Matrix<Real> & grad_u,
                                    Matrix<Real> & epsilon);
  template <UInt dim>
  static inline decltype(auto) gradUToEpsilon(const Matrix<Real> & grad_u);

  template <UInt dim>
  static inline void gradUToE(const Matrix<Real> & grad_u,
                              Matrix<Real> & epsilon);

  template <UInt dim>
  static inline decltype(auto) gradUToE(const Matrix<Real> & grad_u);

  static inline Real stressToVonMises(const Matrix<Real> & stress);

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
  void onNodesAdded(const Array<UInt> & /*unused*/,
                    const NewNodesEvent & /*unused*/) override{};
  void onNodesRemoved(const Array<UInt> & /*unused*/,
                      const Array<UInt> & /*unused*/,
                      const RemovedNodesEvent & /*unused*/) override{};
  void onElementsAdded(const Array<Element> & element_list,
                       const NewElementsEvent & event) override;
  void onElementsRemoved(const Array<Element> & element_list,
                         const ElementTypeMapArray<UInt> & new_numbering,
                         const RemovedElementsEvent & event) override;
  void onElementsChanged(const Array<Element> & /*unused*/,
                         const Array<Element> & /*unused*/,
                         const ElementTypeMapArray<UInt> & /*unused*/,
                         const ChangedElementsEvent & /*unused*/) override{};

  /* ------------------------------------------------------------------------ */
  /* SolidMechanicsModelEventHandler inherited members                        */
  /* ------------------------------------------------------------------------ */
public:
  virtual void beforeSolveStep();
  virtual void afterSolveStep(bool converged = true);

  void onDamageIteration() override;
  void onDamageUpdate() override;
  void onDump() override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Name, name, const std::string &);

  AKANTU_GET_MACRO(Model, model, const SolidMechanicsModel &)

  AKANTU_GET_MACRO(ID, Memory::getID(), const ID &);
  AKANTU_GET_MACRO(Rho, rho, Real);
  AKANTU_SET_MACRO(Rho, rho, Real);

  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /// return the potential energy for the subset of elements contained by the
  /// material
  Real getPotentialEnergy();
  /// return the potential energy for the provided element
  Real getPotentialEnergy(ElementType & type, UInt index);

  /// return the energy (identified by id) for the subset of elements contained
  /// by the material
  virtual Real getEnergy(const std::string & type);
  /// return the energy (identified by id) for the provided element
  virtual Real getEnergy(const std::string & energy_id, ElementType type,
                         UInt index);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementFilter, element_filter, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(GradU, gradu, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Stress, stress, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PotentialEnergy, potential_energy,
                                         Real);
  AKANTU_GET_MACRO(GradU, gradu, const ElementTypeMapArray<Real> &);
  AKANTU_GET_MACRO(Stress, stress, const ElementTypeMapArray<Real> &);
  AKANTU_GET_MACRO(ElementFilter, element_filter,
                   const ElementTypeMapArray<UInt> &);
  AKANTU_GET_MACRO(FEEngine, fem, FEEngine &);

  bool isNonLocal() const { return is_non_local; }

  template <typename T>
  const Array<T> & getArray(const ID & id, ElementType type,
                            GhostType ghost_type = _not_ghost) const;
  template <typename T>
  Array<T> & getArray(const ID & id, ElementType type,
                      GhostType ghost_type = _not_ghost);

  template <typename T>
  const InternalField<T> & getInternal(const ID & id) const;
  template <typename T> InternalField<T> & getInternal(const ID & id);

  template <typename T>
  inline bool isInternal(const ID & id, ElementKind element_kind) const;

  template <typename T>
  ElementTypeMap<UInt> getInternalDataPerElem(const ID & id,
                                              ElementKind element_kind) const;

  bool isFiniteDeformation() const { return finite_deformation; }
  bool isInelasticDeformation() const { return inelastic_deformation; }

  template <typename T> inline void setParam(const ID & param, T value);
  inline const Parameter & getParam(const ID & param) const;

  template <typename T>
  void flattenInternal(const std::string & field_id,
                       ElementTypeMapArray<T> & internal_flat,
                       GhostType ghost_type = _not_ghost,
                       ElementKind element_kind = _ek_not_defined) const;

  /// apply a constant eigengrad_u everywhere in the material
  virtual void applyEigenGradU(const Matrix<Real> & prescribed_eigen_grad_u,
                               GhostType /*ghost_type*/ = _not_ghost);

  bool hasMatrixChanged(const ID & id) {
    if (id == "K") {
      return hasStiffnessMatrixChanged() or finite_deformation;
    }

    return true;
  }

  MatrixType getMatrixType(const ID & id) {
    if (id == "K") {
      return getTangentType();
    }

    if (id == "M") {
      return _symmetric;
    }

    return _mt_not_defined;
  }

  /// specify if the matrix need to be recomputed for this material
  virtual bool hasStiffnessMatrixChanged() { return true; }

  /// specify the type of matrix, if not overloaded the material is not valid
  /// for static or implicit computations
  virtual MatrixType getTangentType() { return _mt_not_defined; }

  /// static method to reteive the material factory
  static MaterialFactory & getFactory();

protected:
  bool isInit() const { return is_init; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// boolean to know if the material has been initialized
  bool is_init;

  std::map<ID, InternalField<Real> *> internal_vectors_real;
  std::map<ID, InternalField<UInt> *> internal_vectors_uint;
  std::map<ID, InternalField<bool> *> internal_vectors_bool;

protected:
  /// Link to the fem object in the model
  FEEngine & fem;

  /// Finite deformation
  bool finite_deformation;

  /// Finite deformation
  bool inelastic_deformation;

  /// material name
  std::string name;

  /// The model to witch the material belong
  SolidMechanicsModel & model;

  /// density : rho
  Real rho;

  /// spatial dimension
  UInt spatial_dimension;

  /// list of element handled by the material
  ElementTypeMapArray<UInt> element_filter;

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

  /// potential energy by element
  InternalField<Real> potential_energy;

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
private:
  /// eigen_grad_u for the parser
  Matrix<Real> eigen_grad_u;
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const Material & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "material_inline_impl.hh"

#include "internal_field_tmpl.hh"
#include "random_internal_field_tmpl.hh"

/* -------------------------------------------------------------------------- */
/* Auto loop                                                                  */
/* -------------------------------------------------------------------------- */
/// This can be used to automatically write the loop on quadrature points in
/// functions such as computeStress. This macro in addition to write the loop
/// provides two tensors (matrices) sigma and grad_u
#define MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type)       \
  auto && grad_u_view =                                                        \
      make_view(this->gradu(el_type, ghost_type), this->spatial_dimension,     \
                this->spatial_dimension);                                      \
                                                                               \
  auto stress_view =                                                           \
      make_view(this->stress(el_type, ghost_type), this->spatial_dimension,    \
                this->spatial_dimension);                                      \
                                                                               \
  if (this->isFiniteDeformation()) {                                           \
    stress_view = make_view(this->piola_kirchhoff_2(el_type, ghost_type),      \
                            this->spatial_dimension, this->spatial_dimension); \
  }                                                                            \
                                                                               \
  for (auto && data : zip(grad_u_view, stress_view)) {                         \
    [[gnu::unused]] Matrix<Real> & grad_u = std::get<0>(data);                 \
    [[gnu::unused]] Matrix<Real> & sigma = std::get<1>(data)

#define MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END }

/// This can be used to automatically write the loop on quadrature points in
/// functions such as computeTangentModuli. This macro in addition to write the
/// loop provides two tensors (matrices) sigma_tensor, grad_u, and a matrix
/// where the elemental tangent moduli should be stored in Voigt Notation
#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_mat)              \
  auto && grad_u_view =                                                        \
      make_view(this->gradu(el_type, ghost_type), this->spatial_dimension,     \
                this->spatial_dimension);                                      \
                                                                               \
  auto && stress_view =                                                        \
      make_view(this->stress(el_type, ghost_type), this->spatial_dimension,    \
                this->spatial_dimension);                                      \
                                                                               \
  auto tangent_size =                                                          \
      this->getTangentStiffnessVoigtSize(this->spatial_dimension);             \
                                                                               \
  auto && tangent_view = make_view(tangent_mat, tangent_size, tangent_size);   \
                                                                               \
  for (auto && data : zip(grad_u_view, stress_view, tangent_view)) {           \
    [[gnu::unused]] Matrix<Real> & grad_u = std::get<0>(data);                 \
    [[gnu::unused]] Matrix<Real> & sigma = std::get<1>(data);                  \
    Matrix<Real> & tangent = std::get<2>(data);

#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END }

/* -------------------------------------------------------------------------- */

#define INSTANTIATE_MATERIAL_ONLY(mat_name)                                    \
  template class mat_name<1>; /* NOLINT */                                     \
  template class mat_name<2>; /* NOLINT */                                     \
  template class mat_name<3>  /* NOLINT */

#define MATERIAL_DEFAULT_PER_DIM_ALLOCATOR(id, mat_name)                       \
  [](UInt dim, const ID &, SolidMechanicsModel & model,                        \
     const ID & id) /* NOLINT */                                               \
      -> std::unique_ptr<                                                      \
          Material> { /* NOLINT */                                             \
                      switch (dim) {                                           \
                      case 1:                                                  \
                        return std::make_unique<mat_name<1>>(/* NOLINT */      \
                                                             model, id);       \
                      case 2:                                                  \
                        return std::make_unique<mat_name<2>>(/* NOLINT */      \
                                                             model, id);       \
                      case 3:                                                  \
                        return std::make_unique<mat_name<3>>(/* NOLINT */      \
                                                             model, id);       \
                      default:                                                 \
                        AKANTU_EXCEPTION(                                      \
                            "The dimension "                                   \
                            << dim                                             \
                            << "is not a valid dimension for the material "    \
                            << #id);                                           \
                      }                                                        \
  }

#define INSTANTIATE_MATERIAL(id, mat_name)                                     \
  INSTANTIATE_MATERIAL_ONLY(mat_name);                                         \
  static bool material_is_alocated_##id [[gnu::unused]] =                      \
      MaterialFactory::getInstance().registerAllocator(                        \
          #id, MATERIAL_DEFAULT_PER_DIM_ALLOCATOR(id, mat_name))

#endif /* AKANTU_MATERIAL_HH_ */
