/**
 * @file   solid_mechanics_model_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  solid mechanics model for IGFEM analysis
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
// #ifndef __AKANTU_SOLID_MECHANICS_MODEL_IGFEM_HH__
// #define __AKANTU_SOLID_MECHANICS_MODEL_IGFEM_HH__

// __BEGIN_AKANTU__

// /* -------------------------------------------------------------------------- */
// /* Solid Mechanics Model for IGFEM analysis                                   */
// /* -------------------------------------------------------------------------- */
// class SolidMechanicsModelIGFEM : public SolidMechanicsModel,
// 				 public SolidMechanicsModelEventHandler{
// public:
//   /* ------------------------------------------------------------------------ */
//   /* Constructors/Destructors                                                 */
//   /* ------------------------------------------------------------------------ */
//   typedef FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_igfem> MyFEEngineIGFEMType;

//   SolidMechanicsModelIGFEM(Mesh & mesh,
// 			   UInt spatial_dimension = _all_dimensions,
// 			   const ID & id = "solid_mechanics_model_igfem",
// 			   const MemoryID & memory_id = 0);

//   virtual ~SolidMechanicsModelIGFEM();

//   /* ------------------------------------------------------------------------ */
//   /* Methods                                                                  */
//   /* ------------------------------------------------------------------------ */
// public:
//   computeRealArray();

//   /* ------------------------------------------------------------------------ */
//   /* Class Members                                                            */
//   /* ------------------------------------------------------------------------ */
// private:

//   /// real displacements array
//   Array <Real> *real_displacement;
//   /// real forces array
//   Array <Real> *real_force;
//   /// real residuals array
//   Array <Real> *real_residual;

// };

// /* -------------------------------------------------------------------------- */
// /* IGFEMMaterialSelector                                                      */
// /* -------------------------------------------------------------------------- */

// class DefaultIGFEMMaterialSelector : public DefaultMaterialSelector {
// public:
//   DefaultMaterialIGFEMMaterialSelector(const SolidMechanicsModelCohesive & model) :
//     DefaultMaterialSelector(model.getMaterialByElement()),
//     igfem_material(model.getIGFEMMaterial()),
//     mesh(model.getMesh()) { }

//   virtual UInt operator()(const Element & element) {
//     if(Mesh::getKind(element.type) == _ek_igfem) 
//       return MaterialSelector::operator(element);
//     else
//       return DefaultMaterialSelector::operator(element);
//   }

// private:
//   const 
// };
// __END_AKANTU__

// #endif /* __AKANTU_SOLID_MECHANICS_MODEL_IGFEM_HH__ */


