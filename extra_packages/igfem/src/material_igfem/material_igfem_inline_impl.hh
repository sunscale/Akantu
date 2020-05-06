/**
 * @file   material_igfem_inline_impl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief Implementation of the inline functions of the parent
 * material for IGFEM
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

} // namespace akantu

#include "igfem_helper.hh"
#include "solid_mechanics_model_igfem.hh"
#include <iostream>

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt MaterialIGFEM::getNbDataForElements(const Array<Element> & elements,
                                                SynchronizationTag tag) const {
  if (tag == SynchronizationTag::_smm_stress) {
    return (this->isFiniteDeformation() ? 3 : 1) * spatial_dimension *
           spatial_dimension * sizeof(Real) *
           this->getModel().getNbIntegrationPoints(elements, "IGFEMFEEngine");
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
inline void MaterialIGFEM::packElementData(CommunicationBuffer & buffer,
                                           const Array<Element> & elements,
                                           SynchronizationTag tag) const {
  if (tag == SynchronizationTag::_smm_stress) {
    if (this->isFiniteDeformation()) {
      packElementDataHelper(piola_kirchhoff_2, buffer, elements,
                            "IGFEMFEEngine");
      packElementDataHelper(gradu, buffer, elements, "IGFEMFEEngine");
    }
    packElementDataHelper(stress, buffer, elements, "IGFEMFEEngine");
  }
}

/* -------------------------------------------------------------------------- */
inline void MaterialIGFEM::unpackElementData(CommunicationBuffer & buffer,
                                             const Array<Element> & elements,
                                             SynchronizationTag tag) {
  if (tag == SynchronizationTag::_smm_stress) {
    if (this->isFiniteDeformation()) {
      unpackElementDataHelper(piola_kirchhoff_2, buffer, elements,
                              "IGFEMFEEngine");
      unpackElementDataHelper(gradu, buffer, elements, "IGFEMFEEngine");
    }
    unpackElementDataHelper(stress, buffer, elements, "IGFEMFEEngine");
  }
}
