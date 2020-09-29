/**
 * @file   igfem_internal_field_tmpl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Implementation of IGFEM internal field
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#ifndef AKANTU_IGFEM_INTERNAL_FIELD_TMPL_HH_
#define AKANTU_IGFEM_INTERNAL_FIELD_TMPL_HH_

namespace akantu {
/* -------------------------------------------------------------------------- */
template <typename T>
IGFEMInternalField<T>::IGFEMInternalField(const ID & id, Material & material)
    : InternalField<T>(
          id, material, material.getModel().getFEEngine("IGFEMFEEngine"),
          dynamic_cast<MaterialIGFEM &>(material).getElementFilter()) {
  this->element_kind = _ek_igfem;
}

/* -------------------------------------------------------------------------- */
template <typename T> IGFEMInternalField<T>::~IGFEMInternalField(){};

} // namespace akantu

#endif /* AKANTU_IGFEM_INTERNAL_FIELD_TMPL_HH_ */
