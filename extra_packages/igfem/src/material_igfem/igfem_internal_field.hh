/**
 * @file   igfem_internal_field.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  IGFEM internal field
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#include "internal_field.hh"

#ifndef __AKANTU_IGFEM_INTERNAL_FIELD_HH__
#define __AKANTU_IGFEM_INTERNAL_FIELD_HH__

__BEGIN_AKANTU__

template <typename T> class IGFEMInternalField : public InternalField<T> {
public:
  IGFEMInternalField(const ID & id, Material & material);
  virtual ~IGFEMInternalField();
};

__END_AKANTU__

#endif /* __AKANTU_IGFEM_INTERNAL_FIELD_HH__ */
