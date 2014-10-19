/**
 * @file   material_random_internal.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  Class describing a material parameter that can be set depending on a
 * random distribution
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_type_map.hh"
#include "aka_random_generator.hh"


#ifndef __AKANTU_MATERIAL_RANDOM_INTERNAL_HH__
#define __AKANTU_MATERIAL_RANDOM_INTERNAL_HH__

__BEGIN_AKANTU__

template<typename Distribution, typename Generator = Rand48Generator<Real> >
class MaterialRandomInternal : public ElementTypeMapArray<Real> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialRandomInternal(const ID & id, const ID & parent_id,
			 const MemoryID & memory_id = 0) :
    ElementTypeMapArray<Real>(id, parent_id, memory_id) {};

  virtual ~MaterialRandomInternal();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  Distribution distribution;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
//#include "material_random_internal_inline_impl.cc"


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_RANDOM_INTERNAL_HH__ */
