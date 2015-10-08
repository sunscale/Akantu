/**
 * @file   material_non_local.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 31 2011
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  Material class that handle the non locality of a law for example damage.
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
#include "aka_common.hh"
#include "material.hh"
#include "fe_engine.hh"


/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_NON_LOCAL_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt dim>
class MaterialNonLocal : public virtual Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialNonLocal(SolidMechanicsModel & model, const ID & id = "");
  virtual ~MaterialNonLocal();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material computed parameter
  virtual void initMaterial();

  virtual void updateResidual(GhostType ghost_type) {};

  /// insert the quadrature points in the neighborhoods of the non-local manager
  virtual void insertQuadsInNeighborhoods(GhostType ghost_type = _not_ghost);

  /// update the values in the non-local internal fields
  void updateNonLocalInternals(ElementTypeMapReal & non_local_flattened, const ID & field_id, 
			       const UInt nb_component);
  /// constitutive law
  virtual void computeNonLocalStresses(GhostType ghost_type = _not_ghost) = 0;

protected:

  /// associate the non-local variables of the material to their neighborhoods
  virtual void nonLocalVariableToNeighborhood() = 0;
  
  virtual inline UInt getNbDataForElements(const Array<Element> & elements,
					   SynchronizationTag tag) const{ return 0;};

  virtual inline void packElementData(CommunicationBuffer & buffer,
				      const Array<Element> & elements,
				      SynchronizationTag tag) const {};

  virtual inline void unpackElementData(CommunicationBuffer & buffer,
					const Array<Element> & elements,
					SynchronizationTag tag){};

  virtual inline void onElementsAdded(const Array<Element> & element_list) {};
  virtual inline void onElementsRemoved(const Array<Element> & element_list,
					const ElementTypeMapArray<UInt> & new_numbering,
					const RemovedElementsEvent & event) {};


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_non_local_inline_impl.cc"


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_NON_LOCAL_HH__ */
