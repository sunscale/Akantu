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
#include "aka_grid_dynamic.hh"
#include "fe_engine.hh"

#include "weight_function.hh"

namespace akantu {
  class GridSynchronizer;
}

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_NON_LOCAL_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt dim, template <UInt> class WeightFunction = BaseWeightFunction>
class MaterialNonLocal : public virtual Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialNonLocal(SolidMechanicsModel & model, const ID & id = "");
  virtual ~MaterialNonLocal();

  template<typename T>
  class PairList : public ElementTypeMap< ElementTypeMapArray<T> > {};
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material computed parameter
  virtual void initMaterial();

  virtual void updateResidual(GhostType ghost_type);

  virtual void computeAllNonLocalStresses(GhostType ghost_type = _not_ghost);

  void savePairs(const std::string & filename) const;
  void neighbourhoodStatistics(const std::string & filename) const;

protected:
  void updatePairList(const ElementTypeMapArray<Real> & quadrature_points_coordinates);

  void computeWeights(const ElementTypeMapArray<Real> & quadrature_points_coordinates);

  void createCellList(ElementTypeMapArray<Real> & quadrature_points_coordinates);

  void cleanupExtraGhostElement(const ElementTypeMap<UInt> & nb_ghost_protected);

  void fillCellList(const ElementTypeMapArray<Real> & quadrature_points_coordinates,
		    const GhostType & ghost_type);

  /// constitutive law
  virtual void computeNonLocalStresses(GhostType ghost_type = _not_ghost) = 0;

  template<typename T>
  void weightedAvergageOnNeighbours(const ElementTypeMapArray<T> & to_accumulate,
				    ElementTypeMapArray<T> & accumulated,
				    UInt nb_degree_of_freedom,
				    GhostType ghost_type2 = _not_ghost) const;

  virtual inline UInt getNbDataForElements(const Array<Element> & elements,
					  SynchronizationTag tag) const;

  virtual inline void packElementData(CommunicationBuffer & buffer,
				      const Array<Element> & elements,
				      SynchronizationTag tag) const;

  virtual inline void unpackElementData(CommunicationBuffer & buffer,
					const Array<Element> & elements,
					SynchronizationTag tag);

  //  virtual inline void onElementsAdded(const Array<Element> & element_list);
  virtual inline void onElementsRemoved(const Array<Element> & element_list,
					const ElementTypeMapArray<UInt> & new_numbering,
					const RemovedElementsEvent & event);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  void registerNonLocalVariable(InternalField<Real> & local,
				InternalField<Real> & non_local,
				UInt nb_degree_of_freedom) {
    ID id = local.getID();
    NonLocalVariable & non_local_variable = non_local_variables[id];

    non_local_variable.local = &local;
    non_local_variable.non_local = &non_local;
    non_local_variable.nb_component = nb_degree_of_freedom;
  }

  AKANTU_GET_MACRO(PairList, pair_list, const PairList<UInt> &)
  Real getRadius() const { return weight_func->getRadius(); }
  AKANTU_GET_MACRO(CellList, *spatial_grid, const SpatialGrid<QuadraturePoint> &)

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// the weight function used
  WeightFunction<dim> * weight_func;

private:
  /// the pairs of quadrature points
  PairList<UInt> pair_list;
  /// the weights associated to the pairs
  PairList<Real> pair_weight;

  /// the regular grid to construct/update the pair lists
  SpatialGrid<QuadraturePoint> * spatial_grid;

  /// the types of the existing pairs
  typedef std::set< std::pair<ElementType, ElementType> > pair_type;
  pair_type existing_pairs[2];

  /// count the number of calls of computeStress
  UInt compute_stress_calls;

  struct NonLocalVariable {
    InternalField<Real> * local;
    InternalField<Real> * non_local;
    UInt nb_component;
  };

  std::map<ID, NonLocalVariable> non_local_variables;

  bool is_creating_grid;

  GridSynchronizer * grid_synchronizer;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_non_local_inline_impl.cc"


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_NON_LOCAL_HH__ */
