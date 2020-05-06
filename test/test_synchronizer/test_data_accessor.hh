/**
 * @file   test_data_accessor.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Thu Apr 11 2013
 * @date last modification: Fri Jan 26 2018
 *
 * @brief  Data Accessor class for testing
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "data_accessor.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

using namespace akantu;
/* -------------------------------------------------------------------------- */

class TestAccessor : public DataAccessor<Element> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  inline TestAccessor(const Mesh & mesh,
                      const ElementTypeMapArray<Real> & barycenters);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Barycenter, barycenters, Real);

  /* ------------------------------------------------------------------------ */
  /* Ghost Synchronizer inherited members                                     */
  /* ------------------------------------------------------------------------ */
protected:
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const;
  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const;
  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  const ElementTypeMapArray<Real> & barycenters;
  const Mesh & mesh;
};

/* -------------------------------------------------------------------------- */
/* TestSynchronizer implementation                                            */
/* -------------------------------------------------------------------------- */
inline TestAccessor::TestAccessor(const Mesh & mesh,
                                  const ElementTypeMapArray<Real> & barycenters)
    : barycenters(barycenters), mesh(mesh) {}

inline UInt TestAccessor::getNbData(const Array<Element> & elements,
                                    const SynchronizationTag &) const {
  if (elements.size())
    // return Mesh::getSpatialDimension(elements(0).type) * sizeof(Real) *
    // elements.size();
    return mesh.getSpatialDimension() * sizeof(Real) * elements.size();
  else
    return 0;
}

inline void TestAccessor::packData(CommunicationBuffer & buffer,
                                   const Array<Element> & elements,
                                   const SynchronizationTag &) const {
  UInt spatial_dimension = mesh.getSpatialDimension();
  Array<Element>::const_iterator<Element> bit = elements.begin();
  Array<Element>::const_iterator<Element> bend = elements.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;

    Vector<Real> bary(
        this->barycenters(element.type, element.ghost_type).storage() +
            element.element * spatial_dimension,
        spatial_dimension);
    buffer << bary;
  }
}

inline void TestAccessor::unpackData(CommunicationBuffer & buffer,
                                     const Array<Element> & elements,
                                     const SynchronizationTag &) {
  UInt spatial_dimension = mesh.getSpatialDimension();

  for (const auto & element : elements) {
    Vector<Real> barycenter_loc(
        this->barycenters(element.type, element.ghost_type).storage() +
            element.element * spatial_dimension,
        spatial_dimension);

    Vector<Real> bary(spatial_dimension);
    buffer >> bary;

    auto dist = (barycenter_loc - bary).template norm<L_inf>();
    EXPECT_NEAR(0, dist, 1e-15);
  }
}
