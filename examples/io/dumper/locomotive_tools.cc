/**
 * @file   locomotive_tools.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Mon Aug 17 2015
 * @date last modification: Mon Jan 18 2016
 *
 * @brief  Common functions for the dumper examples
 *
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "aka_array.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include "locomotive_tools.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
void applyRotation(const Vector<Real> & center, Real angle,
                   const Array<Real> & nodes, Array<Real> & displacement,
                   const Array<UInt> & node_group) {
  auto nodes_it = nodes.begin(nodes.getNbComponent());
  auto disp_it = displacement.begin(center.size());
  Array<UInt>::const_scalar_iterator node_num_it = node_group.begin();
  Array<UInt>::const_scalar_iterator node_num_end = node_group.end();

  Vector<Real> pos_rel(center.size());

  for (; node_num_it != node_num_end; ++node_num_it) {
    const Vector<Real> pos = nodes_it[*node_num_it];
    for (UInt i = 0; i < pos.size(); ++i)
      pos_rel(i) = pos(i);

    Vector<Real> dis = disp_it[*node_num_it];

    pos_rel -= center;
    Real radius = pos_rel.norm();

    if (std::abs(radius) < Math::getTolerance())
      continue;

    Real phi_i = std::acos(pos_rel(_x) / radius);
    if (pos_rel(_y) < 0)
      phi_i *= -1;

    dis(_x) = std::cos(phi_i - angle) * radius - pos_rel(_x);
    dis(_y) = std::sin(phi_i - angle) * radius - pos_rel(_y);
  }
}

/* -------------------------------------------------------------------------- */
void fillColour(const Mesh & mesh, ElementTypeMapArray<UInt> & colour) {
  const ElementTypeMapArray<std::string> & phys_data =
      mesh.getData<std::string>("physical_names");
  const Array<std::string> & txt_colour = phys_data(_triangle_3);
  Array<UInt> & id_colour = colour(_triangle_3);

  for (UInt i = 0; i < txt_colour.size(); ++i) {
    std::string phy_name = txt_colour(i);

    if (phy_name == "red")
      id_colour(i) = 3;
    else if (phy_name == "white" || phy_name == "lwheel_1" ||
             phy_name == "rwheel_1")
      id_colour(i) = 2;
    else
      id_colour(i) = 1;
  }
}
