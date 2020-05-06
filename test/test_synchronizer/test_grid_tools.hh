/**
 * @file   test_grid_tools.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Mar 11 2013
 * @date last modification: Wed Feb 03 2016
 *
 * @brief  Tools to help for the akantu::Grid class tests
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

#include <map>

#include "aka_common.hh"
#include "aka_types.hh"

#define TOLERANCE 1e-7

template <UInt dim> class Point {
public:
  Point() : id(0), tol(TOLERANCE) {
    for (UInt i = 0; i < dim; ++i)
      pos[i] = 0.;
  }

  Point(const Point & pt) : id(pt.id), tol(pt.tol) {
    for (UInt i = 0; i < dim; ++i)
      pos[i] = pt.pos[i];
  }

  Point(const Vector<Real> & pt, UInt id = 0) : id(id), tol(TOLERANCE) {
    for (UInt i = 0; i < dim; ++i)
      pos[i] = pt(i);
  }

  bool operator==(const Point & pt) const {
    for (UInt i = 0; i < dim; ++i) {
      //      std::cout << i << " " << pos[i] << " " << pt.pos[i] << " " <<
      //      std::abs(pos[i] - pt.pos[i]);
      if (std::abs(pos[i] - pt.pos[i]) > tol) {
        //        std::cout << " " << false << std::endl;
        return false;
      } // else
      //        std::cout << " " << true << std::endl;
    }
    return true;
  }

  bool operator<(const Point & pt) const {
    UInt i = 0, j = 0;
    for (; (i < dim) && (j < dim); i++, j++) {
      if (pos[i] < pt.pos[j])
        return true;
      if (pt.pos[j] < pos[i])
        return false;
    }
    return (i == dim) && (j != dim);
  }

  bool operator!=(const Point & pt) const { return !(operator==(pt)); }

  Real & operator()(UInt d) { return pos[d]; }
  const Real & operator()(UInt d) const { return pos[d]; }

  void read(const std::string & str) {
    std::stringstream sstr(str);
    for (UInt i = 0; i < dim; ++i)
      sstr >> pos[i];
  }

  void write(std::ostream & ostr) const {
    for (UInt i = 0; i < dim; ++i) {
      if (i != 0)
        ostr << " ";
      //    ostr << std::setprecision(std::numeric_limits<Real>::digits) <<
      //    pos[i];
      ostr << std::setprecision(9) << pos[i];
    }
  }

private:
  UInt id;
  Real pos[dim];
  double tol;
};

template <UInt dim> struct neighbors_map_t {
  typedef std::map<Point<dim>, std::vector<Point<dim>>> type;
};

template <UInt dim>
inline std::ostream & operator<<(std::ostream & stream,
                                 const Point<dim> & _this) {
  _this.write(stream);
  return stream;
}
