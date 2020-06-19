/**
 * @file   test_.cc
 *
 * @author Clement Roux-Langlois <clement.roux@epfl.ch>
 *
 * @date creation: Fri Mar 13 2015
 * @date last modification: Tue june 16 2015
 *
 * @brief  Tests the interface mesh generation
 *
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Spherical_kernel_3.h>

/* -------------------------------------------------------------------------- */

int main() {

  // typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > NT;
  typedef CGAL::Quotient<CGAL::MP_Float> NT;
  // typedef double NT;

  typedef CGAL::Spherical_kernel_3<CGAL::Simple_cartesian<NT>,
                                   CGAL::Algebraic_kernel_for_spheres_2_3<NT>>
      Spherical;

  // CGAL::Point_3<Spherical> b(0.3070532930016993, 0.3958272423848058, 0.),
  //  c(0.480062080085834, -0.1988482243531368, 0.);
  double xb = 0.3, yb = 0.4, zb = 0., xc = 0.5, yc = -0.2, zc = 0.;
  CGAL::Point_3<Spherical> b(xb, yb, zb), c(xc, yc, zc);

  CGAL::Line_3<Spherical> l2(b, c);

  CGAL::Line_arc_3<Spherical> s2(l2, b, c);

  Spherical::Sphere_3 sphere(Spherical::Point_3(0.4, 0, 0), 0.3 * 0.3);

  typedef std::pair<Spherical::Circular_arc_point_3, int> pair_type;
  typedef boost::variant<pair_type> sk_inter_res;
  std::list<sk_inter_res> s_results;

  CGAL::intersection(s2, sphere, std::back_inserter(s_results));

  if (s_results.size() == 1) { // just one point
    if (pair_type * pair = boost::get<pair_type>(&s_results.front())) {
      if (pair->second == 1) { // not a point tangent to the sphere
        Spherical::Circular_arc_point_3 arc_point = pair->first;
        std::cout << " x = " << arc_point.x() << ", y = " << arc_point.y()
                  << std::endl;
        double xint, yint;
        xint = to_double(arc_point.x());
        yint = to_double(arc_point.y());
        std::cout << "x = " << xint << ", y = " << yint << std::endl;
      }
    }
  }

  return EXIT_SUCCESS;
}
