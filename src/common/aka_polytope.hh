/**
 * @file   aka_polytope.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Tue Jun 17 2014
 *
 * @brief  k-DOP class template
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

#ifndef __AKANTU_POLYTOPE_HH__
#define __AKANTU_POLYTOPE_HH__

#include "aka_point.hh"


__BEGIN_AKANTU__



//! k-DOP (Discrete Orientation Polytope) class template
/*! \tparam k - Number of planes used to define the polytope
 */
template <int k>
class Polytope {
  
  
	Real min_[k / 2];  //! Minimum coordinate value along k/2 directions
	Real max_[k / 2];  //! Maximum coordinate value along k/2 directions

public:
  
  
	//! Parameter constructor that uses a container of points
	template <class point_container>
	Polytope(const point_container& pts) {
		Real inf = std::numeric_limits <Real>::infinity();
		for (size_t i = 0; i < k / 2; ++i) {
			min_[i] = inf;
			max_[i] = -inf;
		}
		initialize(pts);
	}

  //! Check for collision with another polytope
	bool operator& (Polytope & p) {
		for (size_t i = 0; i < k / 2; ++i)
			if (p.min_[i] > max_[i] || p.max_[i] < min_[i])
				return false;
		return true;
	}

private:
  
  //! Private member function used to initialize the polytope
	template <class point_container>
	void initialize(const point_container& pts);
};


template <>
template <class point_container>
void Polytope <14>::initialize(const point_container& pts) {
	typedef typename point_container::value_type point_type;

	Real c = 0.57735026918963;


	for (size_t i = 0; i < pts.size(); ++i) {
		const point_type& p = pts[i];

		// set min/max values along axe (±1,0,0), (0,±1,0), (0,0,±1)
		for (size_t j = 0; j < point_type::dim(); ++j) {
			min_[j] = std::min(min_[j], p[j]);
			max_[j] = std::max(max_[j], p[j]);
		}

		// set min/max values along axe (1,1,1) and (-1,-1,-1) (normalized)
		point_type e1(c, c, c);

		Real f1 = p * e1;
		min_[3] = std::min(min_[3], f1);
		max_[3] = std::max(max_[3], f1);

		// set min/max values along axe (1,1,-1) and (-1,-1,1) (normalized)
		point_type e2(c, c, -c);

		Real f2 = p * e2;
		min_[4] = std::min(min_[4], f2);
		max_[4] = std::max(max_[4], f2);

		// set min/max values along axe (1,-1,1) and (-1,1,-1) (normalized)
		point_type e3(c, -c, c);

		Real f3 = p * e3;
		min_[5] = std::min(min_[5], f3);
		max_[5] = std::max(max_[5], f3);

		// set min/max values along axe (-1,1,1) and (1,-1,-1) (normalized)
		point_type e4(-c, c, c);

		Real f4 = p * e4;
		min_[6] = std::min(min_[6], f4);
		max_[6] = std::max(max_[6], f4);
	}
}

__END_AKANTU__

#endif /* __AKANTU_POLYTOPE_HH__ */
