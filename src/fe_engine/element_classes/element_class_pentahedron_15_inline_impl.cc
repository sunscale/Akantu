/**
 * @file   element_class_pentahedron_15_inline_impl.cc
 *
 * @author Sacha Laffely <sacha.laffely@epfl.ch>
 * @author Damien Scantamburlo <damien.scantamburlo@epfl.ch>
 *
 * @date creation: Wed Mar 19 2015
 *
 * @brief  Specialization of the element_class class for the type _pentahedron_15
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
 * @section DESCRIPTION
 *
 * @verbatim

             /z
             |
             |
	     |  2
	     | /|\
	     |/ | \
	     11 8  7
	    /   |   \
	   /    |    \
	   5    3--9--1
	  | \  /      /
	  |  \12     /
	  14  13    10----------/y
	  | /  \   /
 	  |/    \ /
 	  6--15--4
         /
        /
       /
      \x

       x   y    z
* N1  -1   1    0
* N2  -1   0    1
* N3  -1   0    0
* N4   1   1    0
* N5   1   0    1
* N6   1   0    0
* N7  -1   0.5  0.5
* N8  -1   0    0.5
* N9  -1   0.5  0
* N10  0   1    0
* N11  0   0    1
* N12  0   0    0
* N13  1   0.5  0.5
* N14  1   0    0.5
* N15  1   0.5  0

*/

/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_pentahedron_15,
				     _gt_pentahedron_15,
				     _itp_lagrange_pentahedron_15,
				     _ek_regular,
				     3,
				     _git_pentahedron,
				     2);

AKANTU_DEFINE_SHAPE(_gt_pentahedron_15, _gst_prism);

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline void
InterpolationElement<_itp_lagrange_pentahedron_15>::computeShapes(const vector_type & c,
								 vector_type & N) {
	// Shape Functions, Natural coordinates
	N(0)  =  0.5 * c(1) * (1 - c(0)) * (2 * c(1) - 2 - c(0)); ///N1(q_0)
	N(1)  =  0.5 * c(2) * (1 - c(0)) * (2 * c(2) - 2 - c(0)); ///N2(q_0)
	N(2)  =  0.5 * (c(0) - 1) * (1 - c(1) - c(2)) * (c(0) + 2 * c(1) + 2 * c(2)); ///N3(q_0)
	N(3)  =  0.5 * c(1) * (1 + c(0)) * (2 * c(1) - 2 + c(0)); ///N4(q_0)
	N(4)  =  0.5 * c(2) * (1 + c(0)) * (2 * c(2) - 2 + c(0)); ///N5(q_0)
	N(5)  =  0.5 * (-c(0) - 1) * (1 - c(1) - c(2)) * (-c(0) + 2 * c(1) + 2 * c(2)); ///N6(q_0)
	N(6)  = 2 * c(1) * c(2) * (1 - c(0));                     ///N7(q_0)
	N(7)  = 2 * c(2) * (1 - c(1) - c(2)) * (1 - c(0));        ///N8(q_0)
	N(8)  = 2 * c(1) * (1 - c(0)) * (1 - c(1) - c(2)) * (1 - c(0)); ///N9(q_0)
	N(9)  = c(1) * (1 - std::pow(c(0),2));                    ///N10(q_0)
	N(10) = c(2) * (1 - std::pow(c(0),2));                    ///N11(q_0)
	N(11) = (1 - c(1) - c(2)) * (1 - std::pow(c(0),2));       ///N12(q_0)
	N(12) = 2 * c(1) * c(2) * (1 + c(0));                     ///N13(q_0)
	N(13) = 2 * c(2) * (1 - c(1) - c(2)) * (1 + c(0));        ///N14(q_0)
	N(14) = 2 * c(1) * (1 - c(1) - c(2)) * (1 + c(0));        ///N15(q_0)
}

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void
InterpolationElement<_itp_lagrange_pentahedron_15>::computeDNDS(const vector_type & c,
                                                              matrix_type & dnds) {
	//ddx
	dnds(0,0)  =  0.5 * c(1) * (2 * c(0) - 2 * c(1) + 1);;
	dnds(0,1)  =  0.5 * c(2) * (2 * c(0) - 2 * c(2) + 1);;
	dnds(0,2)  = -0.5 * (2 * c(0) + 2 * c(1) + 2 * c(2) - 1) * (c(1) + c(2) -1);;
	dnds(0,3)  =  0.5 * c(1) * (2 * c(0) + 2 * c(1) - 1);;
	dnds(0,4)  =  0.5 * c(2) * (2 * c(0) + 2 * c(2) - 1);;
	dnds(0,5)  = -0.5 * (c(1) + c(2) - 1) * (2 * c(0) - 2 * c(1) - 2 * c(2) + 1);;
	dnds(0,6)  = -2 * c(1) * c(2);;
	dnds(0,7)  =  2 * c(2) * (c(1) + c(2) - 1);;
	dnds(0,8)  =  2 * c(1) * (c(1) + c(2) - 1);;
	dnds(0,9)  = -2 * c(0) * c(1);;
	dnds(0,10) = -2 * c(0) * c(2);;
	dnds(0,11) =  2 * c(0) * (c(1) + c(2) - 1);;
	dnds(0,12) =  2 * c(1) * c(2);;
	dnds(0,13) = -2 * c(2) * (c(1) + c(2) - 1);;
	dnds(0,14) = -2 * c(1) * (c(1) + c(2) - 1);;

	//ddy
	dnds(1,0)  = -0.5 * (c(0) - 1) * (4 * c(1) - c(0) - 2);;
	dnds(1,1)  =  0;;
	dnds(1,2)  = -0.5 * (c(0) - 1) * (4 * c(1) + c(0) + 2 * (2 * c(2) - 1));;
	dnds(1,3)  =  0.5 * (c(0) + 1) * (4 * c(1) + c(0) - 2);;
	dnds(1,4)  =  0;;
	dnds(1,5)  =  0.5 * (c(0) + 1) * (4 * c(1) - c(0) + 2 * (2 * c(2) - 1));;
	dnds(1,6)  = -2 * (c(0) - 1) * c(2);;
	dnds(1,7)  =  2 * c(2) * (c(0) - 1);;
	dnds(1,8)  =  2 * (2 * c(1) + c(2) - 1) * (c(0) - 1);;
	dnds(1,9)  = -(std::pow(c(0),2) - 1);;
	dnds(1,10) =  0;;
	dnds(1,11) =  (std::pow(c(0),2) - 1);;
	dnds(1,12) =  2 * c(2) * (c(0) + 1);;
	dnds(1,13) = -2 * c(2) * (c(0) + 1);;
	dnds(1,14) = -2 * (2 * c(1) + c(2) - 1) * (c(0) + 1);;

	//ddz
	dnds(2,0)  =  0;;
	dnds(2,1)  = -0.5 * (c(0) - 1) * (4 * c(2) - c(0) - 2);;
	dnds(2,2)  = -0.5 * (c(0) - 1) * (4 * c(2) + c(0) + 2 * (2 * c(1) - 1));;
	dnds(2,3)  =  0;;
	dnds(2,4)  =  0.5 * (c(0) + 1) * (4 * c(2) + c(0) - 2);;
	dnds(2,5)  =  0.5 * (c(0) + 1) * (4 * c(2) - c(0) + 2 * (2 * c(1) - 1));;
	dnds(2,6)  = -2 * (c(0) - 1) * c(1);;
	dnds(2,7)  =  2 * (c(0) - 1) * (2 * c(2) + c(1) - 1);;
	dnds(2,8)  =  2 * c(1) * (c(0) - 1);;
	dnds(2,9)  =  0;;
	dnds(2,10) = -(std::pow(c(0),2) - 1);;
	dnds(2,11) =  (std::pow(c(0),2) - 1);;
	dnds(2,12) =  2 * (c(0) + 1) * c(1);;
	dnds(2,13) = -2 * (c(0) + 1) * (2 * c(2) + c(1) - 1);;
	dnds(2,14) = -2 * (c(0) + 1) * c(1);;
}

//A faire : Diamètre bouboule
