/**
 * @file   element_class_kirchhoff_shell_inline_impl.cc
 *
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Fri Jul 04 2014
 * @date last modification: Fri Jul 04 2014
 *
 * @brief  Element class Kirchhoff Shell
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
AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(_kirchhoff_shell,
						_gt_triangle_3,   
						_itp_kirchhoff_shell,
						_triangle_3,
						_ek_structural,  
						3,		      
						_git_triangle, 2);

/* -------------------------------------------------------------------------- */
//cf. element_class_bernoulli...
//element_class_triangle_3_inline... (pour compute Jacobian)
template <>
inline void
InterpolationElement<_itp_kirchhoff_shell>::computeShapes(const Vector<Real> & natural_coords,
							 Vector<Real> & N,
							 const Matrix<Real> & projected_coord,							 
							  UInt id) {
  // projected_coord (x1 x2 x3) (y1 y2 y3)
  
  // natural coordinate
  Real xi = natural_coords(0);
  Real eta = natural_coords(1);

  Real x21=projected_coord(0,0)-projected_coord(0,1); //x1-x2 
  Real x32=projected_coord(0,1)-projected_coord(0,2);
  Real x13=projected_coord(0,2)-projected_coord(0,0);

  Real y21=projected_coord(1,0)-projected_coord(1,1); //y1-y2 
  Real y32=projected_coord(1,1)-projected_coord(1,2);
  Real y13=projected_coord(1,2)-projected_coord(1,0);

  /* Real x21=projected_coord(0,1)-projected_coord(0,0);
  Real x32=projected_coord(0,2)-projected_coord(0,1);
  Real x13=projected_coord(0,0)-projected_coord(0,2);

  Real y21=projected_coord(1,1)-projected_coord(1,0);
  Real y32=projected_coord(1,2)-projected_coord(1,1);
  Real y13=projected_coord(1,0)-projected_coord(1,2);*/

  // natural triangle side length
  Real L4=sqrt(x21*x21+y21*y21);
  Real L5=sqrt(x32*x32+y32*y32);
  Real L6=sqrt(x13*x13+y13*y13);

  // sinus and cosinus
  Real C4=x21/L4; //1
  Real C5=x32/L5; //-1/sqrt(2);
  Real C6=x13/L6; //0
  Real S4=y21/L4; //0;
  Real S5=y32/L5; //1/sqrt(2);
  Real S6=y13/L6; //-1;

  Real N1 = 1-xi-eta;
  Real N2 = xi;
  Real N3 = eta;

  Real P4 = 4*xi*(1-xi-eta);
  Real P5 = 4*xi*eta;
  Real P6 = 4*eta*(1-xi-eta);

  switch (id) {
  case 0: { // N
    N(0) = N1;
    N(1) = N2;
    N(2) = N3;
    break;
  }
  case 1: { // Nwi2
    N(0) = -(1/8)*P4*L4*C4+(1/8)*P6*L6*C6;
    N(1) = -(1/8)*P5*L5*C5+(1/8)*P4*L4*C4; 
    N(2) = -(1/8)*P6*L6*C6+(1/8)*P5*L5*C5;
    break;
  }
  case 2: { // Nwi3
    N(0) = -(1/8)*P4*L4*S4+(1/8)*P6*L6*S6;
    N(1) = -(1/8)*P5*L5*S5+(1/8)*P4*L4*S4; 
    N(2) = -(1/8)*P6*L6*S6+(1/8)*P5*L5*S5;
    break;
 }
  case 3: { // Nxi1
    N(0) = 3/(2*L4)*P4*C4-3/(2*L6)*P6*C6;
    N(1) = 3/(2*L5)*P5*C5-3/(2*L4)*P4*C4;
    N(2) = 3/(2*L6)*P6*C6-3/(2*L5)*P5*C5;
    break;
  }
  case 4: { // Nxi2
    N(0) = N1-(3./4.)*P4*C4*C4-(3./4.)*P6*C6*C6;
    N(1) = N2-(3./4.)*P5*C5*C5-(3./4.)*P4*C4*C4;
    N(2) = N3-(3./4.)*P6*C6*C6-(3./4.)*P5*C5*C5;
    break;
  }
  case 5: { // Nxi3
    N(0) = -(3./4.)*P4*C4*S4-(3./4.)*P6*C6*S6;
    N(1) = -(3./4.)*P5*C5*S5-(3./4.)*P4*C4*S4;
    N(2) = -(3./4.)*P6*C6*S6-(3./4.)*P5*C5*S5;
    break;
  }
 case 6: { // Nyi1
    N(0) = 3/(2*L4)*P4*S4-3/(2*L6)*P6*S6;
    N(1) = 3/(2*L5)*P5*S5-3/(2*L4)*P4*S4;
    N(2) = 3/(2*L6)*P6*S6-3/(2*L5)*P5*S5;
    break;
  }
 case 7: { // Nyi2
    N(0) = -(3./4.)*P4*C4*S4-(3./4.)*P6*C6*S6;
    N(1) = -(3./4.)*P5*C5*S5-(3./4.)*P4*C4*S4;
    N(2) = -(3./4.)*P6*C6*S6-(3./4.)*P5*C5*S5;
    break;
  }
 case 8: { // Nyi3
    N(0) = N1-(3./4.)*P4*S4*S4-(3./4.)*P6*S6*S6;
    N(1) = N2-(3./4.)*P5*S5*S5-(3./4.)*P4*S4*S4;
    N(2) = N3-(3./4.)*P6*S6*S6-(3./4.)*P5*S5*S5;
    break;
  }
  }
}
 
/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_kirchhoff_shell>::computeDNDS(const Vector<Real> & natural_coords,
						       Matrix<Real> & dnds,
						       const Matrix<Real> & projected_coord,
						       UInt id) {

  
  // natural coordinate
  Real xi  = natural_coords(0);
  Real eta = natural_coords(1);

   // projected_coord (x1 x2 x3) (y1 y2 y3)

// donne juste pour pour le patch test 4_5_5 mais donne quelque changement de signe dans la matrice de rotation

  Real x21=projected_coord(0,0)-projected_coord(0,1); //x1-x2 
  Real x32=projected_coord(0,1)-projected_coord(0,2);
  Real x13=projected_coord(0,2)-projected_coord(0,0);

  Real y21=projected_coord(1,0)-projected_coord(1,1); //y1-y2 
  Real y32=projected_coord(1,1)-projected_coord(1,2);
  Real y13=projected_coord(1,2)-projected_coord(1,0);

  // donne juste pour la matrice de rigidité... mais pas pour le patch test 4_5_5

/* Real x21=projected_coord(0,1)-projected_coord(0,0);
  Real x32=projected_coord(0,2)-projected_coord(0,1);
  Real x13=projected_coord(0,0)-projected_coord(0,2);

  Real y21=projected_coord(1,1)-projected_coord(1,0);
  Real y32=projected_coord(1,2)-projected_coord(1,1);
  Real y13=projected_coord(1,0)-projected_coord(1,2);*/ 

  // natural triangle side length
  Real L4=sqrt(x21*x21+y21*y21);
  Real L5=sqrt(x32*x32+y32*y32);
  Real L6=sqrt(x13*x13+y13*y13);

  // sinus and cosinus
  Real C4=x21/L4; 
  Real C5=x32/L5; 
  Real C6=x13/L6; 
  Real S4=y21/L4; 
  Real S5=y32/L5; 
  Real S6=y13/L6; 
 
 
  Real dN1xi = -1;
  Real dN2xi = 1;
  Real dN3xi = 0;

  Real dN1eta = -1;
  Real dN2eta = 0;
  Real dN3eta = 1;

  Real dP4xi = 4-8*xi-4*eta;
  Real dP5xi = 4*eta;
  Real dP6xi = -4*eta;

  Real dP4eta = -4*xi;
  Real dP5eta = 4*xi;
  Real dP6eta = 4-4*xi-8*eta;

  switch (id) {
  case 0: { // N'xi  N'eta
   dnds(0,0) = dN1xi;
   dnds(0,1) = dN2xi;
   dnds(0,2) = dN3xi;

   dnds(1,0) = dN1eta;
   dnds(1,1) = dN2eta;
   dnds(1,2) = dN3eta;
   break;
  }
  case 1: { // Nxi1'xi    Nxi1'eta
   dnds(0,0) = 3/(2*L4)*dP4xi*C4-3/(2*L6)*dP6xi*C6; 
   dnds(0,1) = 3/(2*L5)*dP5xi*C5-3/(2*L4)*dP4xi*C4; 
   dnds(0,2) = 3/(2*L6)*dP6xi*C6-3/(2*L5)*dP5xi*C5; 

   dnds(1,0) = 3/(2*L4)*dP4eta*C4-3/(2*L6)*dP6eta*C6; 
   dnds(1,1) = 3/(2*L5)*dP5eta*C5-3/(2*L4)*dP4eta*C4;
   dnds(1,2) = 3/(2*L6)*dP6eta*C6-3/(2*L5)*dP5eta*C5;
    break;
  }
  case 2: { // Nxi2'xi    Nxi2'eta
    dnds(0,0) = -1-(3./4.)*dP4xi*C4*C4-(3./4.)*dP6xi*C6*C6;
    dnds(0,1) = 1-(3./4.)*dP5xi*C5*C5-(3./4.)*dP4xi*C4*C4; 
    dnds(0,2) =-(3./4.)*dP6xi*C6*C6-(3./4.)*dP5xi*C5*C5; 

    dnds(1,0) = -1-(3./4.)*dP4eta*C4*C4-(3./4.)*dP6eta*C6*C6;
    dnds(1,1) = -(3./4.)*dP5eta*C5*C5-(3./4.)*dP4eta*C4*C4;
    dnds(1,2) = 1-(3./4.)*dP6eta*C6*C6-(3./4.)*dP5eta*C5*C5;
    break;
  }
  case 3: { // Nxi3'xi    Nxi3'eta
    dnds(0,0) = -(3./4.)*dP4xi*C4*S4-(3./4.)*dP6xi*C6*S6; 
    dnds(0,1) = -(3./4.)*dP5xi*C5*S5-(3./4.)*dP4xi*C4*S4; 
    dnds(0,2) = -(3./4.)*dP6xi*C6*S6-(3./4.)*dP5xi*C5*S5; 

    dnds(1,0) = -(3./4.)*dP4eta*C4*S4-(3./4.)*dP6eta*C6*S6;
    dnds(1,1) = -(3./4.)*dP5eta*C5*S5-(3./4.)*dP4eta*C4*S4;
    dnds(1,2) = -(3./4.)*dP6eta*C6*S6-(3./4.)*dP5eta*C5*S5;
    break;
 }
  case 4: { // Nyi1'xi    Nyi1'eta
    dnds(0,0) = 3/(2*L4)*dP4xi*S4-3/(2*L6)*dP6xi*S6;
    dnds(0,1) = 3/(2*L5)*dP5xi*S5-3/(2*L4)*dP4xi*S4;
    dnds(0,2) = 3/(2*L6)*dP6xi*S6-3/(2*L5)*dP5xi*S5;

    dnds(1,0) = 3/(2*L4)*dP4eta*S4-3/(2*L6)*dP6eta*S6; 
    dnds(1,1) = 3/(2*L5)*dP5eta*S5-3/(2*L4)*dP4eta*S4; 
    dnds(1,2) = 3/(2*L6)*dP6eta*S6-3/(2*L5)*dP5eta*S5; 
    break;
  }
 case 5: { // Nyi2'xi     Nyi2'eta
    dnds(0,0) = -(3./4.)*dP4xi*C4*S4-(3./4.)*dP6xi*C6*S6;
    dnds(0,1) = -(3./4.)*dP5xi*C5*S5-(3./4.)*dP4xi*C4*S4;
    dnds(0,2) = -(3./4.)*dP6xi*C6*S6-(3./4.)*dP5xi*C5*S5;

    dnds(1,0) = -(3./4.)*dP4eta*C4*S4-(3./4.)*dP6eta*C6*S6; 
    dnds(1,1) = -(3./4.)*dP5eta*C5*S5-(3./4.)*dP4eta*C4*S4; 
    dnds(1,2) = -(3./4.)*dP6eta*C6*S6-(3./4.)*dP5eta*C5*S5; 
    break;
  }
 case 6: { // Nyi3'xi     Nyi3'eta
    dnds(0,0) = dN1xi-(3./4.)*dP4xi*S4*S4-(3./4.)*dP6xi*S6*S6;
    dnds(0,1) = dN2xi-(3./4.)*dP5xi*S5*S5-(3./4.)*dP4xi*S4*S4;
    dnds(0,2) = dN3xi-(3./4.)*dP6xi*S6*S6-(3./4.)*dP5xi*S5*S5;

    dnds(1,0) = dN1eta-(3./4.)*dP4eta*S4*S4-(3./4.)*dP6eta*S6*S6; 
    dnds(1,1) = dN2eta-(3./4.)*dP5eta*S5*S5-(3./4.)*dP4eta*S4*S4; 
    dnds(1,2) = dN3eta-(3./4.)*dP6eta*S6*S6-(3./4.)*dP5eta*S5*S5; 
    break;
  }
  }
}
