/**
 * @file   model_manager.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Mon Jan 07 2013
 * @date last modification: Tue May 07 2013
 *
 * @brief  higher order object that deals with collections of models
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

// akantu header files
#include "model_manager.hh"

__BEGIN_AKANTU__


#define SWAPD(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)


uint32_t solve_quadratic(Real a, Real b, Real c,
                         Real *x0, Real *x1 )
{
#ifdef DEBUG_POLYSOLVER
  std::cout << "solve_quadratic( "
  << a << ", "
  << b << ", "
  << c << " )\n";
#endif
  
  // Handle linear case
  if( a == 0.0 ) {
    if( b == 0.0 ) {
#ifdef DEBUG_POLYSOLVER
	    std::cout << "linear with no roots\n";
#endif
	    return( 0 );
    } else {
	    *x0 = -c / b;
#ifdef DEBUG_POLYSOLVER
	    std::cout << "linear with one root at x = " << *x0 << "\n";
#endif
	    return( 1 );
    }
  }
  
  Real disc = b*b - 4.0*a*c;
  if( disc > 0 ) {
    if( b == 0 ) {
	    Real r = fabs( 0.5 * sqrt (disc) / a );
	    *x0 = -r;
	    *x1 =  r;
    } else {
	    Real sgnb = (b > 0 ? 1 : -1);
	    Real temp = -0.5 * (b + sgnb * sqrt (disc));
	    Real r1 = temp / a ;
	    Real r2 = c / temp ;
	    
	    if( r1 < r2 ) {
        *x0 = r1 ;
        *x1 = r2 ;
      } else {
        *x0 = r2 ;
        *x1 = r1 ;
      }
    }
    return 2;
  } else if( disc == 0 ) {
    *x0 = -0.5 * b / a ;
    *x1 = -0.5 * b / a ;
    return 2;
  } else
    return 0;
}


/* Finds the real roots of x^4 + a x^3 + b x^2 + c x + d = 0
 */

  
//template<>
uint32_t Kinematic_traits<Consider_t>::solve_quartic(Real a, Real b, Real c, Real d,
                                                     Real *x0, Real *x1, Real *x2, Real *x3) {
  Real u[3];
  Real aa, pp, qq, rr, rc, sc, tc, mt;
  Real w1r, w1i, w2r, w2i, w3r;
  Real v[3], v1, v2, arg, theta;
  Real disc, h;
  int k1 = 0, k2 = 0;
  Real zarr[4];
  
  /* Deal easily with the cases where the quartic is degenerate. The
   * ordering of solutions is done explicitly. */
  if (0 == b && 0 == c)
  {
    if (0 == d)
    {
	    if (a > 0)
      {
        *x0 = -a;
        *x1 = 0.0;
        *x2 = 0.0;
        *x3 = 0.0;
      }
	    else
      {
        *x0 = 0.0;
        *x1 = 0.0;
        *x2 = 0.0;
        *x3 = -a;
      }
	    return 4;
    }
    else if (0 == a)
    {
	    if (d > 0)
      {
        return 0;
      }
	    else
      {
        *x1 = sqrt (sqrt (-d));
        *x0 = -(*x1);
        return 2;
      }
    }
  }
  
  if (0.0 == c && 0.0 == d)
  {
    *x0=0.0;
    *x1=0.0;
    if( solve_quadratic( 1.0, a, b, x2, x3 ) == 0 ) {
	    mt=3;
    } else {
	    mt=1;
    }
  }
  else
  {
    /* For non-degenerate solutions, proceed by constructing and
     * solving the resolvent cubic */
    aa = a * a;
    pp = b - (3.0/8.0) * aa;
    qq = c - (1.0/2.0) * a * (b - (1.0/4.0) * aa);
    rr = d - (1.0/4.0) * (a * c - (1.0/4.0) * aa * (b - (3.0/16.0) * aa));
    rc = (1.0/2.0) * pp;
    sc = (1.0/4.0) * ((1.0/4.0) * pp * pp - rr);
    tc = -((1.0/8.0) * qq * (1.0/8.0) * qq);
    
    /* This code solves the resolvent cubic in a convenient fashion
     * for this implementation of the quartic. If there are three real
     * roots, then they are placed directly into u[].  If two are
     * complex, then the real root is put into u[0] and the real
     * and imaginary part of the complex roots are placed into
     * u[1] and u[2], respectively. Additionally, this
     * calculates the discriminant of the cubic and puts it into the
     * variable disc. */
    {
	    Real qcub = (rc * rc - 3 * sc);
	    Real rcub = (2 * rc * rc * rc - 9 * rc * sc + 27 * tc);
      
	    Real Q = qcub / 9;
	    Real R = rcub / 54;
      
	    Real Q3 = Q * Q * Q;
	    Real R2 = R * R;
      
	    Real CR2 = 729 * rcub * rcub;
	    Real CQ3 = 2916 * qcub * qcub * qcub;
      
	    disc = (CR2 - CQ3) / 2125764.0;
      
	    if (0 == R && 0 == Q)
	    {
        u[0] = -rc / 3;
        u[1] = -rc / 3;
        u[2] = -rc / 3;
	    }
	    else if (CR2 == CQ3)
	    {
        Real sqrtQ = sqrt (Q);
        if (R > 0)
        {
          u[0] = -2 * sqrtQ - rc / 3;
          u[1] = sqrtQ - rc / 3;
          u[2] = sqrtQ - rc / 3;
        }
        else
        {
          u[0] = -sqrtQ - rc / 3;
          u[1] = -sqrtQ - rc / 3;
          u[2] = 2 * sqrtQ - rc / 3;
        }
	    }
	    else if (CR2 < CQ3)
	    {
        Real sqrtQ = sqrt (Q);
        Real sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
        Real theta = acos (R / sqrtQ3);
        if (R / sqrtQ3 >= 1.0) theta = 0.0;
        {
          Real norm = -2 * sqrtQ;
          
          u[0] = norm * cos (theta / 3) - rc / 3;
          u[1] = norm * cos ((theta + 2.0 * M_PI) / 3) - rc / 3;
          u[2] = norm * cos ((theta - 2.0 * M_PI) / 3) - rc / 3;
        }
	    }
	    else
	    {
        Real sgnR = (R >= 0 ? 1 : -1);
        Real modR = fabs (R);
        Real sqrt_disc = sqrt (R2 - Q3);
        Real A = -sgnR * pow (modR + sqrt_disc, 1.0 / 3.0);
        Real B = Q / A;
        Real mod_diffAB = fabs (A - B);
        
        u[0] = A + B - rc / 3;
        u[1] = -0.5 * (A + B) - rc / 3;
        u[2] = -(sqrt (3.0) / 2.0) * mod_diffAB;
	    }
    }
    /* End of solution to resolvent cubic */
    
    /* Combine the square roots of the roots of the cubic
     * resolvent appropriately. Also, calculate 'mt' which
     * designates the nature of the roots:
     * mt=1 : 4 real roots (disc == 0)
     * mt=2 : 0 real roots (disc < 0)
     * mt=3 : 2 real roots (disc > 0)
     */
    
    if (0.0 == disc)
	    u[2] = u[1];
    
    if (0 >= disc)
    {
	    mt = 2;
      
	    /* One would think that we could return 0 here and exit,
	     * since mt=2. However, this assignment is temporary and
	     * changes to mt=1 under certain conditions below.
	     */
      
	    v[0] = fabs (u[0]);
	    v[1] = fabs (u[1]);
	    v[2] = fabs (u[2]);
      
	    v1 = std::max(std::max(v[0], v[1] ), v[2] );
	    /* Work out which two roots have the largest moduli */
	    k1 = 0, k2 = 0;
	    if (v1 == v[0])
	    {
        k1 = 0;
        v2 = std::max( v[1], v[2] );
	    }
	    else if (v1 == v[1])
	    {
        k1 = 1;
        v2 = std::max( v[0], v[2] );
	    }
	    else
	    {
        k1 = 2;
        v2 = std::max( v[0], v[1] );
	    }
      
	    if (v2 == v[0])
	    {
        k2 = 0;
	    }
	    else if (v2 == v[1])
	    {
        k2 = 1;
	    }
	    else
	    {
        k2 = 2;
	    }
      
	    if (0.0 <= u[k1])
	    {
        w1r=sqrt(u[k1]);
        w1i=0.0;
	    }
	    else
	    {
        w1r=0.0;
        w1i=sqrt(-u[k1]);
	    }
	    if (0.0 <= u[k2])
	    {
        w2r=sqrt(u[k2]);
        w2i=0.0;
	    }
	    else
	    {
        w2r=0.0;
        w2i=sqrt(-u[k2]);
	    }
    }
    else
    {
	    mt = 3;
      
	    if (0.0 == u[1] && 0.0 == u[2])
	    {
        arg = 0.0;
	    }
	    else
	    {
        arg = sqrt(sqrt(u[1] * u[1] + u[2] * u[2]));
	    }
	    theta = atan2(u[2], u[1]);
      
	    w1r = arg * cos(theta / 2.0);
	    w1i = arg * sin(theta / 2.0);
	    w2r = w1r;
	    w2i = -w1i;
    }
    
    /* Solve the quadratic to obtain the roots to the quartic */
    w3r = qq / 8.0 * (w1i * w2i - w1r * w2r) /
    (w1i * w1i + w1r * w1r) / (w2i * w2i + w2r * w2r);
    h = a / 4.0;
    
    zarr[0] = w1r + w2r + w3r - h;
    zarr[1] = -w1r - w2r + w3r - h;
    zarr[2] = -w1r + w2r - w3r - h;
    zarr[3] = w1r - w2r - w3r - h;
    
    /* Arrange the roots into the variables z0, z1, z2, z3 */
    if (2 == mt)
    {
	    if (u[k1] >= 0 && u[k2] >= 0)
      {
        mt = 1;
        *x0 = zarr[0];
        *x1 = zarr[1];
        *x2 = zarr[2];
        *x3 = zarr[3];
      }
	    else
	    {
        return 0;
	    }
    }
    else
    {
	    *x0 = zarr[0];
	    *x1 = zarr[1];
    }
  }
  
  /* Sort the roots as usual */
  if (1 == mt)
  {
    /* Roots are all real, sort them by the real part */
    if (*x0 > *x1)
	    SWAPD (*x0, *x1);
    if (*x0 > *x2)
	    SWAPD (*x0, *x2);
    if (*x0 > *x3)
	    SWAPD (*x0, *x3);
    
    if (*x1 > *x2)
	    SWAPD (*x1, *x2);
    if (*x2 > *x3) {
	    SWAPD (*x2, *x3);
	    if (*x1 > *x2)
        SWAPD (*x1, *x2);
    }
    return 4;
  }
  else
  {
    /* 2 real roots */
    if (*x0 > *x1)
	    SWAPD (*x0, *x1);
  }
  
  return 2;
}


__END_AKANTU__

