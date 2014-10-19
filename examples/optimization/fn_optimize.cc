/**
 * @file   fn_optimize.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date   Thu May 22 14:12:00 2014
 *
 * @brief  File used to show how to use the NLopt optimizator to find the
 * minimum of a function
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include <vector>
#include <math.h>

#include "aka_optimize.hh"


typedef struct {
  double a, b;
} my_constraint_data;


//! Functor used for the evaluation of the function and its gradient
class Functor {
  
  int count_;    //!< Function evaluation counter
  
public:
  
  //! Default constructor
  Functor() : count_() {}
  
  //! Return function evaluation counter
  int count() const
  { return count_; }
  
  double operator()(const std::vector<double> &x, std::vector<double> &grad)
  {
    ++count_;
    if (!grad.empty()) {
      grad[0] = 0.0;
      grad[1] = 0.5 / sqrt(x[1]);
    }
    return sqrt(x[1]);
  }
  
  static double wrap(const std::vector<double> &x, std::vector<double> &grad, void *data) {
    return (*reinterpret_cast<Functor*>(data))(x, grad); }
  
};


double myvconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  my_constraint_data *d = reinterpret_cast<my_constraint_data*>(data);
  double a = d->a, b = d->b;
  if (!grad.empty()) {
    grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
    grad[1] = -1.0;
  }
  return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
}


int main(int argc, char *argv[]) {
  
  my_constraint_data data[2] = { {2,0}, {-1,1} };

  std::vector<double> x(2);
  x[0] = 1.234; x[1] = 5.678;
  
  Functor fn;
  akantu::Optimizator ofn(x, fn);
  
  ofn.add_inequality_constraint(myvconstraint, &data[0], 1e-8);
  ofn.add_inequality_constraint(myvconstraint, &data[1], 1e-8);
  
  ofn.result();
  std::cout<<"\nTotal function evaluations: "<<fn.count()<<std::endl;
 
  return 0;
}
