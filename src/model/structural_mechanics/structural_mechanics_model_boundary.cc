/**
 * @file   structural_mechanics_model_boundary.cc
 *
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Fri Jul 04 2014
 *
 * @brief  Implementation of the boundary conditions for StructuralMechanicsModel
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
#include "model.hh"
#include "structural_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::transferNMatrixToSymVoigtNMatrix<_bernoulli_beam_2>(Array<Real> & N_matrix) {
  AKANTU_DEBUG_IN();
  MyFEEngineType & fem = getFEEngineClass<MyFEEngineType>();
  UInt nb_nodes_per_element = getFEEngine().getMesh().getNbNodesPerElement(_bernoulli_beam_2);

  Array<Real>::const_vector_iterator shape_N0 = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 0).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_M0 = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 1).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_L0 = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 2).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Mp = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 3).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Lp = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 4).begin(nb_nodes_per_element);

  N_matrix.clear();
  Array<Real>::matrix_iterator N_it = N_matrix.begin(nb_degree_of_freedom, nb_degree_of_freedom * nb_nodes_per_element);
  Array<Real>::matrix_iterator N_end = N_matrix.end(nb_degree_of_freedom, nb_degree_of_freedom * nb_nodes_per_element);

  for (;N_it != N_end; ++N_it, ++shape_N0, ++shape_M0, ++shape_L0, ++shape_Mp, ++shape_Lp) {
    Matrix<Real> & N = *N_it;
    const Vector<Real> & N0 = *shape_N0;
    const Vector<Real> & M0 = *shape_M0;
    const Vector<Real> & L0 = *shape_L0;
    const Vector<Real> & Mp = *shape_Mp;
    const Vector<Real> & Lp = *shape_Lp;

    N(0,0) = N0(0);
    N(0,3) = N0(1);

    N(1,1) = M0(0);
    N(1,2) = L0(0);
    N(1,4) = M0(1);
    N(1,5) = L0(1);

    N(2,1) = Mp(0);
    N(2,2) = Lp(0);
    N(2,4) = Mp(1);
    N(2,5) = Lp(1);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::transferNMatrixToSymVoigtNMatrix<_bernoulli_beam_3>(Array<Real> & N_matrix) {
  AKANTU_DEBUG_IN();

  ElementType type = _bernoulli_beam_3;

  MyFEEngineType & fem = getFEEngineClass<MyFEEngineType>();
  UInt nb_nodes_per_element = getFEEngine().getMesh().getNbNodesPerElement(type);

  Array<Real>::const_vector_iterator shape_N0 = fem.getShapeFunctions().getShapes(type, _not_ghost, 0).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_M0 = fem.getShapeFunctions().getShapes(type, _not_ghost, 1).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_L0 = fem.getShapeFunctions().getShapes(type, _not_ghost, 2).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Mp = fem.getShapeFunctions().getShapes(type, _not_ghost, 3).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Lp = fem.getShapeFunctions().getShapes(type, _not_ghost, 4).begin(nb_nodes_per_element);

  N_matrix.clear();
  Array<Real>::matrix_iterator N_it = N_matrix.begin(nb_degree_of_freedom, nb_degree_of_freedom * nb_nodes_per_element);
  Array<Real>::matrix_iterator N_end = N_matrix.end(nb_degree_of_freedom, nb_degree_of_freedom * nb_nodes_per_element);

  for (; N_it != N_end; ++N_it, ++shape_N0, ++shape_M0, ++shape_L0, ++shape_Mp, ++shape_Lp) {
    Matrix<Real> & N = *N_it;
    const Vector<Real> & N0 = *shape_N0;
    const Vector<Real> & M0 = *shape_M0;
    const Vector<Real> & L0 = *shape_L0;
    const Vector<Real> & Mp = *shape_Mp;
    const Vector<Real> & Lp = *shape_Lp;

    N(0,0)  =  N0(0);
    N(0,6)  =  N0(1);

    N(1,1)  =  M0(0);
    N(1,5)  =  L0(0);
    N(1,7)  =  M0(1);
    N(1,11) =  L0(1);

    N(2,2)  =  M0(0);
    N(2,4)  = -L0(0);
    N(2,8)  =  M0(1);
    N(2,10) = -L0(1);

    N(3,3)  =  N0(0);
    N(3,9)  =  N0(1);

    N(4,2)  =  Mp(0);
    N(4,4)  = -Lp(0);
    N(4,8)  =  Mp(1);
    N(4,10) = -Lp(1);

    N(5,1)  =  Mp(0);
    N(5,5)  =  Lp(0);
    N(5,7)  =  Mp(1);
    N(5,11) =  Lp(1);
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */

template<>
void StructuralMechanicsModel::transferNMatrixToSymVoigtNMatrix<_kirchhoff_shell>(Array<Real> & N_matrix) {
  AKANTU_DEBUG_IN();

  ElementType type = _kirchhoff_shell;

  MyFEEngineType & fem = getFEEngineClass<MyFEEngineType>();
  UInt nb_nodes_per_element = getFEEngine().getMesh().getNbNodesPerElement(type);

  Array<Real>::const_vector_iterator shape_N0  = fem.getShapeFunctions().getShapes(type, _not_ghost, 0).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Nw2 = fem.getShapeFunctions().getShapes(type, _not_ghost, 1).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Nw3 = fem.getShapeFunctions().getShapes(type, _not_ghost, 2).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Nx1 = fem.getShapeFunctions().getShapes(type, _not_ghost, 3).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Nx2 = fem.getShapeFunctions().getShapes(type, _not_ghost, 4).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Nx3 = fem.getShapeFunctions().getShapes(type, _not_ghost, 5).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Ny1 = fem.getShapeFunctions().getShapes(type, _not_ghost, 6).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Ny2 = fem.getShapeFunctions().getShapes(type, _not_ghost, 7).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Ny3 = fem.getShapeFunctions().getShapes(type, _not_ghost, 8).begin(nb_nodes_per_element);

  N_matrix.clear();
  Array<Real>::matrix_iterator N_it = N_matrix.begin(nb_degree_of_freedom, nb_degree_of_freedom * nb_nodes_per_element);
  Array<Real>::matrix_iterator N_end = N_matrix.end(nb_degree_of_freedom, nb_degree_of_freedom * nb_nodes_per_element);

  for (; N_it != N_end; ++N_it, ++shape_N0, ++shape_Nw2, ++shape_Nw3, ++shape_Nx1, ++shape_Nx2, ++shape_Nx3, ++shape_Ny1, ++shape_Ny2, ++shape_Ny3) {
    Matrix<Real> & N = *N_it;
    const Vector<Real> & N0 = *shape_N0;
    const Vector<Real> & Nw2 = *shape_Nw2;
    const Vector<Real> & Nw3 = *shape_Nw3;
    const Vector<Real> & Nx1 = *shape_Nx1;
    const Vector<Real> & Nx2 = *shape_Nx2;
    const Vector<Real> & Nx3 = *shape_Nx3;
    const Vector<Real> & Ny1 = *shape_Ny1;
    const Vector<Real> & Ny2 = *shape_Ny2;
    const Vector<Real> & Ny3 = *shape_Ny3;

      N(0,0)  =  N0(0);
      N(0,5)  =  N0(1);  
      N(0,10) =  N0(2);

      N(1,1)  =  N0(0);
      N(1,5)  =  N0(1); 
      N(1,11) =  N0(2);

      N(2,2)  =  N0(0);
      N(2,3)  =  Nw2(0);
      N(2,4)  =  Nw3(0);
      N(2,7)  =  N0(1);
      N(2,8)  =  Nw2(1);
      N(2,9)  =  Nw3(1);
      N(2,12) =  N0(2);
      N(2,13) =  Nw2(2);
      N(2,14) =  Nw3(2);

      N(3,2)  =  Nx1(0);
      N(3,3)  =  Nx2(0);
      N(3,4)  =  Nx3(0);
      N(3,7)  =  Nx1(1);
      N(3,8)  =  Nx2(1);
      N(3,9)  =  Nx3(1);
      N(3,12) =  Nx1(2);
      N(3,13) =  Nx2(2);
      N(3,14) =  Nx3(2);

      N(4,2)  =  Ny1(0);
      N(4,3)  =  Ny2(0);
      N(4,4)  =  Ny3(0);
      N(4,7)  =  Ny1(1);
      N(4,8)  =  Ny2(1);
      N(4,9)  =  Ny3(1);
      N(4,12) =  Ny1(2);
      N(4,13) =  Ny2(2);
      N(4,14) =  Ny3(2);

  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */


__END_AKANTU__
