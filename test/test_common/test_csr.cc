/**
 * @file   test_csr.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jul 30 2012
 * @date last modification: Thu Dec 06 2012
 *
 * @brief  Test the CSR (compressed sparse row) data structure
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include <iostream>


#include "aka_common.hh"
#include "aka_csr.hh"

using namespace akantu;

#define N 1000

int main(int argc, char *argv[]) {

  CSR<UInt> csr;

  std::vector<UInt> nb_cols_per_row;

  csr.resizeRows(N);
  csr.clearRows();

  for (UInt i = 0; i < N; ++i) {
    UInt nb_cols(UInt(rand()*double(N)/(RAND_MAX+1.)));
    nb_cols_per_row.push_back(nb_cols);
    for (UInt j = 0; j < nb_cols; ++j) {
      ++csr.rowOffset(i);
    }
  }

  csr.countToCSR();
  csr.resizeCols();

  csr.beginInsertions();
  for (UInt i = 0; i < N; ++i) {
    UInt nb_cols = nb_cols_per_row[i];
    for (UInt j = 0; j < nb_cols; ++j) {
      csr.insertInRow(i, nb_cols - j);
    }
  }
  csr.endInsertions();

  if(csr.getNbRows() != N) {
    AKANTU_DEBUG_ERROR("The number of rows does not correspond: "
		       << csr.getNbRows()
		       << " != "
		       << N);
  }

  for (UInt i = 0; i < csr.getNbRows(); ++i) {
    CSR<UInt>::iterator it = csr.begin(i);
    CSR<UInt>::iterator end = csr.end(i);
    UInt nb_cols = nb_cols_per_row[i];
    for (; it != end; ++it) {
      if(nb_cols != *it) {
	AKANTU_DEBUG_ERROR("The numbers stored in the row " << i << " are not correct: "
			   << nb_cols 
			   << " != "
			   << *it);
      }
      nb_cols--;
    }

    if(nb_cols != 0) {
      AKANTU_DEBUG_ERROR("Not enough columns in the row " << i << ": "
			 << nb_cols);
    }
  }

  for (UInt i = 0; i < csr.getNbRows(); ++i) {
    CSR<UInt>::iterator it = csr.rbegin(i);
    CSR<UInt>::iterator end = csr.rend(i);

    UInt nb_cols = nb_cols_per_row[i];
    UInt j = nb_cols;

    for (; it != end; --it) {
      if((nb_cols - j + 1) != *it) {
	AKANTU_DEBUG_ERROR("Reverse: The numbers stored in the row " << i << " are not correct: "
			   << (nb_cols - j + 1)
			   << " != "
			   << *it);
      }
      j--;
    }

    if(j != 0) AKANTU_DEBUG_ERROR("Reverse: Not enough columns in the row " << i << ": "
				  << j);
  }

 
  return EXIT_SUCCESS;
}
