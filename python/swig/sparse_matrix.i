

%include "sparse_matrix.hh"

%pythoncode %{
import scipy.sparse
import numpy as _np    
class AkantuSparseMatrix (scipy.sparse.coo_matrix) :

    def __init__(self,aka_sparse):
        
        self.aka_sparse = aka_sparse
        matrix_type = self.aka_sparse.getSparseMatrixType()
        sz = self.aka_sparse.getSize()
        row = self.aka_sparse.getIRN()[:,0] -1
        col = self.aka_sparse.getJCN()[:,0] -1
        data = self.aka_sparse.getA()[:,0]

        row = row.copy()
        col = col.copy()
        data = data.copy()

        if matrix_type == _symmetric:
            non_diags = (row != col)
            row_sup = col[non_diags]
            col_sup = row[non_diags]
            data_sup = data[non_diags]
            col  = _np.concatenate((col,col_sup))
            row  = _np.concatenate((row,row_sup))
            data = _np.concatenate((data,data_sup))
       
        scipy.sparse.coo_matrix.__init__(self,(data, (row,col)),shape=(sz,sz))

%}
