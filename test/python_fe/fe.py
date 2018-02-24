# ------------------------------------------------------------------------------
__author__ = "Nicolas Richart"
__copyright__ = "Copyright (C) 2016-2018, EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__credits__ = ["Nicolas Richart"]
__license__ = "L-GPLv3"
__maintainer__ = "Nicolas Richart"
__email__ = "nicolas.richart@epfl.ch"
# ------------------------------------------------------------------------------

import abc


class Solver(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def solveStep(self):
        pass

    def _zero_rows(self, A, rows):
        for r in rows:
            diag = A[r, r]
            # zeros the lines defined in b
            A.data[A.indptr[r]:A.indptr[r + 1]] = 0.
            A[r, r] = diag


class Model(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def applyDirichletBC(self):
        pass

    @abc.abstractmethod
    def f_int(self):
        pass

    @property
    @abc.abstractmethod
    def f_ext(self):
        pass

    @property
    @abc.abstractmethod
    def K(self):
        pass

    @property
    @abc.abstractmethod
    def M(self):
        pass

    @property
    @abc.abstractmethod
    def C(self):
        pass

    @property
    @abc.abstractmethod
    def blocked(self):
        pass
