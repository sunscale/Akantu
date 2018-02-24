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

import copy
import scipy.sparse.linalg as spl
from . import export
from . import fe


@export
class StaticSolver(fe.Solver):
    def __init__(self, model):
        self._model = model

    def _assembleResidual(self):
        self._r = self._model.f_ext - self._model.f_int

    def _assembleJacobian(self):
        self._J = copy.copy(self._model.K)
        self._model.applyDirichletBC()

        self._zero_rows(self._J, self._model.blocked)

    def solveStep(self):
        self._assembleJacobian()
        self._assembleResidual()
        x = spl.spsolve(self._J, self._r)
        self._model.u += x
