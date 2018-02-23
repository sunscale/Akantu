import copy
import numpy.linalg as npl
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from . import export
from . import fe


@export
class DynamicSolver(fe.Solver):
    def __init__(self, model, **kwargs):
        opt = copy.copy(kwargs)
        self._delta_t = opt.pop("delta_t", 0.001)
        self._alpha = opt.pop("alpha", 1./2.)
        self._beta = opt.pop("beta", 1./2.)
        self._type = opt.pop("type", 'disp')
        self._tolerance = opt.pop("tolerance", 1e-10)
        self._max_nloops = opt.pop("max_iterations", 100)

        self._model = model
        self._J = sp.csr_matrix(self._model.K.shape)
        self._coeff = {'disp': {'disp': 1.,
                                'velo': 1. / (self._alpha * self._delta_t),
                                'acce': 1. / (self._alpha * self._beta * self._delta_t ** 2)},  # NOQA: E501
                       'velo': {'disp': self._alpha * self._delta_t,
                                'velo': 1.,
                                'acce': 1. / (self._beta * self._delta_t)},
                       'acce': {'disp': self._alpha * self._beta * self._delta_t ** 2,  # NOQA: E501
                                'velo': self._beta * self._delta_t,
                                'acce': 1.}}

    def _assembleResidual(self):
        self._r = self._model.f_ext - self._model.f_int - \
                  self._model.M * self._model.a
        C = self._model.C
        if C is not None:
            self._r -= C * self._model.v

    def _predictor(self):
        self._model.u += self._delta_t * self._model.v + \
                         self._delta_t ** 2 / 2. * self._model.a
        self._model.v += self._delta_t * self._model.a

    def _corrector(self, delta_):
        self._model.u += self._coeff[self._type]['disp'] * delta_
        self._model.v += self._coeff[self._type]['velo'] * delta_
        self._model.a += self._coeff[self._type]['acce'] * delta_

    def _assembleJacobian(self):
        K = self._model.K
        e = self._coeff[self._type]['disp']
        self._J = e * K

        C = self._model.C
        if C is not None:
            d = self._coeff[self._type]['velo']
            self._J += d * C

        M = self._model.M
        if M is not None:
            c = self._coeff[self._type]['acce']
            self._J += c * M

        self._model.applyDirichletBC()

        self._zero_rows(self._J, self._model.blocked)

    def solveStep(self):
        self._predictor()
        self._nloops = 0
        converged = False
        while not converged:
            self._assembleJacobian()
            self._assembleResidual()
            delta_ = spl.spsolve(self._J, self._r)
            self._corrector(delta_)

            self._nloops += 1
            error = npl.norm(delta_)

            converged = error < self._tolerance or \
                self._nloops > self._max_nloops

        print("{0} {1} -> {2}".format(error, self._nloops, converged))
        if self._nloops >= self._max_nloops:
            raise ValueError('The solver did not converge')

    @property
    def nloops(self):
        return self._nloops
