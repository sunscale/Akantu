import numpy as np
import copy
import scipy.sparse as sp
from . import export
from . import fe


@export
class TrussFE(fe.Model):
    def __init__(self, **kwargs):
        opt = copy.copy(kwargs)
        self._A = opt.pop("A", 1.)
        self._E = opt.pop("E", 1.)
        self._L = opt.pop("L", 1.)
        self._rho = opt.pop("rho", 1.)

        self._F = opt.pop("F", {-10: [-1]})
        self._blocked = opt.pop("blocked", ([0], [0]))

        self._Ne = opt.pop("Ne", 1)
        self._Le = self._L / self._Ne

        self._F_ext = np.zeros((self._Ne + 1))
        self._F_int = np.zeros((self._Ne + 1))
        self._u = np.zeros((self._Ne + 1))
        self._v = np.zeros((self._Ne + 1))
        self._a = np.zeros((self._Ne + 1))

        self.assembleMass()

    def assembleMass(self):
        Me = np.array([[2., 1.],
                       [1., 2.]])
        Me *= self._rho * self._A * self._Le / 6
        self._M = sp.lil_matrix((self._Ne + 1, self._Ne + 1))
        for i in range(self._Ne):
            self._M[i:i+2, i:i+2] += Me

        self._M = self._M.tocsr()

    def assembleStiffness(self):
        Ke = np.array([[1., -1.],
                       [-1., 1.]])
        Ke *= self._E * self._A / self._Le
        self._K = sp.lil_matrix((self._Ne + 1, self._Ne + 1))
        for i in range(self._Ne):
            self._K[i:i+2, i:i+2] += Ke

        self._K = self._K.tocsr()

    def applyNeumannBC(self):
        for force, nodes in self._F.items():
            self._F_ext[nodes] = force

    def applyDirichletBC(self):
        for i, b in enumerate(self._blocked[0]):
            self._u[b] = self._blocked[1][i]

    def computeInternalForces(self):
        self._F_int = self._K * self._u

    @property
    def K(self):
        self.assembleStiffness()
        return self._K

    @property
    def M(self):
        return self._M

    @property
    def C(self):
        return None

    @property
    def u(self):
        self.applyDirichletBC()
        return self._u

    @u.setter
    def u(self, x):
        self._u = x
        self.applyDirichletBC()

    @property
    def v(self):
        return self._v

    @v.setter
    def v(self, v_):
        self._v = v_

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, a_):
        self._a = a_

    @property
    def f_int(self):
        self.computeInternalForces()
        return self._F_int

    @property
    def f_ext(self):
        self.applyNeumannBC()
        return self._F_ext

    @f_ext.setter
    def f_ext(self, f):
        self._F_ext = f

    @property
    def blocked(self):
        return self._blocked[0]

    @property
    def Le(self):
        return self._Le
