#!/usr/bin/env python3

__all__ = ['Shapes']

import numpy as np
import numpy.polynomial.polynomial as poly
import aka_test

class Shapes(object):
    NATURAL_COORDS = {
        (1, 'quadrangle'): np.array([[-1.], [1.], [0.]]),
        (2, 'quadrangle'): np.array([[-1., -1.], [ 1., -1.], [ 1.,  1.], [-1.,  1.],
                                     [ 0., -1.], [ 1.,  0.], [ 0.,  1.], [-1.,  0.]]),
        (3, 'quadrangle'): np.array([[-1., -1., -1.], [ 1., -1., -1.], [ 1.,  1., -1.], [-1.,  1., -1.],
                                     [-1., -1.,  1.], [ 1., -1.,  1.], [ 1.,  1.,  1.], [-1.,  1.,  1.],
                                     [ 0., -1., -1.], [ 1.,  0., -1.], [ 0.,  1., -1.], [-1.,  0., -1.],
                                     [-1., -1.,  0.], [ 1., -1.,  0.], [ 1.,  1.,  0.], [-1.,  1.,  0.],
                                     [ 0., -1.,  1.], [ 1.,  0.,  1.], [ 0.,  1.,  1.], [-1.,  0.,  1.]]),
        (2, 'triangle'):   np.array([[0., 0.], [1., 0.], [0., 1,], [.5, 0.], [.5, .5], [0., .5]]),
        (3, 'triangle'):   np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.],
                                     [.5, 0., 0.], [.5, .5, 0.], [0., .5, 0.],
                                     [0., 0., .5], [.5, 0., .5], [0., .5, .5]]),
    }

    QUADRATURE_W = {
        (1, 'quadrangle', 1): np.array([2.]),
        (1, 'quadrangle', 2): np.array([1., 1.]),
        (2, 'triangle', 1): np.array([1./2.]),
        (2, 'triangle', 2): np.array([1., 1., 1.])/6.,
        (3, 'triangle', 1): np.array([]),
        (3, 'triangle', 2): np.array([]),
        (2, 'quadrangle', 1): np.array([1., 1., 1., 1.]),
        (2, 'quadrangle', 2): np.array([]),
        (3, 'quadrangle', 1): np.array([1., 1., 1., 1.,
                                        1., 1., 1., 1.]),
        (3, 'quadrangle', 2): np.array([]),

    }
    QUADRATURE_G = {
        (1, 'quadrangle', 1): np.array([[0.]]),
        (1, 'quadrangle', 2): np.array([[-1.], [1.]])/np.sqrt(3.),
        (2, 'triangle', 1): np.array([[1., 1.]])/3.,
        (2, 'triangle', 2): np.array([[1./6., 1./6.], [2./3, 1./6], [1./6., 2./3.]]),
        (2, 'quadrangle', 1): np.array([[-1., -1.], [-1.,  1.],
                                        [ 1.,  1.], [ 1., -1.]])/np.sqrt(3.),
        (2, 'quadrangle', 2): np.array([]),
        (3, 'quadrangle', 1): np.array([[-1., -1., -1.], [-1.,  1., -1.],
                                        [ 1.,  1., -1.], [ 1., -1., -1.],
                                        [-1., -1.,  1.], [-1.,  1.,  1.],
                                        [ 1.,  1.,  1.], [ 1., -1.,  1.]])/np.sqrt(3.),
        (3, 'quadrangle', 2): np.array([]),
    }

    ELEMENT_TYPES = {
        '_segment_2':      ('quadrangle', 1, 'lagrange', 1, 2),
        '_segment_3':      ('quadrangle', 2, 'lagrange', 1, 3),
        '_triangle_3':     ('triangle', 1, 'lagrange', 2, 3),
        '_triangle_6':     ('triangle', 2, 'lagrange', 2, 6),
        '_quadrangle_4':   ('quadrangle', 1, 'serendip', 2, 4),
        '_quadrangle_8':   ('quadrangle', 2, 'serendip', 2, 8),
        '_tetrahedron_4':  ('triangle', 1, 'lagrange', 3, 4),
        '_tetrahedron_10': ('triangle', 2, 'lagrange', 3, 10),
        '_pentahedron_6':  ('pentahedron', 1, 'lagrange', 3, 6),
        '_pentahedron_15': ('pentahedron', 2, 'lagrange', 3, 15),
        '_hexahedron_8':   ('quadrangle', 1, 'serendip', 3, 8),
        '_hexahedron_20':  ('quadrangle', 2, 'serendip', 3, 20),
    }

    MONOMES = {(1, 'quadrangle'): np.array([[0], [1], [2], [3], [4], [5]]),
               (2, 'triangle'): np.array([[0, 0],                     #        1
                                          [1, 0], [0, 1],             #     x     y
                                          [2, 0], [1, 1], [0, 2]]),   # x^2   x.y   y^2
               (2, 'quadrangle'): np.array([[0, 0],
                                            [1, 0], [1, 1], [0, 1],
                                            [2, 0], [2, 1], [1, 2], [0, 2]]),
               (3, 'triangle'): np.array([[0, 0, 0],
                                          [1, 0, 0], [0, 1, 0], [0, 0, 1],
                                          [2, 0, 0], [1, 1, 0], [0, 2, 0], [0, 1, 1], [0, 0, 2],
                                          [1, 0, 1]]),
               (3, 'quadrangle'): np.array([[0, 0, 0],
                                            [1, 0, 0], [0, 1, 0], [0, 0, 1],
                                            [1, 1, 0], [0, 1, 1], [1, 0, 1], [1, 1, 1],
                                            [2, 0, 0], [0, 2, 0], [0, 0, 2],
                                            [2, 1, 0], [2, 0, 1], [2, 1, 1],
                                            [1, 2, 0], [0, 2, 1], [1, 2, 1],
                                            [1, 0, 2], [0, 1, 2], [1, 1, 2]]),
    }

    def __init__(self,  element):
        self._shape, self._order, self._inter_poly, self._dim, self._nnodes = self.ELEMENT_TYPES[element]
        self._ksi = self.NATURAL_COORDS[(self._dim, self._shape)][:self._nnodes]
        self._g = self.QUADRATURE_G[(self._dim, self._shape, self._order)]
        self._w = self.QUADRATURE_W[(self._dim, self._shape, self._order)]
        momo = self.MONOMES[(self._dim, self._shape)][:self._nnodes]

        _shape = list(momo[0])

        for s in range(len(_shape)):
            _shape[s] = max(momo[:,s])+1
        self._poly_shape = tuple(_shape)

        self._monome = []
        for m in momo:
            p = np.zeros(self._poly_shape)
            p[tuple(m)] = 1
            self._monome.append(p)

    def polyval(self, x, p):
        if 1 == self._dim:
            return poly.polyval(x[0], p)
        if 2 == self._dim:
            return poly.polyval2d(x[0], x[1], p)
        if 3 == self._dim:
            return poly.polyval3d(x[0], x[1], x[2], p)

    def precompute(self, **kwargs):
        X = np.array(kwargs["X"],  copy=False)
        nb_element = X.shape[0]
        X = X.reshape(nb_element, self._nnodes, self._dim)

        # evaluate polynomial constant for shapes
        _x = self._ksi
        _xe = np.zeros((self._nnodes, self._nnodes))
        for n in range(self._nnodes):
            _xe[:,n] = [self.polyval(_x[n], m) for m in self._monome]

        _a = np.linalg.inv(_xe)
        _n = np.zeros((self._nnodes,) + self._monome[0].shape)
        # set shapes polynomials
        for n in range(self._nnodes):
            for m in range(len(self._monome)):
                _n[n] += _a[n, m] * self._monome[m]

        # sanity check on shapes
        for n in range(self._nnodes):
            for m in range(self._nnodes):
                v = self.polyval(_x[n], _n[m])
                ve = 1. if n == m else 0.
                test = np.isclose(v, ve)
                if not test:
                    raise("Most probably an error in the shapes evaluation")

        # compute shapes derivatives
        _b = np.zeros((self._dim, self._nnodes,) + self._poly_shape)
        for d in range(self._dim):
            for n in range(self._nnodes):
                _der = poly.polyder(_n[n], axis=d)
                _mshape = np.array(self._poly_shape)
                _mshape[d] = _mshape[d] - _der.shape[d]
                _mshape = tuple(_mshape)
                _comp = np.zeros(_mshape)

                if 1 == self._dim:
                    _bt = np.hstack((_der, _comp))
                else:
                    if 0 == d:
                        _bt = np.vstack((_der, _comp))
                    if 1 == d:
                        _bt = np.hstack((_der, _comp))
                    if 2 == d:
                        _bt = np.dstack((_der, _comp))

                _b[d, n] = _bt

        _nb_quads = len(self._g)
        _nq = np.zeros((_nb_quads, self._nnodes))
        _bq = np.zeros((_nb_quads, self._dim, self._nnodes))

        # evaluate shapes and shapes derivatives on gauss points
        for q in range(_nb_quads):
            _g = self._g[q]
            for n in range(self._nnodes):
                _nq[q, n] = self.polyval(_g, _n[n])

                for d in range(self._dim):
                    _bq[q, d, n] = self.polyval(_g, _b[d, n])

        _j = np.array(kwargs['j'], copy=False).reshape((nb_element, _nb_quads))
        _B = np.array(kwargs['B'], copy=False).reshape((nb_element, _nb_quads, self._dim, self._nnodes))
        _N = np.array(kwargs['N'], copy=False).reshape((nb_element, _nb_quads, self._nnodes))

        for e in range(nb_element):
            for q in range(_nb_quads):
                _J = np.matmul(_bq[q], X[e])

                _N[e, q] = _nq[q]
                _B[e, q] = np.matmul(np.linalg.inv(_J), _bq[q])
                _j[e, q] = np.linalg.det(_J) * self._w[q]
