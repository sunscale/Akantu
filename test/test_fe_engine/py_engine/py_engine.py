#!/usr/bin/env python3

# -*- coding: utf-8 -*-
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

__all__ = ['Shapes']

import numpy as np
import numpy.polynomial.polynomial as poly

# pylint: disable=missing-docstring, invalid-name, too-many-instance-attributes
# flake8: noqa


class Shapes:
    """Python version of the shape functions for test purposes"""
    # pylint: disable=bad-whitespace, line-too-long

    NATURAL_COORDS = {
        (1, 'quadrangle'):  np.array([[-1.], [1.], [0.]]),
        (2, 'quadrangle'):  np.array([[-1., -1.], [ 1., -1.], [ 1.,  1.], [-1.,  1.],
                                      [ 0., -1.], [ 1.,  0.], [ 0.,  1.], [-1.,  0.]]),
        (3, 'quadrangle'):  np.array([[-1., -1., -1.], [ 1., -1., -1.], [ 1.,  1., -1.], [-1.,  1., -1.],
                                      [-1., -1.,  1.], [ 1., -1.,  1.], [ 1.,  1.,  1.], [-1.,  1.,  1.],
                                      [ 0., -1., -1.], [ 1.,  0., -1.], [ 0.,  1., -1.], [-1.,  0., -1.],
                                      [-1., -1.,  0.], [ 1., -1.,  0.], [ 1.,  1.,  0.], [-1.,  1.,  0.],
                                      [ 0., -1.,  1.], [ 1.,  0.,  1.], [ 0.,  1.,  1.], [-1.,  0.,  1.]]),
        (2, 'triangle'):    np.array([[0., 0.], [1., 0.], [0., 1], [.5, 0.], [.5, .5], [0., .5]]),
        (3, 'triangle'):    np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.],
                                      [.5, 0., 0.], [.5, .5, 0.], [0., .5, 0.],
                                      [0., 0., .5], [.5, 0., .5], [0., .5, .5]]),
        (3, 'pentahedron'): np.array([[-1.,  1.,  0.], [-1.,  0.,  1.], [-1.,  0.,  0.],
                                      [ 1.,  1.,  0.], [ 1.,  0.,  1.], [ 1.,  0.,  0.],
                                      [-1.,  .5,  .5], [-1.,  0.,  .5], [-1.,  .5,  0.],
                                      [ 0.,  1.,  0.], [ 0.,  0.,  1.], [ 0.,  0.,  0.],
                                      [ 1.,  .5,  .5], [ 1.,  0.,  .5], [ 1.,  .5,  0.],
                                      [ 0.,  .5,  .5], [ 0.,  0.,  .5], [ 0.,  .5,  0.]]),
    }

    QUADRATURE_W = {
        (1, 'quadrangle', 1): np.array([2.]),
        (1, 'quadrangle', 2): np.array([1., 1.]),
        (2, 'triangle', 1): np.array([1./2.]),
        (2, 'triangle', 2): np.array([1., 1., 1.])/6.,
        (3, 'triangle', 1): np.array([1./6.]),
        (3, 'triangle', 2): np.array([1., 1., 1., 1.])/24.,
        (2, 'quadrangle', 1): np.array([1., 1., 1., 1.]),
        (2, 'quadrangle', 2): np.array([1., 1., 1., 1.]),
        (3, 'quadrangle', 1): np.array([1., 1., 1., 1.,
                                        1., 1., 1., 1.]),
        (3, 'quadrangle', 2): np.array([1., 1., 1., 1.,
                                        1., 1., 1., 1.]),
        (3, 'pentahedron', 1): np.array([1., 1., 1., 1., 1., 1.])/6.,
        (3, 'pentahedron', 2): np.array([1., 1., 1., 1., 1., 1.])/6.,
    }

    _tet_a = (5. - np.sqrt(5.))/20.
    _tet_b = (5. + 3.*np.sqrt(5.))/20.
    QUADRATURE_G = {
        (1, 'quadrangle', 1): np.array([[0.]]),
        (1, 'quadrangle', 2): np.array([[-1.], [1.]])/np.sqrt(3.),
        (2, 'triangle', 1): np.array([[1., 1.]])/3.,
        (2, 'triangle', 2): np.array([[1./6., 1./6.], [2./3, 1./6], [1./6., 2./3.]]),
        (3, 'triangle', 1): np.array([[1., 1., 1.]])/4.,
        (3, 'triangle', 2): np.array([[_tet_a, _tet_a, _tet_a],
                                      [_tet_b, _tet_a, _tet_a],
                                      [_tet_a, _tet_b, _tet_a],
                                      [_tet_a, _tet_a, _tet_b]]),
        (2, 'quadrangle', 1): np.array([[-1., -1.], [ 1., -1.],
                                        [-1.,  1.], [ 1.,  1.]])/np.sqrt(3.),
        (2, 'quadrangle', 2): np.array([[-1., -1.], [ 1., -1.],
                                        [-1.,  1.], [ 1.,  1.]])/np.sqrt(3.),
        (3, 'quadrangle', 1): np.array([[-1., -1., -1.],
                                        [ 1., -1., -1.],
                                        [-1.,  1., -1.],
                                        [ 1.,  1., -1.],
                                        [-1., -1.,  1.],
                                        [ 1., -1.,  1.],
                                        [-1.,  1.,  1.],
                                        [ 1.,  1.,  1.]])/np.sqrt(3.),
        (3, 'quadrangle', 2): np.array([[-1., -1., -1.],
                                        [ 1., -1., -1.],
                                        [-1.,  1., -1.],
                                        [ 1.,  1., -1.],
                                        [-1., -1.,  1.],
                                        [ 1., -1.,  1.],
                                        [-1.,  1.,  1.],
                                        [ 1.,  1.,  1.]])/np.sqrt(3.),
        (3, 'pentahedron', 1): np.array([[-1./np.sqrt(3.), 1./6., 1./6.],
                                         [-1./np.sqrt(3.), 2./3., 1./6.],
                                         [-1./np.sqrt(3.), 1./6., 2./3.],
                                         [ 1./np.sqrt(3.), 1./6., 1./6.],
                                         [ 1./np.sqrt(3.), 2./3., 1./6.],
                                         [ 1./np.sqrt(3.), 1./6., 2./3.]]),
        (3, 'pentahedron', 2): np.array([[-1./np.sqrt(3.), 1./6., 1./6.],
                                         [-1./np.sqrt(3.), 2./3., 1./6.],
                                         [-1./np.sqrt(3.), 1./6., 2./3.],
                                         [ 1./np.sqrt(3.), 1./6., 1./6.],
                                         [ 1./np.sqrt(3.), 2./3., 1./6.],
                                         [ 1./np.sqrt(3.), 1./6., 2./3.]]),
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
                                            [1, 1, 0], [1, 0, 1],
                                            [0, 1, 1], [1, 1, 1],
                                            [2, 0, 0], [0, 2, 0], [0, 0, 2],
                                            [2, 1, 0], [2, 0, 1], [2, 1, 1],
                                            [1, 2, 0], [0, 2, 1], [1, 2, 1],
                                            [1, 0, 2], [0, 1, 2], [1, 1, 2]])}

    SHAPES = {
        (3, 'pentahedron', 1): np.array([
            [[[ 0.,  0.], [ 1.,  0.]], [[ 0.,  0.], [-1.,  0.]]],
            [[[ 0.,  1.], [ 0.,  0.]], [[ 0., -1.], [ 0.,  0.]]],
            [[[ 1., -1.], [-1.,  0.]], [[-1.,  1.], [ 1.,  0.]]],
            [[[ 0.,  0.], [ 1.,  0.]], [[ 0.,  0.], [ 1.,  0.]]],
            [[[ 0.,  1.], [ 0.,  0.]], [[ 0.,  1.], [ 0.,  0.]]],
            [[[ 1., -1.], [-1.,  0.]], [[ 1., -1.], [-1.,  0.]]]
        ])/2.,
        (3, 'pentahedron', 2): np.array([
            # 0
            [[[ 0. ,  0. ,  0. ], [-1. ,  0. ,  0. ], [ 1. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 0.5,  0. ,  0. ], [-1. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 0.5,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],

            # 1
            [[[ 0. , -1. ,  1. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. ,  0.5, -1. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. ,  0.5,  0. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],

            # 2
            [[[ 0. , -1. ,  1. ], [-1. ,  2. ,  0. ], [ 1. ,  0. ,  0. ]],
             [[-0.5,  1.5, -1. ], [ 1.5, -2. ,  0. ], [-1. ,  0. ,  0. ]],
             [[ 0.5, -0.5,  0. ], [-0.5,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],

            # 3
            [[[ 0. ,  0. ,  0. ], [-1. ,  0. ,  0. ], [ 1. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [-0.5,  0. ,  0. ], [ 1. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 0.5,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],

            # 4
            [[[ 0. , -1. ,  1. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. , -0.5,  1. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. ,  0.5,  0. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],

            # 5
            [[[ 0. , -1. ,  1. ], [-1. ,  2. ,  0. ], [ 1. ,  0. ,  0. ]],
             [[ 0.5, -1.5,  1. ], [-1.5,  2. ,  0. ], [ 1. ,  0. ,  0. ]],
             [[ 0.5, -0.5,  0. ], [-0.5,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],

            # 6
            [[[ 0. ,  0. ,  0. ], [ 0. ,  2. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 0. , -2. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],

            # 7
            [[[ 0. ,  2. , -2. ], [ 0. , -2. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. , -2. ,  2. ], [ 0. ,  2. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],

            # 8
            [[[ 0. ,  0. ,  0. ], [ 2. , -2. ,  0. ], [-2. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [-2. ,  2. ,  0. ], [ 2. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],

            # 9
            [[[ 0. ,  0. ,  0. ], [ 1. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [-1. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],

            # 10
            [[[ 0. ,  1. ,  0. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. , -1. ,  0. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],

            # 11
            [[[ 1. , -1. ,  0. ], [-1. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[-1. ,  1. ,  0. ], [ 1. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],

            # 12
            [[[ 0. ,  0. ,  0. ], [ 0. ,  2. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 0. ,  2. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],

            # 13
            [[[ 0. ,  2. , -2. ], [ 0. , -2. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. ,  2. , -2. ], [ 0. , -2. ,  0. ], [ 0. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],

            # 14
            [[[ 0. ,  0. ,  0. ], [ 2. , -2. ,  0. ], [-2. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 2. , -2. ,  0. ], [-2. ,  0. ,  0. ]],
             [[ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ], [ 0. ,  0. ,  0. ]]],
        ])}
    # pylint: enable=bad-whitespace, line-too-long

    def __init__(self, element):
        self._shape, self._order, self._inter_poly, self._dim, \
            self._nnodes = self.ELEMENT_TYPES[element]
        self._ksi = self.NATURAL_COORDS[(self._dim, self._shape)][:self._nnodes]
        self._g = self.QUADRATURE_G[(self._dim, self._shape, self._order)]
        self._w = self.QUADRATURE_W[(self._dim, self._shape, self._order)]
        self._poly_shape = ()
        self._monome = []

    def polyval(self, x, p):
        if self._dim == 1:
            return poly.polyval(x[0], p)
        if self._dim == 2:
            return poly.polyval2d(x[0], x[1], p)
        if self._dim == 3:
            return poly.polyval3d(x[0], x[1], x[2], p)
        return None

    def shape_from_monomes(self):
        momo = self.MONOMES[(self._dim, self._shape)][:self._nnodes]

        _shapes = list(momo[0])

        for s, _ in enumerate(_shapes):
            _shapes[s] = max(momo[:, s])+1
        self._poly_shape = tuple(_shapes)

        self._monome = []
        for m in momo:
            p = np.zeros(self._poly_shape)
            p[tuple(m)] = 1
            self._monome.append(p)

        # evaluate polynomial constant for shapes
        _x = self._ksi
        _xe = np.zeros((self._nnodes, self._nnodes))
        for n in range(self._nnodes):
            _xe[:, n] = [self.polyval(_x[n], m) for m in self._monome]

        _a = np.linalg.inv(_xe)
        _n = np.zeros((self._nnodes,) + self._monome[0].shape)
        # set shapes polynomials
        for n in range(self._nnodes):
            for m in range(len(self._monome)):
                _n[n] += _a[n, m] * self._monome[m]

        return _n

    def compute_shapes(self):
        if (self._dim, self._shape) in self.MONOMES:
            return self.shape_from_monomes()

        _n = self.SHAPES[(self._dim, self._shape, self._order)]
        self._poly_shape = _n[0].shape
        return _n

    # pylint: disable=too-many-locals,too-many-branches
    def precompute(self, **kwargs):
        X = np.array(kwargs["X"], copy=False)
        nb_element = X.shape[0]
        X = X.reshape(nb_element, self._nnodes, self._dim)

        _x = self._ksi
        _n = self.compute_shapes()

        # sanity check on shapes
        for n in range(self._nnodes):
            for m in range(self._nnodes):
                v = self.polyval(_x[n], _n[m])
                ve = 1. if n == m else 0.
                test = np.isclose(v, ve)
                if not test:
                    raise Exception("Most probably an error in the shapes evaluation")

        # compute shapes derivatives
        _b = np.zeros((self._dim, self._nnodes,) + self._poly_shape)
        for d in range(self._dim):
            for n in range(self._nnodes):
                _der = poly.polyder(_n[n], axis=d)
                _mshape = np.array(self._poly_shape)
                _mshape[d] = _mshape[d] - _der.shape[d]
                _mshape = tuple(_mshape)
                _comp = np.zeros(_mshape)

                if self._dim == 1:
                    _bt = np.hstack((_der, _comp))
                else:
                    if d == 0:
                        _bt = np.vstack((_der, _comp))
                    if d == 1:
                        _bt = np.hstack((_der, _comp))
                    if d == 2:
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
        _B = np.array(kwargs['B'], copy=False).reshape((nb_element, _nb_quads,
                                                        self._nnodes, self._dim))
        _N = np.array(kwargs['N'], copy=False).reshape((nb_element, _nb_quads, self._nnodes))
        _Q = kwargs['Q']
        if np.linalg.norm(_Q - self._g.T) > 1e-15:
            raise Exception('Not using the same quadrature'
                            ' points norm({0} - {1}) = {2}'.format(_Q, self._g.T,
                                                                   np.linalg.norm(_Q - self._g.T)))

        for e in range(nb_element):
            for q in range(_nb_quads):
                _J = np.matmul(_bq[q], X[e])

                if np.linalg.norm(_N[e, q] - _nq[q]) > 1e-10:
                    print(f"{e},{q}")
                    print(_N[e, q])
                    print(_nq[q])
                _N[e, q] = _nq[q]
                _tmp = np.matmul(np.linalg.inv(_J), _bq[q])
                _B[e, q] = _tmp.T

                _j[e, q] = np.linalg.det(_J) * self._w[q]
