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
