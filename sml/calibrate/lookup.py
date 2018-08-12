from abc import abstractmethod
from math import sqrt
from numbers import Number

import numpy as np
from scipy.optimize import root, minimize, OptimizeWarning

__all__ = [
    'Ratio',
    'Huang'
]

class Base(object):
    __slots__ = ('fw', 'fh', 'z0')

    def __init__(self, fw, fh, z0, tol=1e-5):
        self.fw = fw
        self.fh = fh
        self.z0 = z0
        self.tol = tol

    @abstractmethod
    def _function(self, z, w, h):
        raise NotImplementedError

    def _find_z(self, w, h):
        sol = root(lambda z: self._function(z, w, h), self.z0, tol=self.tol)
        return (sol.x[0] - self.z0) if sol.success else np.nan

    def __call__(self, w, h):
        if isinstance(w, np.ndarray) and isinstance(h, np.ndarray):
            return np.array([self._find_z(_w, _h) for _w, _h in np.c_[w, h]])
        elif isinstance(w, Number) and isinstance(h, Number):
            return self._find_z(w, h) - self.z0
        else:
            raise ValueError("incompatible input")

class Ratio(Base):
    def _function(self, z, w, h):
        return self.fw(z)/self.fh(z) - w/h

class Huang(Base):
    def _function(self, z, w, h):
        return sqrt((sqrt(w)-sqrt(self.fw(z)))**2 + (sqrt(h)-sqrt(self.fh(z)))**2)

    def _find_z(self, w, h):
        sol = minimize(lambda z: self._function(z, w, h), self.z0, method='Nelder-Mead', tol=self.tol)
        return (sol.x[0] - self.z0) if sol.success else np.nan
