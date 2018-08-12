from abc import ABCMeta, abstractmethod
from math import sqrt

import numpy as np

__all__ = [
    'PolynomialModel',
    'HuangModel'
]

class BaseModel(object, metaclass=ABCMeta):
    __slots__ = ()

    @abstractmethod
    def initial_guess(self, xdata, ydata):
        """
        Initial guess on independent variables.

        Parameters
        ----------
        xdata : np.ndarray
            Independent variables.
        ydata : np.ndarray
            Dependent variables.

        Returns
        -------
        arguments : list
            List of arguments in order specified by `__slots__`.
        bounds : list of tuples
            Boundary for each arguments.
        """
        raise NotImplementedError

    @abstractmethod
    def _function(self, x, *args):
        raise NotImplementedError

    @property
    def arguments(self):
        return [getattr(self, attr) for attr in self.__slots__]

    @arguments.setter
    def arguments(self, args):
        for attr, arg in zip(self.__slots__, args):
            setattr(self, attr, arg)

    def __call__(self, x, *args):
        """
        Apply the model with specified independent variable and model arguments.
        If arguments are not specified, saved arguments are used.

        Parameters
        ----------
        x : scalar or np.ndarray
            Independent variable to apply the model (y=f(x)) to.
        """
        if args:
            return self._function(x, *args)
        else:
            return self._function(x, *self.arguments)

    def __str__(self):
        dump = '('
        dump += ', '.join(self.__slots__)
        dump += ')=('
        dump += ', '.join("{:.4f}".format(x) for x in self.arguments)
        dump += ')'
        return dump

class PolynomialModel(BaseModel):
    """
    Describe the relationship between the axial position of a molecule in its
    imaged widths along two perpendicular axes by a pair of second degree
    polynomials.
    """
    __slots__ = ('w0', 'A', 'B', 'zc')

class HuangModel(BaseModel):
    """
    The relationship between the axial position of a molecule and its imaged
    widths along perpendicular axes is given by Huang et al.
    """
    __slots__ = ('w0', 'A', 'B', 'zc', 's')

    def initial_guess(self, z, w):
        iw_m = w.argmin()

        # center of the peak
        zc = z[iw_m]
        w0 = w[iw_m]

        # finding deltas to estimate sigma^2
        if iw_m == 0:
            dz = z[iw_m+1]-zc
            dw = w[iw_m+1]-w0
        elif iw_m == z.size-1:
            dz = zc-z[iw_m-1]
            dw = w[iw_m-1]-w0
        else:
            dz = ((zc-z[iw_m-1]) + (z[iw_m+1]-zc)) / 2.
            dw = ((w[iw_m-1]-w0) + (w[iw_m+1]-w0)) / 2.
        d = dz / sqrt(2. * dw/w0)

        self.w0, self.A, self.B, self.zc, self.s = w0, .5, .5, zc, d

        return self.arguments, None

    def _function(self, x, w0, A, B, zc, s):
        return w0 * np.sqrt( 1 + ((x-zc)/s)**2 * (1 + A*((x-zc)/s) + B*((x-zc)/s)**2) )
