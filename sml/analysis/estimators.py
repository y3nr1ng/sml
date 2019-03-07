from abc import ABCMeta, abstractmethod

import numpy as np
from scipy import optimize

class Base(metaclass=ABCMeta):
    @abstractmethod
    def __call__(self, frame, candidates):
        pass

class GaussianMLE(Base):
    def __init__(self, radius):
        self.radius = radius
        vec = np.linspace(-radius, radius, 2*radius+1)
        self.gx_offsets, self.gy_offsets = np.meshgrid(vec, vec)

    def __call__(self, frame, candidates):
        r = self.radius

        coords = []
        for candidate in candidates:
            y, x = candidate[-2:]
            roi = [] if frame.ndim == 2 else [candidate[0]]
            roi += [slice(y-r, y+r+1), slice(x-r, x+r+1)]
            roi = tuple(roi)

            patch = frame[roi]

            gx = self.gx_offsets + x
            gy = self.gy_offsets + y



            opt, is_success = optimize.minimize()

    def _function(self, x, y):
        xo = float(xo)
        yo = float(yo)
        a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
        b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
        c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
        return offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))

class Phasor(Base):
    """
    Phasor based single-molecule localization microscopy in 3D (pSMLM-3D): An
    algorithm for MHz localization rates using standard CPUs
    doi: 10.1063/1.5005899
    """
    pass
