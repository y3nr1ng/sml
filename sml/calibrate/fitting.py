from collections import namedtuple
from functools import partial
import logging
from math import sqrt

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from .models import HuangModel

logger = logging.getLogger(__name__)

Parameters = namedtuple('Parameters', ['w0', 'A', 'B', 'zc', 'd'])

def fit_function(x, w0, A, B, zc, d):
    return w0 * np.sqrt( 1 + ((x-zc)/d)**2 * (1 + A*((x-zc)/d) + B*((x-zc)/d)**2) )

def _init_parameters(z, w):
    """TBA

    Parameters
    ----------
    z : np.ndarray
        Input array.
    w : np.ndarray

    """
    i = w.argmin()

    # center of the peak
    zc = z[i]
    w0 = w[i]

    # finding deltas to estimate sigma^2
    if i == 0:
        dz = z[i+1]-zc
        dw = w[i+1]-w0
    elif i == z.size-1:
        dz = zc-z[i-1]
        dw = w[i-1]-w0
    else:
        dz = ((zc-z[i-1]) + (z[i+1]-zc)) / 2.
        dw = ((w[i-1]-w0) + (w[i+1]-w0)) / 2.
    print(dz)
    print(dw)
    print(w0)
    d = dz / sqrt(2. * dw/w0)

    p0 = Parameters(w0, .5, .5, zc, d)
    logger.info(("initialize (w0, A, B, zc, d) = "
                 "({:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f})").format(p0.w0, p0.A, p0.B, p0.zc, p0.d))
    return p0

def fit_curve(z, w, na=5, max_iter=200, tol=1e-6):
    p0 = _init_parameters(z, w)
    print(p0)

    p, _ = curve_fit(fit_function, z, w, p0=p0, check_finite=True, bounds=(-np.inf, np.inf), method='trf')
    print(p)

    f = partial(fit_function, w0=p[0], A=p[1], B=p[2], zc=p[3], d=p[4])
    wo = [f(zi) for zi in z]

    return p, wo
