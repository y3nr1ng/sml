from collections import namedtuple
import logging
from math import sqrt

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

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
    i = w.argmax()
    print(i)

    # center of the peak
    zc = z[i]
    w0 = w[i]

    # finding deltas to estimate sigma^2
    if i == 0:
        dz = z[i+1]-zc
        dw = w0-w[i+1]
    elif i == z.size-1:
        dz = zc-z[i-1]
        dw = w0-w[i-1]
    else:
        dz = ((zc-z[i-1]) + (z[i+1]-zc)) / 2.
        dw = ((w0-w[i-1]) + (w0-w[i+1])) / 2.
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

    p, _ = curve_fit(fit_function, z, w, p0=p0, check_finite=True, bounds=(-np.inf, np.inf), method='trf', verbose=2)

    return p
