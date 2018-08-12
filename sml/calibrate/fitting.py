from collections import namedtuple
from functools import partial
import logging
from math import sqrt

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, root

from .models import *

logger = logging.getLogger(__name__)

model_dispatch = {
    'polynomial': PolynomialModel,
    'huang': HuangModel
}

def fit_depth_curve(z, w, model='huang', tol=1e-5):
    """
    TBA

    Parameters
    ----------
    z : np.ndarray
        Depth.
    w : np.ndarray
        Dependent variable of z.
    model : str
        Model to fit the curve, 'huang', 'polynomial' (TBA).
    tol : float
        Tolerance of z-step variation during fitting.
    """
    model = model_dispatch.get(model)()
    if model is None:
        raise ValueError("invalid model is provided")
    p0, bounds = model.initial_guess(z, w)
    logger.debug("p0: {}".format(model))
    try:
        p, _ = curve_fit(model, z, w, p0=p0, bounds=(-np.inf, np.inf), method='trf', xtol=tol)
    except OptimizeWarning:
        pass
    except RuntimeError:
        logger.error("unable to minimize the result")
    except ValueError:
        logger.error("provided samples contain NaNs")

    model.arguments = p
    logger.debug("popt: {}".format(model))
    return model

def generate_lookup_function(z, w, h, model='huang', tol=1e-5):
    fw = fit_depth_curve(z, w, model, tol)
    fh = fit_depth_curve(z, h, model, tol)

    # find intersection
    sol = root(lambda z: fw(z)-fh(z), (z.min()+z.max())/2., tol=tol)
    if not sol.success:
        raise RuntimeError("unable to find z-center")
    elif sol.x.size > 1:
        raise RuntimeWarning("multiple intersections found {}, use first one".format(sol.x))
    z0 = sol.x[0]

    # different lookup method
    def _find_depth_mean(w, h):
        wsol = root(lambda z: fw(z)-w, z0, tol=tol)
        hsol = root(lambda z: fh(z)-h, z0, tol=tol)
        if (not wsol.success) and (not hsol.success):
            raise RuntimeError("unable to determine z position")
        return (wsol.x[0]+hsol.x[0])/2.-z0

    def _find_depth_ratio(w, h):
        sol = root(lambda z: fw(z)/fh(z)-w/h, z0, tol=tol)
        if not sol.success:
            print(sol)
            #raise RuntimeError("unable to determine z position")
        return sol.x[0]-z0

    def _find_depth_diff(w, h):
        sol = root(lambda z: (fw(z)-fh(z))-(w-h), z0, tol=tol)
        if not sol.success:
            print(sol)
            #raise RuntimeError("unable to determine z position")
        return sol.x[0]-z0

    # generate the lookup function
    def find_depth(w, h):
        _find_depth = _find_depth_diff
        if isinstance(w, np.ndarray) and isinstance(h, np.ndarray):
            return np.array([_find_depth(_w, _h) for _w, _h in np.c_[w, h]])
        else:
            return _find_depth(w, h)
    return find_depth, fw, fh, z0
