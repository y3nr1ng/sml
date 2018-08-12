from collections import namedtuple
from functools import partial
import logging
from math import sqrt

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, root, OptimizeWarning

from . import defocus, lookup

logger = logging.getLogger(__name__)

model_dispatch = {
    'polynomial': defocus.Polynomial,
    'huang': defocus.Huang
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
        model.arguments = p
        logger.debug("popt: {}".format(model))
    except (OptimizeWarning, RuntimeError):
        raise RuntimeError("unable to minimize the result")
    except ValueError:
        raise ValueError("provided samples contain NaNs")
    return model

method_dispatch = {
    'ratio': lookup.Ratio,
    'huang': lookup.Huang
}

def generate_lookup_function(z, w, h, method='ratio', model='huang', tol=1e-5):
    """
    TBA

    Parameters
    ----------
    """
    fw = fit_depth_curve(z, w, model, tol)
    fh = fit_depth_curve(z, h, model, tol)

    # find intersection
    sol = root(lambda z: fw(z)-fh(z), (z.min()+z.max())/2., tol=tol)
    if not sol.success:
        raise RuntimeError("unable to find z-center")
    elif sol.x.size > 1:
        raise RuntimeWarning("multiple intersections found {}, use first one".format(sol.x))
    z0 = sol.x[0]

    method = method_dispatch.get(method)(fw, fh, z0, tol)
    return method
