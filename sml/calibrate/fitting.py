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
    try:
        p, _ = curve_fit(model, z, w, p0=p0, bounds=(-np.inf, np.inf), method='trf', xtol=tol)
    except OptimizeWarning:
        pass
    except RuntimeError:
        logger.error("unable to minimize the result")
    except ValueError:
        logger.error("provided samples contain NaNs")

    model.arguments = p
    return model

def generate_lookup_function(z, w, h, model='huang', tol=1e-5):
    fw = fit_depth_curve(z, w, model, tol)
    fh = fit_depth_curve(z, h, model, tol)

    # find intersection
    def root_func(z):
        return fw(z)-fh(z)
    sol = root(root_func, (z.min()+z.max())/2., tol=tol)
    if not sol.success:
        raise RuntimeError("unable to find z-center")
    elif sol.x.size > 1:
        raise RuntimeWarning("multiple intersections found {}, use first one".format(sol.x))
    return sol.x[0], fw, fh
