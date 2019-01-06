from collections import namedtuple
from functools import partial
import logging
from math import sqrt

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, root, OptimizeWarning

from . import defocus, lookup

logger = logging.getLogger(__name__)

defocus_model_dispatch = {
    'polynomial': defocus.Polynomial,
    'huang': defocus.Huang
}

lookup_method_dispatch = {
    'ratio': lookup.Ratio,
    'huang': lookup.Huang
}

def fit_depth_curve(z, w, model='huang', tol=1e-5):
    """
    Generate the model fitting of z v.s. w.

    Parameters
    ----------
    z : np.ndarray
        Depth (spot z-coordinates).
    w : np.ndarray
        Dependent variable of z (spot widths or heights).
    model : str
        Model to fit the curve, 'huang', 'polynomial' (TBA).
    tol : float
        Tolerance of z-step variation during fitting.
    """
    model = defocus_model_dispatch.get(model)()
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

def generate_lookup_function(z, w, h, method='ratio', model='huang', tol=1e-5):
    """
    Generate the calibration curve from the experimental data.

    Parameters
    ----------
    - z[:]:   z coordinates of the spots.
    - w[:]:   spot widths of correspoinding z coordinate.
    - w[:]:   spot heights of correspoinding z coordinate.
    - method: 
    - model:  model to fit the curve: 'huang', 'polynomial' (TBA)
    - tol:    fitting tolerance.
    """
    fw = fit_depth_curve(z, w, model, tol)
    fh = fit_depth_curve(z, h, model, tol)

    # find intersection
    sol = root(lambda z: fw(z)-fh(z), (z.min()+z.max())/2., tol=tol)
    if not sol.success:
        raise RuntimeError("unable to find z-center")
    elif sol.x.size > 1:
        logger.warning("multiple intersections found, {}"
                       ", using first one instead".format(sol.x))
    z0 = sol.x[0]
    logger.debug("intersection at z={:.4f}".format(z0))

    method = lookup_method_dispatch.get(method)(fw, fh, z0, tol)
    return method
