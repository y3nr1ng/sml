from collections import namedtuple
from functools import partial
import logging
from math import sqrt

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from .models import HuangModel

logger = logging.getLogger(__name__)

def fit_curve(z, w, na=5, max_iter=200, tol=1e-6):
    model = HuangModel()
    p0, _ = model.initial_guess(z, w)
    p, _ = curve_fit(model, z, w, p0=p0, bounds=(-np.inf, np.inf), method='trf')

    print(p)

    model.arguments = p
    wo = [model(zi) for zi in z]

    return p, wo
