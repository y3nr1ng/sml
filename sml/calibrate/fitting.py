from collections import namedtuple
import logging
from math import sqrt

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

logger = logging.getLogger(__name__)

Arguments = namedtuple('Arguments', ['w0', 'A', 'B', 'c', 'd'])

#
# f(x) = w0 * sqrt( 1 + ((x-c)/d)**2 * (1 + A*((x-c)/d) + B*((x-c)/d)**2) )
#
def fit_function():
    pass

def _init_parameters(df):
    i_max = df['w0'].idxmax(axis='index')

    df = nan
    dx = nan
    d = dx / sqrt(2. * df/w0)

    a0 = Arguments(w0, .5, .5, c, d)
    logger.info(("initialize (w0, A, B, c, d) = "
                 "({w0}, {A}, {B}, {c}, {d})").format(**a0))
    return a0

def fit_curve(df, na=5, max_iter=200, tol=1e-6):
    a0 = _init_parameters(data)
