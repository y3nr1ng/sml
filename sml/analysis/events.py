import numpy as np
import pandas as pd

"""
typedef struct {
    double    a[6], da[6];	// fitting parameters.
    double    chisq;		// fitting chisq/ndof.
} fitres_t;

typedef struct {
    int       x, y;		// spot coordinate.
    int       cnt;		// number of accumulates
    int       fID;		// the last frame index
    int      *img;		// pointer to frame image.
    fitres_t *res;		// the fitting results.
} sp_t;
"""

columns = {
    # last observed frame index
    'frame': np.uint32,
    # number of accumulated events
    'occurrences': np.uint32,
    # coordinate in the image
    'x': np.uint16,
    'y': np.uint16
    #TODO fitting results
}

class Events(object):
    def __init__(self):
        self._events = pd.DataFrame(columns=columns.keys())
        self._events = self._events.astype(columns)
