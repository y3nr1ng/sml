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

columns_generic = {
    # last observed frame index
    'frame': np.uint32,
    # number of accumulated events
    'occurrences': np.uint32,
    # fitting results
    'intensity': np.float32
    'background': np.float32
}

columns_2d = {
    # coordinate in the image
    'x': np.float32,
    'y': np.float32,
    # fitting results
    'uncertainty_xy': np.float32
}

columns_3d = {
    # coordinate in the image
    'x': np.float32,
    'y': np.float32,
    'z': np.float32,
    # fitting results
    'uncertainty_xy': np.float32,
    'uncertainty_z': np.float32
}

class Events(object):
    def __init__(self, dim):
        if dim == '2d':
            columns = {**columns_generic, **columns_2d}
        elif dim == '3d':
            columns = {**columns_generic, **columns_3d}
        else:
            raise ValueError("invalid result type")
            
        self._events = pd.DataFrame(columns=columns.keys())
        self._events = self._events.astype(columns)

    def preallocate_by(self, events):
        pass

    def add_event(self, event):
        pass
