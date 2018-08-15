from abc import ABCMeta, abstractmethod

import numpy as np
import numpy.ma as ma

class Base(metaclass=ABCMeta):
    @abstractmethod
    def __call__(self, frame):
        pass

class LocalMaximum(Base):
    def __init__(self, radius, threshold=-np.inf):
        self._prev_frame = None
        self.threshold = threshold
        self.radius = int(radius)

    def __call__(self, frame):
        curr_frame = frame.astype(dtype=np.int32, copy=True)
        # retrieve previous frame
        if self._prev_frame is None:
            self._prev_frame = curr_frame
            return None
        else:
            prev_frame = self._prev_frame

        # use difference imaging technique
        d_frame = curr_frame - prev_frame

        # ignore background
        d_frame[d_frame < self.threshold] = 0

        # use masked array for elimination
        m_frame = d_frame.view(ma.MaskedArray)
        # sort by intensity, ascending
        indices = m_frame.argsort(axis=None)
        coords = []
        for index in indices[::-1]:
            y, x = np.unravel_index(index, frame.shape)
            #NOTE subtract PSF instead of replace patch for multi-emitter
            # detection

            if d_frame[y, x] == 0:
                continue

            r = self.radius
            if (x < r or y < r) or (x > frame.shape[1]-r+1 or y > frame.shape[0]-r+1):
                continue
            if not ma.is_masked(m_frame[y-r:y+r, x-r:x+r]):
                m_frame[y-r:y+r, x-r:x+r] = ma.masked
                coords.append((y, x))

        # replace selected result
        mask = ma.getmaskarray(m_frame)
        frame[mask] = prev_frame[mask]
        self._prev_frame = frame

        return coords
