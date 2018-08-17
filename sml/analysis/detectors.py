from abc import ABCMeta, abstractmethod

import numpy as np
import numpy.ma as ma

class Base(metaclass=ABCMeta):
    @abstractmethod
    def __call__(self, frame):
        pass

    def reset(self):
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

        # use masked array for elimination
        d_frame = d_frame.view(ma.MaskedArray)
        coords = self._find_peaks(d_frame)

        # replace selected result
        mask = ma.getmaskarray(d_frame)
        curr_frame[mask] = prev_frame[mask]
        self._prev_frame = curr_frame

        return coords

    def _find_peaks(self, data):
        """
        TBA

        Parameter
        ---------
        data : np.ndarray
            A single timepoint.
        """
        if data.ndim < 2 or data.ndim > 3:
            raise ValueError("not an image or volume")
        ny, nx = data.shape[-2:]

        indices = data.reshape((-1, nx*ny)).argsort(axis=1)
        coords = []
        for z in range(indices.shape[0]):
            indices[z, ::-1] += z * (nx*ny)
            for index in indices[z, ::-1]:
                coord = np.unravel_index(index, dims=data.shape)

                if data[coord] < self.threshold:
                    # data is sorted, so rest of the pixels are background
                    break

                r = self.radius
                y, x = coord[-2:]
                if (x < r or y < r) or (x > nx-r+1 or y > ny-r+1):
                    # border pixel, ignored
                    continue

                if data.ndim == 2:
                    if not ma.is_masked(data[y-r:y+r+1, x-r:x+r+1]):
                        data[y-r:y+r+1, x-r:x+r+1] = ma.masked
                        coords.append(coord)
                else:
                    if not ma.is_masked(data[z, y-r:y+r+1, x-r:x+r+1]):
                        data[z, y-r:y+r+1, x-r:x+r+1] = ma.masked
                        coords.append(coord)
        return np.array(coords).T

    def reset(self):
        self._prev_frame = None
