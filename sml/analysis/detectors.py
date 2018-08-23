from abc import ABCMeta, abstractmethod

import numpy as np
import numpy.ma as ma

class Base(metaclass=ABCMeta):
    @abstractmethod
    def __call__(self, frame):
        pass

    def reset(self):
        pass

class DifferenceImaging(Base):
    def __init__(self, radius, threshold=-np.inf):
        self.prev_frame = None
        self.threshold = threshold
        self.radius = int(radius)

    def __call__(self, frame):
        curr_frame = frame.astype(dtype=np.int32, copy=True)
        # retrieve previous frame
        if self.prev_frame is None:
            self.prev_frame = curr_frame
            return None

        # use difference imaging technique
        diff_frame = curr_frame - self.prev_frame

        # use masked array for elimination
        coords, mask = self._find_peaks(diff_frame)

        # replace selected result
        curr_frame[mask] = self.prev_frame[mask]
        self.prev_frame = curr_frame

        return coords

    def _find_peaks(self, data):
        """
        TBA

        Parameter
        ---------
        data : np.ndarray
            A single timepoint.
        """
        ny, nx = data.shape[-2:]
        nxy = nx * ny
        data = data.view(ma.MaskedArray)

        if data.ndim > 3 or data.ndim < 2:
            raise ValueError("not a volume or image")

        candidates = []

        # ascending sort, by z
        indices = data.reshape((-1, nxy)).argsort(axis=1)
        for z, _indices in enumerate(indices):
            # convert back to 3d coordinate in reverse
            _indices += nxy * z
            coords = np.unravel_index(_indices[::-1], dims=data.shape)

            for coord in zip(*list(coords)):
                # stop at value lower then threshold
                if data[coord] < self.threshold:
                    break

                # ignore border pixel
                r = self.radius
                y, x = coord[-2:]
                if (x < r or y < r) or (x > nx-r+1 or y > ny-r+1):
                    continue

                roi = [] if data.ndim == 2 else [z]
                roi += [slice(y-r, y+r+1), slice(x-r, x+r+1)]
                roi = tuple(roi)
                if not ma.is_masked(data[roi]):
                    data[roi] = ma.masked
                    candidates.append(coord)

        return candidates, ma.getmaskarray(data)

    def reset(self):
        self._prev_frame = None
