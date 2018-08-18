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
        self.prev_frame = None
        self.threshold = threshold
        self.radius = int(radius)

    def __call__(self, frame):
        curr_frame = frame.astype(dtype=np.int32, copy=True)
        # retrieve previous frame
        if self.prev_frame is None:
            self.prev_frame = curr_frame
            return None, None

        # use difference imaging technique
        diff_frame = curr_frame - self.prev_frame

        # use masked array for elimination
        coords, mask = self._find_peaks(diff_frame)

        # replace selected result
        curr_frame[mask] = self.prev_frame[mask]
        self.prev_frame = curr_frame

        return coords, mask

    def _find_peaks(self, data):
        """
        TBA

        Parameter
        ---------
        data : np.ndarray
            A single timepoint.
        """
        shape = data.shape
        ny, nx = shape[-2:]

        data = data.view(ma.MaskedArray)

        if data.ndim == 2:
            data = data.reshape((1, ny, nx))
        elif data.ndim > 3:
            raise ValueError("not an image or volume")

        candidates = []

        # ascending sort, by z
        indices = data.reshape((-1, nx*ny)).argsort(axis=1)
        for z, _indices in enumerate(indices):
            # convert back to 3d coordinate in reverse
            _indices += (nx*ny) * z
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

                if not ma.is_masked(data[z, y-r:y+r+1, x-r:x+r+1]):
                    data[z, y-r:y+r+1, x-r:x+r+1] = ma.masked
                    candidates.append(coord)

        candidates = np.array(candidates)
        try:
            data = data.reshape(-1, ny, nx)
            data = data.squeeze(axis=0)

            print("{}, squeezed".format(candidates.shape))

            # 2d data
            candidates = candidates[:, :0:-1]
        except ValueError:
            # 3d data cannot be squeezed
            pass
        return candidates, ma.getmaskarray(data)

    def reset(self):
        self._prev_frame = None
