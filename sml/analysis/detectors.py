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
    def __init__(self, frameRange, ROI, radius, thresMax=np.inf,
                 thresMin=-np.inf):
        """
        Parameters
        ----------
         - frameRange: int, (N1,N2), range of frames for analysis
           - N1=0    if N1 == None
           - N2=nL-1 if N2 == None
         - ROI: int, ((x1,y1),(x2,y2)), ROI of each frame.
         - radius: int, r or (rx,ry), radius of each spot.
         - thresMax: float, max threshold of intensity
         - thresMin: float, min threshold of intensity
        """
        self.FID1 = None
        self.FID2 = None
        self.ROI_x1 = None
        self.ROI_y1 = None
        self.ROI_x2 = None
        self.ROI_y2 = None

        if not frameRange == None:
            self.FID1 = int(frameRange[0])
            self.FID2 = int(frameRange[1])
        if not ROI == None:
            self.ROI_x1 = int(ROI[0][0])
            self.ROI_y1 = int(ROI[0][1])
            self.ROI_x2 = int(ROI[1][0])
            self.ROI_y2 = int(ROI[1][1])
        if (type(radius) is tuple):
            self.rx = int(radius[0])
            self.ry = int(radius[1])
        else:
            self.rx = int(radius)
            self.ry = int(radius)
        self.thresMax = float(thresMax)
        self.thresMin = float(thresMin)

    def __call__(self, frames):
        nL, ny, nx = frames.shape
        prev_frame = None
        coords = []
        FID1 = self.FID1
        FID2 = self.FID2

        if (self.FID1 == None):
            FID1 = 0
        if (self.FID2 == None):
            FID2 = nL-1
        if (FID1 < 0 or self.FID2 > nL-1):
            raise ValueError("The specified frames out of range.")

        for i in reversed(range(FID1,FID2+1)):
            curr_frame = frames[i,:,:].astype(dtype=np.int32, copy=True)
            if not self.ROI_x1 == None:
                curr_frame = curr_frame[self.ROI_y1:self.ROI_y2,
                                        self.ROI_x1:self.ROI_x2]
            if prev_frame is None:
                prev_frame = curr_frame
                continue

            # use difference imaging technique
            diff_frame = curr_frame - prev_frame

            # use masked array for elimination
            coord, mask = self._find_peaks(diff_frame)
            coords = coords + coord
            print("layer: %d" % i)
            for c in coord:
                print(c)

            if (i % 10 == 0):
                print("Process frame: %d, N_candidates: %d" % (i,len(coords)))

            # replace selected result
            curr_frame[mask] = prev_frame[mask]
            prev_frame = curr_frame

        return coords

    def _find_peaks(self, data):
        """
        Find the cadidate particles from a differential frame.

        Parameter
        ---------
        data : np.ndarray
            A single differential frame to find candidate particles.
        """
        rx = self.rx
        ry = self.ry
        candidates = []

        if data.ndim > 3 or data.ndim < 2:
            raise ValueError("not a volume or image")
        ny, nx = data.shape[-2:]
        nxy = nx * ny

        # insert a mask array into the array of frames.
        data = data.view(ma.MaskedArray)

        # ascending sort, by z
        indices = data.reshape((-1, nxy)).argsort(axis=1)
        for z, _indices in enumerate(indices):
            # convert back to 3d coordinate in reverse
            _indices += nxy * z
            coords = np.unravel_index(_indices[::-1], dims=data.shape)

            for coord in zip(*list(coords)):
                # stop at value beyond the thresholds
                if data[coord] > self.thresMax:
                    continue
                if data[coord] <= self.thresMin:
                    break

                # ignore border pixel
                y, x = coord[-2:]
                if (x < rx or y < ry) or (x > nx-rx+1 or y > ny-ry+1):
                    continue

                # coordinates of ROI: centered at ([z],x,y) with radius r.
                roi = [] if data.ndim == 2 else [z]
                roi += [slice(y-ry, y+ry+1), slice(x-rx, x+rx+1)]
                roi = tuple(roi)
                if not ma.is_masked(data[roi]):
                    data[roi] = ma.masked
                    candidates.append(coord)

        return candidates, ma.getmaskarray(data)

    def reset(self):
        self._prev_frame = None
