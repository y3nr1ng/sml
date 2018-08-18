from glob import glob
import os

import imageio
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import pandas as pd

from sml.analysis.detectors import LocalMaximum

"""
src_path = "../tests/data/android/android_128_raw.tif"
frames = imageio.mvolread(src_path)
"""
"""
src_dir = "../tests/data/3d_dataset"
paths = glob(os.path.join(src_dir, "*.tif"))
frames = [imageio.volread(path) for path in paths]
frames = np.stack(frames, axis=0)
"""
src_dir = "../tests/data/3d_dataset"
paths = glob(os.path.join(src_dir, "*.tif"))
frames = imageio.volread(paths[0])

print("{} frames".format(frames.shape[0]))
print("size {}".format(frames.shape[1:]))

fig, ax = plt.subplots()

r = 4
det = LocalMaximum(r, 10)

im = ax.imshow(np.zeros(frames[0].shape[-2:]))
m = ax.scatter(0, 0, marker='+', c='r', s=2*r+1, linewidth=1, edgecolor='r')

for t, frame in enumerate(frames):
    coords, mask = det(frame)
    if coords is None:
        continue
    print(coords.shape)

    """
    df = pd.DataFrame(coords, columns=["z [nm]", "y [nm]", "x [nm]"])
    df.to_csv("result.csv", index=False)
    """

    if frame.ndim == 3:
        im.set_data(frame[0, :])
        m.set_offsets(coords[coords[:, 0] == 0][:, :0:-1])
    else:
        im.set_data(frame)
        m.set_offsets(coords)
    im.autoscale()

    plt.draw()
    plt.savefig('{:03d}.png'.format(t), bbox_inches='tight')
