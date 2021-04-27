import numpy as np
from math import pi
from trackgen import Track

# we do not need to optimize this track, as it is already closed.

track = Track(
    crns = np.array([True]),
    lpar = np.array([100.0]),
    delTh = np.array([2 * pi]),
)
track.compTrackXY()

# final track
track.plot()
