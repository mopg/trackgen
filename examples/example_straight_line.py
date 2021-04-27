import numpy as np
from math import pi
from trackgen import Track

# we do not need to optimize this track, as it does not need to be closed

track = Track(
    crns = np.array([False]),
    lpar = np.array([500.0]),
    delTh = np.array([0.0]),
)
track.compTrackXY()

# final track
track.plot()
