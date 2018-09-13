import numpy as np
from math import *
from trackgen import *

crns = np.ones( (7,), dtype=bool )

# Change in angle (needs to be zero for straight)
delTh = np.array( [0,pi/4,0,pi/2,pi,pi/4,0], dtype=float )

# length parameter initial guess (radius for corner, length for straight)
lpar = np.array( [20,10,10,10,10,10,20], dtype=float )

track = Track( length = 500. )
