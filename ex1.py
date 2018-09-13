import numpy as np
from math import *
from trackgen import *

# crns = np.array( [False,True,False,True], dtype=bool )
crns = np.array( [False,True,False,True,True,True,False], dtype=bool )

# Change in angle (needs to be zero for straight)
delTh = np.array( [0,pi/4,0,pi/2,pi,pi/4,0], dtype=float )
# delTh = np.array( [0,pi/2,0,pi/2], dtype=float )

# length parameter initial guess (radius for corner, length for straight)
lpar = np.array( [20,10,10,-10,-10,10,20], dtype=float )
# lpar = np.array( [20,-10,20,10], dtype=float )

track = Track( length = 500. )

xe, ye = compEndpoint( crns, lpar, delTh )

print "End point = (%4.3f, %4.3f)" % (xe, ye)

plotTrack( crns, lpar, delTh, 2. )
