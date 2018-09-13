import numpy as np
from math import *
from trackgen import *

crns = np.array( [False,True,False,True,True,True,False], dtype=bool )

# Change in angle (needs to be zero for straight)
delTh = np.array( [0,pi/4,0,pi/2,pi,pi/4,0], dtype=float )

# length parameter initial guess (radius for corner, length for straight)
lpar = np.array( [20,10,10,-10,-10,10,20], dtype=float )

track = Track( length = 500., left = True )

xe, ye, thcum, dx_dlp, dx_ddth, dy_dlp, dy_ddth, dth_ddth = compEndpoint( crns, lpar, delTh )

print "End point   = (%4.3f, %4.3f)" % (xe, ye)
print "Final angle =  %4.3f" % (thcum)

print dx_ddth

plotTrack( crns, lpar, delTh, 2. )
