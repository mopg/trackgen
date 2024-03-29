import numpy as np
from math import pi
from trackgen import Track

# Where are the corners?
crns = np.array( [False,True,False,True,True,True,False,True,True,False], dtype=bool )

# Change in angle (needs to be zero for straight)
delTh = np.array( [0,pi/2,0,pi/2,pi/2,pi/2,0,pi/4,pi/4,0], dtype=float )

# length parameter initial guess (radius for corner, length for straight)
lpar = np.array( [20,10,20,10,-10,10,200,-10,10,200], dtype=float )

track = Track( length = 250., left = True, crns = crns )

sol = track.solve( lpar, delTh, case = 2 )

xe, ye, thcum = track.endpoint()

print( "End point   = (%4.3f, %4.3f)" % (xe, ye) )
print( "Final angle =  %4.3f" % (thcum) )

# final track
track.plot( cones=True )
