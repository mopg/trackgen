import numpy as np
from math import *
import scipy

class Track( object ):
    '''
    Track object holds all parameters defining the track, as well as the
    constraints under which this track was designed.
    '''

    def __init__( self, length = 500., rmin = 9., rmax = 50., lmax = 80.,
                  width = 3., crns = np.zeros( (0,), dtype=bool ) ):
        self.length = length
        self.rmin   = rmin
        self.rmax   = rmax
        self.lmax   = lmax

        self.width  = width

        self.crns   = crns
        self.lpar   = np.zeros( np.size(crns), dtype=float )
        self.delTh  = np.zeros( np.size(crns), dtype=float )

    def solve( self, lpar_init, delTh_init ):

        # check inputs

        # constrain length
        # - length
        # - end angle
        # end point
        return 0


def compLength( crns, lpar, delTh ):
    '''
    Computes final length of track, defined by corner definition `crns`, length
    parameters `lpar`, and angle changes `delTh`.
    Also computes gradient of length with respect to design variables.
    '''

    trlen = 0.

    # set up arrays for gradient info
    dtrl_dlpar  = np.zeros( np.size(lpar) )
    dtrl_ddelth = np.zeros( np.size(delth) )

    for jj in range(0,len(crns)):
        if crns[jj]:
            trlen += lpar[jj] * delTh[jj]
        else:
            trlen += lpar[jj]

    return trlen

def compEndpoint( crns, lpar, delTh ):
    '''
    Computes end point of track, defined by corner definition `crns`, length
    parameters `lpar`, and angle changes `delTh`.
    Also computes gradient of end point with respect to design variables.
    Eventually the end point should become the origin.
    '''

    xend = 0.
    yend = 0.

    # set up arrays for gradient info
    dxend_dlpar  = np.zeros( np.size(lpar) )
    dxend_ddelth = np.zeros( np.size(delth) )
    dyend_dlpar  = np.zeros( np.size(lpar) )
    dyend_ddelth = np.zeros( np.size(delth) )

    thcum = 0.

    for jj in range(0,len(crns)):
        if crns[jj]:
            delx =        abs(lpar[jj]) * cos( delTh[jj] ) # local coordinate frame
            dely = lpar[jj] - lpar[jj]  * sin( delTh[jj] ) # local coordinate frame

            # map to global coordinate frame
            xend += delx * cos(thcum) - dely * sin(thcum)
            yend += dely * cos(thcum) + delx * sin(thcum)

            # update cumulative angle
            thcum += delTh[jj]
        else:
            xend  += lpar[jj] * cos(thcum)
            yend  += lpar[jj] * sin(thcum)

    return (xend, yend)


def compCurvature( delTh ):
    '''
    Computes track curvature and the gradient of track curvature with respect
    to the design variables.
    '''

    return np.norm(delTh)

def plotTrack( ):

    return 0
