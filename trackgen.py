import numpy as np
from math import *
import matplotlib.pyplot as plt
# from scipy import

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

    def plot( self ):

        plotTrack( self.crns, self.lpar, self.delTh )


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
    dxend_ddelth = np.zeros( np.size(delTh) )
    dyend_dlpar  = np.zeros( np.size(lpar) )
    dyend_ddelth = np.zeros( np.size(delTh) )

    thcum = 0.

    for jj in range(0,len(crns)):
        if crns[jj]:
            delx =        abs(lpar[jj]) * sin( delTh[jj] ) # local coordinate frame
            dely = lpar[jj] - lpar[jj]  * cos( delTh[jj] ) # local coordinate frame

            # map to global coordinate frame
            xend += delx * cos(thcum) - dely * sin(thcum)
            yend += dely * cos(thcum) + delx * sin(thcum)

            # update cumulative angle
            thcum += np.sign(lpar[jj]) * delTh[jj]
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

def plotTrack( crns, lpar, delTh, width ):

    nplot = 50 # number of points used for corners

    nseg  = len(crns)
    ncrns = sum(crns)
    npts  = ncrns * nplot + ( nseg - ncrns ) * 2 + 1

    xmid = np.zeros( (npts,) )
    ymid = np.zeros( (npts,) )

    theta = np.zeros( (npts,) )

    thcum = 0.

    ind = 0

    for jj in range( 0, nseg ):
        if crns[jj]:
            phi = np.linspace( 0., delTh[jj], nplot )

            delx =        abs(lpar[jj]) * np.sin( phi ) # local coordinate frame
            dely = lpar[jj] - lpar[jj]  * np.cos( phi ) # local coordinate frame

            # map to global coordinate frame
            xmid[(ind+1):(ind+nplot+1)] = xmid[ind] + delx * cos(thcum) - dely * sin(thcum)
            ymid[(ind+1):(ind+nplot+1)] = ymid[ind] + dely * cos(thcum) + delx * sin(thcum)

            # update cumulative angle
            thcum += np.sign(lpar[jj]) * delTh[jj]
            theta[(ind+1):(ind+nplot+1)] = theta[ind] + np.sign(lpar[jj]) * phi

            ind += nplot

        else:
            xmid[ind+1] = xmid[ind]
            ymid[ind+1] = ymid[ind]
            xmid[ind+2] = xmid[ind] + lpar[jj] * cos(thcum)
            ymid[ind+2] = ymid[ind] + lpar[jj] * sin(thcum)

            theta[ind+1] = theta[ind]
            theta[ind+2] = theta[ind]

            ind += 2

    xb1 = xmid + width/2 * np.sin( theta )
    yb1 = ymid - width/2 * np.cos( theta )
    xb2 = xmid - width/2 * np.sin( theta )
    yb2 = ymid + width/2 * np.cos( theta )

    print "End point = (%4.3f, %4.3f)" % (xmid[-1], ymid[-1])

    plt.plot(xmid,ymid)
    plt.plot(xb1,yb1)
    plt.plot(xb2,yb2)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
