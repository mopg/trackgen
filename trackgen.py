import numpy as np
from math import *
import matplotlib.pyplot as plt
# from scipy.optimize import minimize

class Track( object ):
    '''
    Track object holds all parameters defining the track, as well as the
    constraints under which this track was designed.
    '''

    def __init__( self, length = 500., rmin = 9.,
                  rmax = 50., lmax = 80., left = True,
                  width = 3., crns = np.zeros( (0,), dtype=bool ) ):

        self.length = length
        self.rmin   = rmin
        self.rmax   = rmax
        self.lmax   = lmax
        self.left   = left # track orientation

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

        plotTrack( self.crns, self.lpar, self.delTh, self.width )

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
            trlen += abs(lpar[jj]) * delTh[jj]
            dtrl_dlpar[jj]  = delTh[jj]
            dtrl_ddelth[jj] = abs(lpar[jj])
        else:
            trlen += lpar[jj]
            dtrl_dlpar[jj] = 1.

    return trlen, dtrl_dlpar, dtrl_ddelth

def compEndpoint( crns, lpar, delTh ):
    '''
    Computes end point of track, defined by corner definition `crns`,
    length parameters `lpar`, and angle changes `delTh`.
    Also computes gradient with respect to design variables
    '''

    xend = 0.
    yend = 0.

    # set up arrays for gradient info
    nseg = len( lpar )
    dxend_dlpar   = np.zeros( (nseg,) )
    dxend_ddelth  = np.zeros( (nseg,) )
    dyend_dlpar   = np.zeros( (nseg,) )
    dyend_ddelth  = np.zeros( (nseg,) )

    dthcum_ddelth = np.zeros( (nseg,) ) # no adjoint necessary for this

    thcum = 0.

    # for gradient computation, need to keep intermediate values
    xend  = np.zeros( (nseg,) )
    yend  = np.zeros( (nseg,) )
    thcum = np.zeros( (nseg,) )

    for jj in range(0,len(crns)):
        jjm1 = max(jj-1,0)
        if crns[jj]:
            delx =        abs(lpar[jj]) * sin( delTh[jj] ) # local coordinate frame
            dely = lpar[jj] - lpar[jj]  * cos( delTh[jj] ) # local coordinate frame

            # map to global coordinate frame
            xend[jj] = xend[jjm1] + delx * cos(thcum[jjm1]) - dely * sin(thcum[jjm1])
            yend[jj] = yend[jjm1] + dely * cos(thcum[jjm1]) + delx * sin(thcum[jjm1])

            # update cumulative angle
            thcum[jj] = thcum[jjm1] + np.sign(lpar[jj]) * delTh[jj]
        else:
            xend[jj]  = xend[jjm1] + lpar[jj] * cos(thcum[jjm1])
            yend[jj]  = yend[jjm1] + lpar[jj] * sin(thcum[jjm1])
            thcum[jj] = thcum[jjm1]

    # compute gradient through adjoint
    # NOTE: these are all separate (three different adjoints)
    dxend  = 1.
    dyend  = 1.
    dthcum_x = 0. # this is only used to track the adjoint of dxend
    dthcum_y = 0. # this is only used to track the adjoint of dyend
    for jj in range(len(crns)-1,0,-1):
        dthcumm1_x = 0.
        dthcumm1_y = 0.
        if crns[jj]:

            ## gradient of dthcum_ddelth (no adjoint used)
            dthcum_ddelth[jj] = np.sign(lpar[jj])

            ## Recompute
            delx =        abs(lpar[jj]) * sin( delTh[jj] ) # local coordinate frame
            dely = lpar[jj] - lpar[jj]  * cos( delTh[jj] ) # local coordinate frame

            ## Start adjoint computation

            # map to global coordinate frame
            # xend[jj] = xend[jjm1] + delx * cos(thcum[jjm1]) - dely * sin(thcum[jjm1])
            # yend[jj] = yend[jjm1] + dely * cos(thcum[jjm1]) + delx * sin(thcum[jjm1])
            ddelx_x =   dxend * cos(thcum[jjm1])
            ddely_x = - dyend * sin(thcum[jjm1])
            ddelx_y =   dxend * sin(thcum[jjm1])
            ddely_y =   dyend * cos(thcum[jjm1])
            dthcumm1_x = - sin(thcum[jjm1]) * delx * dxend - dely * cos(thcum[jjm1]) * dxend
            dthcumm1_y = - sin(thcum[jjm1]) * dely * dyend + delx * cos(thcum[jjm1]) * dyend

            # update cumulative angle
            # thcum[jj] = thcum[jjm1] + np.sign(lpar[jj]) * delTh[jj]
            # dthcum_ddelth[jj] = np.sign(lpar[jj]) * dthcum # NOTE: already taken care of
            dthcumm1_y += dthcum_x
            dthcumm1_y += dthcum_y

            # local updates
            # delx =        abs(lpar[jj]) * sin( delTh[jj] )
            # dely = lpar[jj] - lpar[jj]  * cos( delTh[jj] )
            dxend_ddelth[jj] = cos(delTh[jj]) * abs(lpar[jj]) * ddelx_x - \
                              sin(delTh[jj]) * lpar[jj] * ddely_x
            dyend_ddelth[jj] = cos(delTh[jj]) * abs(lpar[jj]) * ddelx_y - \
                              sin(delTh[jj]) * lpar[jj] * ddely_y
            dxend_dlpar[jj]  = ( sin( delTh[jj] ) * np.sign( lpar[jj] ) ) * ddelx_x + \
                              ( 1. - cos( delTh[jj] ) ) * ddely_x
            dyend_dlpar[jj]  = ( sin( delTh[jj] ) * np.sign( lpar[jj] ) ) * ddelx_y + \
                              ( 1. - cos( delTh[jj] ) ) * ddely_y

        else:
            # thcum[jj] = thcum[jjm1]
            dthcumm1_x += dthcum_x
            dthcumm1_y += dthcum_y

            # xend[jj]  = xend[jjm1] + lpar[jj] * cos(thcum[jjm1])
            dthcumm1_x     += -lpar[jj] * sin(thcum[jjm1]) * dxend
            dxend_dlpar[jj] += cos(thcum[jjm1]) * dxend

            # yend[jj]  = yend[jjm1] + lpar[jj] * sin(thcum[jjm1])
            dthcumm1_x     += lpar[jj] * cos(thcum[jjm1]) * dyend
            dyend_dlpar[jj] += sin(thcum[jjm1]) * dyend

        dthcum_x = dthcumm1_x # remember we're doing this in reverse
        dthcum_y = dthcumm1_y # remember we're doing this in reverse

    return ( xend[-1], yend[-1], thcum[-1],
             dxend_dlpar, dxend_ddelth, dyend_dlpar,
             dyend_ddelth, dthcum_ddelth )

def compCurvature( delTh ):
    '''
    Computes track curvature and the gradient of track curvature with respect
    to the design variables.
    '''

    return np.norm(delTh)**2, 2*delTh

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
    plt.axis('equal')
    plt.show()
