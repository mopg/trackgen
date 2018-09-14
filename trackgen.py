import numpy as np
from math import *
import matplotlib.pyplot as plt
from scipy.optimize import minimize, Bounds

class Track( object ):
    '''
    Track object holds all parameters defining the track, as well as the
    constraints under which this track was designed.

    Attributes:
        length      Desired length of track
        rmin        Minimum corner radius
        rmax        Maximum corner radius
        lmax        Maximum straight length
        lmin        Minimum straight length
        dthmax      Maximum angle change of corner
        dthmin      Minimum angle change of corner
        left        Orientation of track (Left-turning if True)
        width       Track width
        crns        Track lay-out
        lpar        Optimized length parameters
        delTh       Optimized angle changes
        optimized   Has this track been optimized yet?
    '''

    def __init__( self, length = 500., rmin = 9.,
                  rmax = 50., lmax = 80., lmin = 5., left = True,
                  dthmin = pi/6, dthmax = pi,
                  width = 3., crns = np.zeros( (0,), dtype=bool ) ):

        self.length = length
        self.rmin   = rmin
        self.rmax   = rmax
        self.lmax   = lmax
        self.lmin   = lmin
        self.dthmax = dthmax
        self.dthmin = dthmin
        self.left   = left # track orientation

        self.width = width

        self.crns  = crns
        self.lpar  = np.zeros( np.shape(crns), dtype=float )
        self.delTh = np.zeros( np.shape(crns), dtype=float )

        self.optimized = False

    def solve( self, lpar_init, delTh_init, case = 0 ):
        '''
        Solves the optimization problem that ensures the track has the correct
        length, curvature, etc. using an SLSQP algorithm.
        - Case 0: maximizes curvature
        - Case 1: minimizes curvature
        - Case 2: only satisfies constraints
        '''

        nseg = len(lpar_init)
        assert nseg == len( delTh_init )
        assert nseg == len( self.crns )

        if case > 2:
            raise ValueError('Case number higher than 2')

        # Decide on objective function
        if case == 0:
            print("Maximizing curvature of track")
            fobj = objMaxCurv
        elif case == 1:
            print("Minimizing curvature of track")
            fobj = objMinCurv
        elif case == 2:
            print("No minimization")
            fobj = objNone

        x0 = np.hstack( (lpar_init, delTh_init) )

        # check inputs

        # equality constraints
        constr = {}
        constr['type'] = 'eq'
        constr['fun']  = eqConstr
        # constr['jac']  = eqConstrGrad
        constr['args'] = (self.crns,self.length,self.left)

        # bounds
        lb = np.zeros( x0.shape )
        ub = np.zeros( x0.shape )
        for jj in range(0,nseg):
            if self.crns[jj]:
                lb[jj+nseg] = self.dthmin
                ub[jj+nseg] = self.dthmax
                if lpar_init[jj] > 0.:
                    lb[jj] =  self.rmin
                    ub[jj] =  self.rmax
                else:
                    ub[jj] = -self.rmin
                    lb[jj] = -self.rmax
            else:
                lb[jj]      = self.lmin
                lb[jj+nseg] = 0.
                ub[jj]      = self.lmax
                ub[jj+nseg] = 0.

        bnds = Bounds( lb, ub )

        soldict = minimize( fobj, x0, args=(self.lmax,self.dthmax), method='SLSQP',
                            jac=True, bounds=bnds, constraints=constr,
                            tol=None, options=None )

        print(soldict.message)

        self.lpar  = soldict.x[0:nseg]
        self.delTh = soldict.x[nseg:]

        self.optimized = False

        return soldict

    def plot( self ):
        '''
        Plots the track defined in the Track object.
        '''
        if self.optimized:
            plotTrack( self.crns, self.lpar, self.delTh, self.width )
        else:
            print("First optimize the track!")

    def endpoint( self ):
        '''
        Returns endpoint of the track. If optimization is successful, should be the origin.
        '''
        return compEndpoint( self.crns, self.lpar, self.delTh )

def eqConstr( x, crns, leng, left ):
    '''
    Computes the value of the equality constraints for `x`.
    '''
    constr = np.zeros( (4,) )

    nseg = int( len(x)/2 )

    # length constraint
    constr[0] = compLength( crns, x[0:nseg], x[nseg:] )
    constr[0] -= leng

    # end point constraints and angle constraint
    constr[1], constr[2], constr[3] = compEndpoint( crns, x[0:nseg], x[nseg:] )
    constr[3] -= (-1 + left*2)*2*pi
    return constr

def compLength( crns, lpar, delTh ):
    '''
    Computes final length of track, defined by corner definition `crns`, length
    parameters `lpar`, and angle changes `delTh`.
    Also computes gradient of length with respect to design variables.
    '''
    trlen = 0.

    for jj in range(0,len(crns)):
        if crns[jj]:
            trlen += abs(lpar[jj]) * delTh[jj]
        else:
            trlen += lpar[jj]

    return trlen

def compEndpoint( crns, lpar, delTh ):
    '''
    Computes end point of track, defined by corner definition `crns`,
    length parameters `lpar`, and angle changes `delTh`.
    Also computes gradient with respect to design variables.
    '''

    xend = 0.
    yend = 0.

    thcum = 0.

    xend  = 0.
    yend  = 0.
    thcum = 0.

    for jj in range(0,len(crns)):
        jjm1 = max(jj-1,0)
        if crns[jj]:
            delx =        abs(lpar[jj]) * sin( delTh[jj] ) # local coordinate frame
            dely = lpar[jj] - lpar[jj]  * cos( delTh[jj] ) # local coordinate frame

            # map to global coordinate frame
            xend += delx * cos(thcum) - dely * sin(thcum)
            yend += dely * cos(thcum) + delx * sin(thcum)

            # update cumulative angle
            thcum += np.sign(lpar[jj]) * delTh[jj]
        else:
            xend += lpar[jj] * cos(thcum)
            yend += lpar[jj] * sin(thcum)

    return xend, yend, thcum

def compCurvature( lpar, delTh, lmax, dthmax ):
    '''
    Computes track curvature and the gradient of track curvature with respect
    to the design variables.
    '''

    curv = ( np.linalg.norm(lpar) / lmax )**2 + \
           ( np.linalg.norm(delTh) / dthmax )**2

    dcurvdlpar  = 2*lpar   / lmax**2
    dcurvddelth = 2*delTh  / dthmax**2

    return curv, dcurvdlpar, dcurvddelth

def objMaxCurv( x, lmax, dthmax ):
    '''
    Objective function for maximum curvature.
    '''
    nseg = int( len( x )/2 )

    c, dcdlpar, dcddelth = compCurvature( x[0:nseg], x[nseg:], lmax, dthmax )

    return -c, np.hstack( ( -dcdlpar, -dcddelth ) )

def objMinCurv( x, lmax, dthmax ):
    '''
    Objective function for minimum curvature.
    '''
    nseg = int( len( x )/2 )

    c, dcdlpar, dcddelth = compCurvature( x[0:nseg], x[nseg:], lmax, dthmax )

    return c, np.hstack( ( dcdlpar, dcddelth ) )

def objNone( x, lmax, dthmax ):
    '''
    Constant objective function.
    '''
    return 1., np.zeros( (len(x),) )

def plotTrack( crns, lpar, delTh, width ):
    '''
    Plots the track in x,y-space using matplotlib.
    '''

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

    print( "End point = (%4.3f, %4.3f)" % (xmid[-1], ymid[-1]) )

    plt.plot(xmid,ymid)
    plt.plot(xb1,yb1)
    plt.plot(xb2,yb2)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.show()
