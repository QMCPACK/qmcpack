import sys

from numpy import *
from numpy.linalg import *

from developer import unavailable
try:
    from scipy.special import betainc
    from scipy.optimize import fmin
    from scipy.spatial import KDTree,Delaunay
    scipy_unavailable = False
except ImportError:
    betainc = unavailable('scipy.special' ,'betainc')
    fmin    = unavailable('scipy.optimize','fmin')
    KDTree  = unavailable('scipy.special' ,'KDTree')
    Delaunay  = unavailable('scipy.special' ,'Delaunay')
    scipy_unavailable = True
#end try


########################################################################
############   ndgrid
########################################################################

# retrieved from
#   http://www.mailinglistarchive.com/html/matplotlib-users@lists.sourceforge.net/2010-05/msg00055.html

#"""
#n-dimensional gridding like Matlab's NDGRID
#
#Typical usage:
#>>> x, y, z = [0, 1], [2, 3, 4], [5, 6, 7, 8]
#>>> X, Y, Z = ndgrid(x, y, z)
#
#See ?ndgrid for details.
#"""

def ndgrid(*args, **kwargs):
    """
    n-dimensional gridding like Matlab's NDGRID
    
    The input *args are an arbitrary number of numerical sequences, 
    e.g. lists, arrays, or tuples.
    The i-th dimension of the i-th output argument 
    has copies of the i-th input argument.
    
    Optional keyword argument:
    same_dtype : If False (default), the result is an ndarray.
                 If True, the result is a lists of ndarrays, possibly with 
                 different dtype. This can save space if some *args 
                 have a smaller dtype than others.

    Typical usage:
    >>> x, y, z = [0, 1], [2, 3, 4], [5, 6, 7, 8]
    >>> X, Y, Z = ndgrid(x, y, z) # unpacking the returned ndarray into X, Y, Z

    Each of X, Y, Z has shape [len(v) for v in x, y, z].
    >>> X.shape == Y.shape == Z.shape == (2, 3, 4)
    True
    >>> X
    array([[[0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]],
           [[1, 1, 1, 1],
            [1, 1, 1, 1],
            [1, 1, 1, 1]]])
    >>> Y
    array([[[2, 2, 2, 2],
            [3, 3, 3, 3],
            [4, 4, 4, 4]],
           [[2, 2, 2, 2],
            [3, 3, 3, 3],
            [4, 4, 4, 4]]])
    >>> Z
    array([[[5, 6, 7, 8],
            [5, 6, 7, 8],
            [5, 6, 7, 8]],
           [[5, 6, 7, 8],
            [5, 6, 7, 8],
            [5, 6, 7, 8]]])
    
    With an unpacked argument list:
    >>> V = [[0, 1], [2, 3, 4]]
    >>> ndgrid(*V) # an array of two arrays with shape (2, 3)
    array([[[0, 0, 0],
            [1, 1, 1]],
           [[2, 3, 4],
            [2, 3, 4]]])
    
    For input vectors of different data types, same_dtype=False makes ndgrid()
    return a list of arrays with the respective dtype.
    >>> ndgrid([0, 1], [1.0, 1.1, 1.2], same_dtype=False)
    [array([[0, 0, 0], [1, 1, 1]]), 
     array([[ 1. ,  1.1,  1.2], [ 1. ,  1.1,  1.2]])]
    
    Default is to return a single array.
    >>> ndgrid([0, 1], [1.0, 1.1, 1.2])
    array([[[ 0. ,  0. ,  0. ], [ 1. ,  1. ,  1. ]],
           [[ 1. ,  1.1,  1.2], [ 1. ,  1.1,  1.2]]])
    """
    same_dtype = kwargs.get("same_dtype", True)
    V = [array(v) for v in args] # ensure all input vectors are arrays
    shape = [len(v) for v in args] # common shape of the outputs
    result = []
    for i, v in enumerate(V):
        # reshape v so it can broadcast to the common shape
        # http://docs.scipy.org/doc/numpy/user/basics.broadcasting.html
        zero = zeros(shape, dtype=v.dtype)
        thisshape = ones_like(shape)
        thisshape[i] = shape[i]
        result.append(zero + v.reshape(thisshape))
    if same_dtype:
        return array(result) # converts to a common dtype
    else:
        return result # keeps separate dtype for each output

#if __name__ == "__main__":
#    import doctest
#    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)


########################################################################
############   End ndgrid
########################################################################


def simstats(x,dim=None):
    shape = x.shape
    ndim  = len(shape)
    if dim==None:
        dim=ndim-1
    #end if
    permute = dim!=ndim-1
    reshape = ndim>2
    nblocks = shape[dim]
    if permute:
        r = range(ndim)
        r.pop(dim)
        r.append(dim)
        permutation = tuple(r)
        r = range(ndim)
        r.pop(ndim-1)
        r.insert(dim,ndim-1)
        invperm     = tuple(r)
        x=x.transpose(permutation)
        shape = tuple(array(shape)[array(permutation)])
        dim = ndim-1
    #end if
    if reshape:        
        nvars = prod(shape[0:dim])
        x=x.reshape(nvars,nblocks)
        rdim=dim
        dim=1
    else:
        nvars = shape[0]
    #end if

    mean  = x.mean(dim)
    var   = x.var(dim)

    N=nblocks

    if ndim==1:
        i=0          
        tempC=0.5
        kappa=0.0
        mtmp=mean
        if abs(var)<1e-15:
            kappa = 1.0
        else:
            ovar=1.0/var
            while (tempC>0 and i<(N-1)):
                kappa=kappa+2.0*tempC
                i=i+1
                #tempC=corr(i,x,mean,var)
                tempC = ovar/(N-i)*sum((x[0:N-i]-mtmp)*(x[i:N]-mtmp))
            #end while
            if kappa == 0.0:
                kappa = 1.0
            #end if
        #end if
        Neff=(N+0.0)/(kappa+0.0)
        if (Neff == 0.0):
            Neff = 1.0
        #end if
        error=sqrt(var/Neff)
    else:
        error = zeros(mean.shape)
        kappa = zeros(mean.shape)
        for v in xrange(nvars):
            i=0          
            tempC=0.5
            kap=0.0
            vtmp = var[v]
            mtmp = mean[v]
            if abs(vtmp)<1e-15:
                kap = 1.0
            else:
                ovar   = 1.0/vtmp
                while (tempC>0 and i<(N-1)):
                    i += 1
                    kap += 2.0*tempC
                    tempC = ovar/(N-i)*sum((x[v,0:N-i]-mtmp)*(x[v,i:N]-mtmp))
                #end while
                if kap == 0.0:
                    kap = 1.0
                #end if
            #end if
            Neff=(N+0.0)/(kap+0.0)
            if (Neff == 0.0):
                Neff = 1.0
            #end if
            kappa[v]=kap
            error[v]=sqrt(vtmp/Neff)
        #end for    
    #end if

    if reshape:
        x     =     x.reshape(shape)
        mean  =  mean.reshape(shape[0:rdim])
        var   =   var.reshape(shape[0:rdim])
        error = error.reshape(shape[0:rdim])
        kappa = kappa.reshape(shape[0:rdim])
    #end if
    if permute:
        x=x.transpose(invperm)
    #end if

    return (mean,var,error,kappa)
#end def simstats



def simplestats(x,dim=None):
    if dim==None:
        dim=len(x.shape)-1
    #end if
    osqrtN = 1.0/sqrt(1.0*x.shape[dim])
    mean   = x.mean(dim)
    error  = x.var(dim)*osqrtN
    return (mean,error)
#end def simplestats


def equilibration_length(x,tail=.5,plot=False,xlim=None,bounces=2):
    bounces = max(1,bounces)
    eqlen = 0
    nx = len(x)
    xt = x[int((1.-tail)*nx+.5):]
    nxt = len(xt)
    if nxt<10:
        return eqlen
    #end if
    #mean  = xh.mean()
    #sigma = sqrt(xh.var())
    xs = array(xt)
    xs.sort()
    mean  = xs[int(.5*(nxt-1)+.5)]
    sigma = (abs(xs[int((.5-.341)*nxt+.5)]-mean)+abs(xs[int((.5+.341)*nxt+.5)]-mean))/2
    crossings = bounces*[0,0]
    if abs(x[0]-mean)>sigma:
        s = -sign(x[0]-mean)
        ncrossings = 0
        for i in range(nx):
            dist = s*(x[i]-mean) 
            if dist>sigma and dist<5*sigma:
                crossings[ncrossings]=i
                s*=-1
                ncrossings+=1
                if ncrossings==2*bounces:
                    break
                #end if
            #end if
        #end for
        bounce = crossings[-2:]
        bounce[1] = max(bounce[1],bounce[0])
        #print len(x),crossings,crossings[1]-crossings[0]+1
        eqlen = bounce[0]+random.randint(bounce[1]-bounce[0]+1)
    #end if
    if plot:
        xlims = xlim
        del plot,xlim
        from matplotlib.pyplot import plot,figure,show,xlim
        figure()
        ix = arange(nx)
        plot(ix,x,'b.-')
        plot([0,nx],[mean,mean],'k-')
        plot([0,nx],[mean+sigma,mean+sigma],'r-')
        plot([0,nx],[mean-sigma,mean-sigma],'r-')
        plot(ix[crossings],x[crossings],'r.')
        plot(ix[bounce],x[bounce],'ro')
        plot([ix[eqlen],ix[eqlen]],[x.min(),x.max()],'g-')
        plot(ix[eqlen],x[eqlen],'go')
        if xlims!=None:
            xlim(xlims)
        #end if
        show()
    #end if
    return eqlen
#end def equilibration_length


def ttest(m1,e1,n1,m2,e2,n2):
    m1 = float(m1)
    e1 = float(e1)
    m2 = float(m2)
    e2 = float(e2)
    v1 = e1**2
    v2 = e2**2
    t  = (m1-m2)/sqrt(v1+v2)
    nu = (v1+v2)**2/(v1**2/(n1-1)+v2**2/(n2-1))
    x = nu/(nu+t**2)
    p = 1.-betainc(nu/2,.5,x)
    return p
#end def ttest



def surface_normals(x,y,z):
    nu,nv = x.shape
    normals = empty((nu,nv,3))
    mi=nu-1
    mj=nv-1
    v1 = empty((3,))
    v2 = empty((3,))
    v3 = empty((3,))
    dr = empty((3,))
    dr[0] = x[0,0]-x[1,0]
    dr[1] = y[0,0]-y[1,0]
    dr[2] = z[0,0]-z[1,0]
    drtol = 1e-4
    for i in xrange(nu):
        for j in xrange(nv):
            iedge = i==0 or i==mi
            jedge = j==0 or j==mj
            if iedge:
                dr[0] = x[0,j]-x[mi,j]
                dr[1] = y[0,j]-y[mi,j]
                dr[2] = z[0,j]-z[mi,j]
                if norm(dr)<drtol:
                    im = mi-1
                    ip = 1
                elif i==0:
                    im=i
                    ip=i+1
                elif i==mi:
                    im=i-1
                    ip=i
                #end if
            else:
                im=i-1
                ip=i+1
            #end if
            if jedge:
                dr[0] = x[i,0]-x[i,mj]
                dr[1] = y[i,0]-y[i,mj]
                dr[2] = z[i,0]-z[i,mj]
                if norm(dr)<drtol:
                    jm = mj-1
                    jp = 1
                elif j==0:
                    jm=j
                    jp=j+1
                elif j==mj:
                    jm=j-1
                    jp=j
                #end if
            else:
                jm=j-1
                jp=j+1
            #end if
            v1[0] = x[ip,j]-x[im,j]
            v1[1] = y[ip,j]-y[im,j]
            v1[2] = z[ip,j]-z[im,j]
            v2[0] = x[i,jp]-x[i,jm]
            v2[1] = y[i,jp]-y[i,jm]
            v2[2] = z[i,jp]-z[i,jm]
            v3 = cross(v1,v2)
            onorm = 1./norm(v3)
            normals[i,j,:]=v3[:]*onorm
        #end for
    #end for
    return normals
#end def surface_normals


simple_surface_coords = [set(['x','y','z']),set(['r','phi','z']),set(['r','phi','theta'])]
simple_surface_min = {'x':-1.00000000001,'y':-1.00000000001,'z':-1.00000000001,'r':-0.00000000001,'phi':-0.00000000001,'theta':-0.00000000001}
def simple_surface(origin,axes,grid):
    matched=False
    gk = set(grid.keys())
    for c in range(3):
        if gk==simple_surface_coords[c]:
            matched=True
            coord=c
        #end if
    #end for
    if not matched:
        print 'Error in simple_surface: invalid coordinate system provided'
        print '  provided coordinates:',gk
        print '  permitted coordinates:'
        for c in range(3):
            print '               ',simple_surface_coords[c]
        #end for
        exit
    #end if
    for k,v in grid.iteritems():
        if min(v)<simple_surface_min[k]:
            print 'Error in simple surface: '+k+' cannot be less than '+str(simple_surface_min[k])
            print '   actual minimum: '+str(min(v))
            sys.exit()
        #end if
        if max(v)>1.00000000001:
            print 'Error in simple surface: '+k+' cannot be more than 1'
            print '   actual maximum: '+str(max(v))
            sys.exit()
        #end if
    #end if
    u=empty((3,))
    r=empty((3,))
    if coord==0:
        xl = grid['x']
        yl = grid['y']
        zl = grid['z']
        dim = (len(xl),len(yl),len(zl))
        npoints = prod(dim)
        points = empty((npoints,3))
        n=0
        for i in xrange(dim[0]):
            for j in xrange(dim[1]):
                for k in xrange(dim[2]):
                    r[0] = xl[i]
                    r[1] = yl[j]
                    r[2] = zl[k]
                    points[n,:] = dot(axes,r) + origin
                    n+=1
                #end for
            #end for
        #end for
    elif coord==1:
        rl   = grid['r']
        phil = 2.*pi*grid['phi']
        zl   = grid['z']
        dim = (len(rl),len(phil),len(zl))
        npoints = prod(dim)
        points = empty((npoints,3))
        n=0
        for i in xrange(dim[0]):
            for j in xrange(dim[1]):
                for k in xrange(dim[2]):
                    u[0] = rl[i]
                    u[1] = phil[j]
                    u[2] = zl[k]
                    r[0] = u[0]*cos(u[1])
                    r[1] = u[0]*sin(u[1])
                    r[2] = u[2]
                    points[n,:] = dot(axes,r) + origin
                    n+=1
                #end for
            #end for
        #end for
    elif coord==2:
        rl     = grid['r']
        phil   = 2.*pi*grid['phi']
        thetal = pi*grid['theta']
        dim = (len(rl),len(phil),len(thetal))
        if dim[0]==1:
            sgn = -1. #this is to 'fix' surface normals
            #sgn = 1. #this is to 'fix' surface normals
        else:
            sgn = 1.
        #end if
        npoints = prod(dim)
        points = empty((npoints,3))
        n=0
        for i in xrange(dim[0]):
            for j in xrange(dim[1]):
                for k in xrange(dim[2]):
                    u[0] = rl[i]
                    u[1] = phil[j]
                    u[2] = thetal[k]
                    r[0] = sgn*u[0]*sin(u[2])*cos(u[1])
                    r[1] = sgn*u[0]*sin(u[2])*sin(u[1])
                    r[2] = sgn*u[0]*cos(u[2])
                    points[n,:] = dot(axes,r) + origin
                    n+=1
                #end for
            #end for
        #end for
    #end if

    if min(dim)!=1:
        print 'Error in simple_surface: minimum dimension must be 1'
        print '   actual minimum dimension:',str(min(dim))
        sys.exit()
    #end if

    dm = []
    for d in dim:
        if d>1:
            dm.append(d)
        #end if
    #end for
    dm=tuple(dm)
     
    x = points[:,0].reshape(dm)
    y = points[:,1].reshape(dm)
    z = points[:,2].reshape(dm)

    return x,y,z
#end def simple_surface



least_squares = lambda p,x,y,f: ((f(p,x)-y)**2).sum()

def func_fit(x,y,fitting_function,p0,minimizer=least_squares):
    f = fitting_function
    p = fmin(minimizer,p0,args=(x,y,f),maxiter=10000,maxfun=10000)
    return p
#end def func_fit


def distance_table(p1,p2,ordering=0):
    n1 = len(p1)
    n2 = len(p2)
    if not isinstance(p1,ndarray):
        p1=array(p1)
    #end if
    if not isinstance(p2,ndarray):
        p2=array(p2)
    #end if
    dt = zeros((n1,n2))
    for i1 in xrange(n1):
        for i2 in xrange(n2):
            dt[i1,i2] = norm(p1[i1]-p2[i2])
        #end for
    #end for
    if ordering==0:
        return dt
    else:
        if ordering==1:
            n=n1
        elif ordering==2:
            n=n2
            dt=dt.T
        else:
            print 'distance_table Error: ordering must be 1 or 2,\n  you provided '+str(ordering)+'\nexiting.'
            exit()
        #end if
        order = empty(dt.shape,dtype=int)
        for i in xrange(n):
            o = dt[i].argsort()
            order[i] = o
            dt[i,:]  = dt[i,o]
        #end for
        return dt,order
    #end if
#end def distance_table



def nearest_neighbors(n,points,qpoints=None,return_distances=False,slow=False):
    extra = 0
    if qpoints==None:
        qpoints=points
        if len(points)>1:
            extra=1
        elif return_distances:
            return array([]),array([])
        else:
            return array([])
        #end if
    #end if
    if n>len(qpoints)-extra:
        print 'nearest_neighbors Error: requested more than the total number of neighbors\n  maximum is: {0}\n  you requested: {1}\nexiting.'.format(len(qpoints)-extra,n)
        exit()
    #end if
    slow = slow or scipy_unavailable
    if not slow:
        kt = KDTree(points)
        dist,ind = kt.query(qpoints,n+extra)
    else:
        dtable,order = distance_table(points,qpoints,ordering=2)
        dist = dtable[:,0:n+extra]
        ind  = order[:,0:n+extra]
    #end if
    if extra==0 and n==1 and not slow:
        nn = atleast_2d(ind).T
    else:
        nn = ind[:,extra:]
    #end if
    if not return_distances:
        return nn
    else:
        return nn,dist
    #end if
#end def nearest_neighbors


def convex_hull(points,dimension=None,tol=None):
    if dimension is None:
        np,dimension = points.shape
    #end if
    d1 = dimension+1
    tri = Delaunay(points)
    all_inds = empty((d1,),dtype=bool)
    all_inds[:] = True
    verts = []
    have_tol = tol!=None
    for ni in range(len(tri.neighbors)):
        n = tri.neighbors[ni]
        ns = list(n)
        if -1 in ns:
            i = ns.index(-1)
            inds = all_inds.copy()
            inds[i] = False
            v = tri.vertices[ni]
            if have_tol:
                iv = range(d1)
                iv.pop(i)
                c = points[v[iv[1]]]
                a = points[v[i]]-c
                b = points[v[iv[0]]]-c
                bn = norm(b)
                d = norm(a-dot(a,b)/(bn*bn)*b)
                if d<tol:
                    inds[i]=True
                #end if
            #end if
            verts.extend(v[inds])
        #end if
    #end for
    verts = list(set(verts))
    return verts
#end def convex_hull
