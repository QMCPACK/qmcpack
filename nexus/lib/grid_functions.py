
from generic import obj
from developer import DevBase,ci,error,unavailable

try:
    import numpy as np
except:
    np = unavailable('numpy','np')
#end try
try:
    import matplotlib.pyplot as plt
except:
    plt = unavailable('matplotlib','plt')
#end try



## n-dimensional version, not using due to axis order conventions
## https://en.wikipedia.org/wiki/N-sphere
#def spherical_to_cartesian(points):
#    if not isinstance(points,np.ndarray):
#        points = np.array(points)
#    #end if
#    sphere_points = points
#    npoints,dim = sphere_points.shape
#    cart_points = []
#    r = sphere_points[:,0]
#    s = np.sin(sphere_points[:,1:])
#    c = np.cos(sphere_points[:,1:])
#    for d in range(dim):
#        cart = r.copy()
#        if d>0:
#            cart *= s[:,:d].prod(axis=1)
#        #end if
#        if d<dim-1:
#            cart *= c[:,d]
#        #end if
#        cart_points.append(cart)
#    #end for
#    cart_points = np.array(cart_points).T
#    return cart_points
##end def spherical_to_cartesian


def polar_to_cartesian(points):
    if not isinstance(points,np.ndarray):
        points = np.array(points)
    #end if
    npoints,dim = points.shape
    if dim!=2:
        error('dimension of points must be 2','polar_to_cartesian')
    #end if
    r   = points[:,0]
    phi = points[:,1] # range is [0,2*pi)
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    cart_points = np.array((x,y),dtype=points.dtype).T
    return cart_points
#end def polar_to_cartesian


def cartesian_to_polar(points):
    if not isinstance(points,np.ndarray):
        points = np.array(points)
    #end if
    npoints,dim = points.shape
    if dim!=2:
        error('dimension of points must be 2','cartesian_to_polar')
    #end if
    x = points[:,0]
    y = points[:,1]
    r   = np.linalg.norm(points,axis=1)
    phi = np.arctan2(y,x)
    phi[np.abs(phi)<1e-12] = 0.0
    phi[phi<0] += 2*np.pi
    pol_points = np.array((r,phi),dtype=points.dtype).T
    return pol_points
#end def cartesian_to_polar


def spherical_to_cartesian(points):
    if not isinstance(points,np.ndarray):
        points = np.array(points)
    #end if
    npoints,dim = points.shape
    if dim!=3:
        error('dimension of points must be 3','spherical_to_cartesian')
    #end if
    r     = points[:,0]
    theta = points[:,1] # range is [0,pi)
    phi   = points[:,2] # range is [0,2*pi)
    s = np.sin(theta)
    x = r*s*np.cos(phi)
    y = r*s*np.sin(phi)
    z = r*np.cos(theta)
    cart_points = np.array((x,y,z),dtype=points.dtype).T
    return cart_points
#end def spherical_to_cartesian


def cartesian_to_spherical(points):
    if not isinstance(points,np.ndarray):
        points = np.array(points)
    #end if
    npoints,dim = points.shape
    if dim!=3:
        error('dimension of points must be 3','cartesian_to_spherical')
    #end if
    x = points[:,0]
    y = points[:,1]
    z = points[:,2]
    r     = np.linalg.norm(points,axis=1)
    theta = np.arccos(z/r)
    phi   = np.arctan2(y,x)
    phi[np.abs(phi)<1e-12] = 0.0
    phi[phi<0] += 2*np.pi
    sphere_points = np.array((r,theta,phi),dtype=points.dtype).T
    return sphere_points
#end def cartesian_to_spherical




def linear_grid_1d(rmin=0.0,rmax=1.0,n=None,dr=None,centered=False):
    if n is not None and dr is not None:
        error('provide only one of "dr" and "n", not both','linear_grid_1d')
    elif dr is not None:
        n = int(round((rmax-rmin)/dr))
    elif n is None:
        error('either "n" or "dr" must be provided','linear_grid_1d')
    #end if
    r = np.linspace(rmin,rmax,n,endpoint=False)
    if centered:
        r += 0.5/n
    #end if
    return r
#end def linear_grid_1d



def unit_grid_points(shape,centered=False):
    linear_grids = []
    for n in shape:
        lin_grid = np.linspace(0.0,1.0,n,endpoint=False)
        if centered:
            lin_grid += 0.5/n
        #end for
        linear_grids.append(lin_grid)
    #end for
    # equivalent to 
    #   points = ndgrid(*linear_grids)
    points = np.meshgrid(*linear_grids,indexing='ij')
    points = np.array(points)
    # reshape and transpose the points
    points.shape = (len(points),np.array(shape).prod())
    points = points.T
    return points
#end def unit_grid_points



def parallelotope_grid_points(axes,shape=None,dr=None,centered=False,return_shape=False,return_axes=False):
    if not isinstance(axes,np.ndarray):
        axes = np.array(axes)
    #end if
    if shape is not None and dr is not None:
        error('provide only one of "dr" and "shape", not both','parallelotope_grid_points')
    elif dr is not None:
        if not isinstance(dr,np.ndarray):
            dr = np.array(dr,dtype=float)
        #end if
        shape = np.array(np.around(np.linalg.norm(axes,axis=1)/dr),dtype=int)
    elif shape is None:
        error('either "shape" or "dr" must be provided','parallelotope_grid_points')
    #end if
    ugrid = unit_grid_points(shape,centered=centered)
    pgrid = np.dot(ugrid,axes)
    if not return_shape and not return_axes:
        return pgrid
    else:
        res = [pgrid]
        if return_shape:
            res.append(shape)
        #end if
        if return_axes:
            res.append(axes)
        #end if
        return res
    #end if
#end def parallelotope_grid_points



def spheroid_grid_points(axes,shape,centered=False):
    if not isinstance(axes,np.ndarray):
        axes = np.array(axes)
    #end if
    grid_dim,space_dim = axes.shape
    if grid_dim not in (2,3):
        error('spheroid grid generation only supported in 2 or 3 dimensions','spheriod_grid_points')
    #end if
    ugrid = unit_grid_points(shape,centered=centered)
    if grid_dim==2:
        ugrid[:,1] *= 2*np.pi
        sgrid = polar_to_cartesian(ugrid)
    elif grid_dim==3:
        ugrid[:,1] *=   np.pi
        ugrid[:,2] *= 2*np.pi
        sgrid = spherical_to_cartesian(ugrid)
    #end if
    sgrid = np.dot(sgrid,axes) # adds radial range and skew
    return sgrid
#end def spheroid_grid_points



def all_none(*vals):
    all = True
    for v in vals:
        all &= v is None
    #end for
    return all
#end def all_none


def any_none(*vals):
    any = False
    for v in vals:
        any |= v is None
    #end for
    return any
#end def any_none



class PlotHandler(DevBase):
    
    @staticmethod
    def reset():
        PlotHandler.fig = None
        PlotHandler.ax  = None
    #end def reset

    def set_cur_fig(self,fig):
        PlotHandler.fig = fig
    #end def set_cur_fig

    def get_cur_fig(self):
        return PlotHandler.fig
    #end def get_cur_fig

    def set_cur_ax(self,ax):
        PlotHandler.ax = ax
    #end def set_cur_ax

    def get_cur_ax(self):
        return PlotHandler.ax
    #end def get_cur_ax
#end class PlotHandler
PlotHandler.reset()



class GBase(PlotHandler):
    None
#end class GBase



class Grid(GBase):

    persistent_data_types = obj(
        points      = (np.ndarray,None ),
        initialized = (bool      ,False),
        )

    @property
    def r(self): # points must be accessed this way
        return self.points
    #end def r

    @property
    def npoints(self):
        return self.points.shape[0]
    #end def npoints

    @property
    def space_dim(self):
        return self.points.shape[1]
    #end def space_dim

    @property
    def dtype(self):
        return self.points.dtype
    #end def dtype


    def __init__(self,**kwargs):

        self.reset()

        if len(kwargs)>0:
            self.initialize(**kwargs)
        #end if
    #end def __init__


    def reset(self):
        cls = self.__class__
        for name,(dtype,default) in cls.persistent_data_types.items():
            self[name] = default
        #end for
    #end def reset


    def initialize(self,**kwargs):
        # remove check argument
        check = kwargs.pop('check',True)

        # call derived initialize
        self.initialize_local(**kwargs)

        # record initialization action
        self.initialized = True

        # check validity of grid
        if check:
            self.check_valid()
        #end if
    #end def initialize


    def initialize_local(self,points=None,dtype=None,copy=True):
        if points is None:
            self.error('cannot initialize grid, "points" is required')
        #end if
        self.set_points(points,dtype=dtype,copy=copy)
    #end def initialize_local


    def set_points(self,points,dtype=None,copy=True):
        if isinstance(points,(tuple,list)):
            if dtype is None:
                dtype = float
            #end if
            points = np.array(points,dtype=dtype)
        elif isinstance(points,np.ndarray):
            if dtype is not None:
                points = np.array(points,dtype=dtype)
            elif copy:
                points = points.copy()
            #end if
        else:
            self.error('cannot set points from data with type "{}"\nplease use tuple, list, or array for inputted points'.format(points.__class__.__name__))
        #end if
        self.points = points
    #end def set_points


    def copy(self,shallow=False):
        if not shallow:
            c = DevBase.copy(self)
        else:
            points = self.points
            del self.points
            c = DevBase.copy(self)
            c.points = points
            self.points = points
        #end if
        return c
    #end def copy


    def translate(self,shift):
        self.points += shift
    #end def translate


    def validity_checks(self):
        cls = self.__class__
        msgs = []
        # check that the grid has been initialized
        if not self.initialized:
            msgs.append('grid has not been initialized')
            return msgs
        #end if
        valid_key_list = list(cls.persistent_data_types.keys())
        valid_key_set  = set(valid_key_list)
        key_set = set(self.keys())
        # check that all data members are present
        missing = valid_key_set - key_set
        if len(missing)>0:
            msgs.append('some data members are missing: {}'.format(sorted(missing)))
        #end if
        # check that extra data members are not present
        extra = key_set - valid_key_set
        if len(extra)>0:
            msgs.append('unknown data members encountered: {}'.format(sorted(extra)))
        #end if
        # check that all data members are of the correct type
        for name in valid_key_list:
            type = cls.persistent_data_types[name]
            if name in self:
                val = self[name]
                if val is None:
                    msgs.append('data member "{}" has not been initialized, it is "None"'.format(name))
                elif not isinstance(val,type):
                    msgs.append('data member "{}" is not of type "{}", it is of type "{}" instead'.format(name,type.__name__,val.__class__.__name__))
                #end if
            #end if
        #end for
        return msgs
    #end def validity_checks


    def check_valid(self,exit=True):
        msgs = self.validity_checks()
        valid = len(msgs)==0
        if not valid and exit:
            self.log('\n')
            for msg in msgs:
                self.log('  '+msg)
            #end for
            self.error('grid is not valid, see messages above')
        #end if
        return valid
    #end def check_valid


    def valid(self):
        return self.check_valid(exit=False)
    #end def valid


    def plot_points(self,fig=True,show=True):
        if fig:
            fig = plt.figure()
            self.set_cur_fig(fig)
        #end if
        r = self.r.T
        if self.space_dim==3:
            from mpl_toolkits.mplot3d import Axes3D
            if fig:
                ax = fig.add_subplot(111, projection='3d')
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_zlabel('z')
                self.set_cur_ax(ax)
            else:
                ax = self.get_cur_ax()
            #end if
            ax.scatter(*r,marker='.')

        else:
            self.not_implemented()
        #end if
        if show:
            plt.show()
        #end if
    #end def plot_points
#end class Grid



class StructuredGrid(Grid):

    persistent_data_types = obj(
        shape    = (tuple     ,None),
        centered = (bool      ,None),
        bconds   = (np.ndarray,None),
        **Grid.persistent_data_types
        )


    @property
    def grid_dim(self):
        return len(self.shape)
    #end def grid_dim

    @property
    def grid_shape(self):
        return self.shape
    #end def grid_shape

    valid_bconds = set(['o','p'])


    def initialize_local(self,**kwargs):
        shape    = kwargs.pop('shape'   ,None)
        centered = kwargs.pop('centered',False)
        bconds   = kwargs.pop('bconds'  ,None)
        Grid.initialize_local(self,**kwargs)
        self.centered = centered
        if shape is None:
            self.error('cannot initialize grid, "shape" is required')
        #end if
        self.set_shape(shape)
        if bconds is None:
            bconds = len(self.shape)*['o']
        #end if
        self.set_bconds(bconds)
    #end def initialize_local


    def set_shape(self,shape):
        if not isinstance(shape,(tuple,list,np.ndarray)):
            self.error('cannot set shape from data with type "{}"\nplease use tuple, list, or array for inputted shape'.format(shape.__class__.__name__))
        #end if
        self.shape = tuple(shape)
    #end def set_shape


    def set_bconds(self,bconds):
        bconds = np.array(bconds,dtype=object)
        for bc in bconds:
            if bc not in StructuredGrid.valid_bconds:
                self.error('boundary conditions are invalid\nboundary conditions in each dimension must be one of: {}\nboundary conditions provided: {}'.format(sorted(StructuredGrid.valid_bconds),bconds))
            #end if
        #end for
        self.bconds = bconds
    #end def set_bconds


    def unit_points(self,points=None):
        if points is None:
            points = self.r
        #end if
        upoints = self.unit_points_bare(points)
        return upoints
    #end def unit_points


    def unit_indices(self,points=None):
        upoints = self.unit_points(points)
        shape = np.array(self.shape,dtype=int)
        ipoints = upoints*shape
        if not self.centered:
            ipoints += 0.5
        #end if
        ipoints = np.array(np.floor(ipoints),dtype=int)
        dim = self.grid_dim
        cum_shape = np.empty((dim,),dtype=int)
        cum_shape[-1] = 1
        for d in range(1,dim):
            cum_shape[dim-d-1] = cum_shape[dim-d]*shape[dim-d-1]
        #end for
        ipoints = np.dot(ipoints,cum_shape)
        return ipoints
    #end def unit_indices


    def plot_boundary(self,n=200,fig=True,show=True):
        if fig:
            fig = plt.figure()
            self.set_cur_fig(fig)
        #end if
        u = np.linspace(0.,1.,n)
        nlines = self.grid_dim*2**(self.grid_dim-1)
        upoints = np.empty((n*nlines,self.grid_dim),dtype=self.dtype)
                
        if self.grid_dim==3:
            bitset = [[0,0],[0,1],[1,0],[1,1]]
        elif self.grid_dim==2:
            bitset = [[0],[1]]
        #end if
        ni = 0
        for d in range(self.grid_dim):
            for bits in bitset:
                k = 0
                for d2 in range(self.grid_dim):
                    if d2!=d:
                        upoints[ni*n:(ni+1)*n,d2] = bits[k]
                        k+=1
                    else:
                        upoints[ni*n:(ni+1)*n,d2] = u
                    #end if
                #end for
                ni+=1
            #end for
        #end for
        bpoints = self.points_from_unit(upoints)
        bpoints.shape = nlines,n,self.space_dim
        if self.space_dim==3:
            from mpl_toolkits.mplot3d import Axes3D
            if fig:
                ax = fig.add_subplot(111, projection='3d')
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_zlabel('z')
                self.set_cur_ax(ax)
            else:
                ax = self.get_cur_ax()
            #end if

            for bp in bpoints:
                ax.plot(*bp.T,color='k')
            #end for

        else:
            self.not_implemented()
        #end if
        if show:
            plt.show()
        #end if
    #end def plot_boundary


    def unit_points_bare(self,points=None):
        self.not_implemented()
    #end def unit_points_bare


    def points_from_unit(self,upoints):
        self.not_implemented()
    #end def points_from_unit

#end class StructuredGrid



class StructuredGridWithAxes(StructuredGrid):

    persistent_data_types = obj(
        axes   = (np.ndarray,None),
        origin = (np.ndarray,None),
        **StructuredGrid.persistent_data_types
        )
    
    def initialize_local(self,
                         axes   = None,
                         origin = None,
                         **kwargs
                         ):
        if axes is None:
            self.error('cannot initialize grid, "axes" is required')
        #end if

        StructuredGrid.initialize_local(self,**kwargs)

        if origin is None:
            origin = self.space_dim*[0]
        #end if
        self.set_axes(axes)
        self.set_origin(origin)
    #end def initialize_local


    def set_axes(self,axes,shift_points=False):
        if not isinstance(axes,(tuple,list,np.ndarray)):
            self.error('cannot set axes from data with type "{}"\nplease use tuple, list, or array for inputted axes'.format(axes.__class__.__name__))
        #end if
        self.axes = np.array(axes,dtype=self.dtype)
    #end def set_axes


    def set_origin(self,origin):
        if not isinstance(origin,(tuple,list,np.ndarray)):
            self.error('cannot set origin from data with type "{}"\nplease use tuple, list, or array for inputted origin'.format(origin.__class__.__name__))
        #end if
        origin = np.array(origin,dtype=self.dtype)
        if self.origin is not None:
            shift = origin-self.origin
        else:
            shift = origin
        #end if
        self.origin = origin
        self.translate(shift)
    #end def set_origin
#end class StructuredGridWithAxes



class ParallelotopeGrid(StructuredGridWithAxes):

    @property
    def corner(self):
        return self.origin
    #end def corner

    @property
    def center(self):
        return self.corner + self.axes.sum(axis=0)/2
    #end def center


    def initialize_local(self,
                 axes     = None,
                 shape    = None,
                 dr       = None,
                 corner   = None,
                 center   = None,
                 centered = False,
                 **kwargs
                 ):
        if shape is None and dr is None:
            self.error('cannot initialize grid, either "shape" or "dr" is required')
        #end if

        points,shape,axes = parallelotope_grid_points(axes,shape=shape,dr=dr,centered=centered,return_shape=True,return_axes=True)

        if center is not None:
            center = np.array(center)
            corner = center - axes.sum(axis=0)/2
        #end if

        kwargs['axes']     = axes
        kwargs['origin']   = corner
        kwargs['shape']    = shape
        kwargs['centered'] = centered
        kwargs['points']   = points
        StructuredGridWithAxes.initialize_local(self,**kwargs)

    #end def initialize_local


    def unit_points_bare(self,points):
        corner = self.corner
        # invert using pseudo-inverse
        #   this is important for grids embedded in higher dim spaces
        axinv  = np.linalg.pinv(self.axes)
        upoints = np.dot(points-corner,axinv)
        return upoints
    #end def unit_points_bare


    def points_from_unit(self,upoints):
        points = np.dot(upoints,self.axes)+self.corner
        return points
    #end def points_from_unit

#end class ParallelotopeGrid



class SpheroidGrid(StructuredGridWithAxes):

    @property
    def center(self):
        return self.origin
    #end def center

    def initialize_local(self,
                         axes     = None,
                         shape    = None,
                         center   = None,
                         centered = False,
                         **kwargs
                         ):
        if shape is None:
            self.error('cannot initialize grid, "shape" is required')
        #end if

        points = spheroid_grid_points(axes,shape=shape,centered=centered)

        kwargs['axes']     = axes
        kwargs['origin']   = center
        kwargs['shape']    = shape
        kwargs['centered'] = centered
        kwargs['points']   = points
        StructuredGridWithAxes.initialize_local(self,**kwargs)
    #end def initialize_local


    def unit_points_bare(self,points):
        center = self.center
        # invert using pseudo-inverse
        #   this is important for grids embedded in higher dim spaces
        axinv  = np.linalg.pinv(self.axes)
        upoints = np.dot(points-center,axinv)
        # map from unit cartesian to unit spherical
        dim = self.grid_dim
        if dim==2:
            upoints = cartesian_to_polar(upoints)
            upoints[:,1] /= 2*np.pi
        elif dim==3:
            upoints = cartesian_to_spherical(upoints)
            upoints[:,1] /=   np.pi
            upoints[:,2] /= 2*np.pi
        else:
            self.error('unit_points not supported for dim={}'.format(dim))
        #end if
        return upoints
    #end def unit_points_bare


    def points_from_unit(self,upoints):
        dim = self.grid_dim
        upoints = np.array(upoints,dtype=self.dtype)
        if dim==2:
            upoints[:,1] *= 2*np.pi
            upoints = polar_to_cartesian(upoints)
        elif dim==3:
            upoints[:,1] *=   np.pi
            upoints[:,2] *= 2*np.pi
            upoints = spherical_to_cartesian(upoints)
        else:
            self.error('points_from_unit not supported for dim={}'.format(dim))
        #end if
        points = np.dot(upoints,self.axes)+self.center
        return points
    #end def points_from_unit
#end class SpheroidGrid




if __name__=='__main__':

    #p = ParallelotopeGrid(
    #    shape = (10,10),
    #    axes  = [[1,0,0],[1,1,1]],
    #    )

    pe = ParallelotopeGrid()

    pc = ParallelotopeGrid(
        shape    = (10,10),
        axes     = [[1,0,0],[1,1,1]],
        centered = True
        )

    se = SpheroidGrid()

    s = SpheroidGrid(
        shape = (10,10),
        axes  = [[1,0,0],[1,1,1]],
        )
    #sc = SpheroidGrid(
    #    shape    = (10,10),
    #    axes     = [[1,0,0],[1,1,1]],
    #    centered = True,
    #    )

    o = obj(s)
    del o.points
    print repr(o)
    print o
    print s.valid()

    ci()
#end if 
