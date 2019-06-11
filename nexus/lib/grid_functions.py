
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
        error('dimension of points must be 2\ndimension of provided points: {}'.format(dim),'polar_to_cartesian')
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
        error('dimension of points must be 2\ndimension of provided points: {}'.format(dim),'cartesian_to_polar')
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


    def setup_mpl_fig(self,fig=True,dim=None,ax1='x',ax2='y',ax3='z'):
        if dim is None:
            self.error('cannot setup mpl figure, "dim" must be provided')
        #end if
        if not fig:
            fig = self.get_cur_fig()
            ax  = self.get_cur_ax()
        else:
            fig = plt.figure()
            if dim==3:
                from mpl_toolkits.mplot3d import Axes3D
                ax = fig.add_subplot(111, projection='3d')
                ax.set_xlabel(ax1)
                ax.set_ylabel(ax2)
                ax.set_zlabel(ax3)
            elif dim==2 or dim==1:
                ax = fig.add_subplot(111)
                ax.set_xlabel(ax1)
                ax.set_ylabel(ax2)
            else:
                self.not_implemented()
            #end if
            self.set_cur_fig(fig)
            self.set_cur_ax(ax)
        #end if
        return fig,ax
    #end def setup_mpl_fig

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
        if len(msgs)==0:
            self.local_validity_checks(msgs)
        #end if
        return msgs
    #end def validity_checks


    def local_validity_checks(self,msgs):
        shape = self.points.shape
        if len(shape)!=2:
            msgs.append('points array must have two shape entries (number of points, dimension of points)\nnumber of shape entries present: {}\npoints shape: {}'.format(len(shape),shape))
        else:
            if shape[0]!=self.npoints:
                msgs.append('unexpected number of points present\npoints expected: {}\npoints encountered: {}'.format(self.npoints,shape[0]))
            #end if
            if shape[1]!=self.space_dim:
                msgs.append('unexpected number of spatial dimensions encountered\nspatial dimensions expected: {}\nspatial dimensions encountered: {}'.format(self.space_dim,shape[1]))
            #end if
        #end if
        if id(self.r)!=id(self.points):
            msgs.append('property function "r" has been overridden\nthis is a developer error')
        #end if
        return msgs
    #end def local_validity_checks


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


    def check_valid_points(self,points,dim,loc):
        if isinstance(points,(tuple,list)):
            points = np.array(points,dtype=self.dtype)
        elif not isinstance(points,np.ndarray):
            self.error('invalid points provided to function "{}"\nprovided points must be an array\ntype provided: {}'.format(loc,points.__class__.__name__))
        #end if
        if len(points.shape)!=2:
            self.error('invalid points provided to function "{}"\npoints must be a 2D array\nnumber of dimensions in provided points array: {}'.format(loc,len(points.shape)))
        elif points.shape[1]!=dim:
            self.error('invalid points provided to function "{}"\npoints must reside in a {}-dimensional space\nspatial dimensions present: {}'.format(loc,dim,points.shape[1]))
        #end if
        return points
    #end def check_valid_points


    def plot_points(self,points=None,fig=True,show=True,default_marker='.',**kwargs):
        fig,ax = self.setup_mpl_fig(fig=fig,dim=self.space_dim)
        if points is None:
            r = self.r.T
        else:
            points = self.check_valid_points(points,self.space_dim,'plot_points')
            r = points.T
        #end if
        if default_marker is not None and 'marker' not in kwargs:
            kwargs['marker'] = default_marker
        #end if
        if self.space_dim!=1:
            ax.scatter(*r,**kwargs)
        else:
            ax.scatter(r,0*r,**kwargs)
        #end if
        ax.set_aspect('equal','box')
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

    valid_bconds = set(['o','p'])

    bcond_types = obj(
        open     = 'o',
        periodic = 'p',
        )


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
        if not isinstance(bconds,(tuple,list,np.ndarray)):
            self.error('cannot set bconds from data with type "{}"\nplease use tuple, list, or array for inputted bconds'.format(bconds.__class__.__name__))
        #end if
        bconds = np.array(bconds,dtype=object)
        for bc in bconds:
            if bc not in StructuredGrid.valid_bconds:
                self.error('boundary conditions are invalid\nboundary conditions in each dimension must be one of: {}\nboundary conditions provided: {}'.format(sorted(StructuredGrid.valid_bconds),bconds))
            #end if
        #end for
        self.bconds = bconds
    #end def set_bconds


    def local_validity_checks(self,msgs):
        msgs = Grid.local_validity_checks(self,msgs)
        if np.prod(self.shape)!=self.npoints:
            msgs.append('grid shape does not match number of points\nnumber of points: {}\nproduct of shape: {}\nshape: {}'.format(self.npoints,np.prod(self.shape),self.shape))
        #end if
        if len(self.shape)!=self.grid_dim:
            msgs.append('number of entries in grid shape does not match grid dimension\nnumber of entries in grid shape: {}\ngrid dimension: {}'.format(len(self.shape),self.grid_dim))
        #end if
        if len(set(self.bconds)-StructuredGrid.valid_bconds)>0:
            msgs.append('boundary conditions are invalid\nboundary conditions in each dimension must be one of: {}\nboundary conditions provided: {}'.format(sorted(StructuredGrid.valid_bconds),self.bconds))
            #end if
        #end for
        return msgs
    #end def local_validity_checks


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


    def inside(self,points,tol=1e-12):
        points = self.check_valid_points(points,self.space_dim,'inside')
        upoints = self.unit_points(points)
        inside = np.ones((len(points),),dtype=bool)
        for d in range(self.grid_dim):
            if self.bconds[d]!=self.bcond_types.periodic:
                u = upoints[:,d]
                inside &= ( (u > -tol) & (u < 1.+tol))
            #end if
        #end for
        return inside
    #end def inside


    def project(self,points):
        points = self.check_valid_points(points,self.space_dim,'project')
        upoints = self.unit_points(points)
        for d in range(self.grid_dim):
            if self.bconds[d]==self.bcond_types.periodic:
                upoints[:,d] -= np.floor(upoints[:,d])
            #end if
        #end for
        return self.points_from_unit(upoints)
    #end def project


    def get_boundary_lines(self,n=200,unit=False):
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
        if not unit:
            bpoints = self.points_from_unit(upoints)
            bpoints.shape = nlines,n,self.space_dim
        else:
            bpoints = upoints
            bpoints.shape = nlines,n,self.grid_dim
        #end if
        return bpoints
    #end def get_boundary_lines


    def plot_boundary(self,n=200,fig=True,show=True):
        fig,ax = self.setup_mpl_fig(fig=fig,dim=self.space_dim)
        if self.grid_dim!=1:
            bpoints = self.get_boundary_lines(n=n)
            for bp in bpoints:
                ax.plot(*bp.T,color='k')
            #end for
        #end if
        ax.set_aspect('equal','box')
        if show:
            plt.show()
        #end if
    #end def plot_boundary


    def plot_unit_points(self,points=None,fig=True,show=True,default_marker='.',**kwargs):
        fig,ax = self.setup_mpl_fig(fig=fig,dim=self.grid_dim,
                                    ax1='a1',ax2='a2',ax3='a3')
        if points is None:
            u = self.unit_points().T
        else:
            points = self.check_valid_points(points,self.space_dim,'plot_unit_points')
            u = self.unit_points(points=points).T
        #end if
        if default_marker is not None and 'marker' not in kwargs:
            kwargs['marker'] = default_marker
        #end if
        if self.grid_dim!=1:
            ax.scatter(*u,**kwargs)
        else:
            ax.scatter(u,0*u,**kwargs)
        #end if
        ax.set_aspect('equal','box')
        if show:
            plt.show()
        #end if
    #end def plot_unit_points


    def plot_unit_boundary(self,n=200,fig=True,show=True):
        fig,ax = self.setup_mpl_fig(fig=fig,dim=self.grid_dim,
                                    ax1='a1',ax2='a2',ax3='a3')
        if self.grid_dim!=1:
            bpoints = self.get_boundary_lines(n=n,unit=True)
            for bp in bpoints:
                ax.plot(*bp.T,color='k')
            #end for
        #end if
        ax.set_aspect('equal','box')
        if show:
            plt.show()
        #end if
    #end def plot_unit_boundary


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


    def local_validity_checks(self,msgs):
        msgs = StructuredGrid.local_validity_checks(self,msgs)
        shape = self.axes.shape
        if len(shape)!=2:
            msgs.append('axes must be a 2 dimensional array\nnumber of dimensions present: {}\naxes present: {}'.format(len(shape),self.axes))
        else:
            if shape[0]!=self.grid_dim:
                msgs.append('number of axes must be equal to the embedded grid dimension\nembedded grid dimension: {}\nnumber of axes present: {}'.format(self.grid_dim,shape[0]))
            #end if
            if shape[1]!=self.space_dim:
                msgs.append('axis dimension must be equal to the dimension of the space\nspace dimension: {}\naxis dimension: {}\naxes present: {}'.format(self.space_dim,shape[1],self.axes))
            #end if
        #end if
        shape = self.origin.shape
        if len(shape)!=1:
            msgs.append('origin must be a 1 dimensional array (single point)\nnumber of dimensions present: {}\norigin present: {}'.format(len(shape),self.origin))
        elif shape[0]!=self.space_dim:
            msgs.append('origin dimension must be equal to the dimension of the space\nspace dimension: {}\naxis dimension: {}\norigin present: {}'.format(self.space_dim,shape[0],self.origin))
        #end if
        return msgs
    #end def local_validity_checks
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

        # reassign bconds to match sphere constraints
        if self.grid_dim==3:
            self.set_bconds(tuple('oop'))
        elif self.grid_dim==2:
            self.set_bconds(tuple('op'))
        else:
            self.error('one dimensional spheroid is not supported')
        #end if
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


    demos = obj(
        plot_grids      = 0,
        plot_inside     = 0,
        plot_projection = 1,
        )

    shapes = {
        1 : (10,),
        2 : (10,10),
        3 : (10,10,10),
        }

    axes = {
        # 1d grid in 1d space
        (1,1) : [[1]],
        # 1d grid in 2d space
        (1,2) : [[1,1]],
        # 1d grid in 3d space
        (1,3) : [[1,1,1]],
        # 2d grid in 2d space
        (2,2) : [[1,0],[1,1]],
        # 2d grid in 3d space
        (2,3) : [[1,0,0],[1,1,1]],
        # 3d grid in 3d space
        (3,3) : [[1,0,0],[1,1,0],[1,1,1]],
        }

    bconds = {
        1 : 'o p'.split(),
        2 : 'oo op po pp'.split(),
        3 : 'ooo oop opo poo opp pop ppo ppp'.split(),
        }

    grid_types = obj(
        parallelotope = ParallelotopeGrid,
        spheroid      = SpheroidGrid,
        )

    supported = obj(
        parallelotope = obj(dims=set(axes.keys())),
        spheroid      = obj(dims=set([(2,2),(2,3),(3,3)])),
        )

    grids = obj()
    gdict = dict(parallelotope='p',spheroid='s')
    cdict = {False:'',True:'c'}
    for grid_name in sorted(grid_types.keys()):
        label = gdict[grid_name]+'e'
        grids[label] = grid_types[grid_name]()
        for grid_dim,space_dim in sorted(axes.keys()):
            if (grid_dim,space_dim) in supported[grid_name].dims:
                for centered in (False,True):
                    label = gdict[grid_name]+str(grid_dim)+str(space_dim)+cdict[centered]
                    grid_inputs = obj(
                        shape    = shapes[grid_dim],
                        axes     = axes[grid_dim,space_dim],
                        centered = centered,
                        )

                    if grid_name!='parallelotope':
                        g = grid_types[grid_name](**grid_inputs)
                        grids[label] = g
                    else:
                        for bc in bconds[grid_dim]:
                            label_bc = label+'_'+bc
                            grid_inputs_bc = obj(
                                bconds = tuple(bc),
                                **grid_inputs
                                )
                            g = grid_types[grid_name](**grid_inputs_bc)
                            grids[label_bc] = g
                            if 'p' not in bc:
                                grids[label] = g
                            #end if
                        #end for
                    #end if

                #end for
            #end if
        #end for
    #end for

    for label in sorted(grids.keys()):
        g = grids[label]
        print ' {:<16}  {}'.format(label,g.bconds)
    #end for

    if demos.plot_grids:
        grids_plot = 'p11 p12 p13 p23 p23c p33 p33c'.split()
        #grids_plot = 's23 s23c s33 s33c'.split()

        unit = False

        for name in grids_plot:
            grid = grids[name]
            if not unit:
                grid.plot_points(show=0)
                grid.plot_boundary(fig=0,show=0)
            else:
                grid.plot_unit_points(show=0)
                grid.plot_unit_boundary(fig=0,show=0)
            #end if
            ax = grid.get_cur_ax()
            plt.title(name)
        #end for
        plt.show()
    #end if


    if demos.plot_inside:
        #grids_plot = 'p22_oo p22_op p22_pp'.split()
        grids_plot = 's22 s23'.split()

        n = 100

        unit = False

        upoints = dict()
        for d in range(1,3+1):
            upoints[d] = 0.75*np.random.randn(n,d)+0.75
        #end for

        for name in grids_plot:
            g = grids[name]
            dim = int(name[1])
            points = g.points_from_unit(upoints[dim])
            inside = g.inside(points)
            if not unit:
                g.plot_points(points[inside],color='b',show=0)
                g.plot_points(points[~inside],color='r',fig=0,show=0)
                g.plot_boundary(fig=0,show=0)
            else:
                g.plot_unit_points(points[inside],color='b',show=0)
                g.plot_unit_points(points[~inside],color='r',fig=0,show=0)
                g.plot_unit_boundary(fig=0,show=0)
            #end if
            ax = g.get_cur_ax()
            plt.title(name)
        #end for
        plt.show()
    #end if


    if demos.plot_projection:
        grids_plot = 'p22_oo p22_op p22_pp'.split()

        n = 100

        unit = False

        upoints = dict()
        for d in range(1,3+1):
            upoints[d] = 0.75*np.random.randn(n,d)+0.75
        #end for

        for name in grids_plot:
            g = grids[name]
            dim = int(name[1])
            points = g.points_from_unit(upoints[dim])
            proj_points = g.project(points)
            
            inside = g.inside(points)
            if not unit:
                g.plot_points(points,color='r',marker='o',facecolors='none',show=0)
                g.plot_points(proj_points,color='b',fig=0,show=0)
                g.plot_boundary(fig=0,show=0)
            else:
                g.plot_unit_points(points,color='r',marker='o',facecolors='none',show=0)
                g.plot_unit_points(proj_points,color='b',fig=0,show=0)
                g.plot_unit_boundary(fig=0,show=0)
            #end if
            ax = g.get_cur_ax()
            plt.title(name)
        #end for
        plt.show()
    #end if

#end if 
