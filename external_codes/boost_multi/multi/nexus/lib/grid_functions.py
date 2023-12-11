##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################

"""
The :py:mod:`grid_functions` module provides generic capabilities to 
handle function data that are represented on discrete grids of points.  

Spatial grids currently supported by the classes in this module include 
sheared parallelotope (simulation cell), spherical, and spherical surface 
grids.  Grids of arbitrary dimension are supported within parallelotope 
domains while 2D (disk) and 3D (ball) grids are supported for 
spherical domains.  Grids over `N`-sphere surfaces are supported in 1D 
(circles) and 2D (spheres). Grids of lower dimension may reside in 
spaces of higher dimension, for example lines, planar plaquettes, and 
disks may all reside in a 3D space.  Both periodic and open boundary 
conditions are supported for parallelotope grids.

Discrete functions are defined as sets of values over these grids.  
Supported functions may be scalar valued, vector valued, or tensor valued. 
For each of these cases, this module aims to support plotting, 
interpolation, integration, and differentiation.

The main classes intended for instantiation and use are 
:py:class:`ParallelotopeGridFunction`, :py:class:`SpheroidGridFunction`, 
and :py:class:`SpheroidSurfaceGridFunction`.  The point grid classes 
corresponding to these grid functions may also be instantiated and used 
directly.  See :py:class:`ParallelotopeGrid`, :py:class:`SpheroidGrid`, 
and :py:class:`SpheroidSurfaceGrid`. 


List of module contents
-----------------------

Coordinate conversion functions:

* :py:func:`polar_to_cartesian`
* :py:func:`cartesian_to_polar`
* :py:func:`spherical_to_cartesian`
* :py:func:`cartesian_to_spherical`

Grid generation functions:

* :py:func:`unit_grid_points`
* :py:func:`parallelotope_grid_points`
* :py:func:`spheroid_grid_points`
* :py:func:`spheroid_surface_grid_points`

Base classes for :py:class:`Grid` and :py:class:`GridFunction` classes:

* :py:class:`PlotHandler`
* :py:class:`GBase`

Abstract :py:class:`Grid` classes:

* :py:class:`Grid`
* :py:class:`StructuredGrid`
* :py:class:`StructuredGridWithAxes`

Abstract :py:class:`GridFunction` classes:

* :py:class:`GridFunction`
* :py:class:`StructuredGridFunction`
* :py:class:`StructuredGridFunctionWithAxes`

Concrete :py:class:`Grid` classes:

* :py:class:`ParallelotopeGrid` 
* :py:class:`SpheroidGrid`
* :py:class:`SpheroidSurfaceGrid`

Concrete :py:class:`GridFunction` classes:

* :py:class:`ParallelotopeGridFunction` 
* :py:class:`SpheroidGridFunction`
* :py:class:`SpheroidSurfaceGridFunction`


Module contents
---------------
"""

import os

from generic import obj
from developer import DevBase,ci,message,error,unavailable
from fileio import StandardFile,XsfFile

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
try:
    from skimage import measure
except:
    measure = unavailable('skimage','measure')
#end try
try:
    import scipy.ndimage as scipy_ndimage
except:
    scipy_ndimage = unavailable('scipy','ndimage')
#end try




def polar_to_cartesian(points,surface=False):
    """
    Conversion from polar to Cartesian coordinates.

    Parameters
    ----------
    points  : `array_like, float, shape (N,d)`
        Real valued points in polar coordinates :math:`(r,\phi)`. `N` is the 
        number of points and `d` is the dimension of the coordinate system.  
        The inputted points must satisfy :math:`r=\mathrm{points[:,0]}`, 
        :math:`\phi=\mathrm{points[:,1]}`, and :math:`\phi\in[0,2\pi)`.
    surface : `bool, optional, default False`
        Points lie only on the boundary (a circle) or not.  If `False` (the 
        default), the inputted points are two-dimensional (`d=2`) with 
        :math:`(r,\phi)` provided.  If `True`, the inputted points are angular 
        only (`d=1` with :math:`\phi=\mathrm{points[:,0]}`). In this case, 
        :math:`r=1`.

    Returns
    -------
    cart_points : `ndarray, float, shape (N,2)`
        Real valued points in Cartesian coordinates.
    """
    if not isinstance(points,np.ndarray):
        points = np.array(points)
    #end if
    npoints,dim = points.shape
    if dim!=2-int(surface):
        error('dimension of points must be {}\ndimension of provided points: {}'.format(2-int(surface),dim),'polar_to_cartesian')
    #end if
    if surface:
        r   = 1.0
        phi = points[:,0]
    else:
        r   = points[:,0]
        phi = points[:,1] # range is [0,2*pi)
    #end if
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    cart_points = np.array((x,y),dtype=points.dtype).T
    return cart_points
#end def polar_to_cartesian


def cartesian_to_polar(points,surface=False):
    """
    Conversion from Cartesian to polar coordinates.

    Parameters
    ----------
    points  : `array_like, float, shape (N,2)`
        Real valued points in Cartesian coordinates :math:`(x, y)`. `N` is 
        the number of points and :math:`x=\mathrm{points[:,0]}`, 
        :math:`y=\mathrm{points[:,1]}`. 
    surface : `bool, optional, default False`
        Inputted points lie only on a circle or not.  It is the user's 
        responsibility to guarantee the correctness of this assertion.
        If `False` (the default), the outputted points are two-dimensional 
        (`d=2`) with :math:`(r,\phi)` returned.  If `True`, only :math:`\phi` 
        is returned (`d=1`).

    Returns
    -------
    pol_points : `ndarray, float, shape (N,d)`
        Real valued points in polar coordinates.
    """
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
    if surface:
        pol_points = np.array((phi,),dtype=points.dtype).T
    else:
        pol_points = np.array((r,phi),dtype=points.dtype).T
    #end if
    return pol_points
#end def cartesian_to_polar


def spherical_to_cartesian(points,surface=False):
    """
    Conversion from spherical to Cartesian coordinates.

    Parameters
    ----------
    points  : `array_like, float, shape (N,d)`
        Real valued points in spherical coordinates :math:`(r,\\theta,\phi)`. 
        `N` is the number of points, `d` is the dimension of the coordinate 
        system.  The inputted points must satisfy 
        :math:`r=\mathrm{points[:,0]}`, :math:`\\theta=\mathrm{points[:,1]}`, 
        :math:`\phi=\mathrm{points[:,1]}`, with :math:`\\theta\in[0,\pi)`, 
        :math:`\phi\in[0,2\pi)`.
    surface : `bool, optional, default False`
        Points lie only on the boundary (a sphere) or not. If `False` (the 
        default), the inputted points are 3D (`d=3`) with 
        :math:`(r,\\theta,\phi)` provided.  If `True`, the inputted points are 
        angular only (`d=2` with :math:`\\theta=\mathrm{points[:,0]}`, 
        :math:`\phi=\mathrm{points[:,0]}`). In this case, :math:`r=1`.

    Returns
    -------
    cart_points : `ndarray, float, shape (N,3)`
        Real valued points in Cartesian coordinates.
    """
    if not isinstance(points,np.ndarray):
        points = np.array(points)
    #end if
    npoints,dim = points.shape
    if dim!=3-int(surface):
        error('dimension of points must be {}\ndimension of provided points: '.format(req_dim,dim),'spherical_to_cartesian')
    #end if
    if surface:
        r     = 1.0
        theta = points[:,0] # range is [0,pi)
        phi   = points[:,1] # range is [0,2*pi)
    else:
        r     = points[:,0]
        theta = points[:,1] # range is [0,pi)
        phi   = points[:,2] # range is [0,2*pi)
    #end if
    s = np.sin(theta)
    x = r*s*np.cos(phi)
    y = r*s*np.sin(phi)
    z = r*np.cos(theta)
    cart_points = np.array((x,y,z),dtype=points.dtype).T
    return cart_points
#end def spherical_to_cartesian


def cartesian_to_spherical(points,surface=False):
    """
    Conversion from Cartesian to spherical coordinates.

    Parameters
    ----------
    points  : `array_like, float, shape (N,3)`
        Real valued points in Cartesian coordinates :math:`(x,y,z)`. `N` is 
        the number of points and :math:`x=\mathrm{points[:,0]}`, 
        :math:`y=\mathrm{points[:,1]}`, :math:`z=\mathrm{points[:,2]}`. 
    surface : `bool, optional, default False`
        Inputted points lie only on a sphere or not.  It is the user's 
        responsibility to guarantee the correctness of this assertion.
        If `False` (the default), the outputted points are 3D (`d=3`) with 
        :math:`(r,\\theta,\phi)` returned.  If `True`, only 
        :math:`(\\theta,\phi)` is returned (`d=2`).

    Returns
    -------
    sphere_points : `ndarray, float, shape (N,d)`
        Real valued points in spherical coordinates.
    """
    if not isinstance(points,np.ndarray):
        points = np.array(points)
    #end if
    npoints,dim = points.shape
    if dim!=3:
        error('dimension of points must be 3\ndimension of provided points: {}'.format(dim),'cartesian_to_spherical')
    #end if
    x = points[:,0]
    y = points[:,1]
    z = points[:,2]
    r     = np.linalg.norm(points,axis=1)
    theta = np.zeros(r.shape,dtype=points.dtype)
    rfinite = r>1e-12
    theta[rfinite] = np.arccos(z[rfinite]/r[rfinite])
    phi   = np.arctan2(y,x)
    phi[np.abs(phi)<1e-12] = 0.0
    phi[phi<0] += 2*np.pi
    if surface:
        sphere_points = np.array((theta,phi),dtype=points.dtype).T
    else:
        sphere_points = np.array((r,theta,phi),dtype=points.dtype).T
    #end if
    return sphere_points
#end def cartesian_to_spherical



def unit_grid_points(shape,centered=False,endpoint=None):
    """
    Generation of uniform grids in N dimensions.

    Parameters
    ----------
    shape : `array_like, int`
        Number of points in the grid in each dimension.  The dimension of 
        the grid is the number of entries in `shape`.
    centered : `bool, optional, default False`
        Locate grid points at lower cell corners (`False`) or cell centers 
        (`True`).
    endpoint : `array_like, bool, optional`
        If `True` for given dimension, add an endpoint at the upper edge 
        of the grid in that dimension (default `False`).  Applies only to 
        non-centered grids. `shape` and `endpoint` must have the same number 
        of entries.

    Returns
    -------
    points : `ndarray, float, shape (prod(shape),len(shape))`
        Array containing all grid points.
    """
    linear_grids = []
    if endpoint is None:
        endpoint = len(shape)*[False]
    #end if
    for n,ep in zip(shape,endpoint):
        ep &= not centered
        lin_grid = np.linspace(0.0,1.0,n,endpoint=ep)
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



def parallelotope_grid_points(axes,
                              shape        = None,
                              cells        = None,
                              dr           = None,
                              centered     = False,
                              endpoint     = None,
                              return_shape = False,
                              return_axes  = False
                              ):
    """
    Generation of uniform grids within parallelotope volumes.

    Parameters
    ----------
    axes : `array_like, float, shape (d,d)`
        Axis vectors defining the cell (parallelotope).
    shape : `array_like, int, optional`
        Number of points in the grid in each dimension.  The dimension of 
        the grid is the number of entries in `shape`.  Either `shape`,  
        `cells`, or `dr` must be provided.
    cells : `array_like, int, optional`
        Number of cells in the grid in each dimension.  The dimension of the 
        grid is the number of entries in `cells`.  Either `shape`, `cells`, 
        or `dr` must be provided.
    dr : `array_like, float, optional`
        Width of grid cells in each dimension.  The dimension of the grid 
        is the number of entries in `dr`.  Either `shape`, `cells`, or 
        `dr` must be provided.
    centered : `bool, optional, default False`
        Locate grid points at lower cell corners (`False`) or cell centers 
        (`True`).
    endpoint : `array_like, bool, optional`
        If `True` for given dimension, add an endpoint at the upper edge 
        of the grid in that dimension (default `False`).  Applies only to 
        non-centered grids. `shape/cells/dr` and `endpoint` must have the 
        same number of entries.
    return_shape : `bool, optional, default False`
        Additionally return the shape of the grid.
    return_axes : `bool, optional, default False`
        Additionally return the cell axes as an `ndarray`.

    Returns
    -------
    pgrid : `ndarray, float, shape (N,d)`
        Array containing the grid points.  `N` is the number of points, `d` 
        is the dimension of the grid.
    shape : `array_like, shape (d,), optional`
        Number of grid points in each dimension.  Returned only if 
        `return_shape=True`.
    axes : `ndarray, shape (d,d), optional`
        Array of axes vector of the parallelotope cell.  Returned only if 
        `return_axes=True`.
    """
    if not isinstance(axes,np.ndarray):
        axes = np.array(axes)
    #end if
    shape_specifiers = (shape,cells,dr)
    specifier_count = 0
    for s in shape_specifiers:
        specifier_count += int(s is not None)
    #end for
    if specifier_count>1:
        error('provide only one of "shape", "cells" or "dr"','parallelotope_grid_points')
    elif specifier_count==0:
        error('either "shape", "cells" or "dr" must be provided','parallelotope_grid_points')
    #end if
    if shape is None:
        if cells is not None:
            shape = np.array(cells,dtype=int)
        elif dr is not None:
            if not isinstance(dr,np.ndarray):
                dr = np.array(dr,dtype=float)
            #end if
            shape = np.array(np.around(np.linalg.norm(axes,axis=1)/dr),dtype=int)
        #end if
        if not centered and endpoint is not None:
            shape += np.array(endpoint,dtype=int)
        #end if
    #end if
    ugrid = unit_grid_points(shape,centered=centered,endpoint=endpoint)
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



def spheroid_grid_points(axes,shape=None,cells=None,centered=False,endpoint=None,return_shape=False):
    """
    Generation of uniform grids within spheroidal volumes.

    Parameters
    ----------
    shape : `array_like, int, optional`
        Number of points in the grid in each dimension.  The dimension of 
        the grid is the number of entries in `shape`.  Either `shape` or  
        `cells` must be provided.
    cells : `array_like, int, optional`
        Number of cells in the grid in each dimension.  The dimension of the 
        grid is the number of entries in `cells`.  Either `shape` or `cells` 
        must be provided.
    centered : `bool, optional, default False`
        Locate grid points at lower cell corners (`False`) or cell centers 
        (`True`).
    endpoint : `array_like, bool, optional`
        If `True` for given dimension, add an endpoint at the upper edge 
        of the grid in that dimension (default `False`).  Applies only to 
        non-centered grids. `shape/cells/dr` and `endpoint` must have the 
        same number of entries.
    return_shape : `bool, optional, default False`
        Additionally return the shape of the grid.

    Returns
    -------
    sgrid : `ndarray, float, shape (N,d)`
        Array containing the grid points.  `N` is the number of points, `d` 
        is the dimension of the grid.
    shape : `array_like, shape (d,), optional`
        Number of grid points in each dimension.  Returned only if 
        `return_shape=True`.
    """
    if not isinstance(axes,np.ndarray):
        axes = np.array(axes)
    #end if
    shape_specifiers = (shape,cells)
    specifier_count = 0
    for s in shape_specifiers:
        specifier_count += int(s is not None)
    #end for
    if specifier_count>1:
        error('provide only one of "shape" or "cells", not both','spheroid_grid_points')
    elif specifier_count==0:
        error('either "shape" or "cells" must be provided','spheroid_grid_points')
    #end if
    if cells is not None:
        shape = np.array(cells,dtype=int)
        if not centered and endpoint is not None:
            shape += np.array(endpoint,dtype=int)
        #end if
    #end if
    grid_dim,space_dim = axes.shape
    if grid_dim not in (2,3):
        error('spheroid grid generation only supported in 2 or 3 dimensions','spheriod_grid_points')
    #end if
    ugrid = unit_grid_points(shape,centered=centered,endpoint=endpoint)
    if grid_dim==2:
        ugrid[:,1] *= 2*np.pi
        sgrid = polar_to_cartesian(ugrid)
    elif grid_dim==3:
        ugrid[:,1] *=   np.pi
        ugrid[:,2] *= 2*np.pi
        sgrid = spherical_to_cartesian(ugrid)
    #end if
    sgrid = np.dot(sgrid,axes) # adds radial range and skew
    if not return_shape:
        return sgrid
    else:
        return sgrid,shape
    #end if
#end def spheroid_grid_points



def spheroid_surface_grid_points(axes,shape=None,cells=None,centered=False,endpoint=None,return_shape=False):
    """
    Generation of uniform grids over the surfaces of spheroidal volumes.

    Parameters
    ----------
    shape : `array_like, int, optional`
        Number of points in the grid in each dimension.  The dimension of 
        the grid is the number of entries in `shape`.  Either `shape` or  
        `cells` must be provided.
    cells : `array_like, int, optional`
        Number of cells in the grid in each dimension.  The dimension of the 
        grid is the number of entries in `cells`.  Either `shape` or `cells` 
        must be provided.
    centered : `bool, optional, default False`
        Locate grid points at lower cell corners (`False`) or cell centers 
        (`True`).
    endpoint : `array_like, bool, optional`
        If `True` for given dimension, add an endpoint at the upper edge 
        of the grid in that dimension (default `False`).  Applies only to 
        non-centered grids. `shape/cells/dr` and `endpoint` must have the 
        same number of entries.
    return_shape : `bool, optional, default False`
        Additionally return the shape of the grid.

    Returns
    -------
    sgrid : `ndarray, float, shape (N,d)`
        Array containing the grid points.  `N` is the number of points, `d` 
        is the dimension of the grid.
    shape : `array_like, shape (d,), optional`
        Number of grid points in each dimension.  Returned only if 
        `return_shape=True`.
    """
    if not isinstance(axes,np.ndarray):
        axes = np.array(axes)
    #end if
    shape_specifiers = (shape,cells)
    specifier_count = 0
    for s in shape_specifiers:
        specifier_count += int(s is not None)
    #end for
    if specifier_count>1:
        error('provide only one of "shape" or "cells", not both','spheroid_grid_points')
    elif specifier_count==0:
        error('either "shape" or "cells" must be provided','spheroid_grid_points')
    #end if
    if cells is not None:
        shape = np.array(cells,dtype=int)
        if not centered and endpoint is not None:
            shape += np.array(endpoint,dtype=int)
        #end if
    #end if
    grid_dim,space_dim = axes.shape
    grid_dim-=1
    if grid_dim not in (1,2):
        error('spheroid surface grid generation only supported in 1 or 2 dimensions','spheriod_grid_points')
    #end if
    ugrid = unit_grid_points(shape,centered=centered,endpoint=endpoint)
    if grid_dim==1:
        ugrid[:,0] *= 2*np.pi
        sgrid = polar_to_cartesian(ugrid,surface=True)
    elif grid_dim==2:
        ugrid[:,0] *=   np.pi
        ugrid[:,1] *= 2*np.pi
        sgrid = spherical_to_cartesian(ugrid,surface=True)
    #end if
    sgrid = np.dot(sgrid,axes) # adds radial range and skew
    if not return_shape:
        return sgrid
    else:
        return sgrid,shape
    #end if
#end def spheroid_surface_grid_points



class PlotHandler(DevBase):
    """
    Handler class for plotting.

    Currently provides support for 2D and 3D matplotlib plots.


    Attributes
    ----------
    fig : 
        Handle of the current figure.
    ax : 
        Handle of the current figure's axes.
    """
    
    @staticmethod
    def reset():
        """
        Drop handling of the current figure.
        """
        PlotHandler.fig = None
        PlotHandler.ax  = None
    #end def reset

    def set_cur_fig(self,fig):
        """
        Store handle of current figure.
        """
        PlotHandler.fig = fig
    #end def set_cur_fig

    def get_cur_fig(self):
        """
        Return handle of current figure.
        """
        return PlotHandler.fig
    #end def get_cur_fig

    def set_cur_ax(self,ax):
        """
        Store handle of current figure's axes.
        """
        PlotHandler.ax = ax
    #end def set_cur_ax

    def get_cur_ax(self):
        """
        Return handle of current figure's axes.
        """
        return PlotHandler.ax
    #end def get_cur_ax


    def setup_mpl_fig(self,fig=True,dim=None,ax1='x',ax2='y',ax3='z'):
        """
        Setup a matplotlib style plot.

        Parameters
        ----------
        fig : `bool, optional, default True`
            If `True`, create a new figure.  Reuse the current one otherwise.
        dim : `int`
            The dimension of the figure.  Supports 1, 2, and 3 dimensional 
            figures.
        ax1 : `str, optional, default x`
            Label of the first axis.
        ax2 : `str, optional, default y`
            Label of the second axis, if present.
        ax3 : `str, optional, default z`
            Label of the third axis, if present.

        Returns:
        --------
        fig : 
            Handle of the current figure.
        ax : 
            Handle of the current figure's axes.
        """
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
    """
    Base class for Grid and GridFunction classes.

    This class should not be instantiated directly.

    Attributes
    ----------
    initialized : `bool`
        Set to true if the instance has been initialized in non-vacuous fashion.
    """

    descriptor = 'gbase'

    #: (`obj`) Collection of attributes for the class.  Used to check assigned 
    #: members for type conformity and to assign default values.
    persistent_data_types = obj(
        initialized = (bool,False),
        )

    vlogger = None

    @staticmethod
    def reset_vlog():
        GBase.vlogger = None
    #end def reset_vlog

    @staticmethod
    def set_vlog(vlog):
        GBase.vlogger = vlog
    #end def set_vlog

    def vlog(self,*args,**kwargs):
        if GBase.vlogger is not None:
            GBase.vlogger(*args,**kwargs)
        #end if
    #end def vlog


    def __init__(self,*args,**kwargs):
        self.reset()

        if len(args)>0 or len(kwargs)>0:
            self.initialize(**kwargs)
        #end if
    #end def __init__


    def reset(self):
        """
        (`External API`) Reset all attributes to default values.
        """
        cls = self.__class__
        for name,(dtype,default) in cls.persistent_data_types.items():
            self[name] = default
        #end for
    #end def reset


    def initialize(self,*args,**kwargs):
        """
        (`Internal API`) Initialize the instance, starting from the Grid base 
        class and then down the inheritance hierarchy.

        Parameters
        ----------
        check : `bool, optional, default True`
            Check all the assigned attributes for type and shape validity.
        **kwargs : 
            Arbitrary keyword arguments corresponding to instance attributes.
        """
        # remove check argument
        check = kwargs.pop('check',True)

        # call derived initialize
        self.initialize_local(*args,**kwargs)

        # record initialization action
        self.initialized = True

        # check validity of grid
        if check:
            self.check_valid()
        #end if
    #end def initialize

    
    def read(self,filepath,format=None,check=True):
        if isinstance(filepath,StandardFile):
            format = filepath.sftype.lower()
        elif not isinstance(filepath,str):
            self.error('Cannot read file.\nExpected a file path.\nInstead received type: {}\nWith value: {}\nPlease provide a file path and try again.'.format(filepath.__class__.__name__,filepath))
        elif not os.path.exists(filepath):
            self.error('Cannot read file.  File path does not exist.\nFile path provided: {}'.format(filepath))
        elif format is not None:
            if not isinstance(format,str):
                self.error('Cannot read file.\nExpected text (string) for file format.\nInstead received type: {}\nWith value: {}'.format(format.__class__.__name__,format))
            #end if
        else:
            format = filepath.rsplit('.',1)[1].lower()
        #end if
        self.reset()
        self.read_local(filepath,format)
        if check:
            self.check_valid()
        #end if
    #end def read


    def valid(self):
        """
        (`External API`) Return the validity status of a Grid object.

        Returns
        -------
        valid : `bool`
            True if no problems are found.
        """
        return self.check_valid(exit=False)
    #end def valid


    def check_valid(self,exit=True):
        """
        (`External API`) Check the validity of a Grid object.  

        Parameters
        ----------
        exit : `bool, default True`
            If `True` (the default), exit with an error if the object is 
            invalid.

        Returns
        -------
        valid : `bool`
            True if no problems are found.
        """
        cls = self.__class__
        msgs = self.validity_checks()
        valid = len(msgs)==0
        if not valid and exit:
            for msg in msgs:
                self.error(msg,exit=False,trace=False)
            #end for
            self.error('{} is not valid, see error messages above'.format(cls.descriptor))
        #end if
        return valid
    #end def check_valid


    def validity_checks(self):
        """
        (`Internal API`)  Check validity of all assigned attributes, starting 
        from the Grid base class and then down the inheritance hierarchy.
        """
        cls = self.__class__
        msgs = []
        # check that the grid has been initialized
        if not self.initialized:
            msgs.append('{} has not been initialized'.format(cls.descriptor))
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
            type,default = cls.persistent_data_types[name]
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


    def initialize_local(self,*args,**kwargs):
        """
        (`Internal API`) Virtual function used to assign attributes local 
        to the current derived class.
        """
        self.not_implemented()
    #end def initialize_local


    def local_validity_checks(self,msgs):
        """
        (`Internal API`) Virtual function used to check the validity of 
        attributes local to the current derived class.
        """
        self.not_implemented()
    #end def local_validity_checks


    def read_local(self,filepath,format):
        self.not_implemented()
    #end def read_local


    # test needed
    def ensure_array(self,dtype=None,**arrays_in):
        arrays = obj()
        for k,ai in arrays_in.items():
            if isinstance(ai,(tuple,list)):
                if dtype is None:
                    a = np.array(ai)
                else:
                    a = np.array(ai,dtype=dtype)
                #end if
            elif isinstance(ai,np.ndarray):
                a = ai
            else:
                self.error('Cannot ensure array value.\nReceived data with type: {}\nOnly tuple, list, and ndarray are supported.'.format(ai.__class__.__name__))
            #end if
            arrays[k] = a
        #end for
        return arrays
    #end def ensure_array
#end class GBase



class Grid(GBase):
    """
    Base class for `M` dimensional grids embedded within `N` dimensional spaces.

    Universal grid properties are handled/represented at this level.  This 
    includes the fact that a grid contains a set of points and resides within
    a space of a particular dimension.  Derived classes add more specific 
    information and functionality for particular kinds of grids.

    General initialization, checking, and copying of all types of grids is 
    handled at this level.

    This class should not be instantiated directly.

    Parameters
    ----------
    points : `array_like, float, optional`
        Array of grid points.

    Attributes
    ----------
    points : `ndarray, float, shape (N,d)`
        Array containing the grid points.  `N` is the number of grid points, 
        `d` is the dimension of the space.
    r : `ndarray, float, property`
        Array containing the grid points.  User-facing alias for `points`.
    npoints : `int, property`
        Number of grid points.
    space_dim : `int, property`
        Dimension of the space the grid resides in.
    dtype : `property`
        Datatype of the grid point values.
    """

    descriptor = 'grid'

    #: (`obj`) Collection of attributes for the class.  Used to check assigned 
    #: members for type conformity and to assign default values.
    persistent_data_types = obj(
        points      = (np.ndarray,None ),
        **GBase.persistent_data_types
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


    def initialize_local(self,points=None,dtype=None,copy=True):
        """
        (`Internal API`) Sets `points` and `dtype` attributes.

        Parameters
        ----------
        points : `array_like, float, shape (N,d)`
            Array of grid points.  `N` is the number of points, `d` is the 
            dimension of the space (not necessarily the same as the dimension 
            of the grid).
        dtype : `optional`
            Data type of the grid point values.  Should be similar to `float`.
        copy : `bool, optional, default True`
            If `True`, perform a deep copy of the grid point data, othewise
            assign directly (shallow/pointer copy).
        """
        if points is None:
            self.error('cannot initialize grid, "points" is required')
        #end if
        self.set_points(points,dtype=dtype,copy=copy)
    #end def initialize_local


    def set_points(self,points,dtype=None,copy=True):
        """
        (`Internal API`) Set the grid points, checking for type validity.

        Parameters
        ----------
        points : `array_like, float, shape (N,d)`
            Array of grid points.  `N` is the number of points, `d` is the 
            dimension of the space (not necessarily the same as the dimension 
            of the grid).
        dtype : `optional`
            Data type of the grid point values.  Should be similar to `float`.
        copy : `bool, optional, default True`
            If `True`, perform a deep copy of the grid point data, othewise
            assign directly (shallow/pointer copy).
        """
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
        """
        (`External API`) Return a copy of the grid instance.

        Parameters
        ----------
        shallow : `bool, optional, default False`
            If `False` (the default), perform a deep copy of the object.  
            Otherwise, perform a deep copy of all attributes except for 
            `points` which is copied shallowly.
        """
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
        """
        (`External API`)  Translate all points in the grid by a vector.

        Parameters
        ----------
        shift : `array_like, float, shape (d,)`
            Vector displacement used to translate the grid points.
        """
        self.points += shift
    #end def translate


    def local_validity_checks(self,msgs):
        """
        (`Internal API`) Check validity of `shape` and `points`.

        Parameters
        ----------
        msgs : `list, str`
            List containing error messages.  Empty if no problems are found.
        """
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


    def check_valid_points(self,points,dim,loc):
        """
        (`Internal API`) Ensure inputted points are valid.

        This function upcasts inputted points to an `ndarray` and checks that 
        the shape is correct.  It is used by other functions to guard against 
        improper input for `points`.

        Parameters
        ----------
        points : `array_like, float, shape (N,d)`
            Array of `N` points of dimension `d`.
        dim : `int`
            Dimension of the space the points must reside in (`d`).
        loc : `str`
            Name of the calling function.  This is used to format appropriate 
            error messages.
        """
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
        """
        (`External API`)  Make a scatter plot of a set of points in 2D or 3D.

        By default, the set of grid points is plotted.

        Parameters
        ----------
        points : `array_like, float, optional`
            Set of points to plot.  If no points are provided, the grid points 
            are plotted.
        fig : `bool, optional, default True`
            If `True`, make a new figure.  Reuse the current one otherwise.
        show : `bool, optional, default True`
            If `True`, display the plot immediately on the screen.
        default_marker : `str, optional, default` "."
            Default marker symbol for the scatter plot.  Used if "marker" is 
            not provided as a keyword argument.
        **kwargs : 
            Arbitrary keyword arguments passed to `pyplot.scatter`.
        """
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


    def grid_function(self,*args,**kwargs):
        gf = grid_function_from_grid(self)
        return gf(*args,**kwargs)
    #end def grid_function
#end class Grid



class StructuredGrid(Grid):
    """
    Base class for structured grids.

    A structured grid has the property that each grid point in the `M` 
    dimensional space, excepting those on the boundary, are connected to `2M` 
    neighbors.  We further assume that the domain can be described by an `M` 
    dimensional coordinate system that can be mapped onto a unit cube of 
    dimension `M`.  The general structured grid can then be described through  
    a mapping of a uniform grid over the unit M-cube onto the target space. 

    This class represents grids of this type and adds appropriate descriptive 
    properties and functions to the bare Grid class, including the fact that 
    the grid has a particluar shape (number of grid cells in each dimension), 
    that its points can reside at cell centers or edges, and that it has 
    boundaries (as follows from the mapping of the `M`-cube facets) and 
    accompanying boundary conditions (e.g. periodic).  Further, the grid can 
    exist on the surface of another space with little additional complication.

    An advantage of framing an entire class of grids through mappings to the 
    unit `M`-cube is that operations such as interpolation and integration of 
    functions on these grids can be formulated in a centralized fashion.

    The actual mappings back and forth from the unit `M`-cube are left to 
    derived classes.

    This class should not be instantiated directly.

    Parameters
    ----------
    shape : `array_like, int, shape (d,)`
        The number of grid points in each dimension. `d` is the dimension of 
        the grid (`grid_dim`).
    centered : `bool, optional, default False`
        Locate grid points at lower cell corners (`False`) or cell centers 
        (`True`).
    bconds : `array_like, str, shape (d,), {'o','p'}, optional, default d*['o']`
        Boundary conditions for each dimension.  Options are open (`'o'`) 
        and periodic (`'p'`).  `d` is the dimension of the grid (`grid_dim`).
    surface : `bool`
        If `True`, then the grid is known to reside on the surface of another 
        space.  Otherwise, the grid spans the volume of that space.

    Attributes
    ----------
    shape : `tuple, int`
        The number of grid points in each dimension.
    centered : `bool`
        Grid points are located at lower cell corners (`False`) or cell 
        centers (`True`).
    bconds : `ndarray, str`
        Boundary conditions for each dimension.
    surface : `bool`
        The grid is resides on the surface of a space (`True`), otherwise, it 
        spans the volume of the space (`False`).  This attribute is intended 
        to be immutable once set.
    grid_dim : `int, property`
        The dimension of the grid.  Must be less than or equal to `space_dim`.
    grid_shape : `tuple, int, property`
        The number of grid points in each dimension.
    cell_grid_shape : `tuple, int, property`
        The number of grid cells in each dimension.
    ncells : `int, property`
        The total number of grid cells.
    flat_points_shape : `tuple, int, property`
        The shape of the `points` array in its default (flat) representation. 
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `flat_points_shape` is `(N*M*P,D)`.
    full_points_shape : `tuple, int, property`
        The shape of the points array in its full representation.
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `full_points_shape` is `(N,M,P,D)`.
    """

    #: (`obj`) Collection of attributes for the class.  Used to check assigned 
    #: members for type conformity and to assign default values.
    persistent_data_types = obj(
        shape    = (tuple     ,None),
        centered = (bool      ,None),
        bconds   = (np.ndarray,None),
        surface  = (bool      ,None),
        **Grid.persistent_data_types
        )

    #: (`set`)  Set of valid boundary condition types.
    valid_bconds = set(['o','p'])

    #: (`obj`)  Keyword mapping of boundary condition types.
    bcond_types = obj(
        open     = 'o',
        periodic = 'p',
        )


    @property
    def grid_dim(self):
        return len(self.shape)
    #end def grid_dim

    @property
    def grid_shape(self):
        return self.shape
    #end def grid_shape

    @property
    def cell_grid_shape(self):
        shape = np.array(self.shape,dtype=int)
        if not self.centered:
            shape -= np.array(self.has_endpoints(),dtype=int)
        #end if
        return tuple(shape)
    #end def cell_grid_shape

    @property
    def ncells(self):
        return np.prod(self.cell_grid_shape)
    #end def ncells

    @property
    def npoints(self):
        return np.prod(self.shape)
    #end def npoints

    @property
    def flat_points_shape(self):
        space_dim = self.r.shape[-1]
        npoints = np.prod(self.shape)
        return (npoints,space_dim)
    #end def flat_points_shape

    @property
    def full_points_shape(self):
        space_dim = self.r.shape[-1]
        return self.shape+(space_dim,)
    #end def flat_points_shape

    @property
    def periodic(self):
        return (self.bconds==self.bcond_types.periodic).all()
    #end def periodic


    def initialize_local(self,**kwargs):
        """
        (`Internal API`) Sets `shape`, `centered`, `bconds`, and `surface` 
        attributes.

        The `surface` attribute is set to `False`.  It is the responsibility 
        of the derived classes to set this in an appropriate way.
        """
        shape    = kwargs.pop('shape'   ,None)
        centered = kwargs.pop('centered',False)
        bconds   = kwargs.pop('bconds'  ,None)
        Grid.initialize_local(self,**kwargs)
        self.centered = centered
        self.surface  = False
        if shape is None:
            self.error('cannot initialize grid, "shape" is required')
        #end if
        if not isinstance(shape,(tuple,list,np.ndarray)):
            self.error('cannot set shape from data with type "{}"\nplease use tuple, list, or array for inputted shape'.format(shape.__class__.__name__))
        #end if
        shape = np.array(shape,dtype=int)
        if bconds is None:
            bconds = len(shape)*[self.bcond_types.open]
        #end if
        self.set_bconds(bconds)
        self.set_shape(shape)
    #end def initialize_local


    def set_shape(self,shape):
        """
        (`Internal API`) Sets the `shape` attribute in a protected way.
        """
        if not isinstance(shape,(tuple,list,np.ndarray)):
            self.error('cannot set shape from data with type "{}"\nplease use tuple, list, or array for inputted shape'.format(shape.__class__.__name__))
        #end if
        if 'points' in self and self.points is not None:
            npoints = np.array(shape,dtype=int).prod()
            if npoints!=len(self.points):
                self.error('cannot set shape\ngrid shape provided does not match the number of points in the grid\npoints in grid: {}\npoints from shape: {}'.format(len(self.points),npoints))
            #end if
        #end if
        self.shape = tuple(shape)
    #end def set_shape


    def set_bconds(self,bconds):
        """
        (`Internal API`) Sets the `bconds` attribute in a protected way.
        """
        if isinstance(bconds,str):
            bconds = tuple(bconds)
        #end if
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


    def has_endpoints(self,bconds=None,grid_dim=None):
        """
        (`Internal API`) Determine whether a grid should have endpoints at its
        upper boundaries.

        Endpoints should not be present in periodic boundary conditions.
        """
        if grid_dim is not None:
            if bconds is None:
                bconds = grid_dim*[self.bcond_types.open]
            elif isinstance(bconds,str):
                bconds = tuple(bconds)
            #end if
            bconds = np.array(bconds,dtype=object)
        #end if
        if bconds is None:
            bconds = self.bconds
        #end if
        return bconds==self.bcond_types.open
    #end def has_endpoints


    def local_validity_checks(self,msgs):
        """
        (`Internal API`) Check the validity of the `shape` and `bconds` 
        attributes.

        Parameters
        ----------
        msgs : `list, str`
            List containing error messages.  Empty if no problems are found.
        """
        msgs = Grid.local_validity_checks(self,msgs)
        if np.prod(self.shape)!=self.npoints:
            msgs.append('grid shape does not match number of points\nnumber of points: {}\nproduct of shape: {}\nshape: {}'.format(self.npoints,np.prod(self.shape),self.shape))
        #end if
        if len(self.shape)!=self.grid_dim:
            msgs.append('number of entries in grid shape does not match grid dimension\nnumber of entries in grid shape: {}\ngrid dimension: {}'.format(len(self.shape),self.grid_dim))
        #end if
        if len(self.bconds)!=self.grid_dim:
            msgs.append('number of entries in bconds does not match grid dimension\nnumber of entries in bconds: {}\ngrid dimension: {}'.format(len(self.bconds),self.grid_dim))
        if len(set(self.bconds)-StructuredGrid.valid_bconds)>0:
            msgs.append('boundary conditions are invalid\nboundary conditions in each dimension must be one of: {}\nboundary conditions provided: {}'.format(sorted(StructuredGrid.valid_bconds),self.bconds))
            #end if
        #end for
        return msgs
    #end def local_validity_checks


    def reshape_full(self):
        """
        (`Internal API`) Transform the grid points array into full shape.

        This should only be a local and temporary change of state.  It should 
        be reversed by calling `reshape_flat` as soon as possible.
        """
        self.r.shape = self.full_points_shape
    #end def reshape_full


    def reshape_flat(self):
        """
        (`Internal API`) Transform the grid points array into the default flat
        shape.

        This function is meant to reverse the temporary state change induced 
        by `reshape_full`.
        """
        self.r.shape = self.flat_points_shape
    #end def reshape_flat


    # test needed
    def flat_indices(self,full_indices):
        if not isinstance(full_indices,np.ndarray):
            full_indices = np.array(full_indices,dtype=int)
        #end if
        if len(full_indices.shape)!=2:
            self.error('full_indices must have shape (# points)x(# dimensions)\nShape received: {}'.format(full_indices.shape))
        elif full_indices.shape[-1]!=self.grid_dim:
            self.error('full_indices must have same dimension as the grid.\nfull_indices dimension: {}\nGrid dimension: {}'.format(full_indices.shape[-1],self.grid_dim))
        #end if
        grid_shape = self.grid_shape
        D = len(grid_shape)
        grid_shape_prod = np.empty((D,),dtype=int)
        grid_shape_prod[-1] = 1
        for d in range(D-1):
            grid_shape_prod[D-d-2] = grid_shape[D-d-1]*grid_shape_prod[D-d-1]
        #end for
        flat_indices = np.dot(full_indices,grid_shape_prod)
        return flat_indices
    #end def flat_indices


    def unit_points(self,points=None,project=False):
        """
        (`External API`) Map a set of points into the unit cube.

        Parameters
        ----------
        points : `array_like, float, shape(N,d), optional`
            Array of points in the full coordinate space. `N` is the number of
            points and `d` must be equal to `space_dim`.  The grid points are 
            used if no points are provided.

        Returns
        -------
        upoints : `ndarray, float, shape (N,dg)`
            Array of points in the unit coordinate space.  `N` is the number 
            of points and `dg` is equal to `grid_dim`.
        """
        if points is None:
            points = self.r
        else:
            points = self.check_valid_points(points,self.space_dim,'unit_points')
        #end if
        upoints = self.unit_points_bare(points)
        if project:
            for d in range(self.grid_dim):
                if self.bconds[d]==self.bcond_types.periodic:
                    upoints[:,d] -= np.floor(upoints[:,d])
                #end if
            #end for
        #end if
        return upoints
    #end def unit_points


    def unit_metric(self,upoints=None):
        """
        (`External API`)  Compute the integration metric in the unit coordinate 
        space for a set of points defined there.

        Parameters
        ----------
        upoints : `array_like, float, shape (N,d), optional`
            Array of points in the unit coordinate space.  `N` is the number 
            of points and `d` must be equal to `grid_dim`.  The unit 
            representation of the grid points is used if no points are 
            provided.

        Returns
        -------
        umetric : `ndarray, float, shape (N,)`
            Array containing the integration metric in the unit space at the 
            set of points provided.  `N` is the number of points.
        """
        return self.unit_metric_bare(upoints)
    #end def unit_metric


    def cell_indices(self,points=None,project=True):
        """
        (`External API`) Given a set of points, find the index of the grid cell 
        bounding each point.

        Parameters
        ----------
        points : `array_like, float, shape (N,d), optional`
            Array of points in the full coordinate space. `N` is the number of
            points and `d` must be equal to `space_dim`.  The grid points are 
            used if no points are provided.

        Returns
        -------
        ipoints : `ndarray, int, shape (N,)`
            Array of cell indices.  `N` is the number of points.
        """
        upoints = self.unit_points(points,project=project)
        shape = np.array(self.cell_grid_shape,dtype=int)
        ipoints = upoints*shape
        ipoints = np.array(np.floor(ipoints),dtype=int)
        dim = self.grid_dim
        cum_shape = np.empty((dim,),dtype=int)
        cum_shape[-1] = 1
        for d in range(1,dim):
            cum_shape[dim-d-1] = cum_shape[dim-d]*shape[dim-d-1]
        #end for
        ipoints = np.dot(ipoints,cum_shape)
        return ipoints
    #end def cell_indices


    def inside(self,points,tol=1e-12):
        """
        (`External API`)  Given a set of points, determine which are inside the 
        boundary of the grid.

        Parameters
        ----------
        points : `array_like, float, shape (N,d)`
            Array of points in the full coordinate space. `N` is the number of
            points and `d` must be equal to `space_dim`.

        Returns
        -------
        inside : `ndarray, bool, shape (N,)`
            Mask array that is `True` if a point is inside and `False` 
            otherwise.  `N` is the number of points.
        """
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
        """
        (`External API`) Project a set of points into the grid domain, if possible.

        The points are first projected into the unit cube and then, if in 
        periodic boundary conditions, folded back into the unit cube.  Note 
        that following this operation, some points may still fall outside the 
        grid boundary if there are any open boundaries.

        In derived classes, the projection onto the unit cube may have the 
        additional effect of projection points from a higher dimensional 
        embedding space onto the grid domain.

        Parameters
        ----------
        points : `array_like, float, shape (N,d)`
            Array of points in the full coordinate space. `N` is the number of
            points and `d` must be equal to `space_dim`.

        Returns
        -------
        proj_points : `ndarray, float, shape (N,d)`
            Array of points projected onto the grid domain where possible.
        """
        points = self.check_valid_points(points,self.space_dim,'project')
        upoints = self.unit_points(points,project=True)
        return self.points_from_unit(upoints)
    #end def project


    def get_boundary_lines(self,n=200,unit=False):
        """
        (`Internal API`)  Generate a set of points on the edges of the 
        boundary.

        This function is used primarily to make plots of the grid domain 
        boundaries.

        Parameters
        ----------
        n : `int, optional, default 200`
            Number of points to generate along each boundary edge.
        unit : `bool, optional, default False`
            Generate the points in unit coordinates (`True`), or in the full 
            space (`False`).

        Returns
        -------
        bpoints : `ndarray, float, shape (NL,n,d)`
            Array containing the boundary lines. `NL` is the number of lines 
            and `d` is either the dimension of the full space (`space_dim, 
            unit=False`) or the embedded space (`grid_dim, unit=True`).
        """
        u = np.linspace(0.,1.,n)
        nlines = self.grid_dim*2**(self.grid_dim-1)
        upoints = np.empty((n*nlines,self.grid_dim),dtype=self.dtype)
        if self.grid_dim==3:
            bitset = [[0,0],[0,1],[1,0],[1,1]]
        elif self.grid_dim==2:
            bitset = [[0],[1]]
        elif self.grid_dim==1:
            bitset = [[0]]
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
        """
        (`External API`) Plot the edges of the boundary in the full space.

        Parameters
        ----------
        n : `int, optional, default 200`
            Number of points to generate along each boundary edge.
        fig : `bool, optional, default True`
            Create a fresh figure (`True`) or not (`False`).
        show : `bool, optional, default True`
            Display the figure (`True`) or not (`False`).
        """
        fig,ax = self.setup_mpl_fig(fig=fig,dim=self.space_dim)
        if self.grid_dim!=1 or self.surface:
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
        """
        (`External API`) Make a scatter plot of a set of points in unit 
        coordinates.

        The inputted points are projected into the unit cube before plotting.  
        By default, the grid points are plotted.

        Parameters
        ----------
        points : `array_like, float, optional`
            Set of points to plot.  If no points are provided, the grid points 
            are plotted.
        fig : `bool, optional, default True`
            If `True`, make a new figure.  Reuse the current one otherwise.
        show : `bool, optional, default True`
            If `True`, display the plot immediately on the screen.
        default_marker : `str, optional, default` "."
            Default marker symbol for the scatter plot.  Used if "marker" is 
            not provided as a keyword argument.
        **kwargs : 
            Arbitrary keyword arguments passed to `pyplot.scatter`.
        """
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
        """
        (`External API`) Plot the edges of the boundary in the unit space.

        Parameters
        ----------
        n : `int, optional, default 200`
            Number of points to generate along each boundary edge.
        fig : `bool, optional, default True`
            Create a fresh figure (`True`) or not (`False`).
        show : `bool, optional, default True`
            Display the figure (`True`) or not (`False`).
        """
        fig,ax = self.setup_mpl_fig(fig=fig,dim=self.grid_dim,
                                    ax1='a1',ax2='a2',ax3='a3')
        if self.grid_dim!=1 or self.surface:
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
        """
        (`Internal API`)  Derived class function to map points into the unit 
        cube.
        """
        self.not_implemented()
    #end def unit_points_bare


    def points_from_unit(self,upoints):
        """
        (`External API`) Map points from the unit space back into the full space.

        Parameters
        ----------
        upoints : `array_like, float, shape (N,dg), optional`
            Array of points in the unit coordinate space.  `N` is the number 
            of points and `dg` must be equal to `grid_dim`.

        Returns
        -------
        points : `ndarray, float, shape (N,ds)`
            Array of points in the full coordinate space.  `ds` is the 
            dimension of the full space (`space_dim`).
        """
        self.not_implemented()
    #end def points_from_unit


    def unit_metric_bare(self,upoints):
        """
        (`Internal API`) Derived class function that computes the integration 
        metric in the unit coordinate space for a set of points defined there.
        """
        self.not_implemented()
    #end def unit_metric_bare


    def volume(self):
        """
        (`External API`)  Compute the volume of the space bounding the grid.

        Returns
        -------
        volume : `float`
            Volume of the space within the grid boundary.
        """
        self.not_implemented()
    #end def volume


    def cell_volumes(self):
        """
        (`External API`) Compute the volumes of the grid cells.

        Returns
        -------
        cell_vols : `ndarray, shape (N,)`
            Array containing the volume of each grid cell.  `N` is the number 
            of points in the grid.
        """
        self.not_implemented()
    #end def cell_volumes
#end class StructuredGrid



class StructuredGridWithAxes(StructuredGrid):
    """
    Base class for structured grids with linear axes that act as scaffolding 
    for the coordinate system.

    Examples are sheared cartesian and spherical coordinate systems.  Each 
    derived class optionally enacts a coordinate transformation on top of 
    the sheared axes.  With the introduction of the axes, this class also 
    handles the notion of the origin of the coordinate system.

    This class should not be instantiated directly.

    Parameters
    ----------
    axes : `array_like, float, shape (dg,ds)`
        Axes used to form the coordinate system.  `dg` is the dimension of 
        the grid (embedded space, `dg=grid_dim`), `ds` is the dimension of 
        the points in the grid (full space, `ds=space_dim`).
    origin : `array_like, float, shape (ds,), optional, default ds*[0]`
        Origin of the space.  `ds` is the dimension of the full space.
    
    Attributes
    ----------
    axes : `ndarray, float, shape (dg,ds)`
        Axes used to form the coordinate system.  `dg` is the dimension of 
        the grid (embedded space, `dg=grid_dim`), `ds` is the dimension of 
        the points in the grid (full space, `ds=space_dim`).
    origin : `ndarray, float, shape (ds,)`
        Origin of the space.  `ds` is the dimension of the full space.
    """

    #: (`obj`) Collection of attributes for the class.  Used to check assigned 
    #: members for type conformity and to assign default values.
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
        """
        (`Internal API`) Sets `axes` and `origin` attributes.

        The `surface` attribute is set to `False`.  It is the responsibility 
        of the derived classes to set this in an appropriate way.
        """
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


    def set_axes(self,axes):
        """
        (`Internal API`) Sets the `axes` attribute in a protected way.
        """
        if not isinstance(axes,(tuple,list,np.ndarray)):
            self.error('cannot set axes from data with type "{}"\nplease use tuple, list, or array for inputted axes'.format(axes.__class__.__name__))
        #end if
        self.axes = np.array(axes,dtype=self.dtype)
    #end def set_axes


    def set_origin(self,origin):
        """
        (`Internal API`) Sets the `origin` attribute in a protected way.
        """
        if not isinstance(origin,(tuple,list,np.ndarray)):
            self.error('cannot set origin from data with type "{}"\nplease use tuple, list, or array for inputted origin'.format(origin.__class__.__name__))
        #end if
        origin = np.array(origin,dtype=self.dtype)
        if self.origin is not None:
            shift = origin-self.origin
        else:
            shift = origin
            self.origin = 0*origin
        #end if
        self.translate(shift)
    #end def set_origin


    def translate(self,shift):
        """
        (`External API`)  Translate all points in the grid by a vector.

        Parameters
        ----------
        shift : `array_like, float, shape (d,) or (1,d)`
            Vector displacement used to translate the grid points.
        """
        self.origin += shift
        Grid.translate(self,shift)
    #end def translate


    def local_validity_checks(self,msgs):
        """
        (`Internal API`) Check the validity of the `axes` and `origin` 
        attributes.

        Parameters
        ----------
        msgs : `list, str`
            List containing error messages.  Empty if no problems are found.
        """
        msgs = StructuredGrid.local_validity_checks(self,msgs)
        shape = self.axes.shape
        if len(shape)!=2:
            msgs.append('axes must be a 2 dimensional array\nnumber of dimensions present: {}\naxes present: {}'.format(len(shape),self.axes))
        else:
            if not self.surface:
                if shape[0]!=self.grid_dim:
                    msgs.append('number of axes must be equal to the embedded grid dimension\nembedded grid dimension: {}\nnumber of axes present: {}'.format(self.grid_dim,shape[0]))
                #end if
            else:
                if shape[0]!=self.grid_dim+1:
                    msgs.append('number of axes must be equal to the embedded grid dimension + 1 for a surface\nembedded grid dimension + 1: {}\nnumber of axes present: {}'.format(self.grid_dim+1,shape[0]))
                #end if
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


    def axes_volume(self):
        """
        (`Internal API`) Compute the volume enclosed by the axes.

        The volume enclosed by the axes is not in general the volume 
        enclosed by the grid boundaries, but it is used to compute the full 
        grid volume once any additional coordinate transformations are 
        accounted for.

        This function uses singular value decomposition to compute the 
        volume.  This is important for grids embedded in higher dimensional 
        spaces.

        Returns
        -------
        ax_vol : `float`
            Volume enclosed by the axes.  
        """
        return np.abs(np.prod(np.linalg.svd(self.axes,compute_uv=False)))
    #end def axes_volume


    def indices_for_map_coord(self,points):
        if not self.centered:
            ucorner = np.array([0,0,0],dtype=float)
        else:
            ucorner = 0.5/np.array(self.cell_grid_shape)
        #end if
        grid_shape = np.array(self.grid_shape)
        ipoints = (self.unit_points(points)-ucorner)*grid_shape
        return ipoints.T
    #end def indices_for_map_coord

#end class StructuredGridWithAxes



class ParallelotopeGrid(StructuredGridWithAxes):
    """
    A regular structured grid over a sheared cell (parallelotope).

    This type of grid is standard for representing the domain of one, two, 
    and three dimensional functions in scientific applications.

    An `M`-dimensional parallelotope grid may be embedded in an `N`-dimensional 
    space, e.g. a 2D planar plaquette with arbitrary orientation in an open 
    3D space.  Such grids are useful to describe planar or line cuts through 
    higher dimensional spaces.

    Below, `M` and `N` are referred to as `dg` and `ds`, or the dimension of 
    the grid (embedded space) and the full space (embedding space) respectively.
    These spaces may be chosen to be the same.

    This class is intended for direct instantiation and use.

    Parameters
    ----------
    axes : `array_like, float, shape (dg,ds)`
        Axes used to form the coordinate system.  `dg` is the dimension of 
        the grid (embedded space, `dg=grid_dim`), `ds` is the dimension of 
        the points in the grid (full space, `ds=space_dim`).
    bconds : `array_like, str, shape (d,), {'o','p'}, optional, default d*['o']`
        Boundary conditions for each dimension.  Options are open (`'o'`) 
        and periodic (`'p'`).  `d` is the dimension of the grid (`grid_dim`).
    centered : `bool, optional, default False`
        Locate grid points at lower cell corners (`False`) or cell centers 
        (`True`).
    cells : `array_like, int, shape (dg,), optional`
        The number of grid cells in each dimension. An additional grid point 
        will be present at the upper edge of any dimension with open boundary 
        conditions.  `dg` is the dimension of the grid (`grid_dim`).  One of 
        `cells`, `dr`, or `shape` must be provided to generate the grid.
    dr : `array_like, float, shape (dg,), optional`
        The approximate width of grid cells in each dimension. The number of 
        grid cells is determined by adjusting `dr` in each dimension such that 
        an integer number of cells results.  `dg` is the dimension of the grid 
        (`grid_dim`).  One of `cells`, `dr`, or `shape` must be provided to 
        generate the grid.
    shape : `array_like, int, shape (dg,), optional`
        The number of grid points in each dimension. In periodic boundary 
        conditions, the number of grid cells and grid points match in each 
        dimension.  In dimensions with open boundary conditions, there is 
        one fewer grid cell than grid point since a grid point is present 
        at the upper edge of open boundaries, but not periodic ones.  `dg` is 
        the dimension of the grid (`grid_dim`).  One of `cells`, `dr`, or 
        `shape` must be provided to generate the grid.
    corner : `array_like, float, shape (ds,), optional, default ds*[0]`
        The location of the lower corner of the cell.  This point is also 
        the origin of the grid coordinate system.  `ds` is the dimension of the 
        full space (`space_dim`). One of `corner`, `center`, or `origin` must 
        be provided, with `corner` or `center` preferred.
    center : `array_like, float, shape (ds,), optional`
        The location of the middle of the cell.  It is sometimes convenient 
        to define the location of the cell relative to its mid-point rather 
        than its lower corner. `ds` is the dimension of the full space 
        (`space_dim`). One of `corner`, `center`, or `origin` must be provided,
        with `corner` or `center` preferred.
    origin : `array_like, float, shape (ds,), optional, default ds*[0]`
        Origin of the space.  For a parallelotope grid, this point is the 
        lower corner of the cell.  `ds` is the dimension of the full space 
        (`space_dim`). One of `corner`, `center`, or `origin` must be provided,
        with `corner` or `center` preferred.

    Attributes
    ----------
    points : `ndarray, float, shape (N,d)`
        Array containing the grid points.  `N` is the number of grid points, 
        `d` is the dimension of the space.
    axes : `ndarray, float, shape (dg,ds)`
        Axes used to form the coordinate system.  `dg` is the dimension of 
        the grid (embedded space, `dg=grid_dim`), `ds` is the dimension of 
        the points in the grid (full space, `ds=space_dim`).
    bconds : `ndarray, str`
        Boundary conditions for each dimension.
    centered : `bool`
        Grid points are located at lower cell corners (`False`) or cell 
        centers (`True`).
    shape : `tuple, int`
        The number of grid points in each dimension.
    origin : `ndarray, float, shape (ds,)`
        Origin of the space.  `ds` is the dimension of the full space.
    surface : `bool`
        (`Internal`) The grid is resides on the surface of a space (`True`), 
        otherwise, it spans the volume of the space (`False`).  This attribute 
        is intended to be immutable once set.
    initialized : `bool`
        (`Internal`) Set to true if the instance has been initialized in 
        non-vacuous fashion.
    r : `ndarray, float, property`
        Array containing the grid points.  User facing alias for `points`.
    npoints : `int, property`
        Number of grid points.
    space_dim : `int, property`
        Dimension of the space the grid resides in.
    grid_dim : `int, property`
        The dimension of the grid.  Must be less than or equal to `space_dim`.
    grid_shape : `tuple, int, property`
        The number of grid points in each dimension.
    cell_grid_shape : `tuple, int, property`
        The number of grid cells in each dimension.
    dr : `ndarray, float, shape(dg,ds), property`
        Vector displacements between neighboring grid points.
    ncells : `int, property`
        The total number of grid cells.
    flat_points_shape : `tuple, int, property`
        The shape of the `points` array in its default (flat) representation. 
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `flat_points_shape` is `(N*M*P,D)`.
    full_points_shape : `tuple, int, property`
        The shape of the points array in its full representation.
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `full_points_shape` is `(N,M,P,D)`.
    corner : `ndarray, float, shape (ds,), property`
        Location of the lower corner of the cell.  `ds` is the dimension of the
        full space.
    center : `ndarray, float, shape (ds,), property`
        Location of the lower center of the cell.  `ds` is the dimension of the
        full space.
    dtype : `property`
        Datatype of the grid point values.
    """

    @property
    def dr(self):
        dr = np.empty(self.axes.shape,self.dtype)
        cells = self.cell_grid_shape
        for d in range(self.grid_dim):
            dr[d] = self.axes[d]/cells[d]
        #end for
        return dr
    #end def dr

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
                         cells    = None,
                         dr       = None,
                         corner   = None,
                         center   = None,
                         centered = False,
                         **kwargs
                         ):
        """
        (`Internal API`)  Initialize the parallelotope grid points and set the
        origin if not provided.
        """
        if shape is None and cells is None and dr is None:
            self.error('cannot initialize grid, either "shape", "cells", or "dr" is required')
        elif shape is not None:
            grid_dim = len(shape)
        elif cells is not None:
            grid_dim = len(cells)
        elif dr is not None:
            grid_dim = len(dr)
        #end if
        bconds = kwargs.get('bconds',None)

        endpoint = self.has_endpoints(bconds=bconds,grid_dim=grid_dim)

        points,shape,axes = parallelotope_grid_points(
            axes         = axes,
            shape        = shape,
            cells        = cells,
            dr           = dr,
            centered     = centered,
            endpoint     = endpoint,
            return_shape = True,
            return_axes  = True
            )

        if center is not None:
            center = np.array(center)
            corner = center - axes.sum(axis=0)/2
        #end if

        kwargs['axes']     = axes
        if corner is not None:
            kwargs['origin']   = corner
        #end if
        kwargs['shape']    = shape
        kwargs['centered'] = centered
        kwargs['points']   = points
        StructuredGridWithAxes.initialize_local(self,**kwargs)

    #end def initialize_local


    def read_local(self,filepath,format):
        if format=='xsf':
            self.read_xsf(filepath)
        else:
            self.error('Cannot read file.\nUnrecognized file format encountered.\nUnrecognized file format: {}\nValid options are: xsf'.format(format))
        #end if
    #end def read_local


    def read_xsf(self,filepath):
        if isinstance(filepath,XsfFile):
            xsf = filepath
        else:
            xsf = XsfFile(filepath)
        #end if

        d = xsf.get_density()

        cells = d.grid-1
        c = d.cell.sum(axis=0)/cells/2

        centered = False
        corner   = None
        if np.abs(c-d.corner).max()<1e-6:
            centered = True
        else:
            corner = d.corner
        #end if

        self.initialize(
            bconds   = tuple('ppp'),
            axes     = d.cell.copy(),
            cells    = cells,
            centered = centered,
            corner   = corner,
            )
    #end def read_xsf


    def unit_points_bare(self,points):
        """
        (`Internal API`)  Maps points in a parallelotope into the unit cube.
        """
        corner = self.corner
        # invert using pseudo-inverse
        #   this is important for grids embedded in higher dim spaces
        axinv  = np.linalg.pinv(self.axes)
        upoints = np.dot(points-corner,axinv)
        return upoints
    #end def unit_points_bare


    def points_from_unit(self,upoints):
        """
        (`External API`) Map points from the unit space back into the full space.

        Parameters
        ----------
        upoints : `array_like, float, shape (N,dg), optional`
            Array of points in the unit coordinate space.  `N` is the number 
            of points and `dg` must be equal to `grid_dim`.

        Returns
        -------
        points : `ndarray, float, shape (N,ds)`
            Array of points in the full coordinate space.  `ds` is the 
            dimension of the full space (`space_dim`).
        """
        points = np.dot(upoints,self.axes)+self.corner
        return points
    #end def points_from_unit


    def unit_metric_bare(self,upoints):
        """
        (`Internal API`) Compute the parallelotope integration metric in the 
        unit coordinate space for a set of points defined there.
        """
        if upoints is None:
            upoints = self.unit_points()
        #end if
        umetric = np.zeros((len(upoints),),dtype=self.dtype)
        umetric += self.volume()
        return umetric
    #end def unit_metric_bare


    def volume(self):
        """
        (`External API`)  Compute the volume of the parallelotope bounding the 
        grid.

        Returns
        -------
        volume : `float`
            Volume of the space within the grid boundary.
        """
        return self.axes_volume()
    #end def volume


    def cell_volumes(self):
        """
        (`External API`) Compute the volumes of the parallelotope grid cells.

        Returns
        -------
        cell_vols : `ndarray, shape (N,)`
            Array containing the volume of each grid cell.  `N` is the number 
            of points in the grid.
        """
        ncells = self.ncells
        cell_vols = np.zeros((ncells,),dtype=self.dtype)
        cell_vols += self.volume()/ncells
        return cell_vols
    #end def cell_volumes

#end class ParallelotopeGrid



class SpheroidGrid(StructuredGridWithAxes):
    """
    A regular structured grid within a sheared spheroid volume.

    This type of grid is only defined in two or three (possibly embedded) 
    dimensions.  The grid spans the volume within the spheroidal surface.

    An `M`-dimensional spheroid grid may be embedded in an `N`-dimensional 
    space, e.g. a 2D planar disk with arbitrary orientation in an open 
    3D space.  Such grids are useful to describe planar disk cuts through 
    higher dimensional spaces.

    Below, `M` and `N` are referred to as `dg` and `ds`, or the dimension of 
    the grid (embedded space) and the full space (embedding space) respectively.
    These spaces may be chosen to be the same.

    This class is intended for direct instantiation and use.

    Parameters
    ----------
    axes : `array_like, float, shape (dg,ds)`
        Axes used to form the coordinate system.  The axes are used in place 
        of the normal `x`, `y` (and, if present, `z`) axes so that the 
        resultant spheroid may be spatially skewed.  `dg` is the dimension of 
        the grid (embedded space, `dg=grid_dim`), `ds` is the dimension of the 
        points in the grid (full space, `ds=space_dim`).
    centered : `bool, optional, default False`
        Locate grid points at lower cell corners (`False`) or cell centers 
        (`True`).
    cells : `array_like, int, shape (dg,), optional`
        The number of grid cells in each dimension. An additional grid point 
        will be present at the upper edge of the radial dimension and, if in 
        3D spherical coordinates, also the upper edge of the polar angle 
        dimension (:math:`\\theta`).  `dg` is the dimension of the grid 
        (`grid_dim`).  Either `cells` or `shape` must be provided to generate 
        the grid.
    shape : `array_like, int, shape (dg,), optional`
        The number of grid points in each dimension. The number of grid cells 
        in the azimuthal direction (:math:`\phi`) matches the number of grid 
        points.  In the radial (and, if present, 3D polar) dimension, there is 
        one fewer grid cell than grid points as the grid points go all the way 
        to the edge of the boundary in those directions.  `dg` is the dimension 
        of the grid (`grid_dim`).  Either `cells` or `shape` must be provided 
        to generate the grid.
    center : `array_like, float, shape (ds,), optional`
        The location of the center of the spheroid, which is also the origin 
        of the coordinate system for the grid.  `ds` is the dimension of the 
        full space (`space_dim`). Either `center` or `origin` must be provided,
        with `center` preferred.
    origin : `array_like, float, shape (ds,), optional, default ds*[0]`
        Origin of the space.  For a spheroid grid, this point is the 
        center of the spheroid.  `ds` is the dimension of the full space 
        (`space_dim`). Either `center` or `origin` must be provided, with 
        `center` preferred.

    Attributes
    ----------
    points : `ndarray, float, shape (N,d)`
        Array containing the grid points.  `N` is the number of grid points, 
        `d` is the dimension of the space.
    axes : `ndarray, float, shape (dg,ds)`
        Axes used to form the coordinate system.  `dg` is the dimension of 
        the grid (embedded space, `dg=grid_dim`), `ds` is the dimension of 
        the points in the grid (full space, `ds=space_dim`).
    centered : `bool`
        Grid points are located at lower cell corners (`False`) or cell 
        centers (`True`).
    shape : `tuple, int`
        The number of grid points in each dimension.
    origin : `ndarray, float, shape (ds,)`
        Origin of the space.  `ds` is the dimension of the full space.
    surface : `bool`
        (`Internal`) The grid is resides on the surface of a space (`True`), 
        otherwise, it spans the volume of the space (`False`).  This attribute 
        is intended to be immutable once set.
    bconds : `ndarray, str`
        (`Internal`) Boundary conditions for each dimension.
    isotropic : `bool`
        (`Internal`) Whether or not the spheroid is isotropic (regular sphere).
    initialized : `bool`
        (`Internal`) Set to true if the instance has been initialized in 
        non-vacuous fashion.
    r : `ndarray, float, property`
        Array containing the grid points.  User facing alias for `points`.
    npoints : `int, property`
        Number of grid points.
    space_dim : `int, property`
        Dimension of the space the grid resides in.
    grid_dim : `int, property`
        The dimension of the grid.  Must be less than or equal to `space_dim`.
    grid_shape : `tuple, int, property`
        The number of grid points in each dimension.
    cell_grid_shape : `tuple, int, property`
        The number of grid cells in each dimension.
    ncells : `int, property`
        The total number of grid cells.
    flat_points_shape : `tuple, int, property`
        The shape of the `points` array in its default (flat) representation. 
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `flat_points_shape` is `(N*M*P,D)`.
    full_points_shape : `tuple, int, property`
        The shape of the points array in its full representation.
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `full_points_shape` is `(N,M,P,D)`.
    center : `ndarray, float, shape (ds,), property`
        Location of the center of the cell.  `ds` is the dimension of the
        full space.
    dtype : `property`
        Datatype of the grid point values.
    """

    #: (`obj`) Collection of attributes for the class.  Used to check assigned 
    #: members for type conformity and to assign default values.
    persistent_data_types = obj(
        isotropic = (bool,None),
        **StructuredGridWithAxes.persistent_data_types
        )

    @property
    def center(self):
        return self.origin
    #end def center

    def initialize_local(self,
                         axes     = None,
                         shape    = None,
                         cells    = None,
                         center   = None,
                         centered = False,
                         **kwargs
                         ):
        """
        (`Internal API`)  Initialize the spheroid grid points.
        """
        if shape is None and cells is None:
            self.error('cannot initialize grid, either "shape" or "cells" is required')
        elif shape is not None:
            grid_dim = len(shape)
        elif cells is not None:
            grid_dim = len(cells)
        #end if
        if grid_dim not in (2,3):
            self.error('only 2 and 3 dimensional spheroids grids are supported\nrequested dimension: {}'.format(grid_dim))
        #end if

        # bconds match sphere constraints by default
        bconds = kwargs.get('bconds',None)
        if bconds is None:
            if grid_dim==3:
                bconds = tuple('oop')
            elif grid_dim==2:
                bconds = tuple('op')
            #end if
        #end if

        endpoint = self.has_endpoints(bconds=bconds,grid_dim=grid_dim)

        points,shape = spheroid_grid_points(axes,shape=shape,cells=cells,centered=centered,endpoint=endpoint,return_shape=True)

        kwargs['axes']     = axes
        if center is not None:
            kwargs['origin']   = center
        #end if
        kwargs['shape']    = shape
        kwargs['centered'] = centered
        kwargs['bconds']   = bconds
        kwargs['points']   = points
        StructuredGridWithAxes.initialize_local(self,**kwargs)

    #end def initialize_local


    def set_axes(self,axes):
        """
        (`Internal API`) Sets the `axes` attribute in a protected way.
        """
        StructuredGridWithAxes.set_axes(self,axes)
        self.set_isotropic()
    #end def set_axes


    def set_isotropic(self,tol=1e-6):
        """
        (`Internal API`) Determine whether the spheroid is isotropic (constant 
        radius) and internally store the result.

        Parameters
        ----------
        tol : `float, optional, default 1e-6`
            Tolerance used to judge equal axis length and axis orthogonality.
        """
        isotropic = True
        ax_norm = np.linalg.norm(self.axes[0])
        for i,ax1 in enumerate(self.axes):
            axn1 = np.linalg.norm(ax1)
            for j,ax2 in enumerate(self.axes):
                axn2 = np.linalg.norm(ax2)
                if i==j:
                    isotropic &= np.abs(axn1-ax_norm)/ax_norm < tol
                else:
                    isotropic &= np.abs(np.dot(ax1,ax2))/(axn1*axn2) < tol
                #end if
            #end for
        #end for
        self.isotropic = bool(isotropic)
    #end def set_isotropic


    def unit_points_bare(self,points):
        """
        (`Internal API`)  Maps points in a spheroid into the unit cube.
        """
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
        """
        (`External API`) Map points from the unit space back into the full space.

        Parameters
        ----------
        upoints : `array_like, float, shape (N,dg), optional`
            Array of points in the unit coordinate space.  `N` is the number 
            of points and `dg` must be equal to `grid_dim`.

        Returns
        -------
        points : `ndarray, float, shape (N,ds)`
            Array of points in the full coordinate space.  `ds` is the 
            dimension of the full space (`space_dim`).
        """
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


    def unit_metric_bare(self,upoints):
        """
        (`Internal API`) Compute the spheroid integration metric in the 
        unit coordinate space for a set of points defined there.
        """
        if upoints is None:
            upoints = self.unit_points()
        else:
            upoints = np.array(upoints,dtype=self.dtype)
        #end if
        dim = self.grid_dim
        umetric = np.zeros((len(upoints),),dtype=self.dtype)
        r = upoints[:,0]
        if dim==2:
            umetric += 2*np.pi*r # r
        elif dim==3:
            umetric += 2*np.pi**2*r**2*np.sin(np.pi*upoints[:,1]) # r^2 sin(theta)
        else:
            self.error('unit_metric not supported for dim={}'.format(dim))
        #end if
        umetric *= self.axes_volume()
        return umetric
    #end def unit_metric_bare


    def radius(self):
        """
        (`External API`) Return the radius of the spheroid, if isotropic. 
        """
        if not self.isotropic:
            self.error('radius is not supported for anisotropic spheroid surface grids')
        #end if
        return np.linalg.norm(self.axes[0])
    #end def radius


    def volume(self):
        """
        (`External API`) Compute the volume of the spheroid bounding the grid.

        Returns
        -------
        volume : `float`
            Volume of the space within the grid boundary.
        """
        vol_axes = self.axes_volume()
        if self.grid_dim==2:
            vol_spheroid = np.pi # circle of radius 1
        elif self.grid_dim==3:
            vol_spheroid = 4*np.pi/3 # sphere of radius 1
        else:
            self.error('volume not supported for dim={}'.format(self.grid_dim))
        #end if
        return vol_axes*vol_spheroid
    #end def volume


    def cell_volumes(self):
        """
        (`External API`) Compute the volumes of the spheroid grid cells.

        Returns
        -------
        cell_vols : `ndarray, shape (N,)`
            Array containing the volume of each grid cell.  `N` is the number 
            of points in the grid.
        """
        vol_axes = self.axes_volume()
        dim = self.grid_dim
        shape = self.cell_grid_shape
        ncells = self.ncells
        ugrid = unit_grid_points(shape,centered=True)
        du = np.ones((dim,),dtype=self.dtype)/shape
        cell_vols = np.zeros((ncells,),dtype=self.dtype)
        if dim==2:
            ugrid[:,1] *= 2*np.pi
            du[1] *= 2*np.pi
            for i,uc in enumerate(ugrid):
                cv  = uc[0]*du[0]*du[1]
                cell_vols[i] =  cv*vol_axes
            #end for
        elif dim==3:
            ugrid[:,1] *=   np.pi
            ugrid[:,2] *= 2*np.pi
            du[1] *=   np.pi
            du[2] *= 2*np.pi
            for i,uc in enumerate(ugrid):
                cv  = (uc[0]*uc[0]+du[0]*du[0]/12.0)*du[0] # r
                cv *= 2.0*np.sin(uc[1])*np.sin(.5*du[1])   # theta
                cv *= du[2]                                # phi
                cell_vols[i] = cv*vol_axes
            #end for
        else:
            self.error('cell_volumes not supported for dim={}'.format(dim))
        #end if
        return cell_vols
    #end def cell_volumes


    def radii(self):
        self.reshape_full()
        rrad = np.array(self.r[:,0,0,-1].ravel())
        self.reshape_flat()
        return rrad
    #end def radii
#end class SpheroidGrid



class SpheroidSurfaceGrid(StructuredGridWithAxes):
    """
    A regular structured grid over the surface of a sheared spheroid volume.

    This type of grid is only defined in one or two embedded dimensions.  The 
    grid spans the area over the spheroidal surface.

    An `M`-dimensional spheroid surface grid may be embedded in an 
    `N`-dimensional space, e.g. a 2D spherical surface or 1D circular ring 
    with arbitrary orientation in an open 3D space.  Such grids are useful to 
    describe curved surface cuts through higher dimensional spaces.

    Below, `M` and `N` are referred to as `dg` and `ds`, or the dimension of 
    the grid (embedded space) and the full space (embedding space) respectively.
    These spaces may be chosen to be the same.

    This class is intended for direct instantiation and use.

    Parameters
    ----------
    axes : `array_like, float, shape (dg,ds)`
        Axes used to form the coordinate system.  The axes are used in place 
        of the normal `x`, `y` (and, if present, `z`) axes so that the 
        resultant spheroid may be spatially skewed.  `dg` is the dimension of 
        the grid (embedded space, `dg=grid_dim`), `ds` is the dimension of the 
        points in the grid (full space, `ds=space_dim`).
    centered : `bool, optional, default False`
        Locate grid points at lower cell corners (`False`) or cell centers 
        (`True`).
    cells : `array_like, int, shape (dg,), optional`
        The number of grid cells in each (angular) dimension. In 3D spherical 
        coordinates, an additional grid point will be present at the upper edge 
        of the polar angle dimension (:math:`\\theta`).  `dg` is the dimension 
        of the grid (`grid_dim`).  Either `cells` or `shape` must be provided 
        to generate the grid.
    shape : `array_like, int, shape (dg,), optional`
        The number of grid points in each dimension. The number of grid cells 
        in the azimuthal direction (:math:`\phi`) matches the number of grid 
        points.  In 3D spherical coordinates, along the polar dimension there is
        one fewer grid cell than grid points as the grid points go all the way 
        to the edge of the boundary in that direction.  `dg` is the dimension 
        of the grid (`grid_dim`).  Either `cells` or `shape` must be provided 
        to generate the grid.
    center : `array_like, float, shape (ds,), optional`
        The location of the center of the spheroid.  `ds` is the dimension of 
        the full space (`space_dim`). Either `center` or `origin` must be 
        provided, with `center` preferred.
    origin : `array_like, float, shape (ds,), optional, default ds*[0]`
        Origin of the space, which resides at the center of the spheroid. `ds` 
        is the dimension of the full space (`space_dim`). Either `center` or 
        `origin` must be provided, with `center` preferred.

    Attributes
    ----------
    points : `ndarray, float, shape (N,d)`
        Array containing the grid points.  `N` is the number of grid points, 
        `d` is the dimension of the space.
    axes : `ndarray, float, shape (dg,ds)`
        Axes used to form the coordinate system.  `dg` is the dimension of 
        the grid (embedded space, `dg=grid_dim`), `ds` is the dimension of 
        the points in the grid (full space, `ds=space_dim`).
    centered : `bool`
        Grid points are located at lower cell corners (`False`) or cell 
        centers (`True`).
    shape : `tuple, int`
        The number of grid points in each dimension.
    origin : `ndarray, float, shape (ds,)`
        Origin of the space.  `ds` is the dimension of the full space.
    surface : `bool`
        (`Internal`) The grid is resides on the surface of a space (`True`), 
        otherwise, it spans the volume of the space (`False`).  This attribute 
        is intended to be immutable once set.
    bconds : `ndarray, str`
        (`Internal`) Boundary conditions for each dimension.
    isotropic : `bool`
        (`Internal`) Whether or not the spheroid is isotropic (regular sphere).
    initialized : `bool`
        (`Internal`) Set to true if the instance has been initialized in 
        non-vacuous fashion.
    r : `ndarray, float, property`
        Array containing the grid points.  User facing alias for `points`.
    npoints : `int, property`
        Number of grid points.
    space_dim : `int, property`
        Dimension of the space the grid resides in.
    grid_dim : `int, property`
        The dimension of the grid.  Must be less than or equal to `space_dim`.
    grid_shape : `tuple, int, property`
        The number of grid points in each dimension.
    cell_grid_shape : `tuple, int, property`
        The number of grid cells in each dimension.
    ncells : `int, property`
        The total number of grid cells.
    flat_points_shape : `tuple, int, property`
        The shape of the `points` array in its default (flat) representation. 
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `flat_points_shape` is `(N*M*P,D)`.
    full_points_shape : `tuple, int, property`
        The shape of the points array in its full representation.
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `full_points_shape` is `(N,M,P,D)`.
    center : `ndarray, float, shape (ds,), property`
        Location of the center of the cell.  `ds` is the dimension of the
        full space.
    dtype : `property`
        Datatype of the grid point values.
    """

    #: (`obj`) Collection of attributes for the class.  Used to check assigned 
    #: members for type conformity and to assign default values.
    persistent_data_types = obj(
        isotropic = (bool,None),
        **StructuredGridWithAxes.persistent_data_types
        )

    @property
    def center(self):
        return self.origin
    #end def center

    def initialize_local(self,
                         axes     = None,
                         shape    = None,
                         cells    = None,
                         center   = None,
                         centered = False,
                         **kwargs
                         ):
        """
        (`Internal API`)  Initialize the spheroid surface grid points.
        """
        if shape is None and cells is None:
            self.error('cannot initialize grid, either "shape" or "cells" is required')
        elif shape is not None:
            grid_dim = len(shape)
        elif cells is not None:
            grid_dim = len(cells)
        #end if

        # force bconds to match sphere surface constraints
        if grid_dim==2:
            bconds = tuple('op')
        elif grid_dim==1:
            bconds = tuple('p')
        else:
            self.error('only 1 and 2 dimensional spheroid surface grids are supported\nrequested dimension: {}'.format(grid_dim))
        #end if

        endpoint = self.has_endpoints(bconds=bconds,grid_dim=grid_dim)

        points,shape = spheroid_surface_grid_points(axes,shape=shape,cells=cells,centered=centered,endpoint=endpoint,return_shape=True)

        kwargs['axes']     = axes
        if center is not None:
            kwargs['origin']   = center
        #end if
        kwargs['shape']    = shape
        kwargs['centered'] = centered
        kwargs['bconds']   = bconds
        kwargs['points']   = points
        StructuredGridWithAxes.initialize_local(self,**kwargs)

        self.surface = True

    #end def initialize_local


    def set_axes(self,axes):
        """
        (`Internal API`) Sets the `axes` attribute in a protected way.
        """
        StructuredGridWithAxes.set_axes(self,axes)
        self.set_isotropic()
    #end def set_axes


    def set_isotropic(self,tol=1e-6):
        """
        (`Internal API`) Determine whether the spheroid is isotropic (constant 
        radius) and internally store the result.

        Parameters
        ----------
        tol : `float, optional, default 1e-6`
            Tolerance used to judge equal axis length and axis orthogonality.
        """
        isotropic = True
        ax_norm = np.linalg.norm(self.axes[0])
        for i,ax1 in enumerate(self.axes):
            axn1 = np.linalg.norm(ax1)
            for j,ax2 in enumerate(self.axes):
                axn2 = np.linalg.norm(ax2)
                if i==j:
                    isotropic &= np.abs(axn1-ax_norm)/ax_norm < tol
                else:
                    isotropic &= np.abs(np.dot(ax1,ax2))/(axn1*axn2) < tol
                #end if
            #end for
        #end for
        self.isotropic = bool(isotropic)
    #end def set_isotropic


    def unit_points_bare(self,points):
        """
        (`Internal API`)  Maps points on a spheroid surface into the unit cube.
        """
        center = self.center
        # invert using pseudo-inverse
        #   this is important for grids embedded in higher dim spaces
        axinv  = np.linalg.pinv(self.axes)
        upoints = np.dot(points-center,axinv)
        # map from unit cartesian to unit spherical
        dim = self.grid_dim
        if dim==1:
            upoints = cartesian_to_polar(upoints,surface=True)
            upoints[:,0] /= 2*np.pi
        elif dim==2:
            upoints = cartesian_to_spherical(upoints,surface=True)
            upoints[:,0] /=   np.pi
            upoints[:,1] /= 2*np.pi
        else:
            self.error('unit_points not supported for dim={}'.format(dim))
        #end if
        return upoints
    #end def unit_points_bare


    def points_from_unit(self,upoints):
        """
        (`External API`) Map points from the unit space back into the full space.

        Parameters
        ----------
        upoints : `array_like, float, shape (N,dg), optional`
            Array of points in the unit coordinate space.  `N` is the number 
            of points and `dg` must be equal to `grid_dim`.

        Returns
        -------
        points : `ndarray, float, shape (N,ds)`
            Array of points in the full coordinate space.  `ds` is the 
            dimension of the full space (`space_dim`).
        """
        dim = self.grid_dim
        upoints = np.array(upoints,dtype=self.dtype)
        if dim==1:
            upoints[:,0] *= 2*np.pi
            upoints = polar_to_cartesian(upoints,surface=True)
        elif dim==2:
            upoints[:,0] *=   np.pi
            upoints[:,1] *= 2*np.pi
            upoints = spherical_to_cartesian(upoints,surface=True)
        else:
            self.error('points_from_unit not supported for dim={}'.format(dim))
        #end if
        points = np.dot(upoints,self.axes)+self.center
        return points
    #end def points_from_unit


    def unit_metric_bare(self,upoints):
        """
        (`Internal API`) Compute the spheroid surface integration metric in the 
        unit coordinate space for a set of points defined there.
        """
        if not self.isotropic:
            self.error('unit_metric is not supported for anisotropic spheroid surface grids')
        #end if
        if upoints is None:
            upoints = self.unit_points()
        else:
            upoints = np.array(upoints,dtype=self.dtype)
        #end if
        dim = self.grid_dim
        umetric = np.zeros((len(upoints),),dtype=self.dtype)
        r = self.radius()
        if dim==1:
            umetric += 2*np.pi*r
        elif dim==2:
            umetric += 2*np.pi**2*r**2*np.sin(np.pi*upoints[:,1]) # r^2 sin(theta)
        else:
            self.error('unit_metric not supported for dim={}'.format(dim))
        #end if
        return umetric
    #end def unit_metric_bare


    def radius(self):
        """
        (`External API`) Return the radius of the spheroid, if isotropic. 
        """
        if not self.isotropic:
            self.error('radius is not supported for anisotropic spheroid surface grids')
        #end if
        return np.linalg.norm(self.axes[0])
    #end def radius


    def volume(self):
        """
        (`External API`) Compute the area of the spheroid surface containing the 
        grid (isotropic only).

        Returns
        -------
        volume : `float`
            Volume of the space within the grid boundary.
        """
        if not self.isotropic:
            self.error('volume is not supported for anisotropic spheroid surface grids')
        #end if
        dim = self.grid_dim
        r = self.radius()
        if dim==1:
            vol = 2*np.pi*r
        elif dim==2:
            vol = 4*np.pi*r**2
        else:
            self.error('volume is not supported for dim={}'.format(dim))
        #end if
        return vol
    #end def volume


    def cell_volumes(self):
        """
        (`External API`) Compute the areas of the spheroid surface grid cells 
        (isotropic only).

        Returns
        -------
        cell_vols : `ndarray, shape (N,)`
            Array containing the area of each grid cell.  `N` is the number 
            of points in the grid.
        """
        if not self.isotropic:
            self.error('volume is not supported for anisotropic spheroid surface grids')
        #end if
        dim = self.grid_dim
        ncells = self.ncells
        cell_vols = np.zeros((ncells,),dtype=self.dtype)
        if dim==1:
            cell_vols += self.volume()/ncells
        elif dim==2:
            r = self.radius()
            shape = self.cell_grid_shape
            ugrid = unit_grid_points(shape,centered=True)
            du = np.ones((dim,),dtype=self.dtype)/shape
            ugrid[:,0] *=   np.pi
            ugrid[:,1] *= 2*np.pi
            du[0] *=   np.pi
            du[1] *= 2*np.pi
            for i,uc in enumerate(ugrid):
                cv  = r**2                                # r
                cv *= 2.0*np.sin(uc[0])*np.sin(.5*du[0])  # theta
                cv *= du[1]                               # phi
                cell_vols[i] = cv
            #end for
        else:
            self.error('cell_volumes is not supported for dim={}'.format(dim))
        #end if
        return cell_vols
    #end def cell_volumes
#end class SpheroidSurfaceGrid



class GridFunction(GBase):
    """
    Base class for `P` dimensional functions defined over `M` dimensional grids
    that are embedded in `N` dimensional spaces.

    The main aims of this class hierarchy are to provide interpolation, 
    integration, and (where possible) plotting capabilities for functions of 
    this type.  These operations are common in large scale post-processing 
    and analysis of scientific data.

    This class should not be instantiated directly.

    Parameters
    ----------
    grid : `Grid, optional`
        Grid of points in a `d` dimensional space.  If `grid` is not provided, 
        additional parameters must be given to initialize a Grid object.
    values : `array_like, float/complex, shape (N,P), (N,)`
        Array of function values defined at the grid points.  `N` is the number
        of points and `P` is the number of function values.  With `P>1`, the 
        function is vector or tensor valued.  If the array has shape (`N`,), 
        then `P` is set to `1`.
    copy_grid : `bool, optional, default True`
        Copy provided grid (`True`) or not (`False`).
    copy_values : `bool, optional, default True`
        Copy provided values (`True`) or not (`False`).
    copy : `bool, optional`
        Copy provided grid and values (`True`) or not (`False`).
    dtype : `optional`
        Data type for local function values.
    grid_dtype : `optional`
        Data type for grid point locations.
    **kwargs: 
        Arbitrary set of parameters used to create a Grid object.  See 
        documentation for the `Grid` class and its derived classes for allowed 
        inputs.  Used/allowed only if `grid` is not provided.

    Attributes
    ----------
    grid : `Grid`
        Grid of points in a `d` dimensional space.
    values : `ndarray, float/complex, shape (N,P)`
        Array of function values defined at the grid points.  `N` is the number
        of points and `P` is the number of function values.  With `P>1`, the 
        function is vector or tensor valued.
    space_dim : `int, property`
        Dimension of the space the grid resides in.  Referred to as `d` above.
    npoints : `int, property`
        Number of grid points.  Referred to as `N` above.
    nvalues : `int, property`
        Number of function values at each grid point.  Referred to as `P` above.
    r : `ndarray, float, property`
        Array containing the grid points.
    f : `ndarray, float/complex, property`
        Array containing the function values.  User-facing alias for `values`.
    dtype : `property`
        Datatype of the function values.
    """

    #: Descriptive string for class.  Used in the GBase base class when 
    #: printing error messages.
    descriptor = 'grid function'

    #: Grid class type that must be associated (contained by) a particular
    #: `GridFunction` class.  Must be a sub-class of `Grid`.  Required only 
    #: for grid function classes that support direct instantiation.
    grid_class = None

    #: (`obj`) Collection of attributes for the class.  Used to check assigned 
    #: members for type conformity and to assign default values.
    persistent_data_types = obj(
        grid        = (Grid      , None),
        values      = (np.ndarray, None),
        value_shape = (tuple     , None),
        **GBase.persistent_data_types
        )

    @property
    def space_dim(self):
        return self.grid.space_dim
    #end def space_dim

    @property
    def npoints(self):
        return self.grid.npoints
    #end def npoints

    @property
    def value_dim(self):
        return len(self.value_shape)
    #end def value_dim

    @property
    def nvalues(self):
        return np.prod(self.value_shape)
    #end def nvalues

    @property
    def r(self):
        return self.grid.r
    #end def r

    @property 
    def f(self):
        return self.values
    #end def f

    @property
    def dtype(self):
        return self.values.dtype
    #end def dtype


    def initialize_local(self,
                         grid        = None,
                         values      = None,
                         copy        = None,
                         copy_grid   = True,
                         copy_values = True,
                         dtype       = None,
                         grid_dtype  = None,
                         value_shape = None,
                         **kwargs):
        """
        (`Internal API`) Sets `grid` and `values` attributes.

        Parameters
        ----------
        grid : `Grid, optional`
            Grid of points in a `d` dimensional space.  If `grid` is not 
            provided, additional parameters must be given to initialize a 
            `Grid` object.
        values : `array_like, float/complex, shape (N,P), (N,)`
            Array of function values defined at the grid points.  `N` is the 
            number of points and `P` is the number of function values.  With 
            `P>1`, the function is vector or tensor valued.  If the array 
            has shape (`N`,), then `P` is set to `1`.
        copy_grid : `bool, optional, default True`
            Copy provided grid (`True`) or not (`False`).
        copy_values : `bool, optional, default True`
            Copy provided values (`True`) or not (`False`).
        copy : `bool, optional`
            Copy provided grid and values (`True`) or not (`False`).
        dtype : `optional`
            Data type for local function values.
        grid_dtype : `optional`
            Data type for grid point locations.
        **kwargs: 
            Arbitrary set of parameters used to create a `Grid` object.  See 
            documentation for the `Grid` class and its derived classes for 
            allowed inputs.  Used/allowed only if `grid` is not provided.
        """

        if grid_dtype is not None:
            kwargs['dtype'] = grid_dtype
        #end if

        # process copy inputs
        if copy is not None:
            copy_grid   = copy
            copy_values = copy
        #end if
        
        # process grid inputs
        cls = self.__class__
        if grid is None:
            grid = cls.grid_class(**kwargs)
            if not grid.initialized:
                if values is not None:
                    self.error('grid must be initialized to describe the domain of the function values\nplease add inputs to initialize the appropriate grid')
                #end if
                return
            #end if
        elif not isinstance(grid,cls.grid_class):
            self.error('received "grid" input with incorrect type\ntype provided: {}\ntype required: {}'.format(grid.__class__.__name__,cls.grid_class.__name__))
        elif len(kwargs)>0:
            self.error('received both a grid object and parameters intended for grid initialization\nplease remove the following parameters and try again: {}'.format(sorted(kwargs.keys())))
        elif copy_grid:
            grid = grid.copy()
        #end if

        # process values input
        if values is None:
            self.error('function values must be provided at the grid points')
        elif isinstance(values,(tuple,list)):
            if dtype is None:
                dtype = float
            #end if
            values = np.array(values,dtype=dtype)
        elif isinstance(values,np.ndarray):
            if copy_values:
                if dtype is None:
                    dtype = values.dtype
                #end if
                values = np.array(values,dtype=dtype)
            #end if
        else:
            self.error('provided function values are of incorrect type\nvalues must be tuple, list, or ndarray\nyou provided: {}'.format(values.__class__.__name__))
        #end if
        
        # process value_shape input
        if len(values.shape)==1 or values.shape==grid.shape:
            value_shape = (1,)
        elif value_shape is None:
            if len(values)==grid.npoints:
                value_shape = values.shape[1:]
            elif len(value.shape)>len(grid.shape) and value.shape[:len(grid.shape)]==grid.shape:
                value_shape = values.shape[len(grid.shape):]
            #end if
        elif isinstance(value_shape,(list,np.ndarray)):
            value_shape = tuple(value_shape)
        #end if
        if value_shape is not None:
            nvtot = values.size
            nv    = np.prod(value_shape)
            if nvtot%nv!=0 or nvtot//nv!=grid.npoints:
                self.error('value_shape and total number of values are inconsistent.\nTotal number of values: {}\nvalue_shape: {}\nExpected number of values per grid point: {}\nActual number of values per grid point: {}'.format(nvtot,value_shape,nv,nvtot/nv))
            #end if
            values.shape = (grid.npoints,nv)
        #end if

        # assign grid and values
        self.grid        = grid
        self.values      = values
        self.value_shape = value_shape
    #end def initialize_local


    def local_validity_checks(self,msgs):
        """
        (`Internal API`) Check validity of `grid` and `values`.

        Parameters
        ----------
        msgs : `list, str`
            List containing error messages.  Empty if no problems are found.
        """
        cls = self.__class__
        if not isinstance(self.grid,cls.grid_class):
            msgs.append('Grid is not of the required type for current grid function.\nGrid function type: {}\nGrid type required: {}'.format(cls.__name__,self.grid.__class__.__name__))
        #end if
        self.grid.local_validity_checks(msgs)
        if len(self.values)!=self.npoints:
            msgs.append('Number of function values and number of grid points do not match.\nNumber of grid points: {}\nNumber of function values: {}'.format(self.npoints,len(self.values)))
        #end if
        if len(self.values.shape)!=2:
            msgs.append('Function values has incorrect shape.\nExpected shape is (# of points, # of function values at each point)\nShape received: {}'.format(self.values.shape))
        #end if
        if len(self.value_shape)<1:
            self.error('"value_shape" must have at least one entry.')
        #end if
        if np.prod(self.value_shape)!=self.values.size//self.npoints:
            self.error('"value_shape" and "values" are inconsistent.\nNumber of values per point based on "values": {}\nNumber of values per point based on "value_shape": {}'.format(self.values.size/self.npoints,np.prod(self.value_shape)))
        #end if
        if self.values.shape!=(self.npoints,self.nvalues):
            self.error('Function values has incorrect shape.\nExected shape: {}\nShape received: {}'.format((self.npoints,self.nvalues),self.values.shape))
        #end if
    #end def local_validity_checks


    def reshape_values_full(self):
        self.values.shape = (self.npoints,)+self.value_shape
    #end def reshape_values_full


    def reshape_values_flat(self):
        self.values.shape = (self.npoints,self.nvalues)
    #end def reshape_values_flat
#end class GridFunction



class StructuredGridFunction(GridFunction):
    """
    Base class for functions defined on structured grids.

    This class handles plotting functions within the unit coordinate space.  
    It will handle unified interpolation and integration of (potentially 
    multi-valued) discrete functions.

    This class should not be instantiated directly.

    Attributes
    ----------
    grid_dim : `int, property`
        The dimension of the grid.
    grid_shape : `tuple, int, property`
        The number of grid points in each dimension.
    cell_grid_shape : `tuple, int, property`
        The number of grid cells in each dimension.
    ncells : `int, property`
        The total number of grid cells.
    flat_points_shape : `tuple, int, property`
        The shape of the `points` array in its default (flat) representation. 
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `flat_points_shape` is `(N*M*P,D)`.
    full_points_shape : `tuple, int, property`
        The shape of the points array in its full representation.
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `full_points_shape` is `(N,M,P,D)`.
    """

    @property
    def grid_dim(self):
        return self.grid.grid_dim
    #end def grid_dim

    @property
    def grid_shape(self):
        return self.grid.grid_shape
    #end def grid_shape

    @property
    def cell_grid_shape(self):
        return self.grid.cell_grid_shape
    #end def cell_grid_shape

    @property
    def ncells(self):
        return self.grid.ncells
    #end def ncells

    @property
    def flat_points_shape(self):
        return self.grid.flat_points_shape
    #end def flat_points_shape

    @property
    def full_points_shape(self):
        return self.grid.full_points_shape
    #end def full_points_shape

    @property
    def flat_values_shape(self):
        None
    #end def flat_values_shape

    @property
    def periodic(self):
        return self.grid.periodic
    #end def periodic
    

    def reshape_points_full(self):
        self.values.shape = self.grid_shape+(self.nvalues,)
        self.grid.reshape_full()
    #end def reshape_points_full


    def reshape_points_flat(self):
        self.values.shape = (self.npoints,self.nvalues)
        self.grid.reshape_flat()
    #end def reshape_points_flat


    def reshape_full(self):
        self.values.shape = self.grid_shape+self.value_shape
        self.grid.reshape_full()
    #end def reshape_full


    def reshape_flat(self):
        self.values.shape = (self.npoints,self.nvalues)
        self.grid.reshape_flat()
    #end def reshape_flat


    def get_values_with_upper_ghost(self):
        if 'values_with_upper_ghost' in self:
            return self.values_with_upper_ghost
        #end if
        self.reshape_points_full()

        g = np.array(self.grid_shape)
        v = np.empty(tuple(g+1)+(self.nvalues,),dtype=self.values.dtype)
        dim = self.grid_dim
        if dim==1:
            n1 = self.npoints
            v[:-1] = self.values
            v[-1]  = self.values[0]
        elif dim==2:
            v[:-1,:-1] = self.values
            v[ -1,:-1] = self.values[0,:]
            v[:-1, -1] = self.values[:,0]
            v[ -1, -1] = self.values[0,0]
        elif dim==3:
            v[:-1,:-1,:-1] = self.values
            v[ -1,:-1,:-1] = self.values[0,:,:]
            v[:-1, -1,:-1] = self.values[:,0,:]
            v[:-1,:-1, -1] = self.values[:,:,0]
            v[:-1, -1, -1] = self.values[:,0,0]
            v[ -1,:-1, -1] = self.values[0,:,0]
            v[ -1, -1,:-1] = self.values[0,0,:]
            v[ -1, -1, -1] = self.values[0,0,0]
        else:
            self.error('values_with_upper_ghost is not implemented for dimensions greater than 3.\nDimensionality of the current dataset: {}'.format(dim))
        #end if

        self.reshape_points_flat()

        self.values_with_upper_ghost = v

        return v
    #end def values_with_upper_ghost


    def clear_ghost(self):
        ghost_fields = [
            'values_with_upper_ghost',
            ]
        for f in ghost_fields:
            if f in self:
                del self[f]
            #end if
        #end for
    #end def clear_ghost


    def plot_unit_contours(self,boundary=False,fig=True,show=True,**kwargs):
        """
        (`External API`) Make 2D contour plots in the unit coordinate space.

        Parameters
        ----------
        boundary : `bool, optional, default False`
            Draw lines at the grid boundaries (`True`) or not (`False`).
        fig : `bool, optional, default True`
            Create a fresh figure (`True`) or not (`False`).
        show : `bool, optional, default True`
            Display the figure (`True`) or not (`False`).
        **kwargs : `optional`
            Arbitrary keyword argments passed to `pyplot.contour`.
        """
        if self.grid_dim!=2:
            self.error('cannot plot contours in unit coordinates\ngrid must have dimension 2 to make contour plots\ndimension of grid for this function: {}'.format(self.grid_dim))
        #end if
        X,Y = self.grid.unit_points().T
        X.shape = self.grid_shape
        Y.shape = self.grid_shape
        Zm = self.f.T
        for Z in Zm:
            Z.shape = self.grid_shape
            fig,ax = self.setup_mpl_fig(fig=fig,dim=self.grid_dim)
            ax.contour(X,Y,Z,**kwargs)
            if boundary:
                self.grid.plot_unit_boundary(fig=False,show=False)
            #end if
        #end for
        if show:
            plt.show()
        #end if
    #end def plot_unit_contours
    

    def plot_unit_surface(self,fig=True,show=True,**kwargs):
        """
        (`External API`) Make 2D surface plots in the unit coordinate space.

        Parameters
        ----------
        fig : `bool, optional, default True`
            Create a fresh figure (`True`) or not (`False`).
        show : `bool, optional, default True`
            Display the figure (`True`) or not (`False`).
        **kwargs : `optional`
            Arbitrary keyword argments passed to `pyplot.plot_surface`.
        """
        if self.grid_dim!=2:
            self.error('cannot plot contours in unit coordinates\ngrid must have dimension 2 to make contour plots\ndimension of grid for this function: {}'.format(self.grid_dim))
        #end if
        X,Y = self.grid.unit_points().T
        X.shape = self.grid_shape
        Y.shape = self.grid_shape
        Zm = self.f.T
        for Z in Zm:
            Z.shape = self.grid_shape
            fig,ax = self.setup_mpl_fig(fig=fig,dim=self.grid_dim+1)
            ax.plot_surface(X,Y,Z,**kwargs)
        #end for
        if show:
            plt.show()
        #end if
    #end def plot_unit_surface
    

    def plot_unit_isosurface(self,level=None,fig=True,show=True,**kwargs):
        """
        (`External API`) Make 3D isosurface plots in the unit coordinate space.

        Parameters
        ----------
        level : `float, optional`
            Isosurface value to plot.  If not provided, the average of the max 
            and min function values are used.
        fig : `bool, optional, default True`
            Create a fresh figure (`True`) or not (`False`).
        show : `bool, optional, default True`
            Display the figure (`True`) or not (`False`).
        **kwargs : `optional`
            Arbitrary keyword argments passed to `pyplot.plot_trisurf`.
        """
        if self.grid_dim!=3:
            self.error('cannot plot isosurface in unit coordinates\ngrid must have dimension 3 to make contour plots\ndimension of grid for this function: {}'.format(self.grid_dim))
        #end if
        fm = self.f.T
        # unit grid is used implicitly
        spacing = tuple(1./(np.array(self.cell_grid_shape,dtype=float)))
        for f in fm:
            if level is None:
                level = (f.max()+f.min())/2
            #end if
            f.shape = self.grid_shape
            ret = measure.marching_cubes(f,level,spacing=spacing)
            verts = ret[0] 
            faces = ret[1]
            fig,ax = self.setup_mpl_fig(fig=fig,dim=self.grid_dim)
            ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],**kwargs)
        #end for
        if show:
            plt.show()
        #end if
    #end def plot_unit_isosurface

#end class StructuredGridFunction



class StructuredGridFunctionWithAxes(StructuredGridFunction):
    """
    Base class for functions over structured grids with linear axes that act 
    as scaffolding for the coordinate system.

    The associated Grids for these classes can handle the mapping back and 
    forth between the unit and full coordinate spaces.  As such, this class
    handles the plotting of function values in the full space.

    This class should not be instantiated directly.
    """


    def interpolate(self,r,type=None,copy=False,**kw):
        # https://stackoverflow.com/questions/16217995/fast-interpolation-of-regularly-sampled-3d-data-with-different-intervals-in-x-y
        kw = obj(kw)
        grid = None
        if isinstance(r,Grid):
            grid = r
            r = grid.r
        #end if
        if type is None:
            type = 'map_coordinates'
        #end if
        if type=='map_coordinates':
            if 'mode' not in kw:
                if self.periodic:
                    kw.mode = 'wrap'
                else:
                    kw.mode = 'nearest'
                #end if
            #end if
            if 'order' not in kw:
                kw.order = 3
            #end if
            indices = self.grid.indices_for_map_coord(r)
            # needed because of off-by-one error in map_coordinates
            #  see: https://github.com/scipy/scipy/issues/2640
            v = self.get_values_with_upper_ghost()
            if self.nvalues>1:
                self.error('Interpolation is not yet supported for nvalues>1.')
            #end if
            v_shape = v.shape
            v.shape = v_shape[:-1]
            values = scipy_ndimage.map_coordinates(v, indices, **kw)
            v.shape = v_shape
        else:
            self.error('Interpolation of type "{}" is not supported.\nValid options are: map_coordinates'.format(type))
        #end if

        if grid is None:
            return values
        elif grid is not None:
            if copy:
                grid = grid.copy()
            #end if
            gf = grid.grid_function(
                grid        = grid,
                values      = values,
                value_shape = tuple(self.value_shape),
                copy        = False,
                )
            return gf
        #end if
    #end def interpolate


    def plot_contours(self,a1=(1,0),a2=(0,1),boundary=False,fig=True,show=True,**kwargs):
        """
        (`External API`) Make 2D contour plots in the full coordinate space.

        Parameters
        ----------
        boundary : `bool, optional, default False`
            Draw lines at the grid boundaries (`True`) or not (`False`).
        fig : `bool, optional, default True`
            Create a fresh figure (`True`) or not (`False`).
        show : `bool, optional, default True`
            Display the figure (`True`) or not (`False`).
        **kwargs : `optional`
            Arbitrary keyword argments passed to `pyplot.contour`.
        """
        if self.grid_dim!=2:
            self.error('cannot plot contours\ngrid must have dimension 2 to make contour plots\ndimension of grid for this function: {}'.format(self.grid_dim))
        #end if
        #if self.space_dim!=2:
        #    self.error('cannot plot contours\ngrid points must reside in a 2D space to make contour plots\ndimension of the space for this function: {}'.format(self.space_dim))
        ##end if
        if self.space_dim==2:
            X,Y = self.r.T
        else:
            ax     = np.dot(np.array([a1,a2]),self.grid.axes)
            ax[0]  = ax[0]/np.linalg.norm(ax[0])
            ax[1] -= np.dot(ax[1],ax[0])*ax[0]
            ax[1]  = ax[1]/np.linalg.norm(ax[1])
            X,Y    = np.dot(ax,self.r.T)
            ax_trans = ax
        #end if
        X.shape = self.grid_shape
        Y.shape = self.grid_shape
        Zm = self.f.T
        for Z in Zm:
            Z.shape = self.grid_shape
            fig,ax = self.setup_mpl_fig(fig=fig,dim=self.grid_dim)
            ax.contour(X,Y,Z,**kwargs)
            if boundary:
                if self.space_dim==2:
                    self.grid.plot_boundary(fig=False,show=False)
                else:
                    bpoints = self.grid.get_boundary_lines()
                    bpoints = np.inner(ax_trans,bpoints)
                    bpoints = np.transpose(bpoints,(1,2,0))
                    for bp in bpoints:
                        ax.plot(*bp.T,color='k')
                    #end for
                #end if
            #end if
        #end for
        if show:
            plt.show()
        #end if
    #end def plot_contours
    

    def plot_surface(self,fig=True,show=True,**kwargs):
        """
        (`External API`) Make 2D surface plots in the full coordinate space.

        Parameters
        ----------
        fig : `bool, optional, default True`
            Create a fresh figure (`True`) or not (`False`).
        show : `bool, optional, default True`
            Display the figure (`True`) or not (`False`).
        **kwargs : `optional`
            Arbitrary keyword argments passed to `pyplot.plot_surface`.
        """
        if self.grid_dim!=2:
            self.error('cannot plot surface\ngrid must have dimension 2 to make contour plots\ndimension of grid for this function: {}'.format(self.grid_dim))
        #end if
        X,Y = self.r.T
        X.shape = self.grid_shape
        Y.shape = self.grid_shape
        Zm = self.f.T
        for Z in Zm:
            Z.shape = self.grid_shape
            fig,ax = self.setup_mpl_fig(fig=fig,dim=self.grid_dim+1)
            ax.plot_surface(X,Y,Z,**kwargs)
        #end for
        if show:
            plt.show()
        #end if
    #end def plot_surface
    

    def plot_isosurface(self,level=None,fig=True,show=True,**kwargs):
        """
        (`External API`) Make 3D isosurface plots in the full coordinate space.

        Parameters
        ----------
        level : `float, optional`
            Isosurface value to plot.  If not provided, the average of the max 
            and min function values are used.
        fig : `bool, optional, default True`
            Create a fresh figure (`True`) or not (`False`).
        show : `bool, optional, default True`
            Display the figure (`True`) or not (`False`).
        **kwargs : `optional`
            Arbitrary keyword argments passed to `pyplot.plot_trisurf`.
        """
        if self.grid_dim!=3:
            self.error('cannot plot isosurface \ngrid must have dimension 3 to make contour plots\ndimension of grid for this function: {}'.format(self.grid_dim))
        #end if
        if self.space_dim!=3:
            self.error('cannot plot isosurface\ngrid points must reside in a 3D space to make contour plots\ndimension of the space for this function: {}'.format(self.space_dim))
        #end if
        fm = self.f.T
        # unit grid is used implicitly
        spacing = tuple(1./(np.array(self.cell_grid_shape,dtype=float)))
        for f in fm:
            if level is None:
                level = (f.max()+f.min())/2
            #end if
            f.shape = self.grid_shape
            ret = measure.marching_cubes(f,level,spacing=spacing)
            verts = ret[0] 
            faces = ret[1]
            verts = self.grid.points_from_unit(verts)
            fig,ax = self.setup_mpl_fig(fig=fig,dim=self.grid_dim)
            ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],**kwargs)
        #end for
        if show:
            plt.show()
        #end if
    #end def plot_isosurface


#end class StructuredGridFunctionWithAxes



class ParallelotopeGridFunction(StructuredGridFunctionWithAxes):
    """
    Represents functions over parallelotope grids.

    Most of the functionality is enabled by parent classes.  Instances of 
    this class must only own a ParallelotopeGrid.

    This class is intended for direct instantiation and use.

    Parameters
    ----------
    grid : `ParallelotopeGrid, optional`
        Grid of points in a `d` dimensional space.  If `grid` is not provided, 
        additional parameters must be given to initialize a `ParallelotopeGrid`
        object.
    values : `array_like, float/complex, shape (N,P), (N,)`
        Array of function values defined at the grid points.  `N` is the number
        of points and `P` is the number of function values.  With `P>1`, the 
        function is vector or tensor valued.  If the array has shape (`N`,), 
        then `P` is set to `1`.
    copy_grid : `bool, optional, default True`
        Copy provided grid (`True`) or not (`False`).
    copy_values : `bool, optional, default True`
        Copy provided values (`True`) or not (`False`).
    copy : `bool, optional`
        Copy provided grid and values (`True`) or not (`False`).
    dtype : `optional`
        Data type for local function values.
    grid_dtype : `optional`
        Data type for grid point locations.
    **kwargs: 
        Arbitrary set of parameters used to create a `ParallelotopeGrid` 
        object.  See documentation for the `ParallelotopeGrid` class and for 
        allowed inputs.  Used/allowed only if `grid` is not provided.

    Attributes
    ----------
    grid : `ParallelotopeGrid`
        Grid of points in a `d` dimensional space.
    values : `ndarray, float/complex, shape (N,P)`
        Array of function values defined at the grid points.  `N` is the number
        of points and `P` is the number of function values.  With `P>1`, the 
        function is vector or tensor valued.
    space_dim : `int, property`
        Dimension of the space the grid resides in.  Referred to as `d` above.
    grid_dim : `int, property`
        The dimension of the grid.
    npoints : `int, property`
        Number of grid points.  Referred to as `N` above.
    nvalues : `int, property`
        Number of function values at each grid point. Referred to as `P` above.
    grid_shape : `tuple, int, property`
        The number of grid points in each dimension.
    cell_grid_shape : `tuple, int, property`
        The number of grid cells in each dimension.
    ncells : `int, property`
        The total number of grid cells.
    flat_points_shape : `tuple, int, property`
        The shape of the `points` array in its default (flat) representation. 
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `flat_points_shape` is `(N*M*P,D)`.
    full_points_shape : `tuple, int, property`
        The shape of the points array in its full representation.
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `full_points_shape` is `(N,M,P,D)`.
    r : `ndarray, float, property`
        Array containing the grid points.
    f : `ndarray, float/complex, property`
        Array containing the function values.  User-facing alias for `values`.
    dtype : `property`
        Datatype of the function values.
    """
    grid_class = ParallelotopeGrid


    def read_local(self,filepath,format):
        if format=='xsf':
            self.read_xsf(filepath)
        else:
            self.error('Cannot read file.\nUnrecognized file format encountered.\nUnrecognized file format: {}\nValid options are: xsf'.format(format))
        #end if
    #end def read_local


    def read_xsf(self,filepath):
        if isinstance(filepath,XsfFile):
            xsf  = filepath
            copy = True
        else:
            xsf  = XsfFile(filepath)
            copy = False,
        #end if

        grid = self.grid_class()
        grid.read_xsf(xsf)

        xsf.remove_ghost()
        d = xsf.get_density()

        values = d.values_noghost.ravel()
        if copy:
            values = values.copy()
        #end if

        self.initialize(
            grid   = grid,
            values = values,
            )
    #end def read_xsf


    # test needed
    def read_from_points(self,points,values,axes,tol=1e-6,average=False):
        self.vlog('Reading grid function values from scattered data.')

        # check data types and shapes
        d = self.ensure_array(
            points = points,
            values = values,
            axes   = axes,
            dtype  = float,
            )
        points = d.points
        values = d.values
        axes   = d.axes.copy()
        del d
        if len(points)!=len(values):
            self.error('"points" and "values" must have the same length.\nNumber of points: {}\nNumber of values: {}'.format(len(points),len(values)))
        elif len(points.shape)!=2:
            self.error('Shape of "points" array must be (# points) x (# dimensions).\nShape provided: {}'.format(points.shape))
        #end if
        N,D = points.shape
        if axes.shape!=(D,D):
            self.error('"axes" must have shape {}\nShape provided: {} '.format((D,D),axes.shape))
        #end if

        # reshape values (for now GridFunction does not support more structured values)
        values.shape = len(values),values.size//len(values)

        # normalize the axes
        for d in range(D):
            axes[d] /= np.linalg.norm(axes[d])
        #end for

        # make the points rectilinear
        self.vlog('Transforming points to unit coords',n=1,time=True)
        rpoints = np.dot(points,np.linalg.inv(axes))

        # search for layers in each dimension
        def xlayers(xpoints,tol):
            xmin = xpoints.min()
            xmax = xpoints.max()
            nbins = np.uint64(np.round(np.ceil((xmax-xmin+tol)/tol)))
            dx = (xmax-xmin+tol)/nbins
            layers = obj()
            for x in xpoints:
                n = np.uint64(x/dx)
                if n not in layers:
                    layers[n] = obj(xsum=x,nsum=1)
                else:
                    l = layers[n]
                    l.xsum += x
                    l.nsum += 1
                #end if
            #end for
            for l in layers:
                l.xmean = l.xsum/l.nsum
            #end for
            lprev = None
            for n in sorted(layers.keys()):
                l = layers[n]
                if lprev is not None and np.abs(l.xmean-lprev.xmean)<tol:
                    lprev.xsum += l.xsum
                    lprev.nsum += l.nsum
                    lprev.xmean = lprev.xsum/lprev.nsum
                    del layers[n]
                else:
                    lprev = l
                #end if
            #end for
            xlayers = np.empty((len(layers),),dtype=float)
            i = 0
            for n in sorted(layers.keys()):
                l = layers[n]
                xlayers[i] = l.xmean
                i += 1
            #end for
            order  = xlayers.argsort()
            xlayers = xlayers[order]
            return xlayers,xmin,xmax
        #end def xlayers
        def index_by_layer(xpoints,tol):
            xlayer,xmin,xmax = xlayers(xpoints,tol)
            dxlayer = xlayer[1:]-xlayer[:-1]
            dxmin   = dxlayer.min()
            dxmax   = dxlayer.max()
            if np.abs(dxmax-dxmin)>2*tol:
                error('Could not determine layer separation.\nLayers are not evenly spaced.\nMin layer spacing: {}\nMax layer spacing: {}\nSpread   : {}\nTolerance: {}'.format(dxmin,dxmax,dxmax-dxmin,2*tol),'read_from_points')
            #end if
            dx = dxlayer.mean()
            ipoints = np.array(np.around((xpoints-xmin)/dx),dtype=int)
            return ipoints,xmin,xmax
        #end def index_by_layer

        # create a grid consistent with the detected layer separations
        self.vlog('Initializing point index array',n=1,time=True)
        grid_shape  = np.empty((D, ),dtype=int  )
        grid_axes   = np.zeros((D,D),dtype=float)
        grid_corner = np.empty((D, ),dtype=float)
        ipoints     = np.empty((N,D),dtype=int)
        for d in range(D):
            self.vlog('Indexing points along dim {}'.format(d),n=2,time=True)
            ixpoints,xmin,xmax = index_by_layer(rpoints[:,d],tol)
            grid_shape[d]  = ixpoints.max()+1
            grid_axes[d,d] = xmax-xmin
            grid_corner[d] = xmin
            ipoints[:,d]   = ixpoints
        #end for
        grid_axes   = np.dot(grid_axes,axes)
        grid_corner = np.dot(grid_corner,axes)
        grid_bconds = D*('o',) # assumed for now

        self.vlog('Constructing regular bounding grid',n=1,time=True)
        grid = self.grid_class(
            shape    = grid_shape,
            axes     = grid_axes,
            corner   = grid_corner,
            bconds   = grid_bconds,
            centered = False,
            )
        
        self.vlog('Checking grid point mapping',n=1,time=True)
        # check that the generated grid contains the inputted points
        ipflat = grid.flat_indices(ipoints)
        dmax = np.linalg.norm(points-grid.points[ipflat],axis=1).max()
        if dmax>tol:
            self.error('Generated grid points do not match those read in.\nMaximum deviation: {}\nTolerance        : {}'.format(dmax,tol))
        #end if
        # count number of times each grid point is mapped to
        point_counts = np.bincount(ipflat,minlength=grid.npoints)
        # if not averaging, check for one-to-one mapping
        max_count = point_counts.max()
        if not average and max_count>1:
            self.error('Mapping to grid points is not one-to-one.\nMax no. of read points mapped to a grid point: {}'.format(max_count))
        #end if

        # map the inputted values onto the generated grid
        self.vlog('Mapping data values onto grid',n=1,time=True)
        grid_values = np.zeros((grid.npoints,values.shape[1]),dtype=float)
        if not average or max_count==1:
            grid_values[ipflat] = values
        else:
            self.vlog('Averaging multi-valued points',n=2,time=True)
            for i,v in zip(ipflat,values):
                grid_values[i] += v
            #end for
            for i,c in enumerate(point_counts):
                if c>1:
                    grid_values[i] /= c
                #end if
            #end for
        #end if

        # initialize the GridFunction object
        self.vlog('Constructing GridFunction object',n=1,time=True)
        self.reset()
        self.initialize(
            grid   = grid,
            values = grid_values,
            copy   = False,
            )

        self.vlog('Read complete',n=1,time=True)
    #end def read_from_points

#end class ParallelotopeGridFunction



class SpheroidGridFunction(StructuredGridFunctionWithAxes):
    """
    Represents functions over spheroidal grids.

    Most of the functionality is enabled by parent classes.  Instances of 
    this class must only own a SpheroidGrid.

    This class is intended for direct instantiation and use.

    Parameters
    ----------
    grid : `SpheroidGrid, optional`
        Grid of points in a `d` dimensional space.  If `grid` is not provided, 
        additional parameters must be given to initialize a `SpheroidGrid` 
        object.
    values : `array_like, float/complex, shape (N,P), (N,)`
        Array of function values defined at the grid points.  `N` is the number
        of points and `P` is the number of function values.  With `P>1`, the 
        function is vector or tensor valued.  If the array has shape (`N`,), 
        then `P` is set to `1`.
    copy_grid : `bool, optional, default True`
        Copy provided grid (`True`) or not (`False`).
    copy_values : `bool, optional, default True`
        Copy provided values (`True`) or not (`False`).
    copy : `bool, optional`
        Copy provided grid and values (`True`) or not (`False`).
    dtype : `optional`
        Data type for local function values.
    grid_dtype : `optional`
        Data type for grid point locations.
    **kwargs: 
        Arbitrary set of parameters used to create a `SpheroidGrid` object. See
        documentation for the `SpheroidGrid` class for allowed inputs.  
        Used/allowed only if `grid` is not provided.

    Attributes
    ----------
    grid : `SpheroidGrid`
        Grid of points in a `d` dimensional space.
    values : `ndarray, float/complex, shape (N,P)`
        Array of function values defined at the grid points.  `N` is the number
        of points and `P` is the number of function values.  With `P>1`, the 
        function is vector or tensor valued.
    space_dim : `int, property`
        Dimension of the space the grid resides in.  Referred to as `d` above.
    grid_dim : `int, property`
        The dimension of the grid.
    npoints : `int, property`
        Number of grid points.  Referred to as `N` above.
    nvalues : `int, property`
        Number of function values at each grid point. Referred to as `P` above.
    grid_shape : `tuple, int, property`
        The number of grid points in each dimension.
    cell_grid_shape : `tuple, int, property`
        The number of grid cells in each dimension.
    ncells : `int, property`
        The total number of grid cells.
    flat_points_shape : `tuple, int, property`
        The shape of the `points` array in its default (flat) representation. 
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `flat_points_shape` is `(N*M*P,D)`.
    full_points_shape : `tuple, int, property`
        The shape of the points array in its full representation.
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `full_points_shape` is `(N,M,P,D)`.
    r : `ndarray, float, property`
        Array containing the grid points.
    f : `ndarray, float/complex, property`
        Array containing the function values.  User-facing alias for `values`.
    dtype : `property`
        Datatype of the function values.
    """
    grid_class = SpheroidGrid
#end class SpheroidGridFunction



class SpheroidSurfaceGridFunction(StructuredGridFunctionWithAxes):
    """
    Represents functions over spheroidal surface grids.

    Most of the functionality is enabled by parent classes.  Instances of 
    this class must only own a SpheroidalSurfaceGrid.

    This class is intended for direct instantiation and use.

    Parameters
    ----------
    grid : `SpheroidSurfaceGrid, optional`
        Grid of points in a `d` dimensional space.  If `grid` is not provided, 
        additional parameters must be given to initialize a 
        `SpheroidSurfaceGrid` object.
    values : `array_like, float/complex, shape (N,P), (N,)`
        Array of function values defined at the grid points.  `N` is the number
        of points and `P` is the number of function values.  With `P>1`, the 
        function is vector or tensor valued.  If the array has shape (`N`,), 
        then `P` is set to `1`.
    copy_grid : `bool, optional, default True`
        Copy provided grid (`True`) or not (`False`).
    copy_values : `bool, optional, default True`
        Copy provided values (`True`) or not (`False`).
    copy : `bool, optional`
        Copy provided grid and values (`True`) or not (`False`).
    dtype : `optional`
        Data type for local function values.
    grid_dtype : `optional`
        Data type for grid point locations.
    **kwargs: 
        Arbitrary set of parameters used to create a `SpheroidSurfaceGrid` 
        object.  See documentation for the `SpheroidSurfaceGrid` class for 
        allowed inputs.  Used/allowed only if `grid` is not provided.

    Attributes
    ----------
    grid : `SpheroidSurfaceGrid`
        Grid of points in a `d` dimensional space.
    values : `ndarray, float/complex, shape (N,P)`
        Array of function values defined at the grid points.  `N` is the number
        of points and `P` is the number of function values.  With `P>1`, the 
        function is vector or tensor valued.
    space_dim : `int, property`
        Dimension of the space the grid resides in.  Referred to as `d` above.
    grid_dim : `int, property`
        The dimension of the grid.
    npoints : `int, property`
        Number of grid points.  Referred to as `N` above.
    nvalues : `int, property`
        Number of function values at each grid point. Referred to as `P` above.
    grid_shape : `tuple, int, property`
        The number of grid points in each dimension.
    cell_grid_shape : `tuple, int, property`
        The number of grid cells in each dimension.
    ncells : `int, property`
        The total number of grid cells.
    flat_points_shape : `tuple, int, property`
        The shape of the `points` array in its default (flat) representation. 
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `flat_points_shape` is `(N*M*P,D)`.
    full_points_shape : `tuple, int, property`
        The shape of the points array in its full representation.
        If `grid_shape` is `(N,M,P)` and `space_dim` is `D`, then 
        `full_points_shape` is `(N,M,P,D)`.
    r : `ndarray, float, property`
        Array containing the grid points.
    f : `ndarray, float/complex, property`
        Array containing the function values.  User-facing alias for `values`.
    dtype : `property`
        Datatype of the function values.
    """
    grid_class = SpheroidSurfaceGrid
#end class SpheroidSurfaceGridFunction



# test needed
def parallelotope_grid_function(
    loc = 'parallelotope_grid_function',
    **kwargs
    ):
    if 'points' not in kwargs:
        gf = ParallelotopeGridFunction(**kwargs)
    else:
        required = set(('points','values','axes'))
        optional = set(('tol','average'))
        present  = set(kwargs.keys())
        if len(required-present)>0:
            error('Grid function cannot be created.\nWhen "points" is provided, "axes" and "values" must also be given.\nInputs provided: {}'.format(sorted(present)),loc)
        elif len(present-required-optional)>0:
            error('Grid function cannot be created.\nUnrecognized inputs provided.\nUnrecognized inputs: {}\nValid options are: {}'.format(sorted(present-required-optional),sorted(required|optional)))
        #end if
        gf = ParallelotopeGridFunction()
        gf.read_from_points(**kwargs)
    #end if
    return gf
#end def parallelotope_grid_function



# test needed
def grid_function(
    type = 'parallelotope',
    loc  = 'grid_function',
    **kwargs
    ):
    filepath = kwargs.pop('filepath',None)
    if filepath is not None:
        return read_grid_function(filepath,loc=loc)
    #end if
    gf = None
    if type=='parallelotope':
        gf = parallelotope_grid_function(loc=loc,**kwargs)
    elif type=='spheroid':
        gf = SpheroidGridFunction(**kwargs)
    elif type=='spheroid_surface':
        gf = SpheroidSurfaceGridFunction(**kwargs)
    else:
        error('Grid function type "{}" is not recognized.\nValid options are: parallelotope, spheroid, or spheroid_surface'.format(type),loc)
    #end if
    return gf
#end def grid_function


def grid(
    type = 'parallelotope',
    loc  = 'grid',
    **kwargs
    ):
    filepath = kwargs.pop('filepath',None)
    if filepath is not None:
        return read_grid(filepath,loc=loc)
    #end if
    g = None
    if type=='parallelotope':
        g = ParallelotopeGrid(**kwargs)
    elif type=='spheroid':
        g = SpheroidGrid(**kwargs)
    elif type=='spheroid_surface':
        g = SpheroidSurfaceGrid(**kwargs)
    else:
        error('Grid type "{}" is not recognized.\nValid options are: parallelotope, spheroid, or spheroid_surface'.format(type),loc)
    #end if
    return g
#end def grid



#test needed
gf_file_type_map = obj(
    xsf = ParallelotopeGridFunction,
    )
def read_grid_function(filepath,format=None,loc='read_grid_function'):
    filepath,format = process_file_format(filepath,format,loc)
    if format not in gf_file_type_map:
        error('Cannot read file.\nUnrecognized file format for grid function.\nFile format provided: {}\nAllowed formats include: {}'.format(format,sorted(gf_file_type_map.keys())),loc)
    #end if
    gf = gf_file_type_map[format]()
    gf.read(filepath)
    return gf
#end def read_grid_function


# test needed
g_file_type_map = obj(
    xsf = ParallelotopeGrid,
    )
def read_grid(filepath,format=None,loc='read_grid'):
    filepath,format = process_file_format(filepath,format,loc)
    if format not in g_file_type_map:
        error('Cannot read file.\nUnrecognized file format for grid.\nFile format provided: {}\nAllowed formats include: {}'.format(format,sorted(g_file_type_map.keys())),loc)
    #end if
    g = g_file_type_map[format]()
    g.read(filepath)
    return g
#end def read_grid


def process_file_format(filepath,format,loc):
    if isinstance(filepath,StandardFile):
        format = filepath.sftype
    elif not isinstance(filepath,str):
        error('Cannot read file.\nExpected a file path.\nInstead received type: {}\nWith value: {}\nPlease provide a file path and try again.'.format(filepath.__class__.__name__,filepath),loc)
    elif not os.path.exists(filepath):
        error('Cannot read file.  File path does not exist.\nFile path provided: {}'.format(filepath),loc)
    #end if
    if format is None:
        format = filepath.rsplit('.',1)[1].lower()
    #end if
    return filepath,format
#end def process_file_format



gfs = [
    ParallelotopeGridFunction,
    SpheroidGridFunction,
    SpheroidSurfaceGridFunction,
    ]
grid_to_grid_function = obj()
for gf in gfs:
    grid_to_grid_function[gf.grid_class.__name__] = gf
#end for
del gfs

def grid_function_from_grid(grid):
    gname = grid.__class__.__name__
    if gname not in grid_to_grid_function:
        error('Cannot find matching grid function for grid "{}".'.format(gname))
    #end if
    return grid_to_grid_function[gname]
#end def grid_function_from_grid




if __name__=='__main__':


    demos = obj(
        plot_grids         = 0,
        plot_inside        = 0,
        plot_projection    = 0,
        cell_volumes       = 0,
        plot_contours      = 1,
        plot_surface       = 0,
        plot_isosurface    = 0,
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
        #(2,2) : [[1,0],[1,1]],
        (2,2) : [[1,0],[0,1]],
        # 2d grid in 3d space
        #(2,3) : [[1,0,0],[1,1,1]],
        (2,3) : [[1,0,0],[0,1,0]],
        # 3d grid in 3d space
        #(3,3) : [[1,0,0],[1,1,0],[1,1,1]],
        (3,3) : [[1,0,0],[0,1,0],[0,0,1]],
        }

    bconds = {
        1 : 'o p'.split(),
        2 : 'oo op po pp'.split(),
        3 : 'ooo oop opo poo opp pop ppo ppp'.split(),
        }

    grid_types = obj(
        parallelotope    = ParallelotopeGrid,
        spheroid         = SpheroidGrid,
        spheroid_surface = SpheroidSurfaceGrid,
        )

    supported = obj(
        parallelotope    = obj(dims=set(axes.keys())),
        spheroid         = obj(dims=set([(2,2),(2,3),(3,3)])),
        spheroid_surface = obj(dims=set([(1,2),(1,3),(2,3)])),
        )

    gdict = dict(
        parallelotope    = 'p',
        spheroid         = 's',
        spheroid_surface = 'c',
        )

    grids = obj()
    cdict = {False:'',True:'c'}
    for grid_name in sorted(grid_types.keys()):
        label = gdict[grid_name]+'e'
        grids[label] = grid_types[grid_name]()
        for grid_dim,space_dim in sorted(axes.keys()):
            if (grid_dim,space_dim) in supported[grid_name].dims:
                for centered in (False,True):
                    label = gdict[grid_name]+str(grid_dim)+str(space_dim)+cdict[centered]
                    cell_grid_dim = grid_dim
                    axes_grid_dim = grid_dim
                    if grid_name=='spheroid_surface':
                        axes_grid_dim += 1
                    #end if
                    grid_inputs = obj(
                        #shape    = shapes[grid_dim],
                        cells    = shapes[cell_grid_dim],
                        axes     = axes[axes_grid_dim,space_dim],
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
        if not g.initialized:
            continue
        #end if
        print(' {:<16}  {} {}  {}  {}  {}'.format(label,g.grid_dim,g.space_dim,len(g.axes),g.bconds,g.shape))
    #end for

    if demos.plot_grids:
        #grids_plot = 'p23c p23_oo p23_op p23_pp s23c s23'.split()
        #grids_plot = 'p11 p12 p13 p23 p23c p33 p33c'.split()
        #grids_plot = 's23 s23c s33 s33c'.split()
        #grids_plot = 'c12 c12c c13 c13c c23 c23c'.split()
        grids_plot = 'p33 s33'.split()

        unit = False

        for name in grids_plot:
            print(name)
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
        #grids_plot = 's22 s23'.split()
        grids_plot = 'c23'.split()

        n = 100

        unit = False

        upoints = dict()
        for d in range(1,3+1):
            upoints[d] = 0.75*np.random.randn(n,d)+0.75
        #end for

        for name in grids_plot:
            g = grids[name]
            #dim = int(name[1])
            dim = g.grid_dim
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


    if demos.cell_volumes:
        #grids_check = 'p22 p33'.split()
        #grids_check = 's22 s23 s33'.split()
        #grids_check = 'p11 p12 p13 p22 p23 p33 s22 s23 s33'.split()
        grids_check = 'c12 c13 c23'.split()
        for name in grids_check:
            g = grids[name]
            print(name,g.volume(),g.cell_volumes().sum())
        #end for
    #end if


    if demos.plot_contours:

        gp = ParallelotopeGrid(
            axes  = [[1,0],
                     [1,1]],
            bconds = 'pp',
            cells  = (80,80),
            )

        u = gp.unit_points()
        values = np.cos(3*np.pi*u[:,0])**2*np.sin(2*np.pi*u[:,1])

        fp = ParallelotopeGridFunction(
            grid   = gp,
            values = values,
            )

        fp.plot_unit_contours(boundary=True,show=False)

        fp.plot_contours(boundary=True,show=False)



        gs = SpheroidGrid(
            axes  = [[1,0],
                     [1,1]],
            cells  = (80,80),
            )

        u = gs.unit_points()
        values = np.cos(3*np.pi*u[:,0])**2*np.sin(2*np.pi*u[:,1])

        fs = SpheroidGridFunction(
            grid   = gs,
            values = values,
            )

        fs.plot_unit_contours(boundary=True,show=False)

        fs.plot_contours(boundary=True)
        
    #end if


    if demos.plot_surface:

        gp = ParallelotopeGrid( 
           axes  = [[1,0],
                     [1,1]],
            bconds = 'pp',
            cells  = (80,80),
            )

        u = gp.unit_points()
        values = np.cos(3*np.pi*u[:,0])**2*np.sin(2*np.pi*u[:,1])

        fp = ParallelotopeGridFunction(
            grid   = gp,
            values = values,
            )

        fp.plot_unit_surface(show=0)

        fp.plot_surface(show=0)



        gs = SpheroidGrid(
            axes  = [[1,0],
                     [1,1]],
            cells  = (80,80),
            )
        
        u = gs.unit_points()
        values = np.cos(3*np.pi*u[:,0])**2*np.sin(2*np.pi*u[:,1])
        
        fs = SpheroidGridFunction(
            grid   = gs,
            values = values,
            )
        
        fs.plot_unit_surface(show=False)
        
        fs.plot_surface()
        
    #end if


    if demos.plot_isosurface:

        gp = ParallelotopeGrid(
            axes  = [[1,0,0],
                     [1,1,0],
                     [1,1,1]],
            bconds = 'ppp',
            cells  = (30,30,30),
            )

        u = gp.unit_points()
        r = 2*np.pi*(u-1)
        values = np.cos(r[:,0])+ np.cos(r[:,1])+ np.cos(r[:,2])

        fp = ParallelotopeGridFunction(
            grid   = gp,
            values = values,
            )

        fp.plot_unit_isosurface(level=0,show=0)

        fp.plot_isosurface(level=0,show=0)



        gs = SpheroidGrid(
            axes  = [[1,0,0],
                     [1,1,0],
                     [1,1,1]],
            cells  = (30,30,30),
            )

        u = gs.unit_points()
        r = 2*np.pi*(u-1)
        values = np.cos(r[:,0])+ np.cos(r[:,1])+ np.cos(r[:,2])


        fs = SpheroidGridFunction(
            grid   = gs,
            values = values,
            )
        
        fs.plot_unit_isosurface(level=0,show=False)
        
        fs.plot_isosurface(level=0)
        
    #end if
    

#end if 
