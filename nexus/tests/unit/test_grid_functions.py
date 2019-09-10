

def test_imports():
    import numpy
    import testing
    import grid_functions
#end def test_imports


def test_coord_conversion():

    import numpy as np
    from grid_functions import polar_to_cartesian,cartesian_to_polar
    from grid_functions import spherical_to_cartesian,cartesian_to_spherical

    pi = np.pi
    sqrt = np.sqrt

    th_cos_sin = [
             0 , sqrt(4.)/2 , sqrt(0.)/2 ,
          pi/6 , sqrt(3.)/2 , sqrt(1.)/2 ,
          pi/4 , sqrt(2.)/2 , sqrt(2.)/2 ,
          pi/3 , sqrt(1.)/2 , sqrt(3.)/2 ,
          pi/2 , sqrt(0.)/2 , sqrt(4.)/2 ,
        2*pi/3 ,-sqrt(1.)/2 , sqrt(3.)/2 ,
        3*pi/4 ,-sqrt(2.)/2 , sqrt(2.)/2 ,
        5*pi/6 ,-sqrt(3.)/2 , sqrt(1.)/2 ,
            pi ,-sqrt(4.)/2 , sqrt(0.)/2 ,
        7*pi/6 ,-sqrt(3.)/2 ,-sqrt(1.)/2 ,
        5*pi/4 ,-sqrt(2.)/2 ,-sqrt(2.)/2 ,
        4*pi/3 ,-sqrt(1.)/2 ,-sqrt(3.)/2 ,
        3*pi/2 , sqrt(0.)/2 ,-sqrt(4.)/2 ,
        5*pi/3 , sqrt(1.)/2 ,-sqrt(3.)/2 ,
        7*pi/4 , sqrt(2.)/2 ,-sqrt(2.)/2 ,
       11*pi/6 , sqrt(3.)/2 ,-sqrt(1.)/2 ,
             0 , sqrt(4.)/2 , sqrt(0.)/2 ,
        ]

    th_cos_sin = np.array(th_cos_sin)
    th_cos_sin.shape = len(th_cos_sin)/3,3
    th  = th_cos_sin[:,0]
    cos = th_cos_sin[:,1]
    sin = th_cos_sin[:,2]

    x = cos
    y = sin

    polar_ref = np.empty((len(th),2),dtype=float)
    polar_ref[:,0] = 1.0
    polar_ref[:,1] = th

    cart_ref = np.empty((len(th),2),dtype=float)
    cart_ref[:,0] = x
    cart_ref[:,1] = y

    cart_from_polar = polar_to_cartesian(polar_ref)
    diff = np.abs(cart_ref-cart_from_polar).max()
    p2c_works = bool(diff<1e-12)
    assert(p2c_works)

    polar_from_cart = cartesian_to_polar(cart_ref)
    diff =  np.abs(polar_ref-polar_from_cart).max()
    c2p_works = bool(diff<1e-12)
    assert(c2p_works)

    phi = th
    th = th[1:len(th)/2+1]

    sphere_ref = np.empty((len(th)*len(phi),3),dtype=float)
    cart_ref   = np.empty((len(th)*len(phi),3),dtype=float)
    sphere_ref[:,0] = 1.0
    n=0
    for t in th:
        for p in phi:
            sphere_ref[n,1] = t
            sphere_ref[n,2] = p
            cart_ref[n,0] = np.cos(p)*np.sin(t)
            cart_ref[n,1] = np.sin(p)*np.sin(t)
            cart_ref[n,2] = np.cos(t)
            n+=1
        #end for
    #end for
    
    cart_from_sphere = spherical_to_cartesian(sphere_ref)
    diff = np.abs(cart_ref-cart_from_sphere).max()
    assert(diff<1e-12)

    sphere_from_cart = cartesian_to_spherical(cart_ref)
    diff = np.abs(sphere_ref-sphere_from_cart).max()
    assert(diff<1e-12)

#end def test_coord_conversion



def test_unit_grid_points():
    import numpy as np

    from testing import value_eq

    from grid_functions import unit_grid_points

    lin_grid = np.array([0.00,0.25,0.50,0.75])
    lin_grid_centered = lin_grid+0.125
    lin_grid_endpoint = np.array([0.00,0.25,0.50,0.75,1.00])

    lin_gridh = np.array([0.00,0.50])
    lin_gridh_centered = lin_gridh+0.25
    lin_gridh_endpoint = np.array([0.00,0.50,1.00])

    def make_1d(x):
        p = x.copy()
        p.shape = len(p),1
        return p
    #end def make_1d

    def make_2d(x,y):
        p = []
        for xv in x:
            for yv in y:
                p.append((xv,yv))
            #end for
        #ed for
        p = np.array(p)
        return p
    #end def make_2d

    def make_3d(x,y,z):
        p = []
        for xv in x:
            for yv in y:
                for zv in z:
                    p.append((xv,yv,zv))
                #end for
            #end for
        #ed for
        p = np.array(p)
        return p
    #end def make_3d

    # test 1d grids
    ref = make_1d(lin_grid)
    u = unit_grid_points((4,))
    assert(value_eq(u,ref))

    ref = make_1d(lin_grid_centered)
    u = unit_grid_points((4,),centered=True)
    assert(value_eq(u,ref))

    ref = make_1d(lin_grid_endpoint)
    u = unit_grid_points((5,),endpoint=[True])
    assert(value_eq(u,ref))

    # test 2d grids
    ref = make_2d(lin_grid,lin_gridh)
    u = unit_grid_points((4,2))
    assert(value_eq(u,ref))

    ref = make_2d(lin_grid_centered,lin_gridh_centered)
    u = unit_grid_points((4,2),centered=True)
    assert(value_eq(u,ref))

    ref = make_2d(lin_grid_endpoint,lin_gridh_endpoint)
    u = unit_grid_points((5,3),endpoint=[True,True])
    assert(value_eq(u,ref))
    
    # test 3d grids
    ref = make_3d(lin_grid,lin_gridh,lin_grid)
    u = unit_grid_points((4,2,4))
    assert(value_eq(u,ref))

    ref = make_3d(lin_grid_centered,lin_gridh_centered,lin_grid_centered)
    u = unit_grid_points((4,2,4),centered=True)
    assert(value_eq(u,ref))

    ref = make_3d(lin_grid_endpoint,lin_gridh_endpoint,lin_grid_endpoint)
    u = unit_grid_points((5,3,5),endpoint=[True,True,True])
    assert(value_eq(u,ref))

#end def test_unit_grid_points



def test_parallelotope_grid_points():

    import numpy as np
    from testing import value_eq
    from grid_functions import unit_grid_points
    from grid_functions import parallelotope_grid_points

    axes_dim = {
        1 : [[1]],
        2 : [[1,1],
             [0,1]],
        3 : [[1,1,1],
             [0,1,1],
             [0,0,1]],
        }
    shapes_dim = {
        1 : (5,),
        2 : (5,3),
        3 : (5,3,7),
        }

    for d in range(1,4):
        axes = axes_dim[d]
        shape = shapes_dim[d]

        for c in (False,True):
            for e in (False,True):
                ep = d*[e]
                u = unit_grid_points(shape=shape,centered=c,endpoint=ep)
                ref = np.dot(u,axes)
                p = parallelotope_grid_points(axes=axes,
                                              shape=shape,
                                              centered=c,
                                              endpoint=ep)
                assert(value_eq(p,ref))

                cells = np.array(shape,dtype=int)
                if not c:
                    cells-=np.array(ep,dtype=int)
                #end if
                p = parallelotope_grid_points(axes=axes,
                                               cells=cells,
                                               centered=c,
                                               endpoint=ep)
                assert(value_eq(p,ref))

                dr = 1.0/cells
                for i in range(len(dr)):
                    dr[i] *= np.linalg.norm(axes[i])
                #end for
                p = parallelotope_grid_points(axes=axes,
                                               dr=dr,
                                               centered=c,
                                               endpoint=ep)
                assert(value_eq(p,ref))
                
            #end for
        #end for
    #end for
#end def test_parallelotope_grid_points



def test_spheroid_grid_points():
    import numpy as np
    from testing import value_eq
    from grid_functions import cartesian_to_polar
    from grid_functions import cartesian_to_spherical
    from grid_functions import spheroid_grid_points
    from grid_functions import spheroid_surface_grid_points

    n = 5

    def unique(vals,tol=1e-6):
        vals = np.unique(vals)
        if len(vals)>1:
            vc = vals[0]
            uvals = [vc]
            for v in vals[1:]:
                if np.abs(v-vc)>tol:
                    vc = v
                    uvals.append(vc)
                #end if
            #end for
            vals = np.array(uvals)
        #end if
        return vals
    #end def unique

    def check(r=None,theta=None,phi=None):
        if r is not None:
            r = unique(r)
            dr = unique(r[1:]-r[:-1])
            assert(len(r)==n)
            assert(len(dr)==1)
            assert(value_eq(dr[0],1.0/n))
        #end if
        if theta is not None:
            t = unique(theta)
            dt = unique(t[1:]-t[:-1])
            assert(len(t)==n)
            assert(len(dt)==1)
            assert(value_eq(dt[0],np.pi/n))
        #end if
        if phi is not None:
            p = unique(phi)
            dp = unique(p[1:]-p[:-1])
            assert(len(p)==n)
            assert(len(dp)==1)
            assert(value_eq(dp[0],2*np.pi/n))
        #end if
    #end def check

    # grid within unit disk
    p = spheroid_grid_points(axes=np.eye(2),cells=2*[n],centered=True)
    s = cartesian_to_polar(p)
    check(r=s[:,0],phi=s[:,1])

    # grid within unit ball
    p = spheroid_grid_points(axes=np.eye(3),cells=3*[n],centered=True)
    s = cartesian_to_spherical(p)
    check(r=s[:,0],theta=s[:,1],phi=s[:,2])

    # grid on unit circle
    p = spheroid_surface_grid_points(axes=np.eye(2),cells=[n],centered=True)
    s = cartesian_to_polar(p)
    check(phi=s[:,1])

    # grid on unit sphere
    p = spheroid_surface_grid_points(axes=np.eye(3),cells=2*[n],centered=True)
    s = cartesian_to_spherical(p)
    check(theta=s[:,1],phi=s[:,2])

#end def test_spheroid_grid_points




grid_objects = dict()


def get_grids():
    """
    Generate a variety of grids.
    """

    from generic import obj
    from grid_functions import ParallelotopeGrid
    from grid_functions import SpheroidGrid
    from grid_functions import SpheroidSurfaceGrid

    grids = grid_objects

    if len(grids)==0:

        cells = {
            1 : (5,),     # 5x    grid cells
            2 : (5,5),    # 5x5   grid cells
            3 : (5,5,5),  # 5x5x5 grid cells
            }

        axes = {
            # 1d grid in 1d space
            (1,1,False) : [[1]],      # unsheared
            (1,1,True ) : [[1]],      #   sheared
            # 1d grid in 2d space
            (1,2,False) : [[1,1]],    # unsheared
            (1,2,True ) : [[1,1]],    #   sheared
            # 1d grid in 3d space
            (1,3,False) : [[1,1,1]],  # unsheared
            (1,3,True ) : [[1,1,1]],  #   sheared
            # 2d grid in 2d space
            (2,2,False) : [[1,0],[0,1]], # unsheared
            (2,2,True ) : [[1,0],[1,1]], #   sheared
            # 2d grid in 3d space
            (2,3,False) : [[1,0,0],[0,1,0]], # unsheared
            (2,3,True ) : [[1,0,0],[1,1,1]], #   sheared
            # 3d grid in 3d space
            (3,3,False) : [[1,0,0],[0,1,0],[0,0,1]], # unsheared
            (3,3,True ) : [[1,0,0],[1,1,0],[1,1,1]], #   sheared
            }

        bconds = { # combinations of open/periodic bc's
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
            parallelotope    = obj(dims=set([(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)])),
            spheroid         = obj(dims=set([(2,2),(2,3),(3,3)])),
            spheroid_surface = obj(dims=set([(1,2),(1,3),(2,3)])),
            )

        gdict = dict(
            parallelotope    = 'p',
            spheroid         = 's',
            spheroid_surface = 'c',
            )

        # centering label dictionary
        cdict = {False:'',True:'c'}
        # shearing label dictionary
        sdict = {False:'',True:'s'}
        # translation label dictionary
        tdict = {False:'',True:'t'}
        for grid_name in sorted(grid_types.keys()):
            for grid_dim,space_dim,sheared in sorted(axes.keys()):
                if (grid_dim,space_dim) in supported[grid_name].dims:
                    for centered in (False,True):
                        for translated in (False,True):
                            attrib = sdict[sheared]+cdict[centered]+tdict[translated]
                            label = gdict[grid_name]+str(grid_dim)+str(space_dim)+attrib
                            cell_grid_dim = grid_dim
                            axes_grid_dim = grid_dim
                            if grid_name=='spheroid_surface':
                                axes_grid_dim += 1
                            #end if
                            grid_inputs = obj(
                                cells    = cells[cell_grid_dim],
                                axes     = axes[axes_grid_dim,space_dim,sheared],
                                centered = centered,
                                )
                            if translated:
                                grid_inputs.origin = space_dim*(1,)
                            #end if

                            g = grid_types[grid_name](**grid_inputs)
                            grids[label] = g

                            if grid_name=='parallelotope':
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

    #end if
    return obj(grids)
#end def get_grids


def properties_from_name(grid_name):
    from generic import obj
    grid_types = dict(
        p = 'parallelotope',
        s = 'spheroid',
        c = 'spheroid_surface',
        )
    grid_type = grid_types[grid_name[0]]
    grid_dim  = int(grid_name[1])
    space_dim = int(grid_name[2])
    extra_prop = grid_name[3:]
    bconds = None
    if '_' in extra_prop:
        extra_prop,bc = extra_prop.split('_')
        bconds = tuple(bc)
    #end if
    centered   = 'c' in extra_prop
    sheared    = 's' in extra_prop
    translated = 't' in extra_prop
    properties = obj(
        type       = grid_type,
        grid_dim   = grid_dim,
        space_dim  = space_dim,
        centered   = centered,
        sheared    = sheared,
        translated = translated,
        bconds     = bconds,
        cells      = grid_dim*(5,),
        )
    return properties
#end def properties_from_name



def test_grid_initialization():
    import numpy as np
    from testing import object_eq,value_eq
    from grid_functions import Grid,StructuredGrid,StructuredGridWithAxes
    from grid_functions import ParallelotopeGrid
    from grid_functions import SpheroidGrid
    from grid_functions import SpheroidSurfaceGrid

    grids = []

    grid_types = [
        Grid,
        StructuredGrid,
        StructuredGridWithAxes,
        ParallelotopeGrid,
        SpheroidGrid,
        SpheroidSurfaceGrid,
        ]
    for gt in grid_types:
        g = gt()
        assert(not g.valid())
    #end for

    grids = get_grids()

    # check validity
    for g in grids:
        assert(g.valid())
    #end for

    # check properties
    bcs = set(tuple('op'))
    for name in sorted(grids.keys()):
        g = grids[name]
        p = properties_from_name(name)
        assert(g.initialized)
        assert(g.cell_grid_shape==p.cells)
        cells = np.array(p.cells,dtype=int)
        shape = np.array(g.grid_shape,dtype=int)
        assert(np.abs(shape-cells).max()<=1)
        npoints = shape.prod()
        ncells  = cells.prod()
        assert(g.shape==g.grid_shape)
        assert(len(g.points)==npoints)
        assert(g.npoints==npoints)
        assert(g.ncells==ncells)
        assert(id(g.r)==id(g.points))
        assert(g.space_dim==p.space_dim)
        assert(g.grid_dim==p.grid_dim)
        assert(g.centered==p.centered)
        assert(g.surface==(p.type=='spheroid_surface'))
        assert(g.flat_points_shape==(npoints,p.space_dim))
        assert(g.full_points_shape==g.shape+(p.space_dim,))
        assert(len(set(g.bconds)-bcs)==0)
        if not p.translated:
            assert(value_eq(g.origin,np.zeros((p.space_dim,))))
        else:
            assert(value_eq(g.origin,np.ones((p.space_dim,))))
        #end if
        if p.type!='spheroid_surface':
            assert(g.axes.shape==(p.grid_dim,p.space_dim))
        else:
            assert(g.axes.shape==(p.grid_dim+1,p.space_dim))
        #end if
        if p.type=='parallelotope':
            assert(value_eq(g.corner,g.origin))
            center = g.corner + g.axes.sum(axis=0)/2
            assert(value_eq(g.center,center))
        elif 'spheroid' in p.type:
            assert(value_eq(g.center,g.origin))
            if not p.sheared:
                assert(g.isotropic)
            #end if
        #end if
    #end for

#end def test_grid_initialization



def test_parallelotope_grid_initialization():

    from testing import object_eq
    from grid_functions import ParallelotopeGrid

    pgrids = []

    p1 = ParallelotopeGrid(
        shape = (6,),
        axes  = [[1]],
        )
    pgrids.append(p1)

    p2 = ParallelotopeGrid(
        cells = (5,),
        axes  = [[1]],
        )
    pgrids.append(p2)

    p3 = ParallelotopeGrid(
        dr    = (0.2,),
        axes  = [[1]],
        )
    pgrids.append(p2)

    p4 = ParallelotopeGrid(
        shape = (5,),
        axes  = [[1]],
        centered = True,
        )
    pgrids.append(p4)

    p5 = ParallelotopeGrid(
        cells = (5,),
        axes  = [[1]],
        centered = True,
        )
    pgrids.append(p5)

    p6 = ParallelotopeGrid(
        dr    = (0.2,),
        axes  = [[1]],
        centered = True,
        )
    pgrids.append(p6)

    assert(object_eq(p1,p2))
    assert(object_eq(p1,p3))

    assert(object_eq(p4,p5))
    assert(object_eq(p4,p6))
#end def test_parallelotope_grid_initialization



def test_grid_reset():
    from testing import object_eq
    from generic import obj
    from grid_functions import ParallelotopeGrid
    from grid_functions import SpheroidGrid
    from grid_functions import SpheroidSurfaceGrid
    
    grid_inputs = [
        ( ParallelotopeGrid   , obj(cells=(5,6),axes=[[1,0,0],[1,1,0]]) ),
        ( SpheroidGrid        , obj(cells=(5,6),axes=[[1,0,0],[1,1,0]]) ),
        ( SpheroidSurfaceGrid , obj(cells=(5,6),axes=[[1,0,0],[1,1,0],[1,1,1]]) ),
        ]

    for grid_type,inputs in grid_inputs:
        empty_grid = grid_type()
        grid = grid_type(**inputs)
        grid_copy = grid.copy()
        grid.reset()
        assert(object_eq(grid,empty_grid))
        grid.initialize(**inputs)
        assert(object_eq(grid,grid_copy))
    #end for
#end def test_grid_reset



def test_grid_set_operations():
    import numpy as np
    from generic import obj
    from testing import value_eq

    grids = get_grids()

    grids_check = 'p22c s22c c23c'.split()

    props = obj()
    for name in grids_check:
        props[name] = properties_from_name(name)
    #end for

    # set points
    points = np.linspace(0,1,2*10)
    points.shape = 10,2
    for name in grids_check:
        g = grids[name].copy()
        g.set_points(points)
        assert(value_eq(g.points,points))
        assert(id(g.points)!=id(points))
        g.set_points(points,copy=False)
        assert(value_eq(g.points,points))
        assert(id(g.points)==id(points))
    #end for

    # set shape
    for name in grids_check:
        p = props[name]
        g = grids[name].copy()
        g.set_shape(p.cells)
        assert(g.shape==p.cells)
    #end for

    # set bconds
    for name in grids_check:
        p = props[name]
        g = grids[name].copy()
        bcs = tuple('oo')
        g.set_bconds(bcs)
        assert(tuple(g.bconds)==bcs)
        assert(id(g.bconds)!=id(bcs))
        g.set_bconds('oo')
        assert(tuple(g.bconds)==bcs)
    #end for

    # set axes
    for name in grids_check:
        p = props[name]
        g = grids[name].copy()
        dim = p.grid_dim
        if g.surface:
            dim += 1
        #end if
        axes = np.eye(dim,dtype=float)
        g.set_axes(axes)
        assert(value_eq(g.axes,axes))
        assert(id(g.axes)!=id(axes))
    #end for

    # set origin
    for name in grids_check:
        p = props[name]
        g = grids[name].copy()
        go = g.copy()
        origin = np.array(p.space_dim*(4.5,))
        g.set_origin(origin)
        assert(value_eq(g.origin,origin))
        assert(id(g.origin)!=id(origin))
        shift = origin-go.origin
        dpoints = g.points-go.points
        for dp in dpoints:
            assert(value_eq(dp,shift))
        #end for
    #end for

#end def test_grid_set_operations



def test_cell_indices_2d():

    import numpy as np
    from grid_functions import ParallelotopeGrid
    from grid_functions import SpheroidGrid
    from grid_functions import SpheroidSurfaceGrid

    p = ParallelotopeGrid(
        cells = (10,10),
        axes  = [[1,0,0],[1,1,1]],
        )
    pc = ParallelotopeGrid(
        cells = (10,10),
        axes  = [[1,0,0],[1,1,1]],
        centered = True
        )

    s = SpheroidGrid(
        cells = (10,10),
        axes  = [[1,0,0],[1,1,1]],
        )
    sc = SpheroidGrid(
        cells = (10,10),
        axes  = [[1,0,0],[1,1,1]],
        centered = True,
        )

    c = SpheroidSurfaceGrid(
        cells = (10,10),
        axes  = [[1,0,0],[1,1,0],[1,1,1]],
        )
    cc = SpheroidSurfaceGrid(
        cells = (10,10),
        axes  = [[1,0,0],[1,1,0],[1,1,1]],
        centered = True,
        )

    imap = np.arange(100,dtype=int)

    assert((pc.cell_indices()==imap).all())
    assert((p.cell_indices(pc.r)==imap).all())
    assert((sc.cell_indices()==imap).all())
    assert((s.cell_indices(sc.r)==imap).all())
    assert((cc.cell_indices()==imap).all())
    assert((c.cell_indices(cc.r)==imap).all())
#end def test_cell_indices_2d
