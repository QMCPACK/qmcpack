

def test_imports():
    from grid_functions import *
#end def test_imports


def test_coord_conversion():

    import numpy as np
    from grid_functions import polar_to_cartesian,cartesian_to_polar

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
    cart_ref[:,0] = cos
    cart_ref[:,1] = sin

    cart_from_polar = polar_to_cartesian(polar_ref)
    for i in range(len(cart_ref)):
        pr = cart_ref[i]
        pc = cart_from_polar[i]
    #end for
    diff = np.abs(cart_ref-cart_from_polar).max()
    p2c_works = bool(diff<1e-12)
    assert(p2c_works)

    polar_from_cart = cartesian_to_polar(cart_ref)
    for i in range(len(polar_ref)):
        pr = polar_ref[i]
        pc = polar_from_cart[i]
    #end for
    diff =  np.abs(polar_ref-polar_from_cart).max()
    c2p_works = bool(diff<1e-12)
    assert(c2p_works)

#end def test_coord_conversion


def test_grid_initialization():

    from testing import object_eq
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

    p1 = ParallelotopeGrid(
        shape = (6,),
        axes  = [[1]],
        )
    grids.append(p1)

    p2 = ParallelotopeGrid(
        cells = (5,),
        axes  = [[1]],
        )
    grids.append(p2)

    p3 = ParallelotopeGrid(
        dr    = (0.2,),
        axes  = [[1]],
        )
    grids.append(p2)

    p4 = ParallelotopeGrid(
        shape = (5,),
        axes  = [[1]],
        centered = True,
        )
    grids.append(p4)

    p5 = ParallelotopeGrid(
        cells = (5,),
        axes  = [[1]],
        centered = True,
        )
    grids.append(p5)

    p6 = ParallelotopeGrid(
        dr    = (0.2,),
        axes  = [[1]],
        centered = True,
        )
    grids.append(p6)

    assert(object_eq(p1,p2))
    assert(object_eq(p1,p3))

    assert(object_eq(p4,p5))
    assert(object_eq(p4,p6))

    for g in grids:
        assert(g.valid())
    #end for
#end def test_grid_initialization


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
