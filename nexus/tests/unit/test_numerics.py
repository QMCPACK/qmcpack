
from versions import scipy_available


def test_import():
    import numerics
    from numerics import curve_fit
    from numerics import morse,morse_re,morse_a,morse_De,morse_Einf,morse_width
    from numerics import morse_depth,morse_Ee,morse_k,morse_params
    from numerics import morse_reduced_mass,morse_freq,morse_w,morse_wX
    from numerics import morse_En,morse_zero_point,morse_harmfreq
    from numerics import morse_harmonic_potential,morse_spect_fit
    from numerics import morse_rDw_fit,morse_fit,morse_fit_fine
    from numerics import murnaghan,birch,vinet
    from numerics import murnaghan_pressure,birch_pressure,vinet_pressure
    from numerics import eos_param,eos_eval,eos_fit,eos_Einf,eos_V,eos_B,eos_Bp
    from numerics import jackknife,jackknife_aux,check_jackknife_inputs
    from numerics import ndgrid
    from numerics import simstats,simplestats,equilibration_length,ttest
    from numerics import surface_normals,simple_surface
    from numerics import func_fit
    from numerics import distance_table,nearest_neighbors,voronoi_neighbors
    from numerics import convex_hull

#end def test_import



def test_ndgrid():
    import numpy as np
    from testing import value_eq
    from numerics import ndgrid

    x = [0,1,2.]
    y = [3,4.,5,6]
    z = [7,8,9,10,11.]

    points = np.zeros((3,len(x)*len(y)*len(z)),dtype=float)
    n = 0
    for xv in x:
        for yv in y:
            for zv in z:
                points[:,n] = xv,yv,zv
                n += 1
            #end for
        #end for
    #end for
    points = np.array(points)
    points.shape = 3,len(x),len(y),len(z)
    points = points

    grid = ndgrid(x,y,z)

    assert(value_eq(grid,points))
#end def test_ndgrid



def test_distance_table():
    import numpy as np
    from testing import value_eq
    from numerics import distance_table

    points = [
        [0,0,0],
        [1,0,0],
        [0,1,0],
        [0,0,1],
        [1,1,0],
        [1,0,1],
        [0,1,1],
        [1,1,1],
        ]
    points = np.array(points,dtype=float)

    d0 = 0.0
    d1 = 1.0
    d2 = np.sqrt(2.)
    d3 = np.sqrt(3.)

    dtable = [[d0,d1,d1,d1,d2,d2,d2,d3],
              [d1,d0,d2,d2,d1,d1,d3,d2],
              [d1,d2,d0,d2,d1,d3,d1,d2],
              [d1,d2,d2,d0,d3,d1,d1,d2],
              [d2,d1,d1,d3,d0,d2,d2,d1],
              [d2,d1,d3,d1,d2,d0,d2,d1],
              [d2,d3,d1,d1,d2,d2,d0,d1],
              [d3,d2,d2,d2,d1,d1,d1,d0]]
    dtable = np.array(dtable,dtype=float)
    dt = distance_table(points,points)
    assert(value_eq(dt,dtable))

    dtable = [[d0,d1,d1,d1],
              [d1,d0,d2,d2],
              [d1,d2,d0,d2],
              [d1,d2,d2,d0],
              [d2,d1,d1,d3],
              [d2,d1,d3,d1],
              [d2,d3,d1,d1],
              [d3,d2,d2,d2]]
    dtable = np.array(dtable,dtype=float)
    dt = distance_table(points,points[0:4])
    assert(value_eq(dt,dtable))

    points = [
        [0,0,0],
        [1,1,0],
        [1,1,1],
        ]
    points = np.array(points,dtype=float)

    dtable = [[d0,d2,d3],
              [d0,d1,d2],
              [d0,d1,d3]]
    dtable = np.array(dtable,dtype=float)
    dtable_order = [[0,1,2],
                    [1,2,0],
                    [2,1,0]]
    dtable_order = np.array(dtable_order,dtype=int)
    dt,order = distance_table(points,points,ordering=1)
    assert(value_eq(dt,dtable))
    assert(value_eq(order,dtable_order))

    dtable = [[d0,d2,d3],
              [d0,d1,d2],
              [d0,d1,d3]]
    dtable = np.array(dtable,dtype=float)
    dtable_order = [[0,1,2],
                    [1,2,0],
                    [2,1,0]]
    dtable_order = np.array(dtable_order,dtype=int)
    dt,order = distance_table(points,points,ordering=2)
    assert(value_eq(dt,dtable))
    assert(value_eq(order,dtable_order))

#end def test_distance_table
    


def test_nearest_neighbors():
    import numpy as np
    from testing import value_eq
    from numerics import nearest_neighbors

    points = [
        [0,0,0],
        [1,0,0],
        [0,1,0],
        [1,1,0],
        [0,0,1],
        [1,0,1],
        [0,1,1],
        [1,1,1],
        ]
    points = np.array(points,dtype=float)

    plow  = points[:4]
    phigh = points[4:]

    d0 = 0.0
    d1 = 1.0
    d2 = np.sqrt(2.)
    d3 = np.sqrt(3.)

    nn_ref = [[0,1,2],
              [1,0,3],
              [2,0,3],
              [3,1,2]]
    nn_ref = np.array(nn_ref,dtype=int)

    dist_ref = [[d1,d2,d2],
                [d1,d2,d2],
                [d1,d2,d2],
                [d1,d2,d2]]
    dist_ref = np.array(dist_ref,dtype=float)

    def check_nn(nn):
        for n,nr in zip(nn,nn_ref):
            assert(n[0]==nr[0])
            assert(set(n[1:])==set(nr[1:]))
        #end for
    #end def check_nn

    nn = nearest_neighbors(3,plow,phigh,slow=True)
    check_nn(nn)

    nn,dist = nearest_neighbors(3,plow,phigh,return_distances=True,slow=True)
    check_nn(nn)
    assert(value_eq(dist,dist_ref))

    if scipy_available:
        nn = nearest_neighbors(3,plow,phigh)
        check_nn(nn)

        nn,dist = nearest_neighbors(3,plow,phigh,return_distances=True)
        check_nn(nn)
        assert(value_eq(dist,dist_ref))
    #end if

#end def test_nearest_neighbors



if scipy_available:
    def test_voronoi_neighbors():
        import numpy as np
        from testing import value_eq
        from numerics import voronoi_neighbors

        points = [
            [0,0,0],
            [1,0,0],
            [0,1,0],
            [1,1,0],
            [0,0,1],
            [1,0,1],
            [0,1,1],
            [1,1,1],
            ]
        points = np.array(points,dtype=float)

        nn_pairs_ref = [[5,7],
                        [5,4],
                        [5,1],
                        [5,3],
                        [5,6],
                        [2,0],
                        [2,3],
                        [2,6],
                        [2,4],
                        [0,4],
                        [0,3],
                        [1,4],
                        [1,3],
                        [7,3],
                        [7,6],
                        [3,4],
                        [3,6],
                        [6,4]]
        nn_pairs_ref = np.array(nn_pairs_ref,dtype=int)

        nn_pairs = voronoi_neighbors(points)

        print nn_pairs
        return

        assert(isinstance(nn_pairs,np.ndarray))
        nn_pairs = np.array(nn_pairs,dtype=int)
        assert(value_eq(nn_pairs,nn_pairs_ref))

        for i,j in nn_pairs:
            print i,j
            print np.linalg.norm(points[i]-points[j])
        #end for
    #end def test_voronoi_neighbors
#end if
