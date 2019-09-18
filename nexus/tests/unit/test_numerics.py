
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

        joggle = np.array(
            [[0.60801892, 0.68024807, 0.09037058],
             [0.95800898, 0.43112463, 0.52981569],
             [0.08862067, 0.69084511, 0.35177345],
             [0.37363091, 0.57409599, 0.95654043],
             [0.8310818 , 0.17146777, 0.90490215],
             [0.17600223, 0.89772462, 0.75582196],
             [0.7408217 , 0.22768522, 0.64564984],
             [0.71678216, 0.6409734 , 0.53354209]])
        joggle *= 1e-8

        points += joggle

        nn_pairs_ref = [[0,1],
                        [0,2],
                        [0,4],
                        [0,3],
                        [0,6],
                        [0,5],
                        [7,3],
                        [7,5],
                        [7,6],
                        [3,1],
                        [3,2],
                        [3,6],
                        [3,5],
                        [1,5],
                        [5,4],
                        [5,6],
                        [6,4],
                        [6,2]]
        nn_pairs_ref = np.array(nn_pairs_ref,dtype=int)

        d1 = 1.0
        d2 = float(np.sqrt(2.))

        nn_pairs = voronoi_neighbors(points)

        dist_ref = 3*[d1]+3*[d2]+5*[d1]+2*[d2]+2*[d1]+[d2]+2*[d1]

        assert(isinstance(nn_pairs,np.ndarray))
        nn_pairs = np.array(nn_pairs,dtype=int)
        assert(value_eq(nn_pairs,nn_pairs_ref))

        for n,(i,j) in enumerate(nn_pairs):
            d = np.linalg.norm(points[i]-points[j])
            assert(value_eq(float(d),dist_ref[n]))
        #end for

    #end def test_voronoi_neighbors



    def test_convex_hull():
        import numpy as np
        from testing import value_eq
        from numerics import convex_hull

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
        points = np.array(points,dtype=float)-0.5
        
        points = np.append(2*points,points,axis=0)

        hull = convex_hull(points)

        assert(hull==range(8))
    #end def test_convex_hull
#end if



def random_stream(n):
    import numpy as np

    def rng(m=2**32, a=1103515245, c=12345):
        rng.current = (a*rng.current + c) % m
        return float(rng.current)/m
    #end def rng

    rng.current = 1

    return np.array([rng() for i in range(n)],dtype=float)
#end def random_stream


rstream = random_stream(1000)
