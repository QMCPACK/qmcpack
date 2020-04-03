
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

        assert(hull==list(range(8)))
    #end def test_convex_hull
#end if



import numpy as np
def random_stream(n):

    def rng(m=2**32, a=1103515245, c=12345):
        rng.current = (a*rng.current + c) % m
        return float(rng.current)/m
    #end def rng

    rng.current = 1

    return np.array([rng() for i in range(n)],dtype=float)
#end def random_stream


rstream = random_stream(1000)
rstream_wide = np.array(10*[rstream])



def test_simplestats():
    import numpy as np
    from testing import value_eq
    from numerics import simplestats

    m,e = simplestats(rstream)

    m_ref = 0.507154277182
    e_ref = 0.00916655521114

    assert(value_eq(float(m),m_ref))
    assert(value_eq(float(e),e_ref))

    m_ref = np.array(10*[m_ref])
    e_ref = np.array(10*[e_ref])

    m,e = simplestats(rstream_wide)

    assert(value_eq(m,m_ref))
    assert(value_eq(e,e_ref))
#end def test_simplestats



def test_simstats():
    import numpy as np
    from testing import value_eq
    from numerics import simstats

    m,v,e,k = simstats(rstream)

    m_ref = 0.507154277182
    v_ref = 0.0840257344388
    e_ref = 0.00949486225214
    k_ref = 1.07291426596

    assert(value_eq(float(m),m_ref))
    assert(value_eq(float(v),v_ref))
    assert(value_eq(float(e),e_ref))
    assert(value_eq(float(k),k_ref))

    m_ref = np.array(10*[m_ref])
    v_ref = np.array(10*[v_ref])
    e_ref = np.array(10*[e_ref])
    k_ref = np.array(10*[k_ref])
    
    m,v,e,k = simstats(rstream_wide)
    
    assert(value_eq(m,m_ref))
    assert(value_eq(v,v_ref))
    assert(value_eq(e,e_ref))
    assert(value_eq(k,k_ref))
#end def test_simstats



def test_equilibration_length():
    import numpy as np
    from numerics import equilibration_length

    eq = equilibration_length(rstream,random=False)
    assert(eq==0)

    rs = rstream + np.exp(-np.arange(len(rstream),dtype=float)/10)

    eq = equilibration_length(rs,random=False)
    assert(eq==52)

    eq = equilibration_length(rs,random=True)
#end def test_equilibration_length



if scipy_available:
    def test_ttest():
        import numpy as np
        from testing import value_eq
        from numerics import ttest

        p = ttest(0.,1.,100,0.,1.,100)
        assert(value_eq(float(p),1.0))

        p = ttest(0.,1.,10000,1.,1.,10000)
        assert(value_eq(float(p),0.479508361523))
    #end def test_ttest
#end if



def test_jackknife():
    from testing import value_eq
    from numerics import jackknife

    def value(v):
        return v
    #end def value

    jm_ref = np.array(10*[0.50715428],dtype=float)
    je_ref = np.array(10*[0.00917114],dtype=float)

    jm,je = jackknife(rstream_wide.T,value)

    assert(value_eq(jm,jm_ref))
    assert(value_eq(je,je_ref))
#end def test_jackknife



def test_morse():
    from testing import value_eq
    from unit_converter import convert
    from numerics import morse,morse_re,morse_a,morse_De,morse_Einf,morse_width
    from numerics import morse_depth,morse_Ee,morse_k,morse_params
    from numerics import morse_reduced_mass,morse_freq,morse_w,morse_wX
    from numerics import morse_E0,morse_En,morse_zero_point,morse_harmfreq
    from numerics import morse_harmonic_potential,morse_spect_fit
    from numerics import morse_rDw_fit,morse_fit,morse_fit_fine

    rm = morse_reduced_mass('Ti','O')
    assert(value_eq(rm,21862.2266134))

    r_ref,D_ref,w_ref = 1.620,6.87,1009.18
    r_ref_A = r_ref
    p = morse_rDw_fit(r_ref,D_ref,w_ref,'Ti','O',Dunit='eV')

    r_ref = convert(r_ref,'A','B')
    D_ref = convert(D_ref,'eV','Ha')

    r = morse_re(p)
    D = morse_De(p)
    w = morse_w(p,'Ti','O')
    Einf = morse_Einf(p)
    assert(value_eq(float(r),r_ref))
    assert(value_eq(float(D),D_ref))
    assert(value_eq(float(w),w_ref))
    assert(value_eq(float(Einf),0.0))

    width_ref = 1.0451690611
    width = morse_width(p)
    assert(value_eq(float(width),width_ref))

    depth = morse_depth(p)
    assert(value_eq(depth,D_ref))

    a = morse_a(p)
    assert(value_eq(float(a),1/width_ref))

    Ee = morse_Ee(p)
    assert(value_eq(float(Ee),-D_ref))

    k_ref = 0.462235185922
    k = morse_k(p)
    assert(value_eq(float(k),k_ref))

    assert(value_eq(morse_params(r,a,D,Einf),p))

    assert(value_eq(morse_freq(p,'Ti','O'),w))

    #wX_ref = 2.43157675597e-08 # might be buggy
    #wX = morse_wX(p,'Ti','O')
    #assert(value_eq(float(wX),wX_ref))

    E0_ref = -0.250174011529
    E0 = morse_E0(p,'Ti','O')
    assert(value_eq(float(E0),E0_ref))

    E0 = morse_En(p,0,'Ti','O')
    assert(value_eq(float(E0),E0_ref))

    En_ref = [-0.250174011529,
              -0.245617721989,
              -0.241103305298,
              -0.236630761457,
              -0.232200090465]

    for n in range(5):
        En = morse_En(p,n,'Ti','O')
        assert(value_eq(float(En),En_ref[n]))
    #end for

    zp_ref = 0.00229384708885
    zp = morse_zero_point(p,'Ti','O')
    assert(value_eq(float(zp),zp_ref))

    # not consistent w/ morse_wX
    #p_sf = morse_spect_fit(r_ref_A,w_ref,wX_ref,'Ti','O')

    hf_ref = 0.00459816239012
    hf = morse_harmfreq(p,'Ti','O')
    assert(value_eq(float(hf),hf_ref))

    E = morse(p,r_ref)
    assert(value_eq(float(E),-D_ref))
#end def test_morse



if scipy_available:
    def test_morse_fit():
        import numpy as np
        from testing import value_eq
        from unit_converter import convert
        from numerics import morse
        from numerics import morse_rDw_fit,morse_fit,morse_fit_fine

        r_ref,D_ref,w_ref = 1.620,6.87,1009.18
        r_ref_A = r_ref
        p = morse_rDw_fit(r_ref,D_ref,w_ref,'Ti','O',Dunit='eV')

        r_ref = convert(r_ref,'A','B')


        rfine = np.linspace(0.8*r_ref,1.2*r_ref,100)
        Efine = morse(p,rfine)

        pf = morse_fit(rfine,Efine,p)

        pref = tuple(np.array(p,dtype=float))
        pf   = tuple(np.array(pf,dtype=float))

        assert(value_eq(pf,pref))

        pf,Ef = morse_fit_fine(rfine,Efine,p,rfine,both=True)
        pf   = tuple(np.array(pf,dtype=float))
        assert(value_eq(pf,pref))
        assert(value_eq(Ef,Efine))
    #end def test_morse_fit
#end if




def test_eos():
    import numpy as np
    from testing import value_eq
    from unit_converter import convert
    from numerics import eos_fit,eos_eval,eos_param

    data = np.array([
            [0.875, -83.31851261], 
            [0.900, -83.38085214], 
            [0.925, -83.42172843], 
            [0.950, -83.44502216], 
            [0.975, -83.45476035], 
            [1.025, -83.44564229], 
            [1.050, -83.43127254], 
            [1.000, -83.45412846], 
            [1.100, -83.39070714], 
            [1.125, -83.36663810], 
            ])

    a = 5.539 # lattice constant

    V = (a*data[:,0])**3
    E = convert(data[:,1],'Ry','eV')

    # done originally to get params below
    #pf = eos_fit(V,E,'vinet')

    Einf = -1.13547294e+03
    V0   =  1.62708941e+02
    B0   =  1.34467867e-01
    Bp0  =  4.55846963e+00

    pf = Einf,V0,B0,Bp0

    Ef = eos_eval(pf,V,'vinet')

    assert(value_eq(Ef,E,atol=4e-3))

    assert(value_eq(float(eos_param(pf,'Einf','vinet')),Einf))
    assert(value_eq(float(eos_param(pf,'V','vinet')),V0))
    assert(value_eq(float(eos_param(pf,'B','vinet')),B0))
    assert(value_eq(float(eos_param(pf,'Bp','vinet')),Bp0))
#end def test_eos



if scipy_available:
    def test_eos_fit():
        import numpy as np
        from testing import value_eq
        from unit_converter import convert
        from numerics import eos_fit,eos_eval,eos_param

        data = np.array([
                [0.875, -83.31851261], 
                [0.900, -83.38085214], 
                [0.925, -83.42172843], 
                [0.950, -83.44502216], 
                [0.975, -83.45476035], 
                [1.025, -83.44564229], 
                [1.050, -83.43127254], 
                [1.000, -83.45412846], 
                [1.100, -83.39070714], 
                [1.125, -83.36663810], 
                ])

        a = 5.539 # lattice constant

        V = (a*data[:,0])**3
        E = convert(data[:,1],'Ry','eV')

        # done originally to get params below
        #pf = eos_fit(V,E,'vinet')

        Einf = -1.13547294e+03
        V0   =  1.62708941e+02
        B0   =  1.34467867e-01
        Bp0  =  4.55846963e+00

        pf = Einf,V0,B0,Bp0

        Ef = eos_eval(pf,V,'vinet')

        pf2 = eos_fit(V,Ef,'vinet')

        pf  = np.array(pf,dtype=float)
        pf2 = np.array(pf2,dtype=float)

        assert(value_eq(pf,pf2,atol=1e-3))
    #end def test_eos_fit
#end if
