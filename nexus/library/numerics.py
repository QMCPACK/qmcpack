##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  numerics.py                                                       #
#    A collection of useful numerical functions, currently           #
#    including specialized curve fitting, statistical analysis,      #
#    and spatial analysis.                                           #
#                                                                    #
#  Content summary:                                                  #
#    morse_spect_fit                                                 #
#      Return a Morse potential fit from experimental data           #
#      (eqm. radius, vib. frequency, w_eX_e).                        #
#                                                                    #
#    morse_fit                                                       #
#      Perform a Morse potential fit of simulation data.             #
#      If the data are stochastic, return fits with error bars.      #
#                                                                    #
#    For the morse-related functions listed, see documentation below #
#      morse_freq                                                    #
#      morse_w                                                       #
#      morse_wX                                                      #
#      morse_E0                                                      #
#      morse_En                                                      #
#      morse_zero_point                                              #
#      morse_harmfreq                                                #
#      morse_harmonic_potential                                      #
#                                                                    #
#    jackknife                                                       #
#      Perform jack-knife statistical analysis accepting an          #
#      arbitrary function of N-dimensional simulation data.          #
#      Can be used to obtain error bars of fit parameters,           #
#      eigenvalues, and other statistical results that depend on     #
#      the input data in a non-linear fashion.                       #
#                                                                    #
#    ndgrid                                                          #
#      Function to construct an arbitrary N-dimensional grid.        #
#      Similar to ndgrid from MATLAB.                                #
#                                                                    #
#    simstats                                                        #
#      Compute statistics of N-dimensional Monte Carlo simulation    #
#      data, including mean, variance, error, and autocorrelation.   #
#                                                                    #
#    simplestats                                                     #
#      Compute error assuming uncorrelated data.                     #
#                                                                    #
#    equilibration_length                                            #
#      Estimate equilibration point of Monte Carlo time series data  #
#      using a heuristic algorithm.                                  #
#                                                                    #
#    ttest                                                           #
#      Implementation of Student's T-test                            #
#                                                                    #
#    surface_normals                                                 #
#      Compute vectors normal to a parametric surface.               #
#                                                                    #
#    simple_surface                                                  #
#      Create a parametric surface in Cartesian, cylindrical, or     #
#      spherical coordinates.                                        #
#                                                                    #
#    func_fit                                                        #
#      Perform a fit to an arbitrary function using an arbitrary     #
#      cost metric (e.g. least squares).                             #
#                                                                    #
#    distance_table                                                  #
#      Calculate all N^2 pair distances for a set of N points.       #
#                                                                    #
#    nearest_neighbors                                               #
#      Find k nearest neighbors of N points using a fast algorithm.  #
#                                                                    #
#    voronoi_neighbors                                               #
#      Find nearest neighbors in the Voronoi sense, that is for      #
#      each point, find the points whose Voronoi polyhedra contact   #
#      the Voronoi polyhedron of that point.                         #
#                                                                    #
#    convex_hull                                                     #
#      Find the convex hull of a set of points in N dimensions.      #
#                                                                    #        
#====================================================================#


import sys
import inspect
from numpy import array,ndarray,zeros,linspace,pi,exp,sqrt,polyfit,polyval
from numpy import sum,abs,arange,empty,sin,cos,dot,atleast_2d,ogrid
from numpy import ones_like,sign,random,cross,prod
from numpy.linalg import norm
from generic import obj
from developer import unavailable,warn,error
from unit_converter import convert
from periodic_table import pt as ptable
try:
    from scipy.special import betainc
    from scipy.optimize import fmin
    from scipy.spatial import KDTree,Delaunay,Voronoi
    scipy_unavailable = False
except ImportError:
    betainc = unavailable('scipy.special' ,'betainc')
    fmin    = unavailable('scipy.optimize','fmin')
    KDTree,Delaunay,Voronoi  = unavailable('scipy.spatial' ,'KDTree','Delaunay','Voronoi')
    scipy_unavailable = True
#end try


def numerics_error(*args,**kwargs):
    error(*args,**kwargs)
#end def numerics_error
# cost functions
least_squares = lambda p,x,y,f: ((f(p,x)-y)**2).sum()
absmin        = lambda p,x,y,f: abs(f(p,x)-y).sum()
madmin        = lambda p,x,y,f: abs(f(p,x)-y).max()

cost_functions = obj(
    least_squares = least_squares,
    absmin        = absmin,
    madmin        = madmin,
    )
 
# curve fit based on fmin from scipy
def curve_fit(x,y,f,p0,cost='least_squares',optimizer='fmin'):
    if isinstance(cost,str):
        if cost not in cost_functions:
            numerics_error('"{0}" is an invalid cost function\nvalid options are: {1}'.format(cost,sorted(cost_functions.keys())))
        #end if
        cost = cost_functions[cost]
    #end if
    if optimizer=='fmin':
        p = fmin(cost,p0,args=(x,y,f),maxiter=10000,maxfun=10000,disp=0)
    else:
        numerics_error('optimizers other than fmin are not supported yet','curve_fit')
    return p
#end def curve_fit


# morse potential
#               V(r) =  De ( (1-e^-a(r-re))^2 - 1 ) + E_infinity
morse = lambda p,r: p[2]*((1-exp(-(r-p[0])/p[1]))**2-1)+p[3]
morse_re    = lambda p: p[0]         # equilibrium separation
morse_a     = lambda p: 1./p[1]      # 'a' parameter, related to well width
morse_De    = lambda p: p[2]         # 'De' parameter, related to well depth
morse_Einf  = lambda p: p[3]         # potential energy at infinite separation
morse_width = lambda p: p[1]         # well width
morse_depth = lambda p: morse_De(p)  # well depth
morse_Ee    = lambda p: morse_Einf(p)-morse_De(p)   # potential energy at equilibrium
morse_k     = lambda p: 2*morse_De(p)*morse_a(p)**2 # force constant k = d2V/dr2(r=re), Vh=1/2 k r^2
morse_params= lambda re,a,De,E_inf: (re,1./a,De,E_inf)  # return p given standard inputs

# morse_reduced_mass gives the reduced mass in Hartree units
#   m1 and m2 are masses or atomic symbols
def morse_reduced_mass(m1,m2=None):
    if isinstance(m1,str):
        m1 = ptable[m1].atomic_weight.me
    #end if
    if m2 is None:
        m2 = m1
    elif isinstance(m2,str):
        m2 = ptable[m2].atomic_weight.me
    #end if
    m = 1./(1./m1+1./m2) # reduced mass
    return m
#end def morse_reduced_mass    

# morse_freq returns anharmonic frequency in 1/cm if curve is in Hartree units
def morse_freq(p,m1,m2=None): 
    alpha = 7.2973525698e-3           # fine structure constant
    c     = 1./alpha                  # speed of light, hartree units
    m     = morse_reduced_mass(m1,m2) # reduced mass
    lam   = 2*pi*c*sqrt(m/morse_k(p)) # wavelength
    freq  = 1./(convert(lam,'B','m')*100)
    return freq
#end def morse_freq

# w = \omega_e or frequency in 1/cm
def morse_w(p,m1,m2=None):
    return morse_freq(p,m1,m2)
#end def morse_w

# wX = \omega_e\Chi_e spectroscopic constant in 1/cm
def morse_wX(p,m1,m2=None):
    ocm   = 1./(convert(1.0,'B','m')*100) # 1/Bohr to 1/cm
    alpha = 7.2973525698e-3               # fine structure constant
    m     = morse_reduced_mass(m1,m2)     # reduced mass
    wX = (alpha*morse_a(p)**2)/(4*pi*m)
    return wX
#end def morse_wX

# ground state energy (Hartree units in and out, neglects E_infinity)
#   true ground state (or binding) energy is E0-E_infinity
def morse_E0(p,m1,m2=None):
    m  = morse_reduced_mass(m1,m2)
    E0 = .5*sqrt(morse_k(p)/m) - morse_a(p)**2/(8*m) + morse_Ee(p)
    return E0
#end def morse_E0

# energy of nth vibrational level (Hartree units in and out, neglects E_infinity)
def morse_En(p,n,m1,m2=None):
    m  = morse_reduced_mass(m1,m2)
    En = sqrt(morse_k(p)/m)*(n+.5) - morse_a(p)**2/(2*m)*(n+.5)**2 + morse_Ee(p)
    return En
#end def morse_En

# morse_zero_point gives the zero point energy (always positive)
def morse_zero_point(p,m1,m2=None):
    return morse_E0(p,m1,m2)-morse_Ee(p)
#end def morse_zero_point

# morse_harmfreq returns the harmonic frequency (Hartree units in and out)
def morse_harmfreq(p,m1,m2=None): 
    m     = morse_reduced_mass(m1,m2)
    hfreq = sqrt(morse_k(p)/m)
    return hfreq
#end def morse_harmfreq

# morse_harmonic evaluates the harmonic oscillator fit to the morse potential
def morse_harmonic_potential(p,r):
    return .5*morse_k(p)*(r-morse_re(p))**2 - morse_De(p)
#end def morse_harmonic_potential

# morse_spect_fit returns the morse fit starting from spectroscopic parameters
#   spect. params. are equilibrium bond length, vibration frequency, and wX parameter
#   input units are Angstrom for re and 1/cm for w and wX
#   m1 and m2 are masses in Hartree units, only one need be provided
#   outputted fit is in Hartree units
def morse_spect_fit(re,w,wX,m1,m2=None,Einf=0.0):
    alpha = 7.2973525698e-3            # fine structure constant
    m     = morse_reduced_mass(m1,m2)  # reduced mass
    ocm_to_oB = 1./convert(.01,'m','B')# conversion from 1/cm to 1/Bohr
    w    *= ocm_to_oB                  # convert from 1/cm to 1/Bohr
    wX   *= ocm_to_oB                  # convert from 1/cm to 1/Bohr
    re    = convert(re,'A','B')        # convert from Angstrom to Bohr
    a     = sqrt(4*pi*m*wX/alpha)      # get the 'a' parameter
    De    = (pi*w**2)/(2*alpha*wX)     # get the 'De' parameter
    p     = morse_params(re,a,De,Einf) # get the internal fit parameters
    return p
#end def morse_spect_fit


def morse_rDw_fit(re,De,w,m1,m2=None,Einf=0.0,Dunit='eV'):
    alpha = 7.2973525698e-3            # fine structure constant
    m     = morse_reduced_mass(m1,m2)  # reduced mass
    ocm_to_oB = 1./convert(.01,'m','B')# conversion from 1/cm to 1/Bohr
    w    *= ocm_to_oB                  # convert from 1/cm to 1/Bohr
    re    = convert(re,'A','B')        # convert from Angstrom to Bohr
    De    = convert(De,Dunit,'Ha')     # convert from input energy unit
    wX    = (pi*w**2)/(2*alpha*De)     # get wX
    a     = sqrt(4*pi*m*wX/alpha)      # get the 'a' parameter
    p     = morse_params(re,a,De,Einf) # get the internal fit parameters
    return p
#end def morse_rDw_fit


# morse_fit computes a morse potential fit to r,E data
#  fitting through means, E is one dimensional
#    pf    = morse_fit(r,E)                           returns fitted parameters
#  jackknife statistical fits, E is two dimensional with blocks as first dimension
#    pf,pmean,perror = morse_fit(r,E,jackknife=True)  returns jackknife estimates of parameters
def morse_fit(r,E,p0=None,jackknife=False,cost=least_squares,auxfuncs=None,auxres=None,capture=None):
    if isinstance(E,(list,tuple)):
        E = array(E,dtype=float)
    #end if

    Edata = None
    if len(E)!=E.size and len(E.shape)==2:
        Edata = E
        E     = Edata.mean(axis=0)
    #end if

    pp = None
    if p0 is None:
        # make a simple quadratic fit to get initial guess for morse fit
        pp = polyfit(r,E,2)
        r0   = -pp[1]/(2*pp[0])
        E0   = pp[2]+.5*pp[1]*r0
        d2E  = 2*pp[0]
        Einf = E[-1] #0.0
        #  r_eqm, pot_width, E_bind, E_infinity
        p0 = r0,sqrt(2*(Einf-E0)/d2E),Einf-E0,Einf
    #end if
    
    calc_aux = auxfuncs!=None and auxres!=None
    capture_results = capture!=None
    jcapture    = None
    jauxcapture = None
    if capture_results:
        jcapture    = obj()
        jauxcapture = obj()
        capture.r         = r
        capture.E         = E
        capture.p0        = p0
        capture.jackknife = jackknife
        capture.cost = cost
        capture.auxfuncs  = auxfuncs
        capture.auxres    = auxres
        capture.Edata     = Edata
        capture.pp        = pp
    elif calc_aux:
        jcapture = obj()
    #end if

    # get an optimized morse fit of the means
    pf = curve_fit(r,E,morse,p0,cost)

    # obtain jackknife (mean+error) estimates of fitted parameters and/or fitted curves
    pmean  = None
    perror = None
    if jackknife:
        if Edata is None:
            numerics_error('cannot perform jackknife fit because blocked data was not provided (only the means are present)','morse_fit')
        #end if
        pmean,perror = numerics_jackknife(data     = Edata,
                                          function = curve_fit,
                                          args     = [r,None,morse,pf,cost],
                                          position = 1,
                                          capture  = jcapture)
        # compute auxilliary jackknife quantities, if desired (e.g. morse_freq, etc)
        if calc_aux:
            psamples = jcapture.jsamples
            for auxname,auxfunc in auxfuncs.iteritems():
                auxcap = None
                if capture_results:
                    auxcap = obj()
                    jauxcapture[auxname] = auxcap
                #end if
                auxres[auxname] = jackknife_aux(psamples,auxfunc,capture=auxcap)
            #end for
        #end if
    #end if

    if capture_results:
        capture.pmean  = pmean
        capture.perror = perror
        capture.jcapture     = jcapture
        capture.jauxcapture  = jauxcapture
    #end if

    # return desired results
    if not jackknife:
        return pf
    else:
        return pf,pmean,perror
    #end if
#end def morse_fit


# morse_fit_fine: fit data to a morse potential and interpolate on a fine grid
#   compute direct jackknife variations in the fitted curves 
#   by using morse as an auxilliary jackknife function
def morse_fit_fine(r,E,p0=None,rfine=None,both=False,jackknife=False,cost=least_squares,capture=None):  
    if rfine is None:
        rfine = linspace(r.min(),r.max(),400)
    #end if
    auxfuncs = obj(
        Efine = (morse,[None,rfine])
        )
    auxres = obj()

    res = morse_fit(r,E,p0,jackknife,cost,auxfuncs,auxres,capture)

    if not jackknife:
        pf = res
    else:
        pf,pmean,perror = res
    #end if

    Efine = morse(pf,rfine)
    
    if not jackknife:
        if not both:
            return Efine
        else:
            return pf,Efine
        #end if
    else:
        Emean,Eerror = auxres.Efine
        if not both:
            return Efine,Emean,Eerror
        else:
            return pf,pmean,perror,Efine,Emean,Eerror
        #end if
    #end if
#end def morse_fit_fine
 

# equation of state
murnaghan  = lambda p,V: p[0] + p[2]/p[3]*V*((p[1]/V)**p[3]/(p[3]-1)+1)-p[1]*p[2]/(p[3]-1)
birch      = lambda p,V: p[0] + 9*p[1]*p[2]/16*((p[1]/V)**(2./3)-1)**2*( 2 + (p[3]-4)*((p[1]/V)**(2./3)-1) )
vinet      = lambda p,V: p[0] + 2*p[1]*p[2]/(p[3]-1)**2*( 2 - (2+3*(p[3]-1)*((V/p[1])**(1./3)-1))*exp(-1.5*(p[3]-1)*((V/p[1])**(1./3)-1)) ) 

murnaghan_pressure = lambda p,V: p[1]/p[2]*((p[0]/V)**p[2]-1)
birch_pressure     = lambda p,V: 1.5*p[1]*(p[0]/V)**(5./3)*((p[0]/V)**(2./3)-1)*(1.+.75*(p[2]-1)*((p[0]/V)**(2./3)-1))
vinet_pressure     = lambda p,V: 3.*p[1]*(1.-(V/p[0])**(1./3))*(p[0]/V)**(2./3)*exp(1.5*(p[2]-1)*(1.-(V/p[0])**(1./3)))


eos_funcs = obj(
    murnaghan = murnaghan,
    birch     = birch,
    vinet     = vinet,
    )


eos_Einf = lambda p: p[0]  # energy at infinite separation
eos_V    = lambda p: p[1]  # equilibrium volume
eos_B    = lambda p: p[2]  # bulk modulus
eos_Bp   = lambda p: p[3]  # B prime

eos_param_tmp = obj(
    Einf = eos_Einf,
    V    = eos_V,
    B    = eos_B,
    Bp   = eos_Bp,
    )
eos_param_funcs = obj(
    murnaghan = eos_param_tmp,
    birch     = eos_param_tmp,
    vinet     = eos_param_tmp,
    )

def eos_eval(p,V,type='vinet'):
    if type not in eos_funcs:
        numerics_error('"{0}" is not a valid EOS type\nvalid options are: {1}'.format(sorted(eos_funcs.keys())))
    #end if
    return eos_funcs[type](p,V)
#end def eos_eval


def eos_param(p,param,type='vinet'):
    if type not in eos_param_funcs:
        numerics_error('"{0}" is not a valid EOS type\nvalid options are: {1}'.format(sorted(eos_param_funcs.keys())))
    #end if
    eos_pfuncs = eos_param_funcs[type]
    if param not in eos_pfuncs:
        numerics_error('"{0}" is not an available parameter for a {1} fit\navailable parameters are: {2}'.format(param,type,sorted(eos_pfuncs.keys())))
    #end if
    return eos_pfuncs[param](p)
#end def eos_param


def eos_fit(V,E,type='vinet',p0=None,cost='least_squares',optimizer='fmin',jackknife=False):
    if isinstance(V,(list,tuple)):
        V = array(V,dtype=float)
    #end if
    if isinstance(E,(list,tuple)):
        E = array(E,dtype=float)
    #end if
    
    if type not in eos_funcs:
        numerics_error('"{0}" is not a valid EOS type\nvalid options are: {1}'.format(sorted(eos_funcs.keys())))
    #end if
    eos_func = eos_funcs[type]

    if p0 is None:
        pp = polyfit(V,E,2)
        V0   = -pp[1]/(2*pp[0])
        B0   = -pp[1]
        Bp0  = 0.0
        Einf = E[-1] 
        p0 = Einf,V0,B0,Bp0
    #end if

    # get an optimized fit of the means
    pf = curve_fit(V,E,eos_func,p0,cost)

    return pf
#end def eos_fit




#==============#
#  Statistics  #
#==============#

# jackknife
#   data: a multi-dimensional array with blocks as the first dimension
#         nblocks = data.shape[0]
#   function: a function that takes an array for one block (e.g. data[0])
#             and returns an array or a tuple/list of scalars and arrays
#             ret = function(*args,**kwargs)
#   args: a list-like object of positional input arguments
#   kwargs: a dict-like object of keyword input arguments
#   position: location to place input_array in input arguments
#             if integer, will be placed in args:   args[position] = input_array
#             if string , will be placed in kwargs: kwargs[position] = input_array
#   capture: an object that will contain most jackknife info upon exit
def jackknife(data,function,args=None,kwargs=None,position=None,capture=None):
    capture_results = capture!=None
    if capture_results:
        capture.data         = data
        capture.function     = function
        capture.args         = args
        capture.kwargs       = kwargs
        capture.position     = position
        capture.jdata        = []
        capture.jsamples     = []
    #end if
    # check the requested argument position
    argpos,kwargpos,args,kwargs,position = check_jackknife_inputs(args,kwargs,position)

    # obtain sums of the jackknife samples
    nblocks = data.shape[0]
    nb = float(nblocks)
    jnorm   = 1./(nb-1.)
    data_sum = data.sum(axis=0)
    array_return = False
    for b in xrange(nblocks):
        jdata = jnorm*(data_sum-data[b])
        if argpos:
            args[position] = jdata
        elif kwargpos:
            kwargs[position] = jdata
        #end if
        jsample = function(*args,**kwargs)
        if b==0:
            # determine the return type from the first sample
            # and initialize the jackknife sums
            array_return = isinstance(jsample,ndarray)
            if array_return:
                jsum  = jsample.copy()
                jsum2 = jsum**2
            else:
                jsum  = []
                jsum2 = []
                for jval in jsample:
                    jsum.append(jval)
                    jsum2.append(jval**2)
                #end for
            #end if
        else:
            # accumulate the jackknife sums
            if array_return:
                jsum  += jsample
                jsum2 += jsample**2
            else:
                for c in xrange(len(jsample)):
                    jsum[c]  += jsample[c]
                    jsum2[c] += jsample[c]**2
                #end for
            #end if
        #end if
        if capture_results:
            capture.jdata.append(jdata)
            capture.jsamples.append(jsample)
        #end if
    #end for
    # get the jackknife mean and error
    if array_return:
        jmean = jsum/nb
        jerror = sqrt( (nb-1.)/nb*(jsum2-jsum**2/nb) )
    else:
        jmean  = []
        jerror = []
        for c in xrange(len(jsum)):
            jval  = jsum[c]
            jval2 = jsum2[c]
            jm = jval/nb
            je = sqrt( (nb-1.)/nb*(jval2-jval**2/nb) )
            jmean.append(jm)
            jerror.append(je)
        #end for
    #end if
    if capture_results:
        capture.data_sum     = data_sum
        capture.nblocks      = nblocks
        capture.array_return = array_return
        capture.jsum         = jsum
        capture.jsum2        = jsum2
        capture.jmean        = jmean
        capture.jerror       = jerror
    #end if
    return jmean,jerror
#end def jackknife
numerics_jackknife = jackknife


# get jackknife estimate of auxiliary quantities
#   jsamples is a subset of jsamples data computed by jackknife above
#   auxfunc is an additional function to get a jackknife sample of a derived quantity
def jackknife_aux(jsamples,auxfunc,args=None,kwargs=None,position=None,capture=None):
    # unpack the argument list if compressed
    if not inspect.isfunction(auxfunc):
        if len(auxfunc)==1:
            auxfunc = auxfunc[0]
        elif len(auxfunc)==2:
            auxfunc,args = auxfunc
        elif len(auxfunc)==3:
            auxfunc,args,kwargs = auxfunc
        elif len(auxfunc)==4:
            auxfunc,args,kwargs,position = auxfunc
        else:
            numerics_error('between 1 and 4 fields (auxfunc,args,kwargs,position) can be packed into original auxfunc input, received {0}'.format(len(auxfunc)))
        #end if
    #end if

    # check the requested argument position
    argpos,kwargpos,args,kwargs,position = check_jackknife_inputs(args,kwargs,position)

    capture_results = capture!=None
    if capture_results:
        capture.auxfunc  = auxfunc 
        capture.args     = args    
        capture.kwargs   = kwargs  
        capture.position = position
        capture.jdata    = []
        capture.jsamples = []
    #end if

    nblocks = len(jsamples)
    nb      = float(nblocks)
    for b in xrange(nblocks):
        jdata = jsamples[b]
        if argpos:
            args[position] = jdata
        elif kwargpos:
            kwargs[position] = jdata
        #end if
        jsample = auxfunc(*args,**kwargs)
        if b==0:
            jsum  = jsample.copy()
            jsum2 = jsum**2
        else:
            jsum  += jsample
            jsum2 += jsample**2
        #end if
        if capture_results:
            capture.jdata.append(jdata)
            capture.jsamples.append(jsample)
        #end if
    #end for
    jmean = jsum/nb
    jerror = sqrt( (nb-1.)/nb*(jsum2-jsum**2/nb) )

    if capture_results:
        capture.nblocks = nblocks
        capture.jsum    = jsum
        capture.jsum2   = jsum2
        capture.jmean   = jmean
        capture.jerror  = jerror
    #end if

    return jmean,jerror
#end def jackknife_aux


def check_jackknife_inputs(args,kwargs,position):
    argpos   = False
    kwargpos = False
    if position!=None:
        if isinstance(position,int):
            argpos = True
        elif isinstance(position,str):
            kwargpos = True
        else:
            numerics_error('position must be an integer or keyword, received: {0}'.format(position),'jackknife')
        #end if
    elif args is None and kwargs is None:
        args     = [None]
        argpos   = True
        position = 0
    elif kwargs is None and position is None:
        argpos   = True
        position = 0
    else:
        numerics_error('function argument position for input data must be provided','jackknife')
    #end if
    if args is None:
        args = []
    #end if
    if kwargs is None:
        kwargs = dict()
    #end if
    return argpos,kwargpos,args,kwargs,position
#end def check_jackknife_inputs









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
        sys.exit()
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



#least_squares = lambda p,x,y,f: ((f(p,x)-y)**2).sum()

def func_fit(x,y,fitting_function,p0,cost=least_squares):
    f = fitting_function
    p = fmin(cost,p0,args=(x,y,f),maxiter=10000,maxfun=10000)
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


def voronoi_neighbors(points):
    vor = Voronoi(points)
    neighbor_pairs = vor.ridge_points
    return neighbor_pairs
#end def voronoi_neighbors


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
