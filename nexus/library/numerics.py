##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  numerics.py                                                       #
#    A collection of useful numerical functions, currently           #
#    including specialized curve fitting and jack-knife statistical  #
#    analysis.                                                       #
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
#====================================================================#


import types as pytypes
from numpy import array,ndarray,zeros,linspace,pi,exp,sqrt,polyfit,polyval
from generic import obj
from developer import unavailable,warn,error
from unit_converter import convert
from periodic_table import pt as ptable
from debug import *
try:
    from scipy.optimize import fmin
except ImportError:
    fmin = unavailable('scipy.optimize','fmin')
#end try


def numerics_error(*args,**kwargs):
    error(*args,**kwargs)
#end def numerics_error
# minimizers
least_squares = lambda p,x,y,f: ((f(p,x)-y)**2).sum()
absmin        = lambda p,x,y,f: abs(f(p,x)-y).sum()
madmin        = lambda p,x,y,f: abs(f(p,x)-y).max()
 
# curve fit based on fmin from scipy
def curve_fit(x,y,f,p0,minimizer=least_squares,optimizer='fmin'):
    if optimizer=='fmin':
        p = fmin(minimizer,p0,args=(x,y,f),maxiter=10000,maxfun=10000,disp=0)
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
        m = m1
    else:
        if isinstance(m2,str):
            m2 = ptable[m2].atomic_weight.me
        #end if
        m = 1./(1./m1+1./m2) # reduced mass
    #end if
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
def morse_fit(r,E,p0=None,jackknife=False,minimizer=least_squares,auxfuncs=None,auxres=None,capture=None):
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
        capture.minimizer = minimizer
        capture.auxfuncs  = auxfuncs
        capture.auxres    = auxres
        capture.Edata     = Edata
        capture.pp        = pp
    elif calc_aux:
        jcapture = obj()
    #end if

    # get an optimized morse fit of the means
    pf = curve_fit(r,E,morse,p0,minimizer)

    # obtain jackknife (mean+error) estimates of fitted parameters and/or fitted curves
    pmean  = None
    perror = None
    if jackknife:
        if Edata is None:
            numerics_error('cannot perform jackknife fit because blocked data was not provided (only the means are present)','morse_fit')
        #end if
        pmean,perror = numerics_jackknife(data     = Edata,
                                          function = curve_fit,
                                          args     = [r,None,morse,pf,minimizer],
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
def morse_fit_fine(r,E,p0=None,rfine=None,both=False,jackknife=False,minimizer=least_squares,capture=None):  
    if rfine is None:
        rfine = linspace(r.min(),r.max(),400)
    #end if
    auxfuncs = obj(
        Efine = (morse,[None,rfine])
        )
    auxres = obj()

    res = morse_fit(r,E,p0,jackknife,minimizer,auxfuncs,auxres,capture)

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
    if type(auxfunc)!=pytypes.FunctionType:
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
