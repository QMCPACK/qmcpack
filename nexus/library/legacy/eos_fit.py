##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  eos_fit.py                                                        #
#    Support for Equation of State (EOS) curve fitting.  Can perform #
#    Murnaghan, Birch-Murnaghan, or Vinet fits of energy vs. volume  #
#    or pressure vs. volume equation of state data.  To be merged    #
#    with more advanced facilities contained in numerics.py.         #
#                                                                    #
#  Content summary:                                                  #
#    energy_fit                                                      #
#      Perform an energy vs. volume fit.                             #
#                                                                    #
#    energy_eval                                                     #
#      Evaluate an E(V) fit at specified volumes.                    #
#                                                                    #
#    pressure_fit                                                    #
#      Perform a pressure vs. volume fit                             #
#                                                                    #
#    pressure eval                                                   #
#      Evaluate a P(V) fit at specified volumes.                     #
#                                                                    #                                        
#====================================================================#


from numpy import exp
from developer import unavailable
try:
    from scipy.optimize import fmin
except ImportError:
    fmin = unavailable('scipy.optimize','fmin')
#end try


f = None

#minimizers
least_squares = lambda p,x,y,f: ((f(p,x)-y)**2).sum()

#fits
fits = set(['murnaghan','birch','vinet','spanu','casula'])
fstr = ''
for ft in fits:
    fstr += ft+' '
#end for
em  = lambda p,V: p[0] + p[2]/p[3]*V*((p[1]/V)**p[3]/(p[3]-1)+1)-p[1]*p[2]/(p[3]-1)
eb  = lambda p,V: p[0] + 9*p[1]*p[2]/16*((p[1]/V)**(2./3)-1)**2*( 2 + (p[3]-4)*((p[1]/V)**(2./3)-1) )
ev  = lambda p,V: p[0] + 2*p[1]*p[2]/(p[3]-1)**2*( 2 - (2+3*(p[3]-1)*((V/p[1])**(1./3)-1))*exp(-1.5*(p[3]-1)*((V/p[1])**(1./3)-1)) ) 
es  = lambda p,V: p[0] + exp(p[1]*V)-p[2]/V**4
emo = lambda p,V: p[0]+p[2]*( exp(-2*p[3]*(V-p[1])) - 2*exp(-p[3]*(V-p[1])) )
efit = dict(murnaghan=em,birch=eb,vinet=ev,spanu=es,morse=emo)
pm  = lambda p,V: p[1]/p[2]*((p[0]/V)**p[2]-1)
pb  = lambda p,V: 1.5*p[1]*(p[0]/V)**(5./3)*((p[0]/V)**(2./3)-1)*(1.+.75*(p[2]-1)*((p[0]/V)**(2./3)-1))
pv  = lambda p,V: 3.*p[1]*(1.-(V/p[0])**(1./3))*(p[0]/V)**(2./3)*exp(1.5*(p[2]-1)*(1.-(V/p[0])**(1./3)))
pfit = {'murnaghan':pm,'birch':pb,'vinet':pv}

def energy_fit(ftype,E,V=None,a=None,minimizer=least_squares,E0=None,V0=None,B0=None,B0p=None,goodness=0):
    if ftype not in fits:
        print 'Error in energy_fit'
        print '  invalid fit requested:',ftype
        print '  valid types are:'
        print '   ',fstr
    #end if
    if a != None:
        V = a**3
    #end if
    f = efit[ftype]

    if E0==None:
        E0 = E[-1]
    if V0==None:
        V0 = V[E.argmin()]
    if B0==None:
        B0 = 7.
    if B0p==None:
        B0p = 3.5

    p0 = E0,V0,B0,B0p

    p = fmin(minimizer,p0,args=(V,E,f),maxiter=10000,maxfun=10000,disp=0)
    fopt=minimizer(p,V,E,f)

    E0,V0,B0,B0p = p[0],p[1],p[2],p[3]

    if goodness:
        return E0,V0,B0,B0p,fopt
    else:
        return E0,V0,B0,B0p
    #end if
#end def energy_fit


def energy_eval(ftype,E0,V0,B0,B0p,V=None,a=None):
    if ftype not in fits:
        print 'Error in energy_eval'
        print '  invalid fit requested:',ftype
        print '  valid types are:'
        print '   ',fstr
    #end if
    if a != None:
        V = a**3
    #end if
    
    f = efit[ftype]
    p = E0,V0,B0,B0p
    E = f(p,V)

    return E
#end def energy_eval


def pressure_fit(ftype,P,V=None,a=None,minimizer=least_squares):
    if ftype not in fits:
        print 'Error in pressure_fit'
        print '  invalid fit requested:',ftype
        print '  valid types are:'
        print '   ',fstr
    #end if
    if a != None:
        V = a**3
    #end if
    f = pfit[ftype]

    V0 = V[E.argmin()]
    B0 = 7.
    B0p = 3.5

    p0 = V0,B0,B0p

    p = fmin(minimizer,p0,args=(V,P,f),maxiter=10000,maxfun=10000,disp=0)

    V0,B0,B0p = p[0],p[1],p[2]

    return V0,B0,B0p
#end def pressure_fit


def pressure_eval(ftype,V0,B0,B0p,V=None,a=None):
    if ftype not in fits:
        print 'Error in pressure_eval'
        print '  invalid fit requested:',ftype
        print '  valid types are:'
        print '   ',fstr
    #end if
    if a != None:
        V = a**3
    #end if
    
    f = pfit[ftype]
    p = V0,B0,B0p
    P = f(p,V)

    return P
#end def pressure_eval

