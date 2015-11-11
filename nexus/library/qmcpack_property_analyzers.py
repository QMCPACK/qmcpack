##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  qmcpack_property_analyzers.py                                     #
#    Supports analysis of wavefunctions used/made by QMCPACK.        #
#                                                                    #
#  Content summary:                                                  #
#    PropertyAnalyzer                                                #
#      Base class for specific property analyzers.                   #
#      Currently, only the wavefunction falls in this category.      #
#                                                                    #
#    WavefunctionAnalyzer                                            #
#      Class to read and plot wavefunction data.                     #
#      Currently only interfaces with b-spline Jastrows.             #
#                                                                    #
#    RadialJastrow                                                   #
#      Represents a spherically symmetric (radial) jastrow           #
#        correlation (u) function.                                   #
#      Base class for specific radial jastrow types.                 #
#      See Jastrow1B and Jastrow2B                                   #
#                                                                    #
#    Jastrow1B                                                       #
#      Represents a 1-body Jastrow.                                  #
#                                                                    #
#    Jastrow2B                                                       #
#      Represents a 2-body Jastrow.                                  #
#                                                                    #
#    Bspline                                                         #
#      Represents a cubic bspline curve, supports evaluation.        #
#                                                                    #
#====================================================================#



import os
from numpy import loadtxt,array,ones,dot,floor,zeros
from qmcpack_input import QmcpackInput
from qmcpack_analyzer_base import QAobject,QAanalyzer
from developer import unavailable
from debug import *
try:
    from matplotlib.pyplot import plot,show,figure,xlabel,ylabel,title,legend
except (ImportError,RuntimeError):
    plot,show,figure,xlabel,ylabel,title,legend = unavailable('matplotlib.pyplot','plot','show','figure','xlabel','ylabel','title','legend')
#end try



class Bspline(QAobject):
    A = array([-1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
                3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
               -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
                1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0 ])
    dA =array([0.0, -0.5,  1.0, -0.5,
               0.0,  1.5, -2.0,  0.0,
               0.0, -1.5,  1.0,  0.5,
               0.0,  0.5,  0.0,  0.0 ])
    d2A=array([0.0, 0.0, -1.0,  1.0,
               0.0, 0.0,  3.0, -2.0,
               0.0, 0.0, -3.0,  1.0,
               0.0, 0.0,  1.0,  0.0 ])
    A.shape   = 4,4
    dA.shape  = 4,4
    d2A.shape = 4,4                 

    def __init__(self,params,cusp,rcut):
        p = array(params)
        cusp = float(cusp)
        rcut = float(rcut)
        c = zeros((1,len(p)+4))
        nintervals = len(p) + 1
        dr = rcut/nintervals
        odr = 1./dr

        c[0,0] = p[1] - 2.*dr*cusp
        c[0,1] = p[0]
        c[0,2] = p[1]
        for i in range(2,len(p)):
            c[0,i+1] = p[i]
        #end for
           
        self.p      = p      
        self.rcut   = rcut   
        self.cusp   = cusp   
        self.c      = c      
        self.nintervals = nintervals 
        self.dr     = dr     
        self.odr    = odr
        self.default_range = 0.,rcut
    #end def __init__

    def evaluate(self,r):
        tp=zeros((4,1))
        ni  = self.nintervals
        odr = self.odr
        c   = self.c
        A   = self.A
        dA  = self.dA
        d2A = self.d2A
        v   = zeros(r.shape)
        dv  = zeros(r.shape)
        d2v = zeros(r.shape)
        for p in range(len(r)):               
            ri = r[p]*odr
            i = floor(ri)
            if i<ni:
                t = ri - i
                tp[0,0] = t*t*t
                tp[1,0] = t*t
                tp[2,0] = t
                tp[3,0] = 1.
                v[p]   = dot(c[:,i:i+4],  dot(A,tp))
                dv[p]  = dot(c[:,i:i+4], dot(dA,tp))
                d2v[p] = dot(c[:,i:i+4],dot(d2A,tp))
            else:
                v[p] = 0.
                dv[p] = 0.
                d2v[p] = 0.
            #end if
        #end for
        dv*=odr
        d2v*=odr*odr
        return v,dv,d2v
    #end def evaluate
#end class Bspline


from numpy import linspace
class RadialJastrow(QAobject):
    def __init__(self,ftype,coeff,cusp,rcut):
        self.coeff = coeff
        self.cusp  = cusp
        if ftype.lower()=='bspline':
            self.function = Bspline(coeff,cusp,rcut)
            self.rcut = rcut
        #end if
    #end def __init__

    def evaluate(self,r):
        return self.function.evaluate(r)
    #end def evaluate

    def interpolate(self,r1=None,r2=None,n=200):
        if r1 is None:
            r1,r2 = self.function.default_range
        #end if
        r = linspace(r1,r2,n)
        d0,d1,d2 = self.function.evaluate(r)
        return r,d0,d1,d2
    #end def interpolate

    def plot(self,r1=None,r2=None,color='b',ptype=plot):
        r,d0,d1,d2 = self.interpolate(r1,r2)
        c = color
        ptype(r,d0,ls='-' ,c=c,label='value')
        ptype(r,d1,ls='-.',c=c,label='first derivative')
        ptype(r,d2,ls=':' ,c=c,label='second derivative')
        return r,d0,d1,d2
    #end def plot
#end class RadialJastrow

class Jastrow1B(RadialJastrow):
    def __init__(self,ftype,coeff,rcut):
        cusp = 0.0
        RadialJastrow.__init__(self,ftype,coeff,cusp,rcut)
    #end def __init__
#end class Jastrow1B

class Jastrow2B(RadialJastrow):
    def __init__(self,ftype,coeff,species1,species2,rcut):
        if species1==species2:
            cusp = -1./4
        else:
            cusp = -1./2
        #end if
        RadialJastrow.__init__(self,ftype,coeff,cusp,rcut)
    #end def __init__
#end class Jastrow2B


class PropertyAnalyzer(QAanalyzer):
    None
#end class PropertyAnalyzer


class WavefunctionAnalyzer(PropertyAnalyzer):

    jastrow_types = ['J1','J2','J3']

    def __init__(self,arg0=None,load_jastrow=False,nindent=0):
        QAanalyzer.__init__(self,nindent=nindent)
        self.info.load_jastrow = load_jastrow

        if isinstance(arg0,str):
            self.info.filepath = arg0
        else:
            self.info.wfn_xml = arg0
        #end if

        self.info.fail = False
    #end def __init__
            
            
    def load_data_local(self):
        info = self.info
        if info.load_jastrow:
            self.load_jastrow_data()
        elif 'filepath' in info:
            try:
                qxml = QmcpackInput(info.filepath)
                wavefunction = qxml.get('wavefunction')
                wavefunction = wavefunction.get_single('psi0')
                info.wfn_xml = wavefunction
            except:
                info.wfn_xml = None
                info.fail = True
            #end try
        #end if
        if not info.load_jastrow and not info.fail:
            info.wfn_xml.pluralize()
        #end if
    #end def load_data_local


    def analyze_local(self):
        structure = QAanalyzer.run_info.system.structure

        jnames = {'One-Body':'J1','Two-Body':'J2','Three-Body':'J3'}
        jastrows = QAobject()
        for jt,jn in jnames.iteritems():
            jastrows[jn] = QAobject()
        #end for
        del jastrows.J3

        if len(structure.axes)==3:
            rcut_cell = structure.rwigner()
        else:
            rcut_cell = 10
        #end if

        try:
            J1,J2,J3 = self.info.wfn_xml.get(['J1','J2','J3'])
            if J1!=None:
                jname = 'J1'
                func = J1.function.lower()
                if func=='bspline':
                    for jn,corr in J1.correlations.iteritems():
                        if 'rcut' in corr:
                            rcut = corr.rcut
                        else:
                            rcut = rcut_cell
                        #end if
                        coeff = corr.coefficients.coeff
                        jastrows[jname][jn] = Jastrow1B(func,coeff,rcut)
                    #end for
                #end if
            #end if
            if J2!=None:
                jname = 'J2'
                func = J2.function.lower()
                if func=='bspline':
                    for jn,corr in J2.correlations.iteritems():
                        if 'rcut' in corr:
                            rcut = corr.rcut
                        else:
                            rcut = rcut_cell
                        #end if
                        s1 = corr.speciesa
                        s2 = corr.speciesb
                        coeff = corr.coefficients.coeff
                        jastrows[jname][jn] = Jastrow2B(func,coeff,s1,s2,rcut)
                    #end for
                #end if
            #end if
        except:
            self.warn('Jastrow read failed, some data will not be available')
            self.info.fail = True
        #end try
        self._transfer_from(jastrows)
    #end def analyze_local


    def load_jastrow_data(self):
        ext = '.g'+str(self.batch_index).zfill(3)+'.dat'
        for jt in self.jastrow_types:
            for jn,je in self[jt]._iteritems():
                J = self[jt][jn]
                data = loadtxt(os.path.join(self.sourcepath,jt+'.'+jn+ext))
                J.r  = data[:,0]
                J.d0 = data[:,1]
                J.d1 = data[:,2]
                J.d2 = data[:,3]
            #end for
        #end for
    #end def load_jastrow_data


    def plot_jastrow_data(self):
        for jt in self.jastrow_types:
            if len(self[jt])!=0:
                for jn,je in self[jt]._iteritems():
                    r = je.r
                    bs = Bspline(je.coefficients,r.max())
                    d0,d1,d2 = bs.evaluate(r)

                    figure()
                    plot(je.r,je.d0,'b-' ,label='value')
                    plot(je.r,je.d1,'b-.',label='first derivative')
                    plot(je.r,je.d2,'b:' ,label='second derivative')
                    plot(r,d0,'r-' )
                    plot(r,d1,'r-.')
                    plot(r,d2,'r:' )
                    plot(je.r,0*je.r,'k-')
                    xlabel('r (Bohr)')
                    title(jt+' '+jn)
                    legend()
                #end for
            #end if
        #end for
        show()
    #end def plot_jastrow_data


    def plot_jastrows(self,ptype=plot):
        for name,value in self.iteritems():
            if name in self.jastrow_types:
                for label,jastrow in value.iteritems():
                    jtype = jastrow.__class__.__name__
                    figure()
                    jastrow.plot(ptype=ptype)
                    xlabel('r (Bohr)')
                    ylabel(jtype+' ('+label+')')
                    title(jtype+' for '+label)
                    legend()
                #end for
            #end if
        #end for
        show()
    #end def plot_jastrows

#end class WavefunctionAnalyzer







