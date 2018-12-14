##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  sqd_analyzer.py                                                   #
#    Supports analysis of SQD output data.                           #
#                                                                    #
#  Content summary:                                                  #
#    SqdAnalyzer                                                     #
#      SimulationAnalyzer class for the SQD code.                    #
#                                                                    #
#====================================================================#


import os
from numpy import loadtxt,empty,array,abs
from generic import obj
from simulation import SimulationAnalyzer,Simulation
from sqd_input import SqdInput
from plotting import *


class SqdAnalyzer(SimulationAnalyzer):
    channels = {0:'s',1:'p',2:'d',3:'f',4:'g',5:'h'}

    def __init__(self,arg0,qcut=1e-3,analyze=False):
        if isinstance(arg0,Simulation):
            sim = arg0
            infile = os.path.join(sim.resdir,sim.infile)
        elif isinstance(arg0,str):
            infile = arg0
        else:
            self.error('must specify location of input file with a string\n  you provided '+str(arg0))
        #end if
        self.info = obj(
            infile = infile,
            qcut   = qcut
            )
        if analyze:
            self.analyze()
        #end if
    #end def __init__


    def analyze(self):
        infile = self.info.infile
        qcut   = self.info.qcut
        if not os.path.exists(infile):
            self.error('input file {0} does not exist'.format(infile))
        #end if
        basepath,infile_name = os.path.split(infile)
        input = SqdInput(infile)
        self.info.set(
            path        = basepath,
            infile_name = infile_name,
            input       = input,
            atom        = input.simulation.atom.name,
            Z           = input.simulation.atom.hamiltonian.z
            )
        e = array([])
        self.set(
            exchange    = obj(),
            hartree     = obj(),
            r           = e.copy(),
            orbitals    = obj(),
            density     = e.copy()
            )

        outfiles = self.info.input.get_output_info(list=False)
        path = self.info.path

        log_file      = os.path.join(path,outfiles.log)
        exchange_file = os.path.join(path,outfiles.exchange)
        hartree_file  = os.path.join(path,outfiles.hartree)
        orbital_file  = os.path.join(path,outfiles.orb)

        #read log file
        lgcont = open(log_file).read()
        for line in lgcont.splitlines():
            if line.count('=')==1:
                name,value = line.split('=',1)
                name = name.strip()
                value = float(value.strip())
                self[name]=value
            #end if
        #end for
        self.E_potential = self.V_External+self.V_Hartree+self.V_Exchange
        self.E_kinetic = self.E_tot - self.E_potential

        #read exchange file
        exchange = self.exchange
        excont = open(exchange_file).read()
        for line in excont.splitlines()[1:]:
            ls = line.strip()
            if not ls.startswith('#'):
                orb1,orb2,value = ls.split()
                value = float(value)
                exchange[(orb1,orb2)] = value
            #end if
        #end for

        #read hartree file
        hartree = self.hartree
        hacont = open(hartree_file).read()
        for line in hacont.splitlines()[1:]:
            ls = line.strip()
            if not ls.startswith('#'):
                orb1,orb2,value = ls.split()
                value = float(value)
                hartree[(orb1,orb2)] = value
            #end if
        #end for

        #read orbital file
        orbitals = self.orbitals
        orbcont = open(orbital_file).read()
        for line in orbcont.splitlines():
            ls = line.strip()
            if ls.startswith('# n='):
                n,l,value = ls.replace('#','').replace('n=','').replace('l=','').split()
                n = int(n)
                l = int(l)
                value = float(value)
                if not (n,l) in orbitals:
                    orbitals[(n,l)] = obj(
                        degeneracy = 0,
                        E          = None,
                        psi        = None
                        )
                #end if
                orb = orbitals[(n,l)]
                orb.E = value
                orb.degeneracy += 1
            #end if
        #end for
        states = list(orbitals.keys())
        states.sort()
        data = loadtxt(orbital_file)
        self.r = data[:,0]
        data = data[:,1:]
        if len(states)!=data.shape[1]:
            self.error('number of unique eigenstates does not match the number of orbitals in orbital file: '+outfiles.orb)
        #end if
        for i in range(len(states)):
            orbitals[states[i]].psi = data[:,i]
        #end for

        #calculate the density
        density = 0*self.r
        for orbital in self.orbitals:
            density += orbital.degeneracy*orbital.psi**2
        #end for
        self.density = density
        self.density *= self.info.Z/self.sum()

        #find rcut such that int(density,rcut,infinity) = qcut
        self.set_rcut(self.info.qcut)

        #find the first moment r1=int(r*density)/int(density)
        self.r1 = self.moment(n=1)
    #end def analyze


    def sum(self,rmin=None,rmax=None):
        r = self.r
        d = self.density
        imin = 0
        imax = len(r)-1
        sdmin = r[0]*d[0]/2
        if rmin!=None and rmin>r[imin]/2:
            imin = abs(r-rmin).argmin()
            sdmin = 0
        #end if
        if rmax!=None and rmax<(r[imax]+r[imax-1])/2:
            imax = abs(r-rmax).argmin()
        #end if
        ri = r[imin:imax+1]
        di = d[imin:imax+1]
        sd = ((ri[1:]-ri[:-1])*(di[1:]+di[:-1])/2).sum()
        return sdmin + sd        
    #end def sum


    def cumsum(self):
        r = self.r
        d = self.density
        csd = 0*d
        csd[0]  = r[0]*d[0]/2
        csd[1:] = ((r[1:]-r[:-1])*(d[1:]+d[:-1])/2).cumsum()
        return csd
    #end def cumsum


    def moment(self,n):
        r = self.r
        d = r**n*self.density
        m = r[0]*d[0]/2 + ((r[1:]-r[:-1])*(d[1:]+d[:-1])/2).sum()
        return m/self.info.Z
    #end def moment


    def furthest_peak(self):
        r = self.r
        d = self.density
        rpeak = r[-1]
        for i in range(len(d)-2,-1,-1):
            if d[i]<d[i+1]:
                rpeak = r[i+1]
                break
            #end if
        #end for
        return rpeak
    #end def furthest_peak


    def find_rcut(self,qcut):
        r = self.r
        qcum = self.cumsum()
        qtot = self.info.Z
        icut = abs(qtot-qcum-qcut).argmin()        
        return r[icut]
    #end def find_rcut


    def set_rcut(self,qcut=1e-3):
        self.info.qcut = qcut
        self.rcut = self.find_rcut(qcut)
    #end def set_rcut


    def plot_orbitals(self,rmax=None):
        if rmax is None:
            rmax = self.rcut
        #end if
        colors = {0:'k',1:'r',2:'b',3:'g',4:'m',5:'c'}
        r = self.r
        for (n,l),orbital in self.orbitals.iteritems():
            psi = orbital.psi
            plot(r,psi,colors[l])
            imax = psi.argmax()
            text(r[imax],psi[imax],'{0}{1}({2})'.format(n,self.channels[l],orbital.degeneracy))
        #end for
        xlim([r.min(),rmax])
        title('Hartree-Fock orbitals for {0}'.format(self.info.atom))
        xlabel('r (bohr)')
        ylabel(r'$r\Psi$')
    #end def plot_orbitals
               
            
    def plot_density(self,rmax=None):
        if rmax is None:
            rmax = self.rcut
        #end if
        plot(self.r,self.density,'k')
        xlim([self.r.min(),rmax])
        title('Hartree-Fock density for {0}'.format(self.info.atom))
        xlabel('r (bohr)')
        ylabel(r'$r^2\rho$')
    #end def plot_density
#end class SqdAnalyzer
