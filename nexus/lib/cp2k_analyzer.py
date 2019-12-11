##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  cp2k_analyzer.py                                                 #
#    Supports data analysis for CP2K output.  Can handle log file   #
#                                                                    #
#  Content summary:                                                  #
#    Cp2kAnalyzer                                                   #
#      SimulationAnalyzer class for CP2K.                           #
#      Reads log output and converts data to numeric form.           #
#                                                                    #
#====================================================================#


import os
from numpy import array,fromstring,sqrt,dot,max,equal,zeros,min,where,empty
from generic import obj
from unit_converter import convert
from simulation import SimulationAnalyzer,Simulation
#from cp2k_input import Cp2kInput
#from cp2k_data_reader import read_qexml
from fileio import TextFile
from debug import *
import code
import pdb


class Cp2kAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0=None,outfile_name=None,analyze=False):
        if isinstance(arg0,Simulation):
            sim = arg0
            path = sim.locdir
            outfile_name= sim.outfile
        elif arg0 is not None:
            path = arg0
            if not os.path.exists(path):
                self.error('path to data does not exist\npath provided: {}'.format(path))
            #end if
            if os.path.isfile(path):
                filepath = path
                path,filename = os.path.split(filepath)
                outfile_name = filename
#                if filename.endswith('.in'):
#                    infile_name = filename
#                elif filename.endswith('.out'):
#`               outfile_name = filename
#                else:
#                    self.error('could not determine whether file is input or output\nfile provided: {}'.format(filepath))
                #end if
            #end if
#            if outfile_name is None:
#                outfile_name = infile_name.rsplit('.',1)[0]+'.out'
            #end if
#        else:
#            return
        #end if

#        self.infile_name  = infile_name
        self.outfile_name = outfile_name

        self.path = path
        self.abspath = os.path.abspath(path)
        self.analyze()
    #end def __init__

    
    def analyze(self):
        path = self.path
#        infile_name = self.infile_name
        outfile_name = self.outfile_name

        nx=0

        outfile = os.path.join(path,outfile_name)
        # perform MD analysis
        f = TextFile(outfile)
        n = 0
        md_res = []
        while f.seek('STEP NUMBER',1)!=-1:
            f.seek('TIME [fs]',1)
            t = float(f.readtokens()[-1])/1000.
            f.seek('POTENTIAL ENERGY[hartree]',1)
            PE = float(f.readtokens()[-2])
            f.seek('KINETIC ENERGY [hartree]',1)
            KE = float(f.readtokens()[-2])
            E=PE+KE
            f.seek('TEMPERATURE [K]',1)
            T = float(f.readtokens()[-2])
            f.seek('PRESSURE [bar]',1)
            P = float(f.readtokens()[-2])/10000.0
#for NPT we need Volume and Cell Parameters
            f.seek('VOLUME',1)         
            V = float(f.readtokens()[-2])*0.148184584700
            # stress matrix S, note S00 is S11 in normal terms!
#            S11, S12, S13 = f.readtokens()[-3:]
#            S21, S22, S23 = f.readtokens()[-3:]
#            S31, S32, S33 = f.readtokens()[-3:]
            md_res.append((E,P,t,KE,PE,T,V))
            n+=1
        #end while
        md_res = array(md_res,dtype=float).T
        quantities = ('total_energyH','pressureGPa','timeps','kinetic_energyH',
                      'potential_energyH','temperatureK','volumeA3')
        md = obj()
        for i,q in enumerate(quantities):
            md[q] = md_res[i]
        #end for
        self.md_data = md
        self.md_stats = self.md_statistics()
        return
        #end if
    def md_statistics(self,equil=None,autocorr=None):
        import numpy as np
        from numerics import simstats,simplestats
        mds = obj()
        for q,v in self.md_data.items():
            if equil is not None:
                v = v[equil:]
            else:
                equil=0
            if autocorr is None:
                mean,var,error,kappa = simstats(v)
                mds[q] = mean,error,kappa,equil
            else:
                nv = len(v)
                nb = int(np.floor(float(nv)/autocorr))
                nexclude = nv-nb*autocorr
                v = v[nexclude:]
                v.shape = nb,autocorr
                mean,error = simplestats(v.mean(axis=1))
                mds[q] = mean,error,autocorr,nexclude
        return mds


    def md_plots(self,filename,equil,show=True,removeskips=False,startfromone=False):
        import numpy as np
        import matplotlib.pyplot as plt
        md = self.md_data
        plt.figure()
        plt.subplots_adjust(left=0.2)
        plt.subplots_adjust(bottom=0.15)
        plt.suptitle(filename,fontsize=10)
        plt.axis([None, None, -.1, .1])
        plt.subplot(4,1,1)
        plt.ylim(-0.4,0.4)
        e=equil
        dt=md.timeps[1]-md.timeps[0]
        dtl=md.timeps[-1]-md.timeps[-2]
        lmd=len(md.timeps)
        print(" dt (ps)= {0} {1} length= {2}".format(dt,dtl,lmd))
        if(startfromone):
            i=0
            for ii in md.timeps:
                md.timeps[i]=i*dt
                i+=1
        if(removeskips):
            for i in range(2,lmd):
                if (md.timeps[i]-md.timeps[i-1] > 2*dt ):
                    md.timeps[i]=2.0*md.timeps[i-1]-md.timeps[i-2]
        plt.plot(md.timeps[e:],md.total_energyH[e:]-np.mean(md.total_energyH[e:]),label='Etot(H)')
        plt.plot(md.timeps[e:],md.kinetic_energyH[e:]-np.mean(md.kinetic_energyH[e:]),label='Ekin(H)')
        plt.plot(md.timeps[e:],md.potential_energyH[e:]-np.mean(md.potential_energyH[e:]),label='Epot(H)')
        plt.ylabel('E (H)')
        plt.legend()
        plt.axis([None, None, None,None])
        plt.subplot(4,1,2)
        plt.plot(md.timeps[e:],md.temperatureK[e:])
        plt.ylabel('T (K)')
        plt.subplot(4,1,3)
        plt.plot(md.timeps[e:],md.pressureGPa[e:])
        plt.ylabel('P (GPa)')
        plt.xlabel('time (ps)')
        plt.subplot(4,1,4)
        plt.plot(md.timeps[e:],md.volumeA3[e:])
        plt.ylabel('V (A3)')
        plt.xlabel('time (ps)')
        if show:
            plt.show()
        return 
