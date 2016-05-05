#! /usr/bin/env python

import os
import sys
from subprocess import Popen,PIPE
from numpy import array,sqrt,polyfit,polyval,linspace
from generic import obj
from plotting import *

params = {'legend.fontsize':14,'figure.facecolor':'white','figure.subplot.hspace':0.6,
          'axes.labelsize':16,'xtick.labelsize':14,'ytick.labelsize':14}
rcParams.update(params)


#files = sys.argv[1:]

files = 'O.q0.dmc.in.xml  O.q1.dmc.in.xml'.split()


res = obj()
for file in files:
    if not os.path.exists(file):
        print 'file '+file+' does not exist'
        exit()
    #end if
    prefix    = None
    timesteps = []
    indmc = False
    lines = open(file,'r').read().splitlines()
    for line in lines:
        ls = line.strip()
        if ls.startswith('<project'):
            ls = ls.strip('</>').replace('"','')
            tokens = ls.split()
            for t in tokens:
                if 'id=' in t:
                    prefix = t.split('=')[1]
                #end if
            #end for
        elif 'method="dmc"' in ls:
            indmc = True
        elif '</qmc>' in ls:
            indmc = False
        elif '"timestep"' in ls and indmc:
            ptokens = ls.replace('>','>|').replace('<','|<').split('|') 
            timesteps.append(float(ptokens[2]))
        #end if
    #end for
    command = 'qmca -q e -u eV {0}*scalar*'.format(prefix)
    p = Popen(command,stdout=PIPE,shell=True)
    out,err = p.communicate()

    energies = []
    errs     = []
    for line in out.splitlines():
        tokens = line.split()
        if len(tokens)==8:
            energies.append(float(tokens[5]))
            errs.append(float(tokens[7]))
        #end if
    #end for
    ss = len(energies)-len(timesteps)
    energies = energies[ss:]
    errs = errs[ss:]
    r = obj(
        timesteps = array(timesteps),
        energies  = array(energies),
        errs      = array(errs)
        )
    if 'q0' in file:
        res.q0 = r
    elif 'q1' in file:
        res.q1 = r
    #end if
#end for

if len(res)==2:
    q0 = res.q0
    q1 = res.q1

    ts    = q0.timesteps
    ip    = q1.energies-q0.energies
    iperr = sqrt(q1.errs**2+q0.errs**2)

    tfit  = linspace(0,ts.max(),200)

    q0fit = polyval(polyfit(ts,q0.energies,1),tfit)
    q1fit = polyval(polyfit(ts,q1.energies,1),tfit)
    ipfit = polyval(polyfit(ts,ip,1),tfit)


    print
    print
    print 'Total energy (q=0, eV)'
    for i in range(len(ts)):
        print '  {0:<6.4f}  {1:<6.4f} +/- {2:<6.4f}'.format(ts[i],q0.energies[i],q0.errs[i])
    #end for
    print '----------------------------------'
    print '  {0:<6.4f}  {1:<6.4f}'.format(0.00,q0fit[0])

    print
    print
    print 'Total energy (q=1, eV)'
    for i in range(len(ts)):
        print '  {0:<6.4f}  {1:<6.4f} +/- {2:<6.4f}'.format(ts[i],q1.energies[i],q1.errs[i])
    #end for
    print '----------------------------------'
    print '  {0:<6.4f}  {1:<6.4f}'.format(0.00,q1fit[0])

    print
    print
    print '1st ionization potential (eV)'
    for i in range(len(ts)):
        print '  {0:<6.4f}  {1:<6.4f} +/- {2:<6.4f}'.format(ts[i],ip[i],iperr[i])
    #end for
    print '----------------------------------'
    print '  {0:<6.4f}  {1:<6.4f}'.format(0.00,ipfit[0])
    
    print
    print

    figure()
    errorbar(q0.timesteps,q0.energies,q0.errs,fmt='b.',label='q0')
    plot(tfit,q0fit,'k--',lw=2)
    title('Oxygen q=0 total energy vs. timestep')
    ylabel('Total energy (eV)')
    xlabel('timestep (1/Ha)')
    legend()

    figure()
    errorbar(q1.timesteps,q1.energies,q1.errs,fmt='b.',label='q1')
    plot(tfit,q1fit,'k--',lw=2)
    title('Oxygen q=1 total energy vs. timestep')
    ylabel('Total energy (eV)')
    xlabel('timestep (1/Ha)')
    legend()

    figure()
    errorbar(ts,ip,iperr,fmt='r.',label='1st IP')
    plot(tfit,ipfit,'g--',lw=2)
    title('Oxygen 1st IP vs. timestep')
    ylabel('IP (eV)')
    xlabel('timestep (1/Ha)')
    legend()

    show()
#end if
