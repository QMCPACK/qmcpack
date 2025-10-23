#! /usr/bin/env python3

import os
import sys
from subprocess import Popen,PIPE
from numpy import array,sqrt,polyfit,polyval,linspace
from developer import obj
import matplotlib.pyplot as plt

params = {'legend.fontsize':14,'figure.facecolor':'white','figure.subplot.hspace':0.6,
          'axes.labelsize':16,'xtick.labelsize':14,'ytick.labelsize':14}
plt.rcParams.update(params)


#files = sys.argv[1:]

files = 'O.q0.dmc.in.xml  O.q1.dmc.in.xml'.split()


res = obj()
for file in files:
    if not os.path.exists(file):
        print('file '+file+' does not exist')
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


    print()
    print()
    print('Total energy (q=0, eV)')
    for i in range(len(ts)):
        print('  {0:<6.4f}  {1:<6.4f} +/- {2:<6.4f}'.format(ts[i],q0.energies[i],q0.errs[i]))
    #end for
    print('----------------------------------')
    print('  {0:<6.4f}  {1:<6.4f}'.format(0.00,q0fit[0]))

    print()
    print()
    print('Total energy (q=1, eV)')
    for i in range(len(ts)):
        print('  {0:<6.4f}  {1:<6.4f} +/- {2:<6.4f}'.format(ts[i],q1.energies[i],q1.errs[i]))
    #end for
    print('----------------------------------')
    print('  {0:<6.4f}  {1:<6.4f}'.format(0.00,q1fit[0]))

    print()
    print()
    print('1st ionization potential (eV)')
    for i in range(len(ts)):
        print('  {0:<6.4f}  {1:<6.4f} +/- {2:<6.4f}'.format(ts[i],ip[i],iperr[i]))
    #end for
    print('----------------------------------')
    print('  {0:<6.4f}  {1:<6.4f}'.format(0.00,ipfit[0]))
    
    print()
    print()

    fig_1 = plt.figure(num="q=0 Total Energy vs. Time Step")
    ax_1 = fig_1.add_subplot()
    ax_1.errorbar(q0.timesteps,q0.energies,q0.errs,fmt='b.',label='q0')
    ax_1.plot(tfit,q0fit,'k--',lw=2)
    ax_1.title('Oxygen q=0 total energy vs. timestep')
    ax_1.ylabel('Total energy (eV)')
    ax_1.xlabel('timestep (1/Ha)')
    ax_1.legend()

    fig_2 = plt.figure(num="q=1 Total Energy vs. Time Step")
    ax_2 = fig_2.add_subplot()
    ax_2.errorbar(q1.timesteps,q1.energies,q1.errs,fmt='b.',label='q1')
    ax_2.plot(tfit,q1fit,'k--',lw=2)
    ax_2.title('Oxygen q=1 total energy vs. timestep')
    ax_2.ylabel('Total energy (eV)')
    ax_2.xlabel('timestep (1/Ha)')
    ax_2.legend()

    fig_3 = plt.figure(num="1st Ionization Potential vs. Time Step")
    ax_3 = fig_3.add_subplot()
    ax_3.errorbar(ts,ip,iperr,fmt='r.',label='1st IP')
    ax_3.plot(tfit,ipfit,'g--',lw=2)
    ax_3.title('Oxygen 1st IP vs. timestep')
    ax_3.ylabel('IP (eV)')
    ax_3.xlabel('timestep (1/Ha)')
    ax_3.legend()

    plt.show()
#end if
