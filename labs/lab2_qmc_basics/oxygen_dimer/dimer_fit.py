#! /usr/bin/env python

import sys
import os
from numpy import array,polyfit,polyval,linspace,sqrt
from qmcpack_analyzer import QmcpackAnalyzer
from unit_converter import convert
from plotting import *
from debug import *


args = sys.argv[1:]


prefix = 'qmc'
nfit   = 2
Eatom  = None
Eatom_err = 0
if len(args)>0:
    prefix = str(args[0])
#end if
if len(args)>1:
    nfit = int(args[1])
#end if
if len(args)>2:
    Eatom = float(args[2])
#end if
if len(args)>3:
    Eatom_err = float(args[3])
#end if
if Eatom is None:
    Elabel = 'total energy'
    Eatom = 0
else: 
    Elabel = 'binding energy'
#end if


dirs = os.listdir('./')

scales = []
energies = []
errors = []
for dir in dirs:
    if dir.startswith('scale_'):
        scale = float(dir.split('_')[1])
        infile = os.path.join(dir,prefix+'.in.xml')
        outfile = os.path.join(dir,prefix+'.s001.scalar.dat')
        if os.path.exists(outfile):
            qa = QmcpackAnalyzer(infile,analyze=True)
            le = qa.qmc[1].scalars.LocalEnergy
            scales.append(scale)
            E = convert(le.mean,'Ha','eV')
            Eerr = convert(le.error,'Ha','eV')
            E -= 2*Eatom
            Eerr = sqrt(Eerr**2+(2*Eatom_err)**2)
            energies.append(E)
            errors.append(Eerr)
        #end if
    #end if
#end for



se = 1.2074
ee = -5.165

seps = array(scales)*se
order = seps.argsort()

seps     = seps[order]
energies = array(energies)[order]
errors   = array(errors)[order]

p = polyfit(seps,energies,nfit)
sf = linspace(seps.min(),seps.max(),200)
ef = polyval(p,sf)

print
print 'DMC '+Elabel+' vs separation distance'
for i in range(len(seps)):
    print '  {0:6.4f}  {1:6.4f} +/- {2:6.4f}'.format(seps[i],energies[i],errors[i])
#end for

mloc = ef.argmin()
print
print 'DMC bond length, binding energy:',sf[mloc],ef[mloc]
print 'Exp bond length, binding energy:',se,ee


figure()
errorbar(seps,energies,errors,fmt='r.')
plot(sf,ef,'g-',lw=2,label='n={0} fit'.format(nfit))
xlabel('bond length (A)')
ylabel(Elabel+' (eV)')
title('DMC '+Elabel+' vs. bond length')
legend()

show()
