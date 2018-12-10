#! /usr/bin/env python

from matplotlib.pyplot import figure,plot,xlabel,ylabel,title,show,ylim,legend,xlim,rcParams,savefig,xticks,yticks,errorbar,tight_layout

params = {'legend.fontsize':14,'figure.facecolor':'white','figure.subplot.hspace':0.,
          'axes.labelsize':16,'xtick.labelsize':14,'ytick.labelsize':14}
rcParams.update(params)

from numpy import array
from numpy import polyfit,polyval,linspace

# UNR 1993 DMC Emix data (extracted from figure 7 using http://arohatgi.info/WebPlotDigitizer/app/)
# Reference: C. J. Umrigar et al., J. Chem. Phys. 99 2865 (1993)
# DOI:       http://dx.doi.org/10.1063/1.465195
d = array('''
0.250  -14.993938775510204
0.250  -14.993862244897958
0.225  -14.99262244897959
0.225  -14.9925
0.200  -14.99202551020408
0.200  -14.99190306122449
0.175  -14.991244897959183
0.175  -14.99112244897959
0.150  -14.990341836734693
0.150  -14.9902193877551
0.125  -14.989867346938775
0.125  -14.989729591836733
0.100  -14.989377551020407
0.100  -14.989239795918365
0.075  -14.989224489795918
0.075  -14.989086734693876
0.050  -14.98888775510204
0.050  -14.988688775510203
0.025  -14.988780612244897
0.025  -14.988596938775508
0.010  -14.9894693877551
0.010  -14.989117346938775'''.split(),dtype=float)

d.shape = len(d)/2,2

tau    = d[::2,0]         # UNR timesteps
E_unr  = d[::2,1]         # UNR DMC energy means
Ee_unr = d[1::2,1]-E_unr  # UNR DMC energy errorbars

tfit = linspace(0,0.3,300)
p1 = polyfit(tau,E_unr,2)
p2 = polyfit(tau[:-1],E_unr[:-1],2)
p_unr = tuple((array(p1)+array(p2))/2)
Efit = polyval(p_unr,tfit)

# print out UNR Fig 7 means and error bars
print
print 'Extracted Emix data from UNR Fig. 7'
for i in range(len(tau)):
    print '{0:6.4f}  {1:12.6f} +/- {2:12.6f}'.format(tau[i],E_unr[i],Ee_unr[i])
#end for


# recreate UNR Figure 7 (mixed estimator only)
figure(tight_layout=True)
plot(tfit,Efit,'k-')
errorbar(tau,E_unr,Ee_unr,fmt='bs',capsize=5,markerfacecolor='white')
xlim([0,0.25])
xticks([0,0.05,0.1,0.15,0.2,0.25])
ylim([-14.996,-14.988])
yticks([-14.996,-14.994,-14.992,-14.990,-14.988])
xlabel('Time Step $\\tau$ (Hartree$^{-1}$)')
ylabel('Energy (Hartrees)')
title('Reconstruction of UNR Figure 7',fontsize=16)
savefig('UNR_Fig7_reconstructed.pdf')


# QMCPACK data, generated on 29 Sep 2017 with version 3.2.0
#oic5>pwd
#/home/j1k/projects/ecp/01_ctest_development/02_Li_STO
#oic5>qmca -e 40 -q e *timestep_study_very_long*/*scalar*
#
#09_timestep_study_very_long_v320/Li2.STO.textrap series  0 LocalEnergy = -14.948953 +/- 0.000674 
#09_timestep_study_very_long_v320/Li2.STO.textrap series  1 LocalEnergy = -15.027879 +/- 0.000462 
#09_timestep_study_very_long_v320/Li2.STO.textrap series  2 LocalEnergy = -15.015879 +/- 0.000119 
#09_timestep_study_very_long_v320/Li2.STO.textrap series  3 LocalEnergy = -15.008000 +/- 0.000077 
#09_timestep_study_very_long_v320/Li2.STO.textrap series  4 LocalEnergy = -15.002293 +/- 0.000069 
#09_timestep_study_very_long_v320/Li2.STO.textrap series  5 LocalEnergy = -14.997681 +/- 0.000069 
#09_timestep_study_very_long_v320/Li2.STO.textrap series  6 LocalEnergy = -14.994168 +/- 0.000069 
#09_timestep_study_very_long_v320/Li2.STO.textrap series  7 LocalEnergy = -14.991667 +/- 0.000068 
#09_timestep_study_very_long_v320/Li2.STO.textrap series  8 LocalEnergy = -14.989910 +/- 0.000070 
#09_timestep_study_very_long_v320/Li2.STO.textrap series  9 LocalEnergy = -14.989020 +/- 0.000084 
#09_timestep_study_very_long_v320/Li2.STO.textrap series 10 LocalEnergy = -14.988972 +/- 0.000080 
#09_timestep_study_very_long_v320/Li2.STO.textrap series 11 LocalEnergy = -14.989505 +/- 0.000078 
#09_timestep_study_very_long_v320/Li2.STO.textrap series 12 LocalEnergy = -14.989772 +/- 0.000067 


# convert qmcpack data from qmca into numpy array
qd = array('''
0.250   -15.027879  0.000462
0.225   -15.015879  0.000119
0.200   -15.008000  0.000077
0.175   -15.002293  0.000069
0.150   -14.997681  0.000069
0.125   -14.994168  0.000069
0.100   -14.991667  0.000068
0.075   -14.989910  0.000070
0.050   -14.989020  0.000084
0.025   -14.988972  0.000080
0.010   -14.989505  0.000078
0.005   -14.989772  0.000067
'''.split(),dtype=float)
qd.shape = len(qd)/3,3
 
# extract timestep, mean and errorbar
qtau    = qd[:,0]
E_v320  = qd[:,1]
Ee_v320 = qd[:,2]

ifit = qtau<0.16
p_v320     = polyfit(qtau[ifit],E_v320[ifit],2)
Efit_v320  = polyval(p_v320,tfit)


# Plot QMCPACK v3.2.0 data alongside UNR results
figure(tight_layout=True)
plot(tfit,Efit,'k-')
errorbar(tau,E_unr,Ee_unr,fmt='ks',capsize=5,markerfacecolor='white',label='UNR')
plot(tfit,Efit_v320,'g-')
errorbar(qtau,E_v320,Ee_v320,fmt='g^',capsize=5,label='QMCPACK')
xlim([0,0.27])
xlabel('Time Step $\\tau$ (Ha$^{-1}$)')
ylabel('DMC Energy (Ha)')
legend(loc='lower left')
title('UNR vs. QMCPACK',fontsize=16)
savefig('UNR_vs_QMCPACK.pdf')


# Make a zoomed version of the plot
figure(tight_layout=True)
plot(tfit,Efit,'k-')
errorbar(tau,E_unr,Ee_unr,fmt='ks',capsize=5,markerfacecolor='white',label='UNR')
plot(tfit,Efit_v320,'g-')
errorbar(qtau,E_v320,Ee_v320,fmt='g^',capsize=5,label='QMCPACK')
xlim([0,0.15])
ylim([-14.996,-14.988])
xlabel('Time Step $\\tau$ (Ha$^{-1}$)')
ylabel('DMC Energy (Ha)')
legend(loc='lower left')
title('UNR vs. QMCPACK (zoomed)',fontsize=16)
savefig('UNR_vs_QMCPACK_zoomed.pdf')


show()
