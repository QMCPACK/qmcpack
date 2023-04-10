#!/usr/bin/env python
# coding: utf-8

#By Raymond Clay, rclay@sandia.gov

#This is an independent test of the magnetization density.  What
#we do is take two spinors of the form [a b]e^i*k*r in a box.
#By setting a and b, we can set the desired spin state and calculate
#the expectation values.  

#For each X=(R,S), this routine will run through and compute the 
#spin integral going into the magnetization density expectation value.
#The return values for each spin are the check that QMCPACK will employ for 
#the unit test.  

#To verify that the formalism is correct, we perform a monte carlo integration 
#w.r.t. s1 and s2 for the magnetization density using the same determinant 
#trial function as stated previously.  Choosing a=cos(theta/2) b=e^(i*phi)sin(theta/2), 
#we find that the spin expectation value for each spin should be [0.5,0.5,1/sqrt(2)].  

import numpy as np
import math
import copy


#Unit cell specification.
L=np.zeros((3,3))
L[0]=np.array([5.10509515,-3.23993545,0.00000000])
L[1]=np.array([5.10509515, 3.23993545,0.00000000])
L[2]=np.array([-6.49690625,0.00000000,7.08268015])

#For conversion from cartesian to crystal coords, and for kpoints
Linv=np.linalg.inv(L)
G=Linv.transpose()

print("Volume = ",np.linalg.det(L))
print(L)

def convert_to_crystal_coords(r,L):
    Linv=np.linalg.inv(L)
    return np.dot(r,Linv)
#this takes the crystal coords and remaps them to [0,1)
def wrap_crystal_to_unit_cell(r):
    return r-np.floor(r)


#This is a test of the cell wrapping stuff.
x=np.array([5,0,0])
xc=convert_to_crystal_coords(x,L)
xwrap=wrap_crystal_to_unit_cell(xc)

print(x)
print(xc)
print(xwrap)


#This is the electron particleset declaration.

nelec=2
R=np.zeros((nelec,3))
spins=np.zeros(nelec)
#This was the original configuration used to generate the SPO
#reference values.  These will be fed into plane wave spatial orbitals. 
R[0]=np.array([0.0,0,0])
R[1]=np.array([5,0,0])

spins[0]=1.9
spins[1]=2.5410

print("nelec = ",nelec)
print("R = ",R)
print("spins = ",spins)


#This is in real coordinates.  
def eval_pw_spo(r,klist):
    phase=np.dot(klist,r)
    val=np.exp(phase*1j)
    return val
 


#The original k vectors for the two orbitals in our spinor set.  
klist=np.array([[0.25,0.25,0.25],[0.1,0.2,0.3]])
print(klist)


#evaluate a slater determinant with up orbitals determined by kup
#and down orbitals determined by kdn.
def eval_det(R,spins,kup,kdn):
    if(len(kup)!=len(kdn)):
        raise Exception("Spinor SPO's not the same length")
    nelec=len(R)
    norb=len(kup)
    M=np.zeros((nelec,norb),dtype=np.complex_)
    pi=math.pi

    #This section, we deliberately add phase factors corresponding to 
    #a spin orientation of theta=45deg, phi=45deg.  
    theta=pi/4.0
    phi=pi/4.0
    up_phase=np.cos(theta/2.0)
    dn_phase=np.exp(1j*phi)*np.sin(theta/2.0)
    for iat in range(0,nelec):
        spo_up=up_phase*eval_pw_spo(R[iat],kup)
        spo_dn=dn_phase*eval_pw_spo(R[iat],kdn)
        spinor_row=np.exp(1j*spins[iat])*spo_up+np.exp(-1j*spins[iat])*spo_dn
        M[iat]=spinor_row
    return np.linalg.det(M)
        
    
#the spin matrix elements in overcomplete basis.
def evaluate_sx(s1,s2):
    return 2.0*np.cos(s1+s2)
def evaluate_sy(s1,s2):
    return 2.0*np.sin(s1+s2)
def evaluate_sz(s1,s2):
    return -2.0*1j*np.sin(s2-s1)


#Much like what QMCPACK will do, this will split the 
#interval of 0-2pi into ngrid grid points.  Then it will perform the
#spin integral (with trapezoid rule) for electron iat.  

def integrate(R,spins,iat,kup,kdn,ngrid=500):
    psi0=eval_det(R,spins,kup,kdn)
    ratios=np.zeros(ngrid)
    
    sxgrid=np.zeros(ngrid)
    sygrid=np.zeros(ngrid)
    szgrid=np.zeros(ngrid)
    
    twopi=2.0*math.pi
    
    delt=twopi/(ngrid-1.0)
    spins_=copy.deepcopy(spins)
    for i in range(0,ngrid):
        s=delt*i
        spins_[iat]=s
        
        psi_new=eval_det(R,spins_,kup,kdn)
        ratio=psi_new/psi0
        sxgrid[i]=evaluate_sx(spins[iat],s)*ratio
        sygrid[i]=evaluate_sy(spins[iat],s)*ratio
        szgrid[i]=evaluate_sz(spins[iat],s)*ratio
        
    sx=np.trapz(sxgrid,dx=delt)/twopi
    sy=np.trapz(sygrid,dx=delt)/twopi
    sz=np.trapz(szgrid,dx=delt)/twopi
 
    return np.array([sx,sy,sz])


print("Starting spin values")
print(R,spins)
sup=integrate(R,spins,0,klist,klist,ngrid=5000)
sdn=integrate(R,spins,1,klist,klist,ngrid=5000)
print("sup=",sup)
print("sdn=",sdn)
print("sup+sdn=",sup+sdn)


#This is the beginning of the Monte Carlo integration with respect to 
#s1 and s2 for the total spin expectation value.  Recall that analytically,
#we expect the two electrons to each have a spin expectation value
#of [0.5,0.5,1/sqrt(2)] since we deliberately introduced the appropriate phase 
#factors for up and down channels.  

nsamps=200000
report_block=100
mysum_up=np.zeros(3)
mysum_dn=np.zeros(3)
spins_=np.zeros(spins.shape)
spins_prop_=np.zeros(spins.shape)
accept=0
sx=[]
sy=[]
sz=[]

twopi=2.0*math.pi
for i in range(0,nsamps):

    x=np.random.uniform(low=0,high=twopi)
    y=np.random.uniform(low=0,high=twopi)
    #reference
    psi0=eval_det(R,spins_,klist,klist)
    #prposed move
    spins_prop_[0]=x
    spins_prop_[1]=y
    psinew=eval_det(R,spins_prop_,klist,klist)
    ratio=psinew/psi0
    accept_prob=(ratio.conjugate()*ratio)
    xi=np.random.uniform()
    if(xi < min(accept_prob,1.0)):
        spins_=copy.deepcopy(spins_prop_)
        accept+=1

    mysum_up+=integrate(R,spins_,0,klist,klist)
    mysum_dn+=integrate(R,spins_,1,klist,klist)
    if(i%report_block == 0): 
        accept_ratio=accept/float(report_block)
        block_avg=(mysum_up+mysum_dn)/float(report_block)
        accept=0
        mysum_up=0
        mysum_dn=0
        sx.append(block_avg[0])
        sy.append(block_avg[1])
        sz.append(block_avg[2])
        print(i,accept_ratio,spins_[0],spins_[1],block_avg[0],block_avg[1],block_avg[2])

chop=2
print("Sx = ",np.mean(sx[chop:-1]),np.std(sx[chop:-1])/np.sqrt(len(sx[chop:-1])),np.var(sx[chop:-1]))
print("Sy = ",np.mean(sy[chop:-1]),np.std(sy[chop:-1])/np.sqrt(len(sy[chop:-1])),np.var(sy[chop:-1]))
print("Sz = ",np.mean(sz[chop:-1]),np.std(sz[chop:-1])/np.sqrt(len(sz[chop:-1])),np.var(sz[chop:-1]))





