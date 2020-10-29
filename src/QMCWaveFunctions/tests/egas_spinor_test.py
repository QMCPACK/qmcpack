import numpy as np

#This script computes all the reference values for a Slater determinant of spinors.  
#electron gas orbitals are used for easy analytic treatments.  
#
#I compute the value, all gradients w.r.t. position and spin, and the laplacians. 
#
# A trial spin+position move is generated, and the ratio, value, and gradients for the new configuration
#  are computed.  
#
#  Ray Clay, rclay@sandia.gov, Dec 2 2019.

def spinor_val(r,s,k1,k2):
  kph1=complex(0,np.dot(k1,r))
  kph2=complex(0,np.dot(k2,r))
  spinph=complex(0,s);

  return np.exp(spinph)*np.exp(kph1)+np.exp(-spinph)*np.exp(kph2)

def spinor_grad(r,s,k1,k2):
  kph1=complex(0,np.dot(k1,r))
  kph2=complex(0,np.dot(k2,r))
  spinph=complex(0,s)
  eye=complex(0,1)

  return np.exp(spinph)*eye*k1*np.exp(kph1)+eye*k2*np.exp(-spinph)*np.exp(kph2)

def spinor_lapl(r,s,k1,k2):
  kph1=complex(0,np.dot(k1,r))
  kph2=complex(0,np.dot(k2,r))
  spinph=complex(0,s)
  eye=complex(0,1)
  return -1.0*np.dot(k1,k1)*np.exp(spinph)*np.exp(kph1)+-1.0*np.dot(k2,k2)*np.exp(-spinph)*np.exp(kph2)

def spinor_spingrad(r,s,k1,k2):
  kph1=complex(0,np.dot(k1,r))
  kph2=complex(0,np.dot(k2,r))
  spinph=complex(0,s)
  eye=complex(0,1)

  return eye*np.exp(spinph)*np.exp(kph1)-eye*np.exp(-spinph)*np.exp(kph2)


def spinor_matrix(R,s,kup,kdn):
  M=np.zeros((len(R),len(kup)),dtype=complex);

 
  for iat in range(0,len(R)):
    for norb in range(0,len(kup)):
      M[iat][norb]=spinor_val(R[iat],s[iat],kup[norb],kdn[norb])

  return M

def compute_row_spinor_val(r,s,kup,kdn):
  row=np.zeros(len(kup),dtype=complex)
  for norb in range(0,len(kup)):
    row[norb]=spinor_val(r,s,kup[norb],kdn[norb])
  return row

def compute_row_spinor_grad(r,s,kup,kdn):
  rowx=np.zeros(len(kup),dtype=complex)
  rowy=np.zeros(len(kup),dtype=complex)
  rowz=np.zeros(len(kup),dtype=complex)
  for norb in range(0,len(kup)):
    g=spinor_grad(r,s,kup[norb],kdn[norb])
    rowx[norb]=g[0]
    rowy[norb]=g[1]
    rowz[norb]=g[2]
  return rowx,rowy,rowz
    
def compute_row_spinor_lapl(r,s,kup,kdn):
  row=np.zeros(len(kup),dtype=complex)
  for norb in range(0,len(kup)):
    row[norb]=spinor_lapl(r,s,kup[norb],kdn[norb])
  return row

def compute_row_spinor_spingrad(r,s,kup,kdn):
  row=np.zeros(len(kup),dtype=complex)
  for norb in range(0,len(kup)):
    row[norb]=spinor_spingrad(r,s,kup[norb],kdn[norb])
  return row



    
R=np.zeros((3,3))
R[0]=np.array([ 0.1,-0.3,1.0]) 
R[1]=np.array([-0.1, 0.3,1.0]) 
R[2]=np.array([ 0.1, 0.2,0.3]) 

spins=np.array([0.0,0.2,0.4])

kup=np.zeros((3,3))
kdn=np.zeros((3,3))

kup[0]=np.array([0,0,0]);
kup[1]=np.array([0.1,0.2,0.3]);
kup[2]=np.array([0.4,0.5,0.6]);

kdn[0]=np.array([0,0,0]);
kdn[1]=np.array([-0.1,0.2,-0.3]);
kdn[2]=np.array([0.4,-0.5,0.6]);

Mref= spinor_matrix(R,spins,kup,kdn)
Mtmp = np.copy(Mref)
det_ref=np.linalg.det(Mref)
log_ref=np.log(det_ref)

print("############ REFERENCE CONFIG #################")
print(" R = ",R)
print(" spins = ",spins)
print("   -----TOTAL MATRIX QUANTITIES-----")
print("  M = ", Mref)
print(" ")
print("  det(M) = ",det_ref)
print(" ")
print("  log(det) = ",log_ref)

#### Now we're going to compute the 3 gradient components:
G=np.zeros((3,3),dtype=complex)
L=np.zeros(3,dtype=complex)
SG=np.zeros(3,dtype=complex)

for iat in range(0,3):
  r=R[iat]
  s=spins[iat]
  gxr,gyr,gzr=compute_row_spinor_grad(r,s,kup,kdn)
  Mtmp=np.copy(Mref)
  Mtmp[iat] = gxr
  gx=np.linalg.det(Mtmp)/det_ref
  Mtmp[iat] = gyr
  gy=np.linalg.det(Mtmp)/det_ref
  Mtmp[iat] = gzr
  gz=np.linalg.det(Mtmp)/det_ref

  G[iat][0]=gx
  G[iat][1]=gy
  G[iat][2]=gz

  Mtmp[iat]=compute_row_spinor_lapl(r,s,kup,kdn)
  L[iat] = np.linalg.det(Mtmp)/det_ref - np.dot(G[iat],G[iat]) #QMCPACK returns nabla ln(psi), hence this transformation

  Mtmp[iat]=compute_row_spinor_spingrad(r,s,kup,kdn)
  SG[iat]=np.linalg.det(Mtmp)/det_ref

print(" ")
print("  G = ",G)
print(" ")
print("  L = ",L)
print(" ")
print(" SG = ",SG)
print("---------------------------------------")

iel = 1
#Now print results for iat move.
print(" VALUE  =  ", det_ref)
print("  GX = ",G[iel][0])
print("  GY = ",G[iel][1])
print("  GZ = ",G[iel][2])
print("  SG = ",SG[iel])

print(" ")
print("Now we make a particle/spin move for particle ",iel)
#### Now we're going to compute the ratio and gradients at a proposed move location.  
dr=np.array([0.1,-0.05,0.2]);
ds=0.3;

print("  dr = ",dr)
print("  ds = ",ds)

rnew=R[iel]+dr
snew=spins[iel]+ds

print(" ")
print(" r_old = ",R[iel]," r_new = ",rnew)
print(" s_old = ",spins[iel]," s_new = ",snew)

Mtmp=np.copy(Mref)
Mtmp[iel]=compute_row_spinor_val(rnew,snew,kup,kdn)
det_new = np.linalg.det(Mtmp)


##Now to compute the gradients
gxr_new,gyr_new,gzr_new = compute_row_spinor_grad(rnew,snew,kup,kdn)
Mtmp=np.copy(Mref)
Mtmp[iel] = gxr_new
gx_new=np.linalg.det(Mtmp)/det_new
Mtmp[iel] = gyr_new
gy_new=np.linalg.det(Mtmp)/det_new
Mtmp[iel] = gzr_new
gz_new=np.linalg.det(Mtmp)/det_new

#now to compute spin gradient at new position.
Mtmp[iel]=compute_row_spinor_spingrad(rnew,snew,kup,kdn)
sg_new=np.linalg.det(Mtmp)/det_new

print(" NEW VALUE = ",det_new)
print(" NEW LOG(VALUE) = ",np.log(det_new))
print(" RATIO  =  ", det_new/det_ref)
print("  GX = ", gx_new)
print("  GY = ", gy_new)
print("  GZ = ", gz_new)
print("  SG = ", sg_new)
