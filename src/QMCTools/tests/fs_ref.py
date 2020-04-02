import numpy as np

A = np.matrix([[1,0,0],
               [0,1,0],
               [0,0,1]]) #lattice
print("Direct Lattice:\n{}".format(A))

B = 2*np.pi*(np.linalg.inv(A)).H  #recip lattice
print("Recip Lattice:\n{}".format(B))

rc = 0.5*np.min(np.sqrt(np.sum(np.square(A),1)))
lrdimcut = 40
kc = lrdimcut / rc
print("lrdim: {}, rc: {}, kc: {}".format(lrdimcut,rc,kc))

#electronic structure, p. 85
mmax = np.floor(np.sqrt(np.sum(np.square(A),1)) * (kc/(2*np.pi))) + 1
mmax = np.array(mmax,dtype=int).reshape((3,)) #matrix to array 

kpts = [] #translations of recip lattice
kmag = [] #magnitude
for i in range(-mmax[0], mmax[0] + 1):
    for j in range(-mmax[1], mmax[1] + 1):
        for k in range(-mmax[2], mmax[2] + 1):
            if (i == 0) and (j==0) and (k==0):
                continue
            kvec = np.matrix([i,j,k])
            kcart = np.array(np.dot(kvec,B)).reshape((3,))
            if np.linalg.norm(kcart) > kc:
                continue
            kpts.append(np.array(kvec).reshape((3,)))
            kmag.append(np.linalg.norm(kcart))

kpts = np.array(kpts)
kmag = np.array(kmag)

idx = np.argsort(kmag)
kpts = kpts[idx]
kmag = kmag[idx]

# 1-exp(-k^2) and k is unit k
sks = [] 
with open('simple_Sk.dat','w') as f:
    f.write('#  kx  ky  kz  Sk  err\n')
    for i in range(len(kpts)):
        kcart = np.array(np.dot(kpts[i],B)).reshape((3,))
        kunit = kcart / (2*np.pi) 
        k = np.linalg.norm(kunit)
        sk = 1-np.exp(-0.075*k*k)
        sks.append(sk)
        f.write('{0} {1} {2} {3} {4}\n'.format(kcart[0],kcart[1],kcart[2],sk,0.01))


print("Ewald Handler Corrections: ")
#corrections
vol = np.abs(np.linalg.det(A))
sigma2 = 0.5*kc/rc
vsum = 0
for i in range(len(kpts)):
    k2 = kmag[i]*kmag[i]
    vk = 4*np.pi / (k2 * vol) * np.exp(-0.25*k2/sigma2)
    vsum += 0.5*vk*sks[i]
print("  Discrete: {}".format(vsum))

vint = 1.066688342657357 #from analytic mathematica calculation
