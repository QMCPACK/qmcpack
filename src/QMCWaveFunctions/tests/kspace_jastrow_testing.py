#!/usr/bin/env python
import numpy as np

def makeTriclinic(a, b, c, alpha, beta, gamma):
    A = np.array([a,0,0])
    B = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = np.sqrt(c**2 - cx**2 - cy**2)
    C = np.array([cx, cy, cz])

    return np.array([A,B,C])

def getVolume(L):
    V = np.dot(L[0], np.cross(L[1],L[2]))
    return V

def getReciprocal(L):
    Gt = (2 * np.pi) * np.linalg.inv(L)
    return Gt.T

def registerGvecs(gmax, G):
    def checkInclude(i,j,k):
        if i > 0:
            return True
        elif i == 0:
            if j > 0:
                return True
            elif j == 0 and k > 0:
                return True
        else:
            return False
    gvecs = []
    for m1 in range(-5, 6):
        for m2 in range(-5, 6):
            for m3 in range(-5, 6):
                include = checkInclude(m1, m2, m3)
                gvec = m1 * G[0] + m2 * G[1] + m3 * G[2]
                dot = np.dot(gvec, gvec)
                if dot < gmax**2 and dot > 0 and include:
                    gvecs.append(gvec)

    idx = np.argsort([ np.dot(gvecs[i], gvecs[i]) for i in range(len(gvecs)) ])
    return np.array(gvecs)[idx]

def logJk1(r, R, gvecs, coeffs, V):
    # in qmcpack, we factor out the constant term rhoIG. I calculate it here following the definition in the manual
    # but then just set it to 1. We also take only the real part in the code
    val = 0
    for ig in range(len(gvecs)):
        rhomG = 0
        for ie in range(len(r)):
            dot = np.dot(-gvecs[ig], r[ie])
            rhomG += np.cos(dot) + 1j * np.sin(dot)
        rhoIG = 0
        for iI in range(len(R)):
            dot = np.dot(gvecs[ig], R[iI])
            rhoIG += np.cos(dot) + 1j * np.sin(dot)
        rhoIG = 1
        val += coeffs[ig] * rhoIG * rhomG
    return (val / V).real

def logJk2(r, gvecs, coeffs, V):
    val = 0
    for ig in range(len(gvecs)):
        rhoG = 0
        rhomG = 0
        for ie in range(len(r)):
            dot = np.dot(gvecs[ig], r[ie])
            rhoG += np.cos(dot) + 1j * np.sin(dot)
            dot = np.dot(-gvecs[ig], r[ie])
            rhomG += np.cos(dot) + 1j * np.sin(dot)
        val += coeffs[ig] * rhoG * rhomG
    return (val / V).real

def dlogJk1(r, R, gvecs, coeffs, V):
    # in qmcpack, we factor out the constant term rhoIG. I calculate it here following the definition in the manual
    # but then just set it to 1. 
    derivs = np.zeros(len(coeffs), dtype=complex)
    for ig in range(len(gvecs)):
        rhomG = 0
        for ie in range(len(r)):
            dot = np.dot(-gvecs[ig], r[ie])
            rhomG += np.cos(dot) + 1j*np.sin(dot)
        rhoIG = 0
        for iI in range(len(R)):
            dot = np.dot(gvecs[ig], R[iI])
            rhoIG += np.cos(dot) + 1j*np.sin(dot)
        rhoIG = 1
        val = rhoIG * rhomG / V

        #df(z)/da = df/dz dz/da and df(z)/db = df/dz dz/db for z = a + ib
        #so df(z)/da = df/dz and df(z)/db = i df/dz where f(z) is real even though z is complex
        derivs[ig] = np.real(val.real) + 1.j * np.real(1.j * val) 
    return derivs

def dlogJk2(r, gvecs, coeffs, V):
    derivs = np.zeros(len(coeffs))
    for ig in range(len(gvecs)):
        rhoG = 0
        rhomG = 0
        for ie in range(len(r)):
            dot = np.dot(-gvecs[ig], r[ie])
            rhomG += np.cos(dot) + 1j * np.sin(dot)
            dot = np.dot(gvecs[ig], r[ie])
            rhoG += np.cos(dot) + 1j*np.sin(dot)
        deriv = rhoG * rhomG / V
        assert np.abs(deriv.imag) < 1e-12
        derivs[ig] = deriv.real
    return derivs  

def runChecks(): 
    np.random.seed(123454321)
    L = makeTriclinic(1,2,3, np.radians(65), np.radians(84), np.radians(98))
    G = getReciprocal(L)
    kf = 5
    gvecs = registerGvecs(kf, G)
    
    print("Lattice :\n{}" .format(L))
    print("kmax    : {}".format(kf))
    print("ng      : {}".format(len(gvecs)))
    print("gvecs   :\n{}".format(gvecs))

    #elec
    r = np.array([[0.1,0.2,0.3], [0.4,0.5,0.6]])
    #ion
    R = np.array([[0.7, 0.8, 0.9]])

    #randomize coefficients
    jk1c = []
    jk2c = []
    for i in range(len(gvecs)):
        jk1c.append( np.random.random() + 1j*np.random.random() )
        jk2c.append( np.random.random())

    print("Jk1 coeffs: {}".format(jk1c))
    print("Jk2 coeffs: {}".format(jk2c))

    V = getVolume(L)
    val1 = logJk1(r, R, gvecs, jk1c, V)
    val2 = logJk2(r, gvecs, jk2c, V)
    d1 = dlogJk1(r, R, gvecs, jk1c, V)
    d2 = dlogJk2(r, gvecs, jk2c, V)

    print()
    print("Jk1  : {}".format(val1))
    print("JK2  : {}".format(val2))
    print("dratio Jk1: {}".format( d1))
    print("dratio Jk2: {}".format( d2))

runChecks()
