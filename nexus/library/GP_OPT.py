"""
This code uses Gaussian Processes to estimate optimization surfaces of noisy
data.
"""

import numpy as np


def lhs(n, samples=None, iterations=None):
    """
    Generate a latin-hypercube design
    
    Parameters
    ----------
    n : int
        The number of factors to generate samples for
    
    Optional
    --------
    samples : int
        The number of samples to generate for each factor (Default: n)

    iterations : int
        The number of iterations in the maximin and correlations algorithms
        (Default: 5).
    
    Returns
    -------
    H : 2d-array
        An n-by-samples design matrix that has been normalized so factor 
        values are uniformly spaced between zero and one.
        
    
    Example
    -------
             
    A 4-factor design with 5 samples and 10 iterations::
    
        >>> lhs(4, samples=5, iterations=10)
    
    """
    H = None
    
    if samples is None:
        samples = n
    
    if iterations is None:
        iterations = 5
    
    H = _lhsmaximin(n, samples, iterations)
    
    
    return H

###############################################################################

def _lhsclassic(n, samples):
    # Generate the intervals
    cut = np.linspace(0, 1, samples + 1)    
    
    # Fill points uniformly in each interval
    u = np.random.rand(samples, n)
    a = cut[:samples]
    b = cut[1:samples + 1]
    rdpoints = np.zeros_like(u)
    for j in range(n):
        rdpoints[:, j] = u[:, j]*(b-a) + a
    
    # Make the random pairings
    H = np.zeros_like(rdpoints)
    for j in range(n):
        order = np.random.permutation(range(samples))
        H[:, j] = rdpoints[order, j]
    
    return H
    
    
###############################################################################

def _lhsmaximin(n, samples, iterations):
    maxdist = 0
    
    # Maximize the minimum distance between points
    for i in range(iterations):
        Hcandidate = _lhsclassic(n, samples)

        d = _pdist(Hcandidate)
        if maxdist<np.min(d):
            maxdist = np.min(d)
            H = Hcandidate.copy()
    
    return H
    
###############################################################################

def _pdist(x):
    """
    Calculate the pair-wise point distances of a matrix
    
    Parameters
    ----------
    x : 2d-array
        An m-by-n array of scalars, where there are m points in n dimensions.
    
    Returns
    -------
    d : array
        A 1-by-b array of scalars, where b = m*(m - 1)/2. This array contains
        all the pair-wise point distances, arranged in the order (1, 0), 
        (2, 0), ..., (m-1, 0), (2, 1), ..., (m-1, 1), ..., (m-1, m-2).
    
    Examples
    --------
    ::
    
        >>> x = np.array([[0.1629447, 0.8616334],
        ...               [0.5811584, 0.3826752],
        ...               [0.2270954, 0.4442068],
        ...               [0.7670017, 0.7264718],
        ...               [0.8253975, 0.1937736]])
        >>> _pdist(x)
        array([ 0.6358488,  0.4223272,  0.6189940,  0.9406808,  0.3593699,
                0.3908118,  0.3087661,  0.6092392,  0.6486001,  0.5358894])
              
    """
    
    x = np.atleast_2d(x)
    assert len(x.shape)==2, 'Input array must be 2d-dimensional'
    
    m, n = x.shape
    if m<2:
        return []
    
    d = []
    for i in range(m - 1):
        for j in range(i + 1, m):
            d.append((sum((x[j, :] - x[i, :])**2))**0.5)
    
    return np.array(d)

def ackleyf(x):

    """
    Calculate the d-dimensional Ackley function on [-1,1]^d
    
    see: Floudas, C. A., & Pardalos, P. M. (1990). A collection of test 
    problems for constrained global optimization algorithms. In G. Goods & 
    J. Hartmanis (Eds.) Lecture notes in computer science: Vol. 455. 
    Berlin: Springer
    
    Parameters
    ----------
    x : 1 x d-array, x \in [-1,1]^d
    
    Returns
    -------
    y : value of Ackely function
    
    """
    d = x.shape[1]
    x = 0.5*x
    y = 20.0+np.exp(1.0)-20.0*np.exp(-0.2*np.sqrt(np.sum(x**2)/float(d)))- \
        np.exp(np.sum(np.cos(2.0*np.pi*x))/float(d))
    
    
    return y

def cfgp(r,ep):

    """
    Calculate correlation function for Gaussian Process
    
    
    Parameters
    ----------
    r : val r>=0
    ep: val ep>0
    
    Returns
    -------
    y : value of correlation function
    
    """

    y = np.exp(-ep*r**2)
    
    return y


def cmgp(S1,S2,ep):

    """
    Calculate correlation matrix for Gaussian Process
    
    
    Parameters
    ----------
    S1: n x d array of points
    S2: m x d array of points 
    ep: GP kernel parameter ep>0
    
    Returns
    -------
    K : n x m correlation matrix
    
    """
    
    n = S1.shape[0]
    m = S2.shape[0]
    K = np.zeros([n,m])
    
    for jn in range(n):
        K[jn,:] = np.sqrt(np.sum((np.repeat(np.reshape(S1[jn,:],[1,d]),\
             m,axis=0)-S2)**2,axis=1))
    
    K = cfgp(K,ep)
    
    return K


def sgd(hyper_cube,co_GP,mol_pos,ep,tol=None,verbose=None):

    """
    Calculate minumums of Gaussian Process using stochastic gradient descent
    
    
    Parameters
    ----------
    co_GP:      GP coefficients
    mol_pos:    GP centers 
    ep:         GP kernel parameter ep>0
    tol:        Tolerance for termination
    
    Returns
    -------
    min_loc : Local minimums
    min_E   : Energy Estimation of local minimums
    """
    
    
    if tol is None:
        tol = 1e-4
    
    if verbose is None:
        verbose=0
    
    npoints = mol_pos.shape[0]
    d = mol_pos.shape[1]
    epm = np.min(_pdist(mol_pos))
    
    
    best_pos = np.copy(mol_pos)
    delta_pos = epm*np.ones([npoints,d])/10.0
    best_E = np.matmul(cmgp(best_pos,mol_pos,ep),co_GP)
    
    cm = np.max(delta_pos)
    use_ind = np.where(np.max(delta_pos,axis=1)>epm*tol)
    while cm>epm*tol:
        for jd in range(d):
            did_move = np.zeros([npoints,1])
            t_pos =  np.copy(best_pos)
            t_pos[:,jd] = t_pos[:,jd] + delta_pos[:,jd]
            for jdim in range(d):
                t_pos[np.where(t_pos[:,jdim]>hyper_cube[jdim,1]),jdim]= \
                    hyper_cube[jdim,0]
                    
            c_E = np.matmul(cmgp(t_pos,mol_pos,ep),co_GP)
            for jp in range(npoints):
                if c_E[jp,0]<best_E[jp,0]:
                    did_move[jp,0] = 1
                    best_E[jp,0] = c_E[jp,0]
                    best_pos[jp,:] = t_pos[jp,:]
                    
            t_pos =  np.copy(best_pos)
            t_pos[:,jd] = t_pos[:,jd] - delta_pos[:,jd]
            for jdim in range(d):
                t_pos[np.where(t_pos[:,jdim]<hyper_cube[jdim,0]),jdim]= \
                    hyper_cube[jdim,0]
                    
            c_E = np.matmul(cmgp(t_pos,mol_pos,ep),co_GP)
            for jp in range(npoints):
                if c_E[jp,0]<best_E[jp,0]:
                    did_move[jp,0] = 1
                    best_E[jp,0] = c_E[jp,0]
                    best_pos[jp,:] = t_pos[jp,:]
                if did_move[jp,0] == 1:
                    delta_pos[jp,jd] = 1.25*delta_pos[jp,jd]
                else:
                    delta_pos[jp,jd] = 0.75*delta_pos[jp,jd]
        
        cm = np.max(delta_pos)
        if verbose == 1:
            use_ind = np.where(np.max(delta_pos,axis=1)>epm*tol)
            print cm,np.shape(use_ind)[1]
      
              

    c_class = 0
    ind_marker = -np.ones([npoints,1])
    while np.min(ind_marker)<0:
        indc = np.where(ind_marker<0)
        cdata = best_pos[indc[0],:]
        cvec = np.reshape(cdata[0,:],[1,d])
        m = cdata.shape[0]
        dis_v = np.sqrt(np.sum((np.repeat(cvec,m,axis=0)-cdata)**2,axis=1))
        class_ind = indc[0][np.where(dis_v<10*epm*tol)]
        ind_marker[class_ind,:] = c_class
        c_class = c_class+1
        
    min_loc = np.zeros([c_class,d])
    min_E = np.zeros([c_class,1])  
    
    for jc in range(c_class):
        indc = np.where(ind_marker==jc)
        c_E = best_E[indc[0],:]
        t_pos = best_pos[indc[0],:]
        indb = np.argmin(c_E)
        min_E[jc,:] = c_E[indb,:]
        min_loc[jc,:] = t_pos[indb,:]
        
        
    return min_loc,min_E
   
    
## Parameters ##
d = 22                           # Dimension of problem
nits = 6                        # Number of iterations
sigma = 0                       # Uncertainty in data
hyper_cube = np.ones([d,2])     # Upper and lower bound for each dimension
hyper_cube[:,0] = -hyper_cube[:,0]

## Suggested Parameters ##
npoints = (d+2)*(d+1)+1    
nlhc_its = 10

## Generate random sample from LHS
mol_pos = np.repeat(np.reshape(hyper_cube[:,1]-hyper_cube[:,0],[1,d]) \
    ,npoints,axis=0)*lhs(d, samples=npoints,iterations=nlhc_its)+np.repeat( \
    np.reshape(hyper_cube[:,0],[1,d]),npoints,axis=0)

hyper_cube[:,0] = np.reshape(np.min(mol_pos,axis=0),d)
hyper_cube[:,1] = np.reshape(np.max(mol_pos,axis=0),d)

E = np.zeros([npoints,1])
## Randomly pick center for Ackley function
x_true = 0.05*(2.0*np.random.rand(1,d)-1.0)

## Determine Energies using Ackley function
for jp in range(npoints):
    E[jp,0] = ackleyf(mol_pos[jp,:]-x_true)

## Estimate Optimal Choice for ep
ep  = 1.0/np.max(_pdist(mol_pos))

## Build Gaussian Process
K = cmgp(mol_pos,mol_pos,ep) 
Ks = K + (sigma**2)*np.identity(npoints)
co_GP = np.linalg.solve(Ks,E)

print "Starting Stochastic Grandient Descent"
## Use stochastic gradient descent to determine minimums
min_loc,min_E = sgd(hyper_cube,co_GP,mol_pos,ep)

indb = np.argmin(min_E)
bmin_loc = np.reshape(min_loc[indb,:],[1,d])

print "Error in min position",np.sqrt(np.sum((bmin_loc-x_true)**2)),"its=0"

if d == 1:
    import matplotlib.pyplot as plt
    npts = 100
    plot_points = np.reshape(np.linspace(hyper_cube[0,0],hyper_cube[0,1], \
        num=npts),[npts,1])
    gp_E = np.matmul(cmgp(plot_points,mol_pos,ep),co_GP)
    Et = np.zeros([npts,1])
    for jp in range(npts):
        Et[jp,0] = ackleyf(plot_points[jp,:]-x_true)
    plt.plot(plot_points,gp_E,'b-',plot_points,Et,'g-',mol_pos,E,'r.')
    plt.show()

chyper_cube = np.ones([d,2])
chyper_cube[:] = hyper_cube[:]
for jits in range(nits-1):
    chyper_cube = chyper_cube+np.repeat(np.reshape(bmin_loc,[d,1]),2,axis=1)
    chyper_cube = chyper_cube/2.0
    ## Generate random sample from LHS
    mol_pos = np.repeat(np.reshape(chyper_cube[:,1]-chyper_cube[:,0],[1,d]) \
        ,npoints,axis=0)*lhs(d, samples=npoints,iterations=nlhc_its)+np.repeat( \
        np.reshape(chyper_cube[:,0],[1,d]),npoints,axis=0)
    chyper_cube[:,0] = np.reshape(np.min(mol_pos,axis=0),d)
    chyper_cube[:,1] = np.reshape(np.max(mol_pos,axis=0),d)
    ## Determine Energies using Ackley function
    for jp in range(npoints):
        E[jp,0] = ackleyf(mol_pos[jp,:]-x_true)
    
    ## Estimate Optimal Choice for ep
    ep  = 1.0/np.max(_pdist(mol_pos))
    
    ## Build Gaussian Process
    K = cmgp(mol_pos,mol_pos,ep) 
    Ks = K + (sigma**2)*np.identity(npoints)
    co_GP = np.linalg.solve(Ks,E)
    print "Starting Stochastic Grandient Descent"
    ## Use stochastic gradient descent to determine minimums
    min_loc,min_E = sgd(chyper_cube,co_GP,mol_pos,ep)
    indb = np.argmin(min_E)
    bmin_loc = np.reshape(min_loc[indb,:],[1,d])

    print "Error in min position",np.sqrt(np.sum((bmin_loc-x_true)**2)),\
        "its=",jits+1
    
    if d == 1:
        plot_points = np.reshape(np.linspace(chyper_cube[0,0], \
            chyper_cube[0,1],num=npts),[npts,1])
        gp_E = np.matmul(cmgp(plot_points,mol_pos,ep),co_GP)
        Et = np.zeros([npts,1])
        for jp in range(npts):
            Et[jp,0] = ackleyf(plot_points[jp,:]-x_true)
        plt.plot(plot_points,gp_E,'b-',plot_points,Et,'g-',mol_pos,E,'r.')
        plt.show()
    
    
    
    
    

# Code to test out LOOCV for optimal ep
#Kinv = np.linalg.pinv(K)
#loocv_est = np.zeros([npoints,1])
#
#for jsub in range(npoints):
#    loocv_est[jsub,0] = co_GP[jsub,0]/Kinv[jsub,jsub]
#
#print sum(abs(loocv_est))

#loocv_true = np.zeros([npoints,1])   
#for jsub in range(npoints):   
#    mol_pos_sub = np.zeros([npoints-1,d])
#    E_sub = np.zeros([npoints-1,1])
#    cp = 0
#    for jp in range(npoints):
#        if jp != jsub:
#            mol_pos_sub[cp,:] = mol_pos[jp,:]
#            E_sub[cp] = E[jp]
#            cp = cp+1
#            
#    K_sub = cmgp(mol_pos_sub,mol_pos_sub,ep)
#    Kinv_sub = np.linalg.pinv(K_sub)
#    co_GP_sub = np.matmul(Kinv_sub,E_sub)
#    S1 = np.zeros([1,d])
#    S1[0,:] = mol_pos[jsub,:]
#    loocv_true[jsub] = E[jsub]-np.matmul(cmgp(S1,mol_pos_sub,ep),co_GP_sub)
#    



