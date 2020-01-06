##################################################################
##  (c) Copyright 2018-  by Richard K. Archibald                ##
##                       and Jaron T. Krogel                    ##
##################################################################


#====================================================================#
#  gaussian_process.py                                               #
#    Uses Gaussian Processes to estimate optimization surfaces of    #
#    noisy data.                                                     #
#                                                                    #
#  Content summary:                                                  #
#    lhs                                                             #
#      Function to generate a latin-hypercube design.                #
#                                                                    #
#    _lhsmaximin                                                     #
#      Function to maximize the minimum distance between points.     #
#      Subroutine of lhs.                                            #
#                                                                    #
#    _lhsclassic                                                     #
#      Subroutine of _lhsmaximin                                     #
#                                                                    #
#    _pdist                                                          #
#      Function to calculate the pair-wise point distances of a      #
#      matrix.  Subroutine of lhsmaximin.                            #
#                                                                    #
#    _fdist                                                          #
#      Function to calculate distances between all points.           #
#      Subroutine of find_hyperparameters, find_min, and             #
#      find_new_points.                                              #
#                                                                    #
#    cfgp                                                            #
#      Function to calculate correlation function for Gaussian       #
#      Process.  Subroutine of cmgp.                                 #
#                                                                    #
#    cmgp                                                            #
#      Function to calculate correlation matrix for Gaussian         #
#      Process.  Subroutine of find_min.                             #
#                                                                    #
#    find_hyperparameters                                            #
#      Function to find best Gaussian Process hyperparameters.       #
#      Subroutine of gp_opt.                                         #
#                                                                    #
#    find_min                                                        #
#      Function to find local minimum of Gaussian Process surface.   #
#      Subroutine of gp_opt.                                         #
#                                                                    #
#    find_new_points                                                 #
#      Function to estimate optimal new points for Gaussian Process. #
#      Subroutine of gp_opt.                                         #
#                                                                    #
#    multiscale_E                                                    #
#      Function to perform multiscale normalization of E.            #
#      Subroutine of gp_opt.                                         #
#                                                                    #
#    inv_multiscale_E                                                #
#      Function to reverse multiscale normalization of E.            #
#      Subroutine of gp_opt.                                         #
#                                                                    #
#    norm_pos                                                        #
#      Function to normalize molecular positions to unit hyper cube. #
#      Subroutine of gp_opt.                                         #
#                                                                    #
#    inv_norm_pos                                                    #
#      Function to reverse normalization of molecular positions.     #
#      Subroutine of gp_opt.                                         #
#                                                                    #
#    gp_opt_setup                                                    #
#      Function that sets up parameters and variables for GP.        #
#                                                                    #
#    gp_opt                                                          #
#      Function that performs one iteration of Gaussian Process      #
#      optimization.  This is the main interface to Gaussian Process #
#      iterations.                                                   #
#                                                                    #
#    opt_surface_f                                                   #
#      Example function to optimize.                                 #
#                                                                    #
#    gp_opt_example                                                  #
#      Example using Gaussian Process optimization on opt_surface_f. #
#                                                                    #
#====================================================================#

import os
from numpy import array,zeros,append
from generic import obj
from developer import unavailable,available,DevBase,log,warn,error,ci

import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    plt = unavailable('matplotlib.pyplot')
#end try
try:
    from scipy.spatial import Delaunay
except:
    Delaunay = unavailable('scipy.spatial','Delaunay')
#end try
try:
    from scipy.linalg import logm
except:
    logm = unavailable('scipy.linalg','logm')
#end try


class GParmFlags:
    """
    Define an empty class to store GP parameters and information
    """
    pass

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
    
    d = np.zeros((m**2-m)/2)
    c_v = 0
    for i in range(m - 1):
        d[0+c_v:m-i-1+c_v]=np.sqrt(np.sum((np.repeat(np.reshape(x[i, :],\
                        [1,n]),m-i-1,axis=0)-x[i+1:m+1, :])**2,axis=1))
        c_v = c_v+m-i-1
    
    return d

###############################################################################
    
def _fdist(x):

    """
    Calculate Distance between all points
    
    
    Parameters
    ----------
    x: n x d array of points

    Returns
    -------
    Dm : n x m correlation matrix
    
    """
    
    [n,d] = x.shape
    Dm = np.zeros([n,n])
    
    for jn in range(n):
        Dm[jn,:] = np.sqrt(np.sum((np.repeat(np.reshape(x[jn,:],[1,d]),\
             n,axis=0)-x)**2,axis=1))
    
    
    return Dm
###############################################################################
    
def opt_surface_f(x,opt_fun):

    """
    Optimization Surface Function 
        
    Parameters
    ----------
    x:          n x d-array, x \in Hyper cube
    opt_fun:    Optimization Function case
    
    Returns
    -------
    y : value of optimization surface function
    
    """
    [n,d] = x.shape
      
    if opt_fun == 1:
        y = np.zeros([n,1])
        for jp in range(n):
            y[jp,0] = lennard_jones(x[jp,:])
    elif opt_fun == 2:
        y = np.zeros([n,1])
        for jp in range(n):
            y[jp,0] = beale(x[jp,:])
    elif opt_fun == 3:
        y = np.zeros([n,1])
        for jp in range(n):
            y[jp,0] = rosenbrock(x[jp,:])
    elif opt_fun == 4:
        y = np.zeros([n,1])
        for jp in range(n):
            y[jp,0] = goldstein_price(x[jp,:])       
    else:
        y = np.reshape(2.0*np.sin(np.pi*np.sum(x**2,1)/2)-1,[n,1])
            

    return y


###############################################################################
def lennard_jones(x):
    # parameters, 1+2+3*nx
    x = x.ravel()
    nparams = len(x)
    if nparams!=1 and nparams%3!=0 or nparams==0:
        print 'invalid number of parameters for lennard jones'
        exit()
    #end if
    r = [[   0,0,0],
         [x[0],0,0]]
    if nparams>1:
        r.append([x[1],x[2],0])
    #end if
    if nparams>3:
        i=3
        while i+2<nparams:
            r.append([x[i],x[i+1],x[i+2]])
            i+=3
        #end while
    #end if
    r = np.array(r,dtype=float)
    npoints = len(r)
    E1 = 1.0
    rm = 1.0
    E = 0.0
    for i in range(npoints):
        for j in range(i+1,npoints):
            d = np.linalg.norm(r[i]-r[j])
            ELJ = E1*((rm/d)**12-2*(rm/d)**6)
            E += ELJ
        #end for
    #end for
    return E
#end def lennard_jones

###############################################################################
def rosenbrock(x):
    """
    Global minimum in dim 3-7 at (1,1,...,1), local min at (-1,1,...,1)
    """
    x = x.ravel()
    return ( (1-x[0:-1])**2+100*(x[1:]-x[0:-1]**2)**2 ).sum()

###############################################################################
def beale(x):
    """
    2D function, global minimum at (3,.5)
    """
    x,y = x.ravel()
    return (1.5-x+x*y)**2+(2.25-x+x*y**2)**2+(2.625-x+x*y**3)**2


###############################################################################

def goldstein_price(x):
    """
    2D function, global minimum at (0,-1)
    """
    x,y = x.ravel()
    return (1.+(x+y+1.)**2*(19.-14*x+3*x**2-14*y+6*x*y+3*y**2))*(30+(2*x-3*y)**2*(18.-32*x+12*x**2+48*y-36*x*y+27*y**2))

###############################################################################

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

###############################################################################
    
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
    d = S1.shape[1]
    m = S2.shape[0]
    K = np.zeros([n,m])
    
    for jn in range(n):
        K[jn,:] = np.sqrt(np.sum((np.repeat(np.reshape(S1[jn,:],[1,d]),\
             m,axis=0)-S2)**2,axis=1))
    
    K = cfgp(K,ep)
    
    return K

###############################################################################

def find_hyperparameters(mol_pos,E,dE):

    """
    Finds best Gaussian Process hyperparameters
    
    
    Parameters
    ----------
    mol_pos:    GP centers 
    E:          Energy for GP centers
    dE:         Uncertainty in energy
    
    Returns
    -------
    hyp_parm:   GP hyperparameters
    Co_GP:      GP coefficents
    """
    [Nc,d] = mol_pos.shape
    
    W = np.diag(dE[:,0])
    
    D =_fdist(mol_pos)
    
    hyp_parm = np.zeros(2)
    
    hyp_parm[0] = np.log(np.max(np.min(D+np.max(D)*np.eye(Nc),axis=1)))/2.0
    dhyp_parm = np.zeros(2)
    dhyp_parm[:] = 1e-1*abs(hyp_parm[0])
    tol = 1e-2*abs(hyp_parm[0])
    
    K = np.exp(2.0*hyp_parm[1])*cfgp(D,np.exp(-2.0*hyp_parm[0]))
    Kw_inv = np.linalg.pinv(K+W)
    Co_GP = np.matmul(Kw_inv,E)
    
    b_loocv_err = np.matmul(np.transpose(E),Co_GP)+np.trace(logm(K+W))
     
    while max(dhyp_parm)>tol:
        for jp in range(hyp_parm.shape[0]):
            did_move = 0;
            chyp_parm = np.copy(hyp_parm)
            chyp_parm[jp] = chyp_parm[jp]+dhyp_parm[jp]
            K = np.exp(2.0*chyp_parm[1])*cfgp(D,np.exp(-2.0*chyp_parm[0]))
            Kw_inv = np.linalg.pinv(K+W)
            Co_GP = np.matmul(Kw_inv,E)
            
            c_loocv_err = np.matmul(np.transpose(E),Co_GP)+np.trace(logm(K+W))
                
            if c_loocv_err<b_loocv_err:
                b_loocv_err = c_loocv_err
                hyp_parm = np.copy(chyp_parm)
                did_move = 1
            
            chyp_parm = np.copy(hyp_parm)
            chyp_parm[jp] = chyp_parm[jp]-dhyp_parm[jp]
            K = np.exp(2.0*chyp_parm[1])*cfgp(D,np.exp(-2.0*chyp_parm[0]))
            Kw_inv = np.linalg.pinv(K+W)
            Co_GP = np.matmul(Kw_inv,E)
            
            c_loocv_err = np.matmul(np.transpose(E),Co_GP)+np.trace(logm(K+W))
                
            if c_loocv_err<b_loocv_err:
                b_loocv_err = c_loocv_err
                hyp_parm = np.copy(chyp_parm)
                did_move = 1

            
            if did_move == 1:
                dhyp_parm[jp] = 1.25*dhyp_parm[jp]
            else:
                dhyp_parm[jp] = .75*dhyp_parm[jp]
            
        
        
        
    K = np.exp(2.0*hyp_parm[1])*cfgp(D,np.exp(-2.0*hyp_parm[0]))
    Kw_inv = np.linalg.pinv(K+W)
    Co_GP = np.matmul(Kw_inv,E)
        
    return hyp_parm,Co_GP 


###############################################################################

def find_min(Co_GP,mol_pos,E,dE,hyp_parm):

    """
    Finds local minimum of Gaussian Process
    
    
    Parameters
    ----------
    mol_pos:    GP centers 
    E:          Energy for GP centers
    dE:         Uncertainty in energy for GP centers
    Co_GP:      GP coefficents
    hyp_parm:   GP hyperparameters
    
    Returns
    -------
    best_pos:       Estimated local minimum
    best_E:         Estimated local energies
    neig_data_ind:  Neighborhood Data indexes 
    """
    [Nc,d] = mol_pos.shape
    
    if d>1:
        T = Delaunay(mol_pos)
        us_ind = np.setdiff1d(np.arange(Nc),np.unique(T.convex_hull))
        if us_ind.shape[0] == 0:
            us_ind = np.zeros(1, dtype=int)
            us_ind[0] = np.argmin(E)
    else:
        us_ind = np.arange(Nc)
        
    D =_fdist(mol_pos)
    best_pos = mol_pos
    ref_dis = np.min(np.min(D+np.max(D)*np.eye(Nc),axis=1))/10
    K = np.exp(2.0*hyp_parm[1])*cfgp(D,np.exp(-2.0*hyp_parm[0]))
    W = np.diag(dE[:,0])
    Kw_inv = np.linalg.pinv(K+W)

    cnum =  np.linalg.cond(K+W)
    if cnum>1e16:
        print 'Warning high condition number results may suffer'
 
    delta_pos = ref_dis+np.zeros([Nc,d])
    tol_pos = ref_dis/10
    best_E = np.matmul(K,Co_GP)
    Nc_search = us_ind.shape[0]
    
    
    delta_pos = delta_pos[us_ind,:]
    best_E = best_E[us_ind,:]
    best_pos = best_pos[us_ind,:]
    
    
    while np.max(delta_pos)>tol_pos:
        for jp in range(Nc_search):
            for jd in range(d):
                did_move = 0
                c_pos = np.copy(best_pos)
                c_pos[jp,jd]= min(best_pos[jp,jd]+delta_pos[jp,jd],1.0)
                    
                Kc = np.exp(2.0*hyp_parm[1])*cmgp(np.reshape(c_pos[jp,:],[1,d]),mol_pos,np.exp(-2.0*hyp_parm[0]))
                c_E = np.matmul(Kc,Co_GP)
                if c_E<best_E[jp]:
                    best_E[jp] = c_E
                    best_pos = np.copy(c_pos)
                    did_move = 1
                    
                c_pos = np.copy(best_pos)
                c_pos[jp,jd]=max(best_pos[jp,jd]-delta_pos[jp,jd],0.0)
                
                        
                Kc = np.exp(2.0*hyp_parm[1])*cmgp(np.reshape(c_pos[jp,:],[1,d]),mol_pos,np.exp(-2.0*hyp_parm[0]))
                c_E = np.matmul(Kc,Co_GP)
                if c_E<best_E[jp]:
                    best_E[jp] = c_E
                    best_pos = np.copy(c_pos)
                    did_move = 1
                    
                if did_move == 1:
                    delta_pos[jp,jd] = 1.25*delta_pos[jp,jd]
                else:
                    delta_pos[jp,jd] = 0.75*delta_pos[jp,jd]
                
    ## Remove Points on Convex Haul
    
    if d>1:
        us_ind = np.argwhere(T.find_simplex(best_pos)>-1)
        if us_ind.shape[0] == 0:
            us_ind = np.argwhere(T.simplices==np.argmin(E))
            c_E = np.zeros([us_ind.shape[0],1])
            for jp in range(us_ind.shape[0]): 
                c_E[jp] = np.mean(E[T.simplices[us_ind[jp,0],:],:])
            
            best_ind = np.argmin(c_E[:,0])
            best_E = np.reshape(c_E[best_ind,:],[1,1])
            best_pos = np.mean(mol_pos[T.simplices[us_ind[best_ind,0],:],:],axis=0)
        else:
            best_E = best_E[us_ind[:,0],:]
            best_pos = best_pos[us_ind[:,0],:]      
    else:
        us_ind = np.union1d(np.argwhere(best_pos[:,0]>min(mol_pos)),np.argwhere(best_pos[:,0]<max(mol_pos)))
        us_ind = np.reshape(us_ind,[us_ind.shape[0],1])
        best_E = best_E[us_ind[:,0],:]
        best_pos = best_pos[us_ind[:,0],:]
        
    best_E = np.reshape(best_E,[best_E.shape[0],1])    
    best_pos = np.reshape(best_pos,[best_E.shape[0],d]) 
    
    ## Remove in same simplex ##
    if d>1:
        jp = 0
        ind_T = T.find_simplex(best_pos)
        ind_T = np.reshape(ind_T,[ind_T.shape[0],1])
        while jp<best_pos.shape[0]:
            us_ind = np.argwhere(ind_T==ind_T[jp])
            if us_ind.shape[0]>1:
                min_ind = np.argmin(best_E[us_ind[:,0],0])
                best_E[us_ind[0,0],0] = best_E[us_ind[min_ind,0],0]
                best_pos[us_ind[0,0],:] = best_pos[us_ind[min_ind,0],:]
                best_E = np.delete(best_E,us_ind[1:us_ind.shape[0],0],0)
                best_pos = np.delete(best_pos,us_ind[1:us_ind.shape[0],0],0)
                ind_T = np.delete(ind_T,us_ind[1:us_ind.shape[0],0],0)
                
            jp += 1
            
        jp = 0
        while jp<best_pos.shape[0]:
            c_dis = np.sqrt(np.sum((np.repeat(np.reshape(best_pos[jp,:],[1,d]),\
                 best_pos.shape[0],axis=0)-best_pos)**2,axis=1))
            us_ind = np.argwhere(c_dis<ref_dis)
            if us_ind.shape[0]>1:
                min_ind = np.argmin(best_E[us_ind[:,0],0])
                best_E[us_ind[0,0],0] = best_E[us_ind[min_ind,0],0]
                best_pos[us_ind[0,0],:] = best_pos[us_ind[min_ind,0],:]
                best_E = np.delete(best_E,us_ind[1:us_ind.shape[0],0],0)
                best_pos = np.delete(best_pos,us_ind[1:us_ind.shape[0],0],0)
                
            jp += 1
        
        ind_T = T.find_simplex(best_pos)
        ind_T = np.reshape(ind_T,[ind_T.shape[0],1])
        
        if best_E.shape[0] > 1:
            ind_sort = np.argsort(best_E[:,0])
            ind_rm = np.argwhere(best_E[ind_sort,0]>np.mean(E)-np.std(E))
            if ind_rm.shape[0]>0:
                if ind_rm[0,0] == 0:
                    ind_rm = 1
                else:
                    ind_rm = min(ind_rm[0,0],d+1)
                    
                ind_sort = ind_sort[0:ind_rm]
            
            best_E = best_E[ind_sort,:]
            best_pos = best_pos[ind_sort,:]
            ind_T = ind_T[ind_sort,:]
                
        neig_data_ind = T.simplices[ind_T[:,0],:]
    else:
        
        jp = 0
        while jp<best_pos.shape[0]:
            c_dis = np.sqrt(np.sum((np.repeat(np.reshape(best_pos[jp,:],[1,d]),\
                 best_pos.shape[0],axis=0)-best_pos)**2,axis=1))
            us_ind = np.argwhere(c_dis<10*ref_dis)
            if us_ind.shape[0]>1:
                min_ind = np.argmin(best_E[us_ind[:,0],0])
                best_E[us_ind[0,0],0] = best_E[us_ind[min_ind,0],0]
                best_pos[us_ind[0,0],:] = best_pos[us_ind[min_ind,0],:]
                best_E = np.delete(best_E,us_ind[1:us_ind.shape[0],0],0)
                best_pos = np.delete(best_pos,us_ind[1:us_ind.shape[0],0],0)
                
            jp += 1
            
        if best_E.shape[0] > 1:
            ind_sort = np.argsort(best_E[:,0])
            ind_rm = np.argwhere(best_E[ind_sort,0]>np.mean(E)-np.std(E))
            if ind_rm.shape[0]>0:
                if ind_rm[0,0] == 0:
                    ind_rm = 1
                else:
                    ind_rm = ind_rm[0,0]
                    
                ind_sort = ind_sort[0:ind_rm]
            
            best_E = best_E[ind_sort,:]
            best_pos = best_pos[ind_sort,:]    
        
        neig_data_ind = np.zeros([best_pos.shape[0],2], dtype=int)
        for jp in range(best_pos.shape[0]):
            c_dis = np.sqrt(np.sum((np.repeat(np.reshape(best_pos[jp,:],[1,d]),\
                 mol_pos.shape[0],axis=0)-mol_pos)**2,axis=1))
            sort_ind = np.argsort(c_dis)
            neig_data_ind[jp,0]=sort_ind[0]
            neig_data_ind[jp,1]=sort_ind[1]
            
    ## Estimate Error
    K_res = np.exp(2.0*hyp_parm[1])*cmgp(best_pos,mol_pos,np.exp(-2.0*hyp_parm[0]))
    best_dE = np.reshape(2*abs(np.exp(2.0*hyp_parm[1])-\
        np.diag(np.matmul(K_res,np.matmul(Kw_inv,np.transpose(K_res))))),[best_pos.shape[0],1])
          
    return best_pos,best_E,best_dE,neig_data_ind


###############################################################################
def find_new_points(mol_pos,best_pos,delta_r,Nc0):

    """
    Estimate optimal new points for Gaussian Process
    
    
    Parameters
    ----------
    best_pos:   Estimated local minimum
    best_E:     Estimated local energies
    delta_r:    Maximum distance of new search space
    Nc0:        Number of points to add
    
    Returns
    -------
    new_points: New sample points for GP
    """
    
    [Nc,d] = mol_pos.shape
    
    if best_pos.shape[0]>Nc0:
        best_pos = best_pos[0:Nc0,:]
    
    new_points = np.zeros([Nc0,d])
    n_clusters = best_pos.shape[0]
    ind_clusters = np.round(np.linspace(0,n_clusters,n_clusters+1)*Nc0/(n_clusters))
    ind_clusters = ind_clusters.astype(int)  
    shyper_cube = np.zeros([d,2])
    n_min_d = 5
    
    for jc in range(n_clusters):
        n_pnts = ind_clusters[jc+1]-ind_clusters[jc]
        shyper_cube[:,0] = np.maximum(best_pos[jc,:]-delta_r,0)
        shyper_cube[:,1] = np.minimum(best_pos[jc,:]+delta_r,1)
        c_dis = min(np.sqrt(np.sum((np.repeat(np.reshape(best_pos[jc,:],[1,d]),\
             mol_pos.shape[0],axis=0)-mol_pos)**2,axis=1)))
        
        test_points = np.repeat(np.reshape(shyper_cube[:,1]-shyper_cube[:,0],[1,d]) \
            ,n_pnts,axis=0)*lhs(d,n_pnts)+np.repeat(np.reshape(shyper_cube[:,0],[1,d]),n_pnts,axis=0)
        
        D =_fdist(np.append(mol_pos,new_points[0:ind_clusters[jc]+n_pnts,:],axis=0))   
        b_mes_v = np.min(D+np.max(D)*np.eye(Nc+ind_clusters[jc]+n_pnts),axis=1)
        
        ## Only includes best_pos if it is significantly different from mol_pos ##
        if c_dis>min(b_mes_v):
            test_points[0,:] = best_pos[jc,:]
            b_mes_v[0] = c_dis
            
        b_mes = np.sum(b_mes_v)   
        new_points[ind_clusters[jc]:ind_clusters[jc]+n_pnts,:] = test_points
        
        for jit in range(n_min_d):
            cnew_points = np.copy(new_points)
            test_points = np.repeat(np.reshape(shyper_cube[:,1]-shyper_cube[:,0],[1,d]) \
                ,n_pnts,axis=0)*lhs(d,n_pnts)+ \
                np.repeat(np.reshape(shyper_cube[:,0],[1,d]),n_pnts,axis=0)
                                
            D =_fdist(np.append(mol_pos,cnew_points[0:ind_clusters[jc]+n_pnts,:],axis=0))    
            c_mes_v = np.min(D+np.max(D)*np.eye(Nc+ind_clusters[jc]+n_pnts),axis=1)
            ## Only includes best_pos if it is significantly different from mol_pos ##
            if c_dis>min(c_mes_v):
                test_points[0,:] = best_pos[jc,:]
                c_mes_v[0] = c_dis
                
            cnew_points[ind_clusters[jc]:ind_clusters[jc]+n_pnts,:] = test_points
            c_mes = np.sum(c_mes_v)
            
            if c_mes>b_mes:
                new_points = np.copy(cnew_points)
                b_mes = c_mes
            
    
    return new_points
    
###############################################################################

def multiscale_E(E,dE):

    """
    Multiscale normalization of E
    
    
    Parameters
    ----------
    E:              Energy
    dE:             Uncertainty in Energy
    
    Returns
    -------
    E_parm:    Scale parameters
    E_scaled:       Scaled energies
    dE_scaled:      Scaled uncertainties in energies
    """
    
    E_parm = np.zeros(3)
    E_parm[0] = min(E)
    us_ind = np.argsort(E, axis=0)
    E_parm[1] = E[us_ind[np.int(np.round(E.shape[0]/4.0)),0],0]
    E_parm[1] = min(max(E_parm[1],E_parm[0]+100*np.mean(dE[np.argwhere(E<=E_parm[1])])),max(E))
    if E_parm[1] == E_parm[0]:
        E_parm[1] = min(E_parm[0]+1,max(E))
        
    E_scaled = (np.copy(E)-E_parm[0])/(E_parm[1]-E_parm[0])
    dE_scaled = np.copy(dE)/(E_parm[1]-E_parm[0])
    us_ind = np.argwhere(E_scaled[:,0]>1.0)
    
    if us_ind.shape[0]>0:
        dE_scaled[us_ind[:,0],:] = np.log(1+dE_scaled[us_ind[:,0],:]/E_scaled[us_ind[:,0],:])
        E_scaled[us_ind[:,0],:] = np.log(E_scaled[us_ind[:,0],:])+1
        E_parm[2] = max(E_scaled)
        E_scaled[us_ind[:,0],:] = 9.0*(E_scaled[us_ind[:,0],:]-1)/(E_parm[2]-1) +1
        dE_scaled[us_ind[:,0],:] = 9.0*dE_scaled[us_ind[:,0],:]/(E_parm[2]-1)
        dE_scaled[us_ind[:,0],:] = np.maximum(dE_scaled[us_ind[:,0],:],np.mean(dE_scaled[np.argwhere(E<=E_parm[1])]))


        
    return E_parm,E_scaled,dE_scaled

###############################################################################

def inv_multiscale_E(E_parm,E_scaled,dE_scaled):

    """
    Reverses Multiscale normalization of E
    
    
    Parameters
    ----------
    E_parm:         Scale parameters
    E_scaled:       Scaled energies
    
    Returns
    -------
    E:          Energy
    """
    E = np.copy(E_scaled)
    dE = np.copy(dE_scaled)
    if E_parm[2]>0:
        us_ind = np.argwhere(E[:,0]>1.0)
        if us_ind.shape[0]>0:
            E[us_ind[:,0],:] = np.exp((E[us_ind[:,0],:]-1)*(E_parm[2]-1)/9.0)
    
    E = E*(E_parm[1]-E_parm[0])+E_parm[0]
    dE = dE*(E_parm[1]-E_parm[0])
        
    return E,dE

###############################################################################
    
def inv_norm_pos(mol_pos_norm,hyper_cube):

    """
    Reverses normalization of molecular positions to unit hyper cube
    
    
    Parameters
    ----------
    hyper_cube:     Hyper cube
    mol_pos_norm:   Normalized GP centers
    
    
    Returns
    -------
    mol_pos:        GP centers
    """
    
    [n,d] = mol_pos_norm.shape
    mol_pos = np.repeat(np.reshape(hyper_cube[:,1]-hyper_cube[:,0],[1,d]) \
                ,n,axis=0)*mol_pos_norm+ \
                np.repeat(np.reshape(hyper_cube[:,0],[1,d]),n,axis=0)

        
    return mol_pos


###############################################################################
    
def norm_pos(mol_pos,hyper_cube):

    """
    Reverses normalization of molecular positions to unit hyper cube
    
    
    Parameters
    ----------
    hyper_cube:     Hyper cube
    mol_pos_norm:   Normalized GP centers
    
    
    Returns
    -------
    mol_pos:        GP centers
    """
    
    [n,d] = mol_pos.shape
    mol_pos_norm = (mol_pos-np.repeat(np.reshape(hyper_cube[:,0],[1,d]),n,axis=0))/ \
        np.repeat(np.reshape(hyper_cube[:,1]-hyper_cube[:,0],[1,d]),n,axis=0)

        
    return mol_pos_norm

###############################################################################


def gp_opt_setup(GP_info):
    """
    Setup parameters and variables for GP
    
    
    Parameters
    ----------
    GP_info:    GP information 
    E_b:        Scale parameters
    E_scaled:   Scaled energies
    
    
    Returns
    -------
    E:   
    """
    
    ## Define Empty arrays for energies, positions, and uncertainties
    E = []
    dE = []
    mol_pos = []

    best_pos =[]
    best_E = []
    
    d = GP_info.hyper_cube.shape[0]
    try:
        GP_info.n_its
    except:
        GP_info.n_its = 5               # Number of maximum iterations
        
    try:
        GP_info.do_vis
    except:
        GP_info.do_vis = 0              # Visualize Error
        
        
    try:
        GP_info.Nc0
    except:
        GP_info.Nc0 = (d+2)*(d+1)       # Number of points to add per iteration  
        
    try:
        GP_info.delta_r
    except:
        GP_info.delta_r = 0.5           # Initial Radius of Convergence     
    
    try:
        GP_info.delta_r_s
    except:
        GP_info.delta_r_s = 1.25         # Radius shrink factor
    
    
    try:
        GP_info.n_min_d
    except:
        GP_info.n_min_d = 10             # Number of times to resample LHC 
    
    return E,dE,mol_pos,best_pos,best_E,GP_info

###############################################################################

def gp_opt(E,dE,mol_pos,GP_info):
    """
    Setup parameters and variables for GP
    
    
    Parameters
    ----------
    best_pos:   Estimated local minimum
    best_E:     Estimated local energies
    mol_pos:    GP centers
    GP_info:    GP information 
    
    E:          Energy for GP centers
    dE:         Uncertainty in energy for GP centers
    mol_pos:    GP coefficents
    GP_info:    GP information 
    
    Returns
    -------
    new_points:     New points to sample
    best_pos:       Estimated local minimum
    best_E:         Estimated local energies
    GP_info:        GP information    
    """
    
    ## Define Empty arrays for energies, positions, and uncertainties
    d = GP_info.hyper_cube.shape[0]
    
    if GP_info.jits == 0:
        mol_pos_norm = lhs(d,GP_info.Nc0,GP_info.n_min_d)
        new_points = inv_norm_pos(mol_pos_norm,GP_info.hyper_cube)
        best_pos = []
        best_E =[]
        best_dE = []
        neig_data = []
    else:
        GP_info.delta_r = GP_info.delta_r/GP_info.delta_r_s
        mol_pos_norm = norm_pos(mol_pos,GP_info.hyper_cube)
        
            
        E_parm,E_scaled,dE_scaled = multiscale_E(E,dE)
        ## Best Approximation of Optimization Surface ##
        hyp_parm,Co_GP = find_hyperparameters(mol_pos_norm,E_scaled,dE_scaled)
        ## Find Global Min ##
        best_pos_norm,best_E_scaled,best_dE_scaled,neig_data_ind = find_min(Co_GP,mol_pos_norm,E_scaled,dE_scaled,hyp_parm)
        best_E,best_dE = inv_multiscale_E(E_parm,best_E_scaled,best_dE_scaled) 
        best_pos = inv_norm_pos(best_pos_norm,GP_info.hyper_cube)
        
        new_points = find_new_points(mol_pos_norm,best_pos_norm,GP_info.delta_r,GP_info.Nc0)
        new_points = inv_norm_pos(new_points,GP_info.hyper_cube)
        
        neig_data = E[neig_data_ind,0]
        neig_data = np.append(neig_data,dE[neig_data_ind,0],axis=0)
        
        if GP_info.do_vis == 1:
            
            
            neig_bool = neig_data[0:best_E.shape[0],:]- \
                2*neig_data[best_E.shape[0]:2*best_E.shape[0],:]
            neig_bool = neig_bool<np.repeat(best_E,d+1,axis=1)
            neig_percent = np.zeros(neig_bool.shape)
            for jx in range(neig_bool.shape[0]):
                for jy in range(neig_bool.shape[1]):
                    neig_percent[jx,jy] = float(neig_bool[jx,jy])
                
            neig_percent = 100*np.mean(neig_percent,axis=1)
            x = np.zeros((d+1)*best_E.shape[0])
            y = np.zeros((d+1)*best_E.shape[0])
            dy = np.zeros((d+1)*best_E.shape[0])
            Ex = np.zeros(best_E.shape[0])
            for jx in range(best_E.shape[0]):
                x[jx*(d+1):(jx+1)*(d+1)] = jx+1
                y[jx*(d+1):(jx+1)*(d+1)] = neig_data[jx,:]
                dy[jx*(d+1):(jx+1)*(d+1)]= neig_data[jx+best_E.shape[0],:]
                Ex[jx] = jx+1
                
            plt.figure()
            plt.errorbar(x, y, yerr=dy, fmt='o',label='sampled point with uncertainty')
            plt.errorbar(Ex,best_E[:,0], yerr=best_dE[:,0],fmt='+',label='estimated minumum')
            plt.legend()
            plt.title('Neigborhood Information')
            plt.xlabel('Local Minimum')
            plt.ylabel('E')
            plt.show()
            if d == 1:
                npts = 100
                plot_points = np.reshape(np.linspace(GP_info.hyper_cube[0,0],\
                    GP_info.hyper_cube[0,1],num=npts),[npts,1])
                #Et = np.reshape(opt_surface_f(plot_points,GP_info.opt_fun),[npts,1])
                
                plot_points_norm = norm_pos(plot_points,GP_info.hyper_cube)
                
                gp_E = np.matmul(np.exp(2.0*hyp_parm[1])*cmgp(plot_points_norm,\
                    mol_pos_norm,np.exp(-2.0*hyp_parm[0])),Co_GP)
                gp_E_inv,gp_E_dE = inv_multiscale_E(E_parm,gp_E,0*gp_E) 
                plt.figure()
                #plt.plot(plot_points,gp_E_inv,'y-',plot_points,Et,'g-',mol_pos,E,'b.',best_pos, best_E,'ro')
                plt.plot(plot_points,gp_E_inv,'y-',mol_pos,E,'b.',best_pos, best_E,'ro')
                plt.show()
            
            
         
        
    return new_points,best_pos,best_E,best_dE,GP_info,neig_data

###############################################################################
    

def gp_opt_example():
    """
    Example using Gaussian Process optimization for noisy data and multiple mins

    """
    plt.close('all')

    ## Required to define is hypercube  ##


    d = 3                                       # Dimension of problem
    opt_fun = 5                                 # Define type of optimization function

    d = 2
    opt_fun = 4

    ## Parameters Used for Optimization Function Defining Range of Errors in Energy##
    sigma_a = 1e-4;
    sigma_b = 5e-2;


    ## Set Hypercube for each type of optimization function and check for proper dim ##
    GP_info = GParmFlags()

    if opt_fun == 1:     
        ## Do LJ Potential ##
        GP_info.hyper_cube = 2*np.ones([d,2])       # Upper and lower bound for each dimension
        GP_info.hyper_cube[:,0] = .2*GP_info.hyper_cube[:,0]
        if d>1:
            GP_info.hyper_cube[1:d,0] = 0
            GP_info.n_its = 5

            if d == 3:
                print 'True min: f(x) = -3, x = [1 .5 sqrt(3)/2]'
        else:
            print 'True min: f(x) = -1, x = 1'

    elif opt_fun == 2: 
        d = 2
        GP_info.hyper_cube = np.ones([d,2])       # Upper and lower bound for each dimension
        GP_info.hyper_cube[0,0] = 2
        GP_info.hyper_cube[1,0] = -.5
        GP_info.hyper_cube[0,1] = 4
        GP_info.hyper_cube[1,1] = 1

        print 'Only defined for d=2'
        print '2D function, global minimum at (3,.5)'   

    elif opt_fun == 3: 

        if d>7:
            d = 2
        elif d<3:
            d = 3


        GP_info.Nc0 = 2*(d+2)*(d+1)    
        GP_info.hyper_cube = 1.5*np.ones([d,2])       # Upper and lower bound for each dimension
        GP_info.hyper_cube[:,0] = 0.5     

        print 'Only defined for d=3-7'
        print 'Global minimum at (1,1,...,1)'

    elif opt_fun == 4: 
        d = 2
        GP_info.hyper_cube = np.ones([d,2])       # Upper and lower bound for each dimension
        GP_info.hyper_cube[0,0] = -1
        GP_info.hyper_cube[1,0] = -2
        GP_info.hyper_cube[0,1] = 1
        GP_info.hyper_cube[1,1] = 0

        print 'Only defined for d=2'
        print '2D function, global minimum at (0,-1)'  

    else:
        GP_info.hyper_cube = np.ones([d,2])/2       # Upper and lower bound for each dimension
        GP_info.hyper_cube[:,0] = -GP_info.hyper_cube[:,0]

        GP_info.n_its = 4

        print 'True min: f(x) = -1, x = 0'

    ## Setup Gausian Process Structure ##
    try:
        GP_info.hyper_cube
        E,dE,mol_pos,best_pos,best_E,GP_info = gp_opt_setup(GP_info)
    except:
        print 'Please Define GP_info.hyper_cube'
        exit

    ## Do Gausian Process Iteration ##

    for jits in range(GP_info.n_its):
        GP_info.jits = jits
        new_points,best_pos,best_E,best_dE,GP_info,neig_data = gp_opt(E,dE,mol_pos,GP_info) 
        if jits < GP_info.n_its-1:
            ## Simulate random uncertainty in optimization function
            new_dE = (sigma_b-sigma_a)*np.random.rand(GP_info.Nc0,1)+sigma_a
            new_E = opt_surface_f(new_points,opt_fun)+new_dE*np.random.randn(GP_info.Nc0, 1)
            if jits == 0:
                mol_pos = np.copy(new_points)
                E = np.copy(new_E)
                dE = np.copy(new_dE)
            else:
                mol_pos = np.append(mol_pos,new_points,axis=0)
                E = np.append(E,new_E,axis=0)
                dE = np.append(dE,new_dE,axis=0)


    print 'Best Guess in position: '
    print best_pos
    print 'Best Guess in energy: '
    print best_E
#end def gp_opt_example        


###############################################################################




class GPState(DevBase):
    """
    Represents internal state of a Gaussian Process iteration.

    The GPState class is the primary interface to the procedural GP optimization
    approach represented by the gp_opt function and sub-functions.  This class 
    represents the state of the GP at each iteration and allows for the state 
    to be stored and loaded from disk.  This allows for interruptible GP 
    optimizations, which is useful when interfacing to target functions that 
    are expensive to evaluate, such as a QMC energy surface that is evaluated 
    on an HPC machine with Nexus driven workflows.

    GPState is not intended for public instantiation. GaussianProcessOptimizer
    class is the sole intended owner.  Read only access (intended) is granted 
    by GaussianProcessOptimizer to the user when a GPState instance is passed in 
    to the user-provided energy_function.  See GaussianProcessOptimizer.optimize.
    """

    def __init__(self,
                 param_lower = None,
                 param_upper = None,
                 niterations = None,
                 save_path   = './',
                 ):
        """
        Setup initial GP data.

        Parameters
        ----------
        param_lower : Lower bounds of the parameter search domain (list-like).
        param_upper : Upper bounds of the parameter search domain (list-like).
        niterations : Number of GP iterations to perform.
        save_path   : (Optional) directory where state snapshots will be saved.
        """
        if isinstance(param_lower,GPState):
            o = param_lower.copy()
            self.__dict__ = o.__dict__
            return
        #end if
        dim     = -1
        Pdomain = None
        E       = None
        dE      = None
        P       = None
        Popt    = None
        Eopt    = None
        GP_info = None
        if param_lower is not None:
            dim = len(param_lower)
            param_lower = array(param_lower)
            param_upper = array(param_upper)
            hyper_cube = zeros([dim,2],dtype=float)
            hyper_cube[:,0] = param_lower
            hyper_cube[:,1] = param_upper
            Pdomain = hyper_cube.copy()
            # Setup Gaussian Process data structure 
            GP_info = GParmFlags()
            GP_info.n_its = niterations
            GP_info.Nc0 = 2*(dim+2)*(dim+1)
            GP_info.hyper_cube = hyper_cube
            E,dE,P,Popt,Eopt,GP_info = gp_opt_setup(GP_info)
        #end if
        self.save_path = save_path
        self.dim       = dim
        self.Pdomain   = Pdomain
        self.E         = E
        self.dE        = dE
        self.P         = P
        self.Popt      = Popt
        self.Eopt      = Eopt
        self.GP_info   = GP_info
        self.iteration = -1
        self.Pnew      = None
        self.neig_data = None
        self.Popt_hist = []
        self.Eopt_hist = []
    #end def __init__


    def gp_inputs(self):
        """
        Returns inputs required by the gp_opt function.
        """
        return self.list('E dE P GP_info'.split())
    #end def gp_inputs


    def update_gp(self,gp_opt_outputs):
        """
        Stores outputs from the gp_opt_function.
        """
        self.Pnew      = gp_opt_outputs[0]
        self.Popt      = gp_opt_outputs[1]
        self.Eopt      = gp_opt_outputs[2]
        self.dEopt     = gp_opt_outputs[3]
        self.GP_info   = gp_opt_outputs[4]
        self.neig_data = gp_opt_outputs[5]
    #end def update_gp
    

    def parameters(self):
        """
        Returns parameters generated/sampled in the current iteration.
        """
        return self.Pnew.copy()
    #end def parameters


    def sample_parameters(self):
        """
        Calls the gp_opt function to obtain LHC parameter samples and predicted
        optimal parameters and energies.
        """
        self.GP_info.jits = self.iteration
        self.update_gp(gp_opt(*self.gp_inputs()))
        return self.parameters()
    #end def sample_parameters


    def update_energies(self,energies,errors):
        """
        Appends parameter and energy+error data to internal state arrays.
        """

        if self.iteration==0:
            self.P  = self.Pnew.copy()
            self.E  = energies.copy()
            self.dE = errors.copy()
        else:
            self.P  = append(self.P ,self.Pnew,axis=0)
            self.E  = append(self.E ,energies ,axis=0)
            self.dE = append(self.dE,errors   ,axis=0)
        #end if
    #end def update_energies


    def optimal_parameters(self):
        """
        Returns predicted optimal parameters for the current iteration.
        """
        return self.Popt.copy()
    #end def optimal_parameters


    def optimal_energies(self):
        """
        Returns predicted optimal energies for the current iteration.
        """
        return self.Eopt.copy()
    #end def optimal_energies


    def update_optimal_trajectory(self):
        """
        Stores a history of predicted optimal parameters/energies across all 
        iterations.
        """
        Popt = self.optimal_parameters()
        Eopt = self.optimal_energies()
        self.Popt_hist.append(Popt)
        self.Eopt_hist.append(Eopt)
        return Popt,Eopt
    #end def update_optimal_trajectory


    def optimal_trajectory(self):
        """
        Returns the predicted optimal parameters/energies across all iterations.
        """
        return self.Popt_hist,self.Eopt_hist
    #end def optimal_trajectory


    def state_filepath(self,label=None,path=None,prefix='gp_state',filepath=None):
        """
        Returns the filepath used to store the state of the current iteration.
        """
        if filepath is not None:
            return filepath
        #end if
        if path is None:
            path = self.save_path
        #end if
        filename = prefix+'_iteration_{0}'.format(self.iteration)
        if label is not None:
            filename += '_'+label
        #end if
        filename += '.p'
        filepath = os.path.join(path,filename)
        return filepath
    #end def state_filepath


    def save(self,*args,**kwargs):
        """
        Saves the state of the current iteration to disk.
        """
        filepath = self.state_filepath(*args,**kwargs)
        path,file = os.path.split(filepath)
        if not os.path.exists(path):
            os.makedirs(path)
        #end if
        s = obj(self)
        s.save(filepath)
    #end def save


    def load(self,*args,**kwargs):
        """
        Loads the state of the current iteration from disk.
        """
        filepath = self.state_filepath(*args,**kwargs)
        s = obj()
        s.load(filepath)
        self.transfer_from(s)
    #end def load


    def attempt_load(self,*args,**kwargs):
        """
        Attempts to load current iteration state from disk.
        """
        loaded = False
        filepath = self.state_filepath(*args,**kwargs)
        if os.path.exists(filepath):
            self.load(filepath=filepath)
            loaded = True
        #end if
        return loaded
    #end def attempt_load
#end class GPState



class GaussianProcessOptimizer(DevBase):
    """
    Offers the main user interface to Gaussian Process optimization.

    Apart from instantiation, the user should only need to call the "optimize" 
    member function.  

    Examples
    --------
    ::

        # minimize a 2D parabola over the domain -1 < x,y < 1
        def x2(x,s):
            E = (x**2).sum(axis=1,keepdims=1)
            dE = 0*E
            return E,dE
        #end def x2    
        gpo = GaussianProcessOptimizer(5,[-1,-1],[1,1],x2)
        P0,E0 = gpo.optimize()
        print 'optimal parameters:',P0[-1]
        print 'optimal energy    :',E0[-1]
    """

    allowed_modes = 'stateless stateful'.split()

    def __init__(self,
                 niterations     = None,
                 param_lower     = None,
                 param_upper     = None,
                 energy_function = None,
                 state_path      = './opt_gp_states',
                 output_path     = './opt_gp_output',
                 mode            = 'stateful',
                 verbose         = True,
                 exit_on_save    = False,
                 ):
        """
        Parameters
        ----------
        niterations     : Number of GP iterations to perform.
        param_lower     : Lower bounds of the parameter search domain (list-like).
        param_upper     : Upper bounds of the parameter search domain (list-like).
        energy_function : Target function for optimization.  Must accept parameters
                          and a GPState instance as arguments and return energies 
                          and errorbars for the energies (function-like).
        state_path      : (Optional) directory where state snapshots will be saved.
        output_path     : (Optional) directory where parameter/energy data will be 
                          saved for each iteration.
        mode            : (Optional) selects whether or not to save state during 
                          the GP optimization.  Can be "stateless" or "stateful",
                          default is "stateful".
        verbose         : (Optional) selects whether detailed information is 
                          printed (bool).
        exit_on_save    : (Optional) testing parameter that selects whether to 
                          exit each time state is saved to disk.  This enables 
                          the user to check that the optimizer can suspend and 
                          resume properly.
        """
        if niterations is None:
            self.error('niterations is required')
        #end if
        if param_lower is None:
            self.error('param_lower is required')
        #end if
        if param_upper is None:
            self.error('param_upper is required')
        #end if
        if energy_function is None:
            self.error('energy_function is required')
        #end if
        self.niterations     = niterations
        self.state_path      = state_path
        self.output_path     = output_path
        self.state           = GPState(param_lower,param_upper,niterations,state_path)
        self.energy_function = energy_function
        self.mode            = mode
        self.verbose         = verbose
        self.exit_on_save    = exit_on_save
        if self.mode not in self.allowed_modes:
            self.error('mode {0} is not allowed\nallowed modes are: {1}'.format(self.mode,self.allowed_modes))
        #end if
        self.state_hist = obj()
    #end def __init__


    def optimize(self):
        """
        Runs the main optimization cycle.
        """
        self.vlog('Beginning Gaussian Process optimization')
        res = None
        if self.mode=='stateless':
            res = self.optimize_stateless()
        elif self.mode=='stateful':
            res = self.optimize_stateful()
        else:
            self.error('cannot optimize, mode "{0}" has not been implemented'.format(self.mode))
        #end if
        self.vlog('Gaussian Process optimization finished')
        return res
    #end def optimize


    def optimize_stateless(self):
        """
        Optimizes energy_function without saving state.
        """
        state = self.state
        for i in range(self.niterations+1):
            state.iteration+=1
            self.vlog('iteration {0}'.format(state.iteration),n=1)
            self.vlog('sampling parameters',n=2)
            params = state.sample_parameters()
            if i>0:
                self.vlog('updating predicted optimal params/energies',n=2)       
                state.update_optimal_trajectory()
            #end if
            if i<self.niterations:
                self.vlog('calculating new energies',n=2)
                energies,errors = self.energy_function(params,state)
                self.vlog('incorporating new energies',n=2)
                state.update_energies(energies,errors)
            #end if
            self.state_hist[i] = state.copy()
        #end for
        return state.optimal_trajectory()
    #end def optimize_stateless


    def optimize_stateful(self):
        """
        Optimizes energy_function, saving state each iteration.

        All load*/save*/write* functions are ultimately called by this function.
        """
        state = self.state
        for i in range(self.niterations+1):
            state.iteration+=1
            self.vlog('iteration {0}'.format(state.iteration),n=1)
            if not self.load_start():
                self.save_start()
            #end if
            if not self.load_parameters():
                params = state.sample_parameters()
                if i>0:
                    Popt,Eopt = state.update_optimal_trajectory()
                    self.save_optimal(Popt,Eopt)
                #end if
                self.save_parameters(params)
            else:
                params = state.parameters()
            #end if
            if i<self.niterations:
                if not self.load_energies():
                    energies,errors = self.energy_function(params,state)
                    state.update_energies(energies,errors)
                    self.save_energies(energies,errors)
                #else: # wrong, temporary, only for completed nexus runs
                #    energies,errors = self.energy_function(params,state)
                #end if
            #end if
            self.state_hist[i] = state.copy()            
        #end for
        return state.optimal_trajectory()
    #end def optimize_stateful


    def vlog(self,msg,n=0,indent='  '):
        """
        Prints a message if verbose=True.
        """
        if self.verbose:
            self.log(msg,indent=n*indent)
        #end if
    #end def vlog


    def make_output_dir(self,iter='cur'):
        """
        Creates a directory for text file parameter/energy output for the 
        current iteration.
        """
        i = self.state.iteration
        if iter=='cur':
            None
        elif iter=='prev':
            i-=1
        else:
            self.error('invalid iteration requested for save: {0}'.format(iter))
        #end if
        path = os.path.join(self.output_path,'iteration_{0}'.format(i))
        if not os.path.exists(path):
            os.makedirs(path)
        #end if
        return path
    #end def make_output_dir


    def write_param_file(self,param_file,Pset,iter='cur'):
        """
        Writes a parameter file to the current iteration's output directory.
        """
        path = self.make_output_dir(iter)
        filepath = os.path.join(path,param_file)
        self.vlog('writing parameters to {0}'.format(filepath),n=3)
        nsamples,nparams = Pset.shape
        s = ''
        for i in range(nsamples):
            for j in range(nparams):
                s += ' {0: 16.12f}'.format(Pset[i,j])
            #end for
            s += '\n'
        #end for
        f = open(filepath,'w')
        f.write(s)
        f.close()
    #end def write_param_file


    def write_energy_file(self,energy_file,energies,errors=None,iter='cur'):
        """
        Writes an energy file to the current iteration's output directory.
        """
        path = self.make_output_dir(iter)
        filepath = os.path.join(path,energy_file)
        self.vlog('writing energies to {0}'.format(filepath),n=3)
        energies = energies.ravel()
        s = ''
        if errors is None:
            for e in energies:
                s += '{0: 16.12f}\n'.format(e)
            #end for
        else:
            errors = errors.ravel()
            for e,ee in zip(energies,errors):
                s += '{0: 16.12f}  {1: 16.12f}\n'.format(e,ee)
            #end for
        #end if
        f = open(filepath,'w')
        f.write(s)
        f.close()
    #end ef write_energy_file


    def load_state(self,label):
        """
        Loads state information for the current iteration from a file matching 
        "label".
        """
        savefile = self.state.state_filepath(label)
        self.vlog('attempting to load state file {0}'.format(savefile),n=3)
        loaded = self.state.attempt_load(label)
        if loaded:
            self.vlog('state file load successful',n=3)
        else:
            self.vlog('state file is not present',n=3)
        #end if
        return loaded
    #end def load_state


    def save_state(self,label):
        """
        Saves state information for the current iteration to a file matching 
        "label".
        """
        savefile = self.state.state_filepath(label)
        self.vlog('saving state file {0}'.format(savefile),n=3)
        self.state.save(label)
        if self.exit_on_save:
            self.vlog('exiting',n=3)
            exit()
        #end if
    #end def save_state


    def load_start(self):
        """
        Attempt to load state for startup portion of current iteration.
        """
        self.vlog('startup',n=2)
        return self.load_state('start')
    #end def load_start


    def save_start(self):
        """
        Save state for startup portion of current iteration.
        """
        self.save_state('start')
    #end def save_start


    def load_parameters(self):
        """
        Attempt to load state for parameter sampling portion of current iteration.
        """
        self.vlog('sampling parameters',n=2)
        return self.load_state('sample_params')
    #end def load_parameters


    def save_parameters(self,params):
        """
        Save state for parameter sampling portion of current iteration.
        Also write sampled parameters to a text file.
        """
        self.vlog('parameter sampling finished',n=3)
        self.save_state('sample_params')
        self.vlog('saving sampled parameters',n=3)
        self.write_param_file('parameters.dat',params)
    #end def save_parameters


    def save_optimal(self,opt_params,opt_energies):
        """
        Write predicted optimal parameters and energies to a text file.
        """
        self.vlog('saving optimal parameters and energies',n=3)
        self.write_param_file('optimal_parameters.dat',opt_params,iter='prev')
        self.write_energy_file('optimal_energies.dat',opt_energies,iter='prev')
    #end def save_optimal


    def load_energies(self):
        """
        Attempt to load state for energy evaluation portion of current iteration.
        """
        self.vlog('calculating new energies',n=2)
        return self.load_state('update_energies')
    #end def load_energies


    def save_energies(self,energies,errors):
        """
        Save state for energy evaluation portion of current iteration.
        Also write calculated energies to a text file.
        """
        self.vlog('energy calculation finished',n=3)
        self.save_state('update_energies')
        self.vlog('saving calculated energies',n=3)
        self.write_energy_file('energies.dat',energies,errors)
    #end def save_energies

#end class GaussianProcessOptimizer



def rP(P):
    """
    Writes formatted parameter output for GPTestFunction.
    """
    P = P[0]
    s = '('
    for p in P:
        s += '{0: 6.4f},'.format(p)
    #end for
    return s[:-1]+')'
#end def rP


def rE(E):
    """
    Writes formatted energy output for GPTestFunction.
    """
    return '{0: 6.4f}'.format(E[0,0])
#end def rE


class GPTestFunction(DevBase):
    """
    Allows convenient test-driving of GaussianProcessOptimizer on a variety 
    of target functions.
    """

    def __init__(self,
                 name,
                 niterations,
                 param_lower,
                 param_upper,
                 function,
                 Pmin,
                 deterministic = False,
                 sigma         = 1e-6,
                 plot          = False,
                 mode          = 'stateless',
                 verbose       = False,
                 exit_on_save  = False
                 ):
        """
        Parameters
        ----------
        name          : Name of the current test.
        niterations   : Number of GP iterations to perform.
        param_lower   : Lower bounds of the parameter search domain (list-like).
        param_upper   : Upper bounds of the parameter search domain (list-like).
        function      : Target function to optimize.  The simple target function 
                        is wrapped by GPTestFunction.__call__.
        Pmin          : Parameters for the known minimum.
        deterministic : (Optional) flag to mark function as deterministic or not. 
                        If not deterministic, then the target function is expected 
                        to return statistical error bars (bool).
        sigma         : (Optional) "errorbar" to apply to deterministic functions 
                        (float).
        plot          : (Optional) selects whether to plot the target function 
                        along a straight line from the true minimum to the 
                        predicted minimum.
        mode          : (Optional) same as GaussianProcessOptimizer.mode.
        verbose       : (Optional) same as GaussianProcessOptimizer.verbose.
        exit_on_save  : (Optional) same as GaussianProcessOptimizer.exit_on_save.
        """
        Pmin = array(Pmin,dtype=float).ravel()
        Pmin.shape = (1,len(Pmin))
        self.name          = name
        self.niterations   = niterations
        self.param_lower   = param_lower
        self.param_upper   = param_upper
        self.function      = function
        self.Pmin          = Pmin
        self.deterministic = deterministic
        self.sigma         = sigma
        self.plot          = plot
        self.mode          = mode
        self.verbose       = verbose
        self.exit_on_save  = exit_on_save
    #end def __init__


    def __call__(self,P,state=None):
        """
        Calls the simple target function and stores evaluations in arrays 
        of the appropriate shape for gp_opt.  This enables GPTestFunction 
        to act as a functor which is directly used as the energy_function
        of GaussianProcessOptimizer.  Future integrations with Nexus workflows 
        will take this same approach.
        """
        npoints = len(P)
        E  = zeros([npoints,1],dtype=float)
        dE = zeros([npoints,1],dtype=float)
        if self.deterministic:
            for i in range(npoints):
                Pi = P[i].ravel()
                Pi.shape = (1,len(Pi))
                E[i,0]  = self.function(Pi)
                dE[i,0] = self.sigma
            #end for
        else:
            for i in range(npoints):
                Pi = P[i].ravel()
                Pi.shape = (1,len(Pi))
                E[i,0],dE[i,0] = self.function(Pi)
            #end for
        #end if
        return E,dE
    #end def __call__


    def test(self):
        """
        Perform optimization test and print the results of each iteration.
        """
        print
        print 'testing:',self.name
        GPO = GaussianProcessOptimizer(
            niterations     = self.niterations,
            param_lower     = self.param_lower,
            param_upper     = self.param_upper,
            energy_function = self,
            mode            = self.mode,
            verbose         = self.verbose,
            exit_on_save    = self.exit_on_save,
            )

        self.gp = GPO

        P0,E0 = GPO.optimize()
        Eref = []
        Eref_err = []
        for P in P0:
            E,Ee = self(P)
            Eref.append(E)
            Eref_err.append(Ee)
        #end for
        Emin,Emin_err = self(self.Pmin)

        for i in range(len(P0)):
            print rP(P0[i]),rE(E0[i]),rE(Eref[i]),rE(E0[i]-Eref[i]),rE(E0[i]-Emin)
        #end for
        print 20*'-'
        print rP(self.Pmin),rE(Emin)

        if self.plot:
            n=200
            fline = []
            Eline = []
            for i in range(n):
                f = float(i)/(n-1)
                Pf = (1-f)*self.Pmin+f*P0[-1]
                Ef = self(Pf).ravel()[0]
                fline.append(f)
                Eline.append(Ef)
                #print f,Ef
            #end for
            from plotting import figure,plot,show
            figure()
            plot(fline,Eline)
            show()
        #end if
    #end def test
#end class GPTestFunction



if __name__=='__main__':

    #
    # Test drive GaussianProcessOptimizer using various test functions
    #
    #   Most are well-known optimization test functions taken from:
    #     https://en.wikipedia.org/wiki/Test_functions_for_optimization
    #
    #   A Lennard-Jones test function is also provided with known dumbell, 
    #   triangle, and tetrahedra minima presented as 1, 3, and 6 dimensional
    #   optimization problems.
    #
    #   The current algorithm is challenged by functions with large variations
    #   in scale, such as those represented by the Rosenbrock, Beale, 
    #   Goldstein-Price, and Lennard-Jones functions below.
    #
    
    from numpy import abs,sin,cos,exp,sqrt,pi,log
    from numpy.linalg import norm

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
        y = 20.0+exp(1.0)-20.0*exp(-0.2*sqrt((x**2).sum()/float(d)))- \
            exp((cos(2.0*pi*x)).sum()/float(d))

        return y
    #end def ackleyf


    def sphere(x):
        x = x.ravel()
        return (x**2).sum()
    #end def sphere


    def rosenbrock(x):
        x = x.ravel()
        return ( (1-x[0:-1])**2+100*(x[1:]-x[0:-1]**2)**2 ).sum()
    #end def rosenbrock


    def beale(x):
        x,y = x.ravel()
        return (1.5-x+x*y)**2+(2.25-x+x*y**2)**2+(2.625-x+x*y**3)**2
    #end def beale


    def goldstein_price(x):
        x,y = x.ravel()
        return (1.+(x+y+1.)**2*(19.-14*x+3*x**2-14*y+6*x*y+3*y**2))*(30+(2*x-3*y)**2*(18.-32*x+12*x**2+48*y-36*x*y+27*y**2))
    #end def goldstein_price


    def himmelblau(x):
        x,y = x.ravel()
        return (x**2+y-11)**2+(x+y**2-7)**2
    #end def himmelblau


    def holder_table(x):
        x,y = x.ravel()
        return -abs(sin(x)*cos(y)*exp(abs(1.0-sqrt(x**2+y**2)/pi)))
    #end def holder_table

    
    def lennard_jones(x):
        # parameters, 1+2+3*n
        x = x.ravel()
        nparams = len(x)
        if nparams!=1 and nparams%3!=0 or nparams==0:
            error('invalid number of parameters for lennard jones')
        #end if
        r = [[   0,0,0],
             [x[0],0,0]]
        if nparams>1:
            r.append([x[1],x[2],0])
        #end if
        if nparams>3:
            i=3
            while i+2<nparams:
                r.append([x[i],x[i+1],x[i+2]])
                i+=3
            #end while
        #end if
        r = array(r,dtype=float)
        npoints = len(r)
        E1 = 1.0
        rm = 1.0
        E = 0.0
        for i in range(npoints):
            for j in range(i+1,npoints):
                d = norm(r[i]-r[j])
                ELJ = E1*((rm/d)**12-2*(rm/d)**6)
                E += ELJ
            #end for
        #end for
        return E
    #end def lennard_jones
        

    def trunc_lennard_jones(x):
        E = lennard_jones(x)
        E = min(E,10)
        return E
    #end def trunc_lennard_jones
        

    def trunc_log_lennard_jones(x):
        E = lennard_jones(x)
        if E>1:
            E = log(E)+1
        #end if
        return E
    #end def trunc_log_lennard_jones
        

    def trunc_lennard_jones(x):
        E = lennard_jones(x)
        E = min(E,10)
        return E
    #end def trunc_lennard_jones


    sphere2 = GPTestFunction(
        name          = 'sphere2',
        niterations   = 6,
        param_lower   = [-1,-1],
        param_upper   = [ 1, 1],
        function      = sphere,
        Pmin          = [0,0],
        deterministic = True,
        sigma         = 1e-6,
        )

    ackley2 = GPTestFunction(
        name          = 'ackley2',
        niterations   = 5,
        param_lower   = [-1,-1],
        param_upper   = [ 1, 1],
        function      = ackleyf,
        Pmin          = [0,0],
        deterministic = True,
        sigma         = 1e-6,
        )

    ackley4 = GPTestFunction(
        name          = 'ackley4',
        niterations   = 10,
        param_lower   = [-1,-1,-1,-1],
        param_upper   = [ 1, 1, 1, 1],
        function      = ackleyf,
        Pmin          = [0,0,0,0],
        deterministic = True,
        sigma         = 1e-6,
        )

    rosenbrock2 = GPTestFunction(
        name          = 'rosenbrock2',
        niterations   = 8,
        param_lower   = [-2,-2],
        param_upper   = [ 2, 2],
        function      = rosenbrock,
        Pmin          = [1,1],
        deterministic = True,
        sigma         = 1e-6,
        )

    beale2 = GPTestFunction(
        name          = 'beale2',
        niterations   = 8,
        param_lower   = [-4.5,-4.5],
        param_upper   = [ 4.5, 4.5],
        function      = beale,
        Pmin          = [3,0.5],
        deterministic = True,
        sigma         = 1e-6,
        )

    goldstein_price2 = GPTestFunction(
        name          = 'goldstein_price2',
        niterations   = 8,
        param_lower   = [-2,-2],
        param_upper   = [ 2, 2],
        function      = goldstein_price,
        Pmin          = [0,-1],
        deterministic = True,
        sigma         = 1e-6,
        )

    himmelblau2 = GPTestFunction(
        name          = 'himmelblau2',
        niterations   = 10,
        param_lower   = [-5,-5],
        param_upper   = [ 5, 5],
        function      = himmelblau,
        Pmin          = [-3.779310,-3.283186],
        deterministic = True,
        sigma         = 1e-6,
        )

    holdertable2 = GPTestFunction(
        name          = 'holdertable2',
        niterations   = 10,
        param_lower   = [-10,-10],
        param_upper   = [ 10, 10],
        function      = holder_table,
        Pmin          = [8.05502,-9.66459],
        deterministic = True,
        sigma         = 1e-6,
        )

    lennardjones1 = GPTestFunction(
        name          = 'lennardjones1',
        niterations   = 10,
        #param_lower   = [ 0 ],
        param_lower   = [ 0.4 ],
        param_upper   = [ 2 ],
        function      = lennard_jones,
        Pmin          = [1.0],
        deterministic = True,
        sigma         = 1e-6,
        #plot          = True,
        )

    lennardjones3 = GPTestFunction(
        name          = 'lennardjones3',
        niterations   = 5,
        #param_lower   = [ 0, 0, 0],
        param_lower   = [ 0.4, 0.4, 0.4],
        param_upper   = [ 2, 2, 2],
        function      = lennard_jones,
        #function      = trunc_log_lennard_jones,
        Pmin          = [1.0,0.5,sqrt(3.)/2],
        deterministic = True,
        sigma         = 1e-6,
        #plot          = True,
        )

    lennardjones6 = GPTestFunction(
        name          = 'lennardjones6',
        niterations   = 10,
        param_lower   = [ 0, 0, 0, 0, 0, 0],
        param_upper   = [ 2, 2, 2, 2, 2, 2],
        function      = lennard_jones,
        Pmin          = [1.0,0.5,sqrt(3.)/2,0.5,1./(2*sqrt(3.)),sqrt(2./3)],
        deterministic = True,
        sigma         = 1e-6,
        )


    sphere2_save_state = GPTestFunction(
        name          = 'sphere2',
        niterations   = 6,
        param_lower   = [-1,-1],
        param_upper   = [ 1, 1],
        function      = sphere,
        Pmin          = [0,0],
        deterministic = True,
        sigma         = 1e-6,
        mode          = 'stateful',
        verbose       = True,
        exit_on_save  = False,
        )

    #sphere2.test()
    #ackley2.test()
    #ackley4.test()
    #rosenbrock2.test()
    #beale2.test()
    #goldstein_price2.test()
    #himmelblau2.test()
    #holdertable2.test()
    lennardjones1.test()
    #lennardjones3.test()
    #lennardjones6.test()

    #sphere2_save_state.test()

#end if
