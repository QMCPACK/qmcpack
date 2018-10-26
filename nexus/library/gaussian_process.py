##################################################################
##  (c) Copyright 2018-  by Richard K. Archibald                ##
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
#    cfgp                                                            #
#      Function to calculate correlation function for Gaussian       #
#      Process.  Subroutine of cmgp.                                 #
#                                                                    #
#    cmgp                                                            #
#      Function to calculate correlation matrix for Gaussian         #
#      Process.  Subroutine of sgd and gp_sgd.                       #
#                                                                    #
#    sgd                                                             #
#      Function to calculate minimums of Gaussian Process using      #
#      stochastic gradient descent.  Subroutine of gp_sgd.           #
#                                                                    #
#    gp_sgd                                                          #
#      Function that performs one iteration of Gaussian Process      #
#      optimization using stochastic gradient descent.  This is the  #
#      main interface to Gaussian Process iterations.                #
#                                                                    #
#    ackleyf                                                         #
#      Function to calculate the d-dimensional Ackley function on    #
#      [-1,1]^d.                                                     #
#                                                                    #
#    ackley_example                                                  #
#      Example using Gaussian Process optimization on Ackley         #
#      function.                                                     #
#                                                                    #
#====================================================================#


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
    #end if
    
    if iterations is None:
        iterations = 5
    #end if

    H = _lhsmaximin(n, samples, iterations)
    
    return H
#end def lhs

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
    #end for

    # Make the random pairings
    H = np.zeros_like(rdpoints)
    for j in range(n):
        order = np.random.permutation(range(samples))
        H[:, j] = rdpoints[order, j]
    #end for

    return H
#end def _lhsclassic
    
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
        #end if
    #end for
    return H
#end def _lhsmaximin

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
    #end if
    
    d = np.zeros((m**2-m)/2)
    c_v = 0
    for i in range(m - 1):
        d[0+c_v:m-i-1+c_v]=np.sqrt(np.sum((np.repeat(np.reshape(x[i, :],\
                        [1,n]),m-i-1,axis=0)-x[i+1:m+1, :])**2,axis=1))
        c_v = c_v+m-i-1
    #end for

    return d
#end def _pdist

###############################################################################
    
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
#end def ackleyf

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
#end def cfgp

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
    #end for
    
    K = cfgp(K,ep)
    
    return K
#end def cmgp

###############################################################################

def sgd(hyper_cube,co_GP,mol_pos,ep,tol=None,verbose=None):
    """
    Calculate minimums of Gaussian Process using stochastic gradient descent
    
    
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
    #end if
    
    if verbose is None:
        verbose=0
    #end if

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
            #end for

            c_E = np.matmul(cmgp(t_pos,mol_pos,ep),co_GP)
            for jp in range(npoints):
                if c_E[jp,0]<best_E[jp,0]:
                    did_move[jp,0] = 1
                    best_E[jp,0] = c_E[jp,0]
                    best_pos[jp,:] = t_pos[jp,:]
                #end if
            #end for

            t_pos =  np.copy(best_pos)
            t_pos[:,jd] = t_pos[:,jd] - delta_pos[:,jd]
            for jdim in range(d):
                t_pos[np.where(t_pos[:,jdim]<hyper_cube[jdim,0]),jdim]= \
                    hyper_cube[jdim,0]
            #end for

            c_E = np.matmul(cmgp(t_pos,mol_pos,ep),co_GP)
            for jp in range(npoints):
                if c_E[jp,0]<best_E[jp,0]:
                    did_move[jp,0] = 1
                    best_E[jp,0] = c_E[jp,0]
                    best_pos[jp,:] = t_pos[jp,:]
                #end if
                if did_move[jp,0] == 1:
                    delta_pos[jp,jd] = 1.25*delta_pos[jp,jd]
                else:
                    delta_pos[jp,jd] = 0.75*delta_pos[jp,jd]
                #end if
            #end for
        #end for

        cm = np.max(delta_pos)
        if verbose == 1:
            use_ind = np.where(np.max(delta_pos,axis=1)>epm*tol)
            print cm,np.shape(use_ind)[1]
        #end if
    #end while

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
    #end while

    min_loc = np.zeros([c_class,d])
    min_E = np.zeros([c_class,1])  
    
    for jc in range(c_class):
        indc = np.where(ind_marker==jc)
        c_E = best_E[indc[0],:]
        t_pos = best_pos[indc[0],:]
        indb = np.argmin(c_E)
        min_E[jc,:] = c_E[indb,:]
        min_loc[jc,:] = t_pos[indb,:]
    #end for
        
    return min_loc,min_E
#end def sgd

###############################################################################
def gp_sgd(hyper_cube_F,E_F,P_F,E0_F,P0_F,sigma,npoints=None,nlhc_its=None,\
           ep=None,tol=None,verbose=None,gen_new=None):
    """
    Performs one iteration of Gaussian Process optimization using 
    stochastic gradient descent
    
    
    Parameters
    ----------
    
    hyper_cube_F:   Hypercube for all levels
    E_F:            Energies for all levels
    P_F:            Sample points for all levels
    E0_F:           Estimated Minimum Energy
    P0_F:           Estimated Position of Minimum
    sigma:          Error measure
    npoints:        Number of points to create optimization surface 
    nlhc_its:       The number of iterations in the lhc algorithm
    ep:             GP kernel parameter ep>0
    tol:            Tolerance for termination
    
    Returns following for next iteration
    -------
    hyper_cube_F:   Hypercube for all levels
    E_F:            Energies for all levels
    P_F:            Sample points for all levels
    E0_F:           Estimated Minimum Energy
    P0_F:           Estimated Position of Minimum
    """
    
    d = hyper_cube_F[0].shape[0]
    
    if tol is None:
        tol = 1e-4
    #end if
    if verbose is None:
        verbose=0 
    #end if
    if npoints is None:
        npoints = (d+2)*(d+1)+1 
    #end if
    if nlhc_its is None:
        nlhc_its = 10
    #end if
    if gen_new is None:
        gen_new = 1
    #end if

    if len(P_F)==0:
        hyper_cube = hyper_cube_F[0]
        
        ## Generate random sample from LHS
        mol_pos = np.repeat(np.reshape(hyper_cube[:,1]-hyper_cube[:,0],[1,d]) \
            ,npoints,axis=0)*lhs(d, samples=npoints,iterations=nlhc_its)+ \
            np.repeat(np.reshape(hyper_cube[:,0],[1,d]),npoints,axis=0)
        
        hyper_cube[:,0] = np.reshape(np.min(mol_pos,axis=0),d)
        hyper_cube[:,1] = np.reshape(np.max(mol_pos,axis=0),d)
        
        hyper_cube_F[0] = hyper_cube.copy()
        
        P_F.append(mol_pos)

    else:
        jits = len(P_F)-1
        
        mol_pos = P_F[jits]
        E = E_F[jits]
        hyper_cube = hyper_cube_F[jits]
        ## Estimate Optimal Choice for ep
        ep  = 1.0/np.max(_pdist(mol_pos))

        ## Build Gaussian Process
        K = cmgp(mol_pos,mol_pos,ep) 
        Ks = K + (sigma**2)*np.identity(npoints)
        co_GP = np.linalg.solve(Ks,E)
        
        if verbose != 0:
            print "Starting Stochastic Grandient Descent"
        #end if

        ## Use stochastic gradient descent to determine minimums
        min_loc,min_E = sgd(hyper_cube,co_GP,mol_pos,ep)
        
        indb = np.argmin(min_E)
        bmin_loc = np.reshape(min_loc[indb,:],[1,d])
        
        P0_F.append(bmin_loc)
        E0_F.append(min_E[indb])
        
        if gen_new == 1:
            hyper_cube_F.append((hyper_cube_F[jits].copy()+\
                        np.repeat(np.reshape(bmin_loc,[d,1]),2,axis=1))/2)
            
    
            hyper_cube_F[jits+1][:,0] = np.maximum(hyper_cube_F[jits+1][:,0],\
                        hyper_cube_F[jits][:,0])
            hyper_cube_F[jits+1][:,1] = np.minimum(hyper_cube_F[jits+1][:,1],\
                        hyper_cube_F[jits][:,1])
            cLHC = hyper_cube_F[jits+1]
            mol_pos = np.repeat(np.reshape(cLHC[:,1]-cLHC[:,0],[1,d])\
                ,npoints,axis=0)*lhs(d, samples=npoints,iterations=nlhc_its)+\
                np.repeat(np.reshape(cLHC[:,0],[1,d]),npoints,axis=0)
            
            P_F.append(mol_pos)
            
            hyper_cube_F[jits+1][:,0] = np.reshape(np.min(mol_pos,axis=0),d)
            hyper_cube_F[jits+1][:,1] = np.reshape(np.max(mol_pos,axis=0),d)
        #end if

        if verbose != 0:
            print "Error in min position",\
                np.sqrt(np.sum((bmin_loc-x_true)**2)),"its=",jits
                
            if d == 1:
                import matplotlib.pyplot as plt
                npts = 100
                mol_pos = P_F[jits]
                plot_points = np.reshape(np.linspace(hyper_cube[0,0],\
                    hyper_cube[0,1],num=npts),[npts,1])
                gp_E = np.matmul(cmgp(plot_points,mol_pos,ep),co_GP)
                Et = np.zeros([npts,1])
                for jp in range(npts):
                    Et[jp,0] = ackleyf(plot_points[jp,:]-x_true)
                plt.plot(plot_points,gp_E,'b-',plot_points,Et,'g-',\
                         mol_pos,E,'r.')
                plt.show()
            #end if
        #end if
    #end if

    return hyper_cube_F,E_F,P_F,E0_F,P0_F
#end def gp_sgd

###############################################################################





def ackley_example():
    """
    Example using Gaussian Process optimization on Ackley function

    """

    ## Parameters ##
    d = 2                               # Dimension of problem
    nits = 3                            # Number of iterations
    sigma = 1e-6                        # Uncertainty in data
    hyper_cube = np.ones([d,2])         # Upper and lower bound for each dimension
    hyper_cube[:,0] = -hyper_cube[:,0]

    ## Set up gp_sgd structure ##
    hyper_cube_F =[]
    hyper_cube_F.append(hyper_cube)
    E_F = []
    P_F = []
    E0_F = []
    P0_F = []

    ## Begin heirarchical optimization ##
    for jits in range(nits+1):
        if jits < nits:
            hyper_cube_F,E_F,P_F,E0_F,P0_F = gp_sgd(hyper_cube_F,E_F,P_F,E0_F,P0_F,sigma)

            # Determine Energy for new set of points
            npoints = P_F[jits].shape[0]
            E = np.zeros([npoints,1])  
            if jits == 0:
                ## Randomly pick center for Ackley function
                x_true = 0.05*(2.0*np.random.rand(1,d)-1.0)
            #end if

            ## Determine Energies on current points using Ackley function
            for jp in range(npoints):
                E[jp,0] = ackleyf(P_F[jits][jp,:]-x_true)
            #end for
            E_F.append(E)
        else:
            # Final iteration where only minimum is found and no new points generated
            hyper_cube_F,E_F,P_F,E0_F,P0_F = gp_sgd(hyper_cube_F,E_F,P_F,E0_F,P0_F,sigma,gen_new=0)
        #end if
    #end for
#end def ackley_example

