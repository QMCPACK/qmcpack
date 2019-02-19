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


from developer import unavailable,available,DevBase,log,warn,error

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


class GParmFlags(DevBase):
    """
    Define an empty class to store GP parameters and information
    """
    pass
#end class GParmFlags



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
    #end for
    
    return Dm
#end def _fdist

###############################################################################
    
def opt_surface_f(x):

    """
    Optimization Surface Function 
        
    Parameters
    ----------
    x : 1 x d-array, x \in [-1,1]^d
    
    Returns
    -------
    y : value of optimization surface function
    
    """
    [n,d] = x.shape
    y = np.reshape(2.0*np.sin(np.pi*np.sum(x**2,1)/2)-1,[n,1])
    #y = np.reshape(np.prod(np.sin(2*np.pi*x**2),axis=1),[n,1])

    return y
#end def opt_surface

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
    ep,sigma:   GP hyperparameters
    Co_GP:      GP coefficents
    """
    [Nc,d] = mol_pos.shape
    
    W = np.diag(dE[:,0])
    
    D =_fdist(mol_pos)
    ep = .5/np.max(np.min(D+np.max(D)*np.eye(Nc),axis=1))**2
    ep_ceil = 10*ep
    d_ep = 1e-2*ep
    tol = 1e-3*ep
    
    sigma = 1.0
    d_sigma = 1e-2*sigma
    tol_sigma = 1e-3*sigma
    ep_c = ep
    
    K = cfgp(D,ep_c)
    Kw_inv = np.linalg.pinv(K+ sigma*W)
    Co_GP = np.matmul(Kw_inv,E)
    loocv = Co_GP/np.reshape(np.diag(Kw_inv),[Nc,1])
    
    if sigma == 0:
        b_loocv_err = np.sum(-(loocv**2))
        d_sigma = 0
    else:
        sigma_calc = 1-np.diag(np.matmul(K,np.matmul(Kw_inv,np.transpose(K))))
        sigma_calc = np.reshape(np.maximum(sigma_calc,(sigma/10)**2),[Nc,1])
        b_loocv_err = np.sum(-(loocv**2)/sigma_calc-np.log(sigma_calc)/2)- \
            Nc*np.sum(np.abs(np.abs(E-np.matmul(K,Co_GP))-dE)/dE)
    #end if
            
    while d_ep>tol and d_sigma>tol_sigma:
        did_move = 0;
        ep_c = min(ep+d_ep,ep_ceil)
        K = cfgp(D,ep_c)
        Kw_inv = np.linalg.pinv(K+ sigma*W)
        Co_GP = np.matmul(Kw_inv,E)
        loocv = Co_GP/np.reshape(np.diag(Kw_inv),[Nc,1])
        if sigma == 0:
            sigma_calc = 1;
            c_loocv_err = np.sum(-(loocv**2))
        else:
            sigma_calc = 1-np.diag(np.matmul(K,np.matmul(Kw_inv,np.transpose(K))))
            sigma_calc = np.reshape(np.maximum(sigma_calc,(sigma/10)**2),[Nc,1])
            c_loocv_err = np.sum(-(loocv**2)/sigma_calc-np.log(sigma_calc)/2)- \
                Nc*np.sum(np.abs(np.abs(E-np.matmul(K,Co_GP))-dE)/dE)
        #end if 
            
        if c_loocv_err>b_loocv_err:
            b_loocv_err = c_loocv_err
            ep = ep_c
            did_move = 1
        #end if

        ep_c = ep-d_ep
        if ep_c<=0:
            ep_c = ep
        #end if

        K = cfgp(D,ep_c)
        Kw_inv = np.linalg.pinv(K+ sigma*W)
        Co_GP = np.matmul(Kw_inv,E)
        loocv = Co_GP/np.reshape(np.diag(Kw_inv),[Nc,1])
        if sigma == 0:
            sigma_calc = 1;
            c_loocv_err = np.sum(-(loocv**2))
        else:
            sigma_calc = 1-np.diag(np.matmul(K,np.matmul(Kw_inv,np.transpose(K))))
            sigma_calc = np.reshape(np.maximum(sigma_calc,(sigma/10)**2),[Nc,1])
            c_loocv_err = np.sum(-(loocv**2)/sigma_calc-np.log(sigma_calc)/2)- \
                Nc*np.sum(np.abs(np.abs(E-np.matmul(K,Co_GP))-dE)/dE)
        #end if
            
        if c_loocv_err>b_loocv_err:
            b_loocv_err = c_loocv_err
            ep = ep_c
            did_move = 1
        #end if

        if did_move == 1:
            d_ep = 1.25*d_ep
        else:
            d_ep = .75*d_ep
        #end if

        if sigma>0:
            did_move = 0;
            sigma_c = sigma+d_sigma
            K = cfgp(D,ep)
            Kw_inv = np.linalg.pinv(K+ sigma_c*W)
            Co_GP = np.matmul(Kw_inv,E)
            loocv = Co_GP/np.reshape(np.diag(Kw_inv),[Nc,1])

            sigma_calc = 1-np.diag(np.matmul(K,np.matmul(Kw_inv,np.transpose(K))))
            sigma_calc = np.reshape(np.maximum(sigma_calc,(sigma/10)**2),[Nc,1])
            c_loocv_err = np.sum(-(loocv**2)/sigma_calc-np.log(sigma_calc)/2)- \
                Nc*np.sum(np.abs(np.abs(E-np.matmul(K,Co_GP))-dE)/dE)
                
                
            if c_loocv_err>b_loocv_err:
                b_loocv_err = c_loocv_err
                sigma = sigma_c
                did_move = 1
            #end if
 
            sigma_c = sigma-d_sigma
            if sigma_c<=0:
                sigma_c = sigma
            #end if

            Kw_inv = np.linalg.pinv(K+ sigma_c*W)
            Co_GP = np.matmul(Kw_inv,E)
            loocv = Co_GP/np.reshape(np.diag(Kw_inv),[Nc,1])

            sigma_calc = 1-np.diag(np.matmul(K,np.matmul(Kw_inv,np.transpose(K))))
            sigma_calc = np.reshape(np.maximum(sigma_calc,(sigma/10)**2),[Nc,1])
            c_loocv_err = np.sum(-(loocv**2)/sigma_calc-np.log(sigma_calc)/2)- \
                Nc*np.sum(np.abs(np.abs(E-np.matmul(K,Co_GP))-dE)/dE)
                
                
            if c_loocv_err>b_loocv_err:
                b_loocv_err = c_loocv_err
                sigma = sigma_c
                did_move = 1
            #end if

            if did_move == 1:
                d_sigma = 1.25*d_sigma
            else:
                d_sigma = .75*d_sigma
            #end if
        #end if
        K = cfgp(D,ep)
        Kw_inv = np.linalg.pinv(K+ sigma*W)
        Co_GP = np.matmul(Kw_inv,E)
    #end while

    return ep,sigma,Co_GP 
#end def find_hyperparameters

###############################################################################

def find_min(Co_GP,mol_pos,E,dE,ep):
    """
    Finds local minimum of Gaussian Process
    
    
    Parameters
    ----------
    mol_pos:    GP centers 
    E:          Energy for GP centers
    dE:         Uncertainty in energy for GP centers
    Co_GP:      GP coefficents
    ep:         GP exponential parameter
    
    Returns
    -------
    best_pos:   Estimated local minimum
    best_E:     Estimated local energies
    neig_data:  Neighborhood Data
    """
    [Nc,d] = mol_pos.shape
    
    T = Delaunay(mol_pos)
    D =_fdist(mol_pos)
    best_pos = mol_pos
    ref_dis = np.min(np.min(D+np.max(D)*np.eye(Nc),axis=1))/10
    K = cfgp(D,ep)
    
    delta_pos = ref_dis+np.zeros([Nc,d])
    tol_pos = ref_dis/10
    best_E = np.matmul(K,Co_GP)
    
    while np.max(delta_pos)>tol_pos:
        for jp in range(Nc):
            for jd in range(d):
                did_move = 0
                c_pos = np.copy(best_pos)
                c_pos[jp,jd]= best_pos[jp,jd]+delta_pos[jp,jd]
                ind_T = T.find_simplex(c_pos[jp,:])
                if ind_T == -1:
                    c_pos[jp,jd]= best_pos[jp,jd]
                #end if

                Kc = cmgp(np.reshape(c_pos[jp,:],[1,d]),mol_pos,ep)
                c_E = np.matmul(Kc,Co_GP)
                if c_E<best_E[jp]:
                    best_E[jp] = c_E
                    best_pos = np.copy(c_pos)
                    did_move = 1
                #end if

                c_pos = np.copy(best_pos)
                c_pos[jp,jd]=best_pos[jp,jd]-delta_pos[jp,jd]
                ind_T = T.find_simplex(c_pos[jp,:])
                if ind_T == -1:
                    c_pos[jp,jd]= best_pos[jp,jd]
                #end if
                Kc = cmgp(np.reshape(c_pos[jp,:],[1,d]),mol_pos,ep)
                c_E = np.matmul(Kc,Co_GP)
                if c_E<best_E[jp]:
                    best_E[jp] = c_E
                    best_pos = np.copy(c_pos)
                    did_move = 1
                #end if

                if did_move == 1:
                    delta_pos[jp,jd] = 1.25*delta_pos[jp,jd]
                else:
                    delta_pos[jp,jd] = 0.75*delta_pos[jp,jd]
                #end if
            #end for
        #end for
    #end while
     
    ## Remove in same simplex ##
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
        #end if
        jp += 1
    #end while
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
        #end if
        jp += 1
    #end while
    ind_T = T.find_simplex(best_pos)
    ind_T = np.reshape(ind_T,[ind_T.shape[0],1])
    
    if best_E.shape[0] > 1:
        ind_sort = np.argsort(best_E[:,0])
        ind_rm = np.argwhere(best_E[ind_sort,0]>np.mean(E)-np.std(E))
        if ind_rm.shape[0]>0:
            if ind_rm[0,0] == 0:
                ind_rm = 1
            else:
                ind_rm = ind_rm[0,0]
            #end if  
            ind_sort = ind_sort[0:ind_rm]
        #end if
        best_E = best_E[ind_sort,:]
        best_pos = best_pos[ind_sort,:]
        ind_T = ind_T[ind_sort,:]
    #end if
    neig_data = E[T.simplices[ind_T[:,0],:],0]
    neig_data = np.append(neig_data,dE[T.simplices[ind_T[:,0],:],0],axis=0)
        
    return best_pos,best_E,neig_data
#end def find_min

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
        test_points = np.repeat(np.reshape(shyper_cube[:,1]-shyper_cube[:,0],[1,d]) \
            ,n_pnts,axis=0)*lhs(d,n_pnts)+ \
            np.repeat(np.reshape(shyper_cube[:,0],[1,d]),n_pnts,axis=0)
        test_points[0,:] = best_pos[jc,:]   
        new_points[ind_clusters[jc]:ind_clusters[jc]+n_pnts,:] = test_points
        
        D =_fdist(np.append(mol_pos,new_points[0:ind_clusters[jc]+n_pnts,:],axis=0))
        b_mes = np.sum(np.min(D+np.max(D)*np.eye(Nc+ind_clusters[jc]+n_pnts),axis=1))
        
        for jit in range(n_min_d):
            cnew_points = np.copy(new_points)
            test_points = np.repeat(np.reshape(shyper_cube[:,1]-shyper_cube[:,0],[1,d]) \
                ,n_pnts,axis=0)*lhs(d,n_pnts)+ \
                np.repeat(np.reshape(shyper_cube[:,0],[1,d]),n_pnts,axis=0)
            test_points[0,:] = best_pos[jc,:]      
            cnew_points[ind_clusters[jc]:ind_clusters[jc]+n_pnts,:] = test_points
            
            D =_fdist(np.append(mol_pos,cnew_points[0:ind_clusters[jc]+n_pnts,:],axis=0))
            c_mes = np.sum(np.min(D+np.max(D)*np.eye(Nc+ind_clusters[jc]+n_pnts),axis=1))
            
            if c_mes>b_mes:
                new_points = np.copy(cnew_points)
                b_mes = c_mes
            #end if
        #end for
    #end for
    
    return new_points
#end def find_new_points
    
###############################################################################

def multiscale_E(E,dE,scale_thres):
    """
    Multiscale normalization of E
    
    
    Parameters
    ----------
    E:              Energy
    dE:             Uncertainty in Energy
    scale_thres:    Scale parameter
    
    Returns
    -------
    E_a:        Scale parameters
    E_b:        Scale parameters
    E_scaled:   Scaled energies
    dE_scaled:  Scaled uncertainties in energies
    """
    
    E_a = min(E)
    E_scaled = (E-E_a)/scale_thres
    dE_scaled = dE/scale_thres
    us_ind = np.argwhere(E_scaled[:,0]>1)
    if us_ind.shape[0]>0:   
        E_scaled[us_ind[:,0],0] = np.log(E_scaled[us_ind[:,0],0])+1
        dE_scaled[us_ind[:,0],0] = max(dE_scaled)
        E_b = max(E_scaled[us_ind[:,0],0])
        if E_b>2*scale_thres:
            E_scaled[us_ind[:,0],0] = scale_thres*(E_scaled[us_ind[:,0],0]-1)/(E_b-1)+1
        else:
            E_b= 0
        #end if
    else:
        E_scaled = E
        dE_scaled = dE
        E_a = 0
        E_b= 0
    #end if
        
    return E_a,E_b,E_scaled,dE_scaled
#end def multiscale_E

###############################################################################

def inv_multiscale_E(E_a,E_b,E_scaled,scale_thres):
    """
    Reverses Multiscale normalization of E
    
    
    Parameters
    ----------
    E_a:            Scale parameters
    E_b:            Scale parameters
    E_scaled:       Scaled energies
    scale_thres:    Scale parameter
    
    Returns
    -------
    E:          Energy
    """
    if E_a == 0:
        E = E_scaled
    else:
        us_ind = np.argwhere(E_scaled[:,0]>1)
        if us_ind.shape[0]>0:
            if E_b!=0:
                E_scaled[us_ind[:,0],0] = (E_scaled[us_ind[:,0],0]-1)*(E_b-1)/scale_thres+1
            #end if
            E_scaled[us_ind[:,0],0] = np.exp(E_scaled[us_ind[:,0],0]-1)
        #end if
        E = scale_thres*E_scaled+E_a
    #end if
        
    return E
#end def inv_multiscale_E

###############################################################################
    
def norm_pos(mol_pos,hyper_cube):
    """
    Normalizes molecular positions to unit hyper cube
    
    
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
#end def norm_pos

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
#end def inv_norm_pos

###############################################################################


def gp_sdg_setup(GP_info):
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
    
    if 'hyper_cube' not in GP_info:
        error('Please Define GP_info.hyper_cube')
    #end if

    d = GP_info.hyper_cube.shape[0]

    GP_info.set_optional(
        n_its           = 5,          # Number of maximum iterations
        jits            = 0,          # Current iteration index
        do_vis          = 0,          # Visualize Error
        do_scale_energy = 0,          # Scale energy
        scale_thres     = 0.,         # Energy scale parameter
        Nc0             = (d+2)*(d+1),# Number of points to add per iteration
        delta_r         = 1.0,        # Initial Radius of Convergence
        delta_r_s       = 2.0,        # Radius shrink factor
        n_min_d         = 5,          # Number of times to resample LHC
        )

    return E,dE,mol_pos,best_pos,best_E,GP_info
#end def gp_sdg_setup

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
    else:
        GP_info.delta_r = GP_info.delta_r/GP_info.delta_r_s
        mol_pos_norm = norm_pos(mol_pos,GP_info.hyper_cube)
        
        if GP_info.do_scale_energy == 1:
            E_a,E_b,E_scaled,dE_scaled = multiscale_E(E,dE,GP_info.scale_thres)
            ## Best Approximation of Optimization Surface ##
            ep,sigma,Co_GP = find_hyperparameters(mol_pos_norm,E_scaled,dE_scaled)
            ## Find Global Min ##
            best_pos_norm,best_E_scaled,neig_data = find_min(Co_GP,mol_pos_norm,E_scaled,dE_scaled,ep)
            best_E = inv_multiscale_E(E_a,E_b,best_E_scaled,GP_info.scale_thres)
            
        else:
            ep,sigma,Co_GP = find_hyperparameters(mol_pos_norm,E,dE)
            best_pos_norm,best_E,neig_data = find_min(Co_GP,mol_pos_norm,E,dE,ep)
        #end if

        best_pos = inv_norm_pos(best_pos_norm,GP_info.hyper_cube)
        
        new_points = find_new_points(mol_pos_norm,best_pos_norm,GP_info.delta_r,GP_info.Nc0)
        
        new_points = inv_norm_pos(new_points,GP_info.hyper_cube)
        
        if GP_info.do_vis == 1 and available(plt):
            neig_bool = neig_data[0:best_E.shape[0],:]- \
                2*neig_data[best_E.shape[0]:2*best_E.shape[0],:]
            neig_bool = neig_bool<np.repeat(best_E,d+1,axis=1)
            neig_percent = np.zeros(neig_bool.shape)
            for jx in range(neig_bool.shape[0]):
                for jy in range(neig_bool.shape[1]):
                    neig_percent[jx,jy] = float(neig_bool[jx,jy])
                #end for
            #end for
            neig_percent = 100*np.mean(neig_percent,axis=1)
            x = np.zeros((d+1)*best_E.shape[0])
            y = np.zeros((d+1)*best_E.shape[0])
            dy = np.zeros((d+1)*best_E.shape[0])
            Ex = np.zeros(best_E.shape[0])
            for jx in range(best_E.shape[0]):
                x[jx*(d+1):(jx+1)*(d+1)] = jx
                y[jx*(d+1):(jx+1)*(d+1)] = neig_data[jx,:]
                dy[jx*(d+1):(jx+1)*(d+1)]= neig_data[jx+best_E.shape[0],:]
                Ex[jx] = jx
            #end for

            plt.figure()
            plt.errorbar(x, y, yerr=dy, fmt='o',label='sampled point with uncertainty')
            plt.plot(Ex,best_E[:,0], 'r+',label='estimated minumum')
            plt.legend()
            plt.title('Neigborhood Information')
            plt.xlabel('Local Minimum')
            plt.ylabel('E')
            plt.show()
        elif GP_info.do_vis == 1:
            warn('visualization requested but matplotlib is not available\nplotting suppressed')
        #end if
    #end if
        
    return new_points,best_pos,best_E,GP_info
#end def gp_opt

###############################################################################
    

def gp_opt_example():
    """
    Example using Gaussian Process optimization for noisy data and multiple mins

    """

    ## Required to define is hypercube  ##
    d = 4                                       # Dimension of problem
    GP_info = GParmFlags()
    GP_info.hyper_cube = np.ones([d,2])/2       # Upper and lower bound for each dimension
    GP_info.hyper_cube[:,0] = -GP_info.hyper_cube[:,0]/2


    ## Parameters Used for Optimization Function Defining Range of Errors in Energy##
    sigma_a = 1e-3;
    sigma_b = 1e-1;


    ## Setup Gausian Process Structure ##
    E,dE,mol_pos,best_pos,best_E,GP_info = gp_sdg_setup(GP_info)

    ## Do Gaussian Process Iteration ##
    for jits in range(GP_info.n_its):
        GP_info.jits = jits
        new_points,best_pos,best_E,GP_info = gp_opt(E,dE,mol_pos,GP_info)
        print jits,'best_P',best_pos
        print jits,'best_E',best_E
        if jits < GP_info.n_its-1:
            ## Simulate random uncertainty in optimization function
            new_dE = (sigma_b-sigma_a)*np.random.rand(GP_info.Nc0,1)+sigma_a
            new_E = opt_surface_f(new_points)+new_dE*np.random.randn(GP_info.Nc0, 1)
            if jits == 0:
                mol_pos = np.copy(new_points)
                E = np.copy(new_E)
                dE = np.copy(new_dE)
            else:
                mol_pos = np.append(mol_pos,new_points,axis=0)
                E = np.append(E,new_E,axis=0)
                dE = np.append(dE,new_dE,axis=0)
            #end if
        #end if
    #end for
#end def gp_opt_example



if __name__=='__main__':
    gp_opt_example()
#end if
