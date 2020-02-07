"""Simple extraction of afqmc rdms."""
import h5py
import numpy
from math import sqrt
import scipy.stats
from afqmctools.analysis.extraction import (
        get_metadata,
        extract_observable
        )


# enumish
WALKER_TYPE = ['undefined', 'closed', 'collinear', 'non_collinear']


def average_one_rdm(filename, estimator='back_propagated', eqlb=1, skip=1, ix=None):
    """Average AFQMC 1RDM.

    Returns P_{sij} = <c_{is}^+ c_{js}^> as a (nspin, M, M) dimensional array.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    eqlb : int
        Number of blocks for equilibration. Default 1.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).

    Returns
    -------
    one_rdm : :class:`numpy.ndarray`
        Averaged 1RDM.
    one_rdm_err : :class:`numpy.ndarray`
        Error bars for 1RDM elements.
    """
    md = get_metadata(filename)
    mean, err = average_observable(filename, 'one_rdm', eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    nbasis = md['NMO']
    wt = md['WalkerType']
    try:
        walker = WALKER_TYPE[wt]
    except IndexError:
        print('Unknown walker type {}'.format(wt))

    if walker == 'closed':
        return mean.reshape(1,nbasis,nbasis), err.reshape(1,nbasis, nbasis)
    elif walker == 'collinear':
        return mean.reshape((2,nbasis,nbasis)), err.reshape((2, nbasis, nbasis))
    elif walker == 'non_collinear':
        return mean.reshape((1,2*nbasis,2*nbasis)), err.reshape((1,2*nbasis, 2*nbasis))
    else:
        print('Unknown walker type.')
        return None

def average_two_rdm(filename, estimator='back_propagated', eqlb=1, skip=1, ix=None):
    """Average AFQMC 2RDM.

    Returns a list of 2RDMS, where 
      2RDM[s1s2,i,k,j,l] = <c_{i}^+ c_{j}^+ c_{l} c_{k}>.
      For closed shell systems, returns [(a,a,a,a),(a,a,b,b)] 
      For collinear systems, returns [(a,a,a,a),(a,a,b,b),(b,b,b,b)] 

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    eqlb : int
        Number of blocks for equilibration. Default 1.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).

    Returns
    -------
    two_rdm : :class:`numpy.ndarray`
        List of averaged 2RDM.
    two_rdm_err : :class:`numpy.ndarray`
        List of error bars for 2RDM elements.
    """
    md = get_metadata(filename)
    mean, err = average_observable(filename, 'two_rdm', eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    nbasis = md['NMO']
    wt = md['WalkerType']
    try:
        walker = WALKER_TYPE[wt]
    except IndexError:
        print('Unknown walker type {}'.format(wt))

    if walker == 'closed':
        return mean.reshape(2,nbasis,nbasis,nbasis,nbasis), err.reshape(2,nbasis,nbasis,nbasis,nbasis)
    elif walker == 'collinear':
        return mean.reshape(3,nbasis,nbasis,nbasis,nbasis), err.reshape(3,nbasis,nbasis,nbasis,nbasis)
    elif walker == 'non_collinear':
        return mean.reshape(2*nbasis,2*nbasis,2*nbasis,2*nbasis), err.reshape(2*nbasis,2*nbasis,2*nbasis, 2*nbasis)
    else:
        print('Unknown walker type.')
        return None

def average_diag_two_rdm(filename, estimator='back_propagated', eqlb=1, skip=1, ix=None):
    """Average diagonal part of 2RDM.

    Returns <c_{is}^+ c_{jt}^+ c_{jt} c_{is}> as a (2M,2M) dimensional array.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    eqlb : int
        Number of blocks for equilibration. Default 1.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).

    Returns
    -------
    two_rdm : :class:`numpy.ndarray`
        Averaged diagonal 2RDM.
    two_rdm_err : :class:`numpy.ndarray`
        Error bars for diagonal 2RDM elements.
    """
    md = get_metadata(filename)
    mean, err = average_observable(filename, 'diag_two_rdm', eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    nbasis = md['NMO']
    wt = md['WalkerType']
    try:
        walker = WALKER_TYPE[wt]
    except IndexError:
        print('Unknown walker type {}'.format(wt))

    if walker == 'closed':
        dm_size = nbasis*(2*nbasis-1) - nbasis*(nbasis-1) // 2
        assert mean.shape == dm_size
        two_rdm = numpy.zeros((2*nbasis, 2*nbasis), dtype=mean.dtype)
        two_rdm_err = numpy.zeros((2*nbasis, 2*nbasis), dtype=mean.dtype)
        ij = 0
        for i in range(nbasis):
            for j in range(i+1, 2*nbasis):
                two_rdm[i,j] = mean[ij]
                two_rdm_err[i,j] = err[ij]
                two_rdm[j,i] = mean[ij].conj()
                two_rdm_err[j,i] = err[ij].conj()
                ij += 1
        two_rdm[nbasis:,nbasis:] = two_rdm[:nbasis,:nbasis].copy()
    elif walker == 'collinear':
        dm_size = nbasis*(2*nbasis-1)
        assert mean.shape == dm_size
        two_rdm = numpy.zeros((2*nbasis, 2*nbasis), dtype=mean.dtype)
        two_rdm_err = numpy.zeros((2*nbasis, 2*nbasis), dtype=mean.dtype)
        ij = 0
        for i in range(2*nbasis):
            for j in range(i+1, 2*nbasis):
                two_rdm[i,j] = mean[ij]
                two_rdm_err[i,j] = err[ij]
                two_rdm[j,i] = mean[ij].conj()
                two_rdm_err[j,i] = err[ij].conj()
                ij += 1
    elif walker == 'non_collinear':
        print("Non-collinear wavefunction not supported.")
        return None
    else:
        print('Unknown walker type.')
        return None
    # Diagonal is zero
    return two_rdm, two_rdm_err

def average_on_top_pdm(filename, estimator='back_propagated', eqlb=1, skip=1, ix=None):
    """Average on-top pair density matrix.

    Returns n2(r,r) for a given real space grid.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    eqlb : int
        Number of blocks for equilibration. Default 1.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).

    Returns
    -------
    opdm : :class:`numpy.ndarray`
        Averaged diagonal on-top pair density matrix.
    opdm_err : :class:`numpy.ndarray`
        Error bars for diagonal on-top pair density matrix elements.
    """
    md = get_metadata(filename)
    mean, err = average_observable(filename, 'on_top_pdm', eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    # TODO: Update appropriately.
    return mean, err

def average_realspace_correlations(filename, estimator='back_propagated', eqlb=1, skip=1, ix=None):
    """Average on-top pair density matrix.

    Returns <C(r1)C(r2)> and <S(r1)S(r2)> for a given set of points in real space.
    C = (nup + ndown), S = (nup - ndown)

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    eqlb : int
        Number of blocks for equilibration. Default 1.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).

    Returns
    -------
    cc_mean : :class:`numpy.ndarray`
        Averaged charge-charge correlation function. 
    cc_err : :class:`numpy.ndarray`
        Error bars for charge-charge correlation function. 
    ss_mean : :class:`numpy.ndarray`
        Averaged spin-spin correlation function. 
    ss_err : :class:`numpy.ndarray`
        Error bars for spin-spin correlation function.  
    """
    md = get_metadata(filename)
    cc_mean, cc_err = average_observable(filename, 'cc_realspace_correlation', eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    ss_mean, ss_err = average_observable(filename, 'ss_realspace_correlation', eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    # TODO: get npts from metadata local to cc_correlation
    npts = int(sqrt(cc_mean.shape[0]))
    return cc_mean.reshape((npts,-1)), cc_err.reshape((npts,-1)), \
                ss_mean.reshape((npts,-1)), ss_err.reshape((npts,-1))

def average_atom_correlations(filename, estimator='back_propagated', eqlb=1, skip=1, ix=None):
    """Average atom centered correlations.

    Returns <C(I)>, <S(I)>, <C(I)C(J)> and <S(I)S(J)> for a given set of atomic sites.
    C = (nup + ndown), S = (nup - ndown)

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    eqlb : int
        Number of blocks for equilibration. Default 1.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).

    Returns
    -------
    c_mean : :class:`numpy.ndarray`
        Averaged charge. 
    c_err : :class:`numpy.ndarray`
        Error bars for charge. 
    s_mean : :class:`numpy.ndarray`
        Averaged spin. 
    s_err : :class:`numpy.ndarray`
        Error bars for spin.  
    cc_mean : :class:`numpy.ndarray`
        Averaged charge-charge correlation function. 
    cc_err : :class:`numpy.ndarray`
        Error bars for charge-charge correlation function. 
    ss_mean : :class:`numpy.ndarray`
        Averaged spin-spin correlation function. 
    ss_err : :class:`numpy.ndarray`
        Error bars for spin-spin correlation function.  
    """
    md = get_metadata(filename)
    c_mean, c_err = average_observable(filename, 'c_atom_correlation', eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    s_mean, s_err = average_observable(filename, 's_atom_correlation', eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    m_mean, m_err = average_observable(filename, 'm_atom_correlation', eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    cc_mean, cc_err = average_observable(filename, 'cc_atom_correlation', eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    ss_mean, ss_err = average_observable(filename, 'ss_atom_correlation', eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    # TODO: get npts from metadata local to cc_correlation
    npts = int(sqrt(cc_mean.shape[0]))
    return c_mean, c_err, s_mean, s_err, m_mean, m_err, cc_mean.reshape((npts,-1)), cc_err.reshape((npts,-1)), \
                ss_mean.reshape((npts,-1)), ss_err.reshape((npts,-1))

def average_observable(filename, name, eqlb=1, estimator='back_propagated',
                       ix=None, skip=1):
    """Compute mean and error bar for AFQMC HDF5 observable.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    name : string
        Name of observable (see estimates.py for list).
    eqlb : int
        Number of blocks for equilibration. Default 1.
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).

    Returns
    -------
    mean : :class:`numpy.ndarray`
        Averaged quantity.
    err : :class:`numpy.ndarray`
        Error bars for quantity.
    """
    md = get_metadata(filename)
    free_proj = md['FreeProjection']
    if free_proj:
        mean = None
        err = None
        print("# Error analysis for free projection not implemented.")
    else:
        data = extract_observable(filename, name=name, estimator=estimator, ix=ix)
        mean = numpy.mean(data[eqlb:len(data):skip], axis=0)
        err = scipy.stats.sem(data[eqlb:len(data):skip].real, axis=0)
    return mean, err

def get_noons(filename, estimator='back_propagated', eqlb=1, skip=1, ix=None,
              nsamp=20, screen_factor=1):
    """Get NOONs from averaged AFQMC RDM.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    eqlb : int
        Number of blocks for equilibration. Default 1.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).
    nsamp : int
        Number of perturbed RDMs to construct to estimate of error bar. Optional.
        Default: 20.
    screen_factor : int
        Zero RDM elements with abs(P[i,j]) < screen_factor*Perr[i,j]. Optional
        Default: 1.

    Returns
    -------
    noons : :class:`numpy.ndarray`
        NOONS.
    noons_err : :class:`numpy.ndarray`
        Estimate of error bar on NOONs.
    """
    P, Perr = average_one_rdm(filename, estimator='back_propagated', eqlb=1,
                              skip=1, ix=ix)
    # Sum over spin.
    P = numpy.sum(P, axis=0)
    if len(Perr.shape[0]) == 2:
        # Collinear
        Perr = numpy.sqrt((Perr[0]**2 + Perr[1]**2))
    else:
        # Non-collinear / Closed
        Perr = Perr[0]
    P = 0.5 * (P + P.conj().T)
    Perr = 0.5 * (Perr + Perr.T)
    P[numpy.abs(P) < screen_factor*Perr] = 0.0
    noons = numpy.zeros((nsamp, P.shape[-1]))
    for s in range(nsamp):
        # Dangerous. Discards imaginary part.
        Ppert = [numpy.random.normal(loc=a, scale=err, size=1) for a, err in
                 zip(P.ravel(), Perr.ravel())]
        Ppert = numpy.array(Ppert).reshape(P.shape)
        Ppert = 0.5*(Ppert + Ppert.conj().T)
        e, ev = numpy.linalg.eigh(Ppert)
        noons[s] = e[::-1]

    e, ev = numpy.linalg.eigh(P)
    return e[::-1], numpy.std(noons, axis=0, ddof=1)
