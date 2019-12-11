"""Simple extraction of afqmc rdms."""
import h5py
import numpy
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
        return mean.reshape(nbasis,nbasis), err.reshape(nbasis, nbasis)
    elif walker == 'collinear':
        return mean.reshape((2,nbasis,nbasis)), err.reshape((2, nbasis, nbasis))
    elif walker == 'non_collinear':
        return mean.reshape((2*nbasis,2*nbasis)), err.reshape((2*nbasis, 2*nbasis))
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
                ij += 1
        two_rdm[nbasis:,nbasis:] = two_rdm[:nbasis,:nbasis].copy()
        two_rdm_err[nbasis:,nbasis:] = two_rdm_err[:nbasis,:nbasis].copy()
    elif walker == 'collinear':
        dm_size = nbasis*(2*nbasis-1)
        two_rdm[numpy.triu_indices] = mean
        two_rdm_err[numpy.triu_indices] = err
    elif walker == 'non_collinear':
        print("Non-collinear wavefunction not supported.")
        return None
    else:
        print('Unknown walker type.')
        return None
    # Diagonal is zero
    two_rdm = 0.5 * (two_rdm + two_rdm.conj().T)
    two_rdm_err = 0.5 * (two_rdm_err + two_rdm_err.T)
    return two_rdm, two_rdm_err

def average_on_top_pdm(filename, estimator='back_propagated', eqlb=1, skip=1, ix=None):
    """Average on-top pair density matrix.

    Returns n(r,r) for a given real space grid.

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

def average_observable(filename, name, eqlb=1, estimator='back_propagated',
                       ix=None, skip=1):
    """Compute mean and error bar for AFQMC HDF5 observable.

    Returns n(r,r) for a given real space grid.

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
