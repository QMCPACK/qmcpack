"""Simple extraction of afqmc rdms."""
import h5py
import numpy
import scipy.stats

def extract_qmc_dm(filename, dm_name='Mixed'):
    """Extract AFQMC 1RDM from file.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    dm_name: string
        Estimator type. Mixed or BackPropagated. Default: Mixed.

    Returns
    -------
    dm : :class:`numpy.ndarray`
        QMC 1RDM. shape = (ntau, 2, norb, norb).
    weights : :class:`numpy.ndarray`
        Walker weights. shape = (ntau,).
    """
    md = get_metadata(filename)
    walker_type = md['walker_type']
    free_proj = md['free_proj']
    nmo = md['nmo']
    with h5py.File(filename, 'r') as fh5:
        dms = []
        weights = []
        data = fh5[dm_name]
        for G in data.keys():
            if 'denominator' in G:
                weight = data[G][:].view(numpy.complex128).ravel()
                weights.append(weight[0])
                ix = G.split('_')[-1]
                dm = data['full_one_rdm_'+ix][:].view(numpy.complex128).ravel()
                # QMCPACK enum:
                # walker_type == 3: NONCOLLINEAR
                # walker_type == 2: COLLINEAR
                # walker_type == 1: CLOSED
                # walker_type == 0: UNKNOWN
                if walker_type == 2:
                    dm = dm.reshape(2,nmo,nmo)
                    if not free_proj:
                        dm = dm / weight[:,None,None]
                elif walker_type == 1:
                    dm = dm.reshape(nmo,nmo)
                    dm = numpy.array([dm,dm])
                    if not free_proj:
                        dm = dm / weight[:,None,None]
                else:
                    print("Unknown walker type.")
                    return
                dms.append(dm)
    return (numpy.array(dms), numpy.array(weights))

def get_metadata(filename):
    """Extract QMC estimator metadata from h5 file.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).

    Returns
    -------
    md : dict
        Estimator metadata.
    """
    md = {}
    with h5py.File(filename, 'r') as fh5:
        md['nmo'] = fh5['Metadata/NMO'][:][0]
        md['free_proj'] = bool(fh5['Metadata/FreeProjection'][:][0])
        md['walker_type'] = fh5['Metadata/WalkerType'][:][0]
        md['nalpha'] = fh5['Metadata/NAEA'][:][0]
        md['nbeta'] = fh5['Metadata/NAEB'][:][0]
        try:
            md['num_bp']= fh5['Metadata/NumBackProp'][:][0]
            md['num_av'] = fh5['Metadata/NumAverages'][:][0]
            md['num_ref'] = fh5['Metadata/NumReferences'][:][0]
        except KeyError:
            md['num_bp'] = None
        md['dt'] = fh5['Metadata/Timestep'][:][0]
    return md

def get_rdm_len(filename, dm_name='Mixed'):
    """Get number of samples of 1RDM from file.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).

    Returns
    -------
    len_rdm : int
        Number of samples of 1RDM.
    """
    with h5py.File(filename, 'r') as fh5:
        dms = fh5[dm_name].keys()
        # Block size of RDM is not necessarily known so just be dumb and count.
        num_dm = len([n for n in dms if 'denominator' in n])
    return num_dm

def extract_rdm_single(filename, indx, dm_name='Mixed'):
    """Extract single sample of QMC 1RDM.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    indx : int
        Index of sample of RDM to extract.
    dm_name: string
        Estimator type. Mixed or BackPropagated. Default: Mixed.

    Returns
    -------
    dm : :class:`numpy.ndarray`
        QMC 1RDM. shape = (2, norb, norb).
    weight : :class:`numpy.ndarray`
        Sum of walker weights at this sample.
    """
    num_dm = get_rdm_len(filename, dm_name)
    md = get_metadata(filename)
    walker_type = md['walker_type']
    free_proj = md['free_proj']
    nmo = md['nmo']
    with h5py.File(filename, 'r') as fh5:
        keys = list(fh5[dm_name].keys())
        numer = keys[indx]
        denom = keys[indx+num_dm]
        data = fh5[dm_name]
        weight = data[denom][:].view(numpy.complex128).ravel()
        dm = data[numer][:].view(numpy.complex128).ravel()
        # QMCPACK enum:
        # walker_type == 3: NONCOLLINEAR
        # walker_type == 2: COLLINEAR
        # walker_type == 1: CLOSED
        # walker_type == 0: UNKNOWN
        if walker_type == 2:
            if free_proj:
                dm = dm.reshape(2,nmo,nmo)
            else:
                dm = dm.reshape(2,nmo,nmo)/weight[:,None,None]
        elif walker_type == 1:
            dm = numpy.array([dm,dm])
        else:
            print("Unknown walker type.")
            return
    return dm, weight[0]

def get_one_rdm_av(filename, skip, name="Mixed"):
    """Helper routine to compute averaged qmc density matrix.

    TODO: May want to use extract single if DM gets too big.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.scalar.out file).
    skip : float
        Percentage of data to skip at the beginning of the simulation to account
        for equilibration.
    name : string
        Estimator type. Mixed or BackPropagated. Default: Mixed.

    Returns
    -------
    dm : :class:`numpy.ndarray`
        Averaged QMC density matrix.
    dm_err : :class:`numpy.array`
        Associated matrix of estimates for errors.
        .. Warning:
            This is a very simple estimate of the error and does not include
            temporal correlations.
    """
    md = get_metadata(filename)
    dm, weights = extract_qmc_dm(filename, name)
    dm_av = dm[skip:].mean(axis=0)
    if md['free_proj']:
        dm_av /= w.mean()
        dm_err = numpy.zeros(dm_av.shape)
    else:
        dm_err = scipy.stats.sem(numpy.real(dm[skip:]))
    return dm_av, dm_err

def get_qmc_dm_trace(filename, dm_name):
    """Get trace of spin up and down 1RDM as function tau. Useful for debugging.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    dm_name: string
        Estimator type. Mixed or BackPropagated. Default: Mixed.

    Returns
    -------
    dmu : :class:`numpy.ndarray`
        Trace of spin up QMC 1RDM as function of sample block.
    dmd : :class:`numpy.ndarray`
        Trace of spin down QMC 1RDM as function of sample block.
    """
    (rdm, weights) = extract_qmc_dm(filename, dm_name)
    dmu = []
    dmd = []
    for dm in rdm:
        dmu.append(dm[0].trace())
        dmd.append(dm[1].trace())
    return (dmu, dmd)
