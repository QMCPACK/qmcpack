import h5py
import numpy
from afqmctools.utils.io import from_qmcpack_complex

# Map user names to QMCPACK names.
MAP = {
        'one_rdm': {
            'group': 'FullOneRDM',
            'numer': 'one_rdm'
        },
        'two_rdm': {
            'group': 'FullTwoRDM',
            'numer': 'two_rdm'
        },
        'diag_two_rdm': {
            'group': 'DiagTwoRDM',
            'numer': 'diag_two_rdm'
        },
        'on_top_pdm': {
            'group': 'N2R',
            'numer': 'n2r'
        },
        'cc_realspace_correlation': {
            'group': 'OD2RDM/CC',
            'numer': 'od2rdm_CC'
        },
        'ss_realspace_correlation': {
            'group': 'OD2RDM/SS',
            'numer': 'od2rdm_SS'
        },
        'cc_atom_correlation': {
            'group': 'ATOM_CORRELATORS/CC',
            'numer': 'correlator2D_CC'
        },
        'ss_atom_correlation': {
            'group': 'ATOM_CORRELATORS/SS',
            'numer': 'correlator2D_SS'
        },
        'c_atom_correlation': {
            'group': 'ATOM_CORRELATORS/C',
            'numer': 'correlator1D_C'
        },
        's_atom_correlation': {
            'group': 'ATOM_CORRELATORS/S',
            'numer': 'correlator1D_S'
        },
        'm_atom_correlation': {
            'group': 'ATOM_CORRELATORS/M',
            'numer': 'correlator1D_M'
        },
        'gen_fock_plus': {
            'group': 'GenFockPlus',
            'numer': 'gfockp'
        },
        'gen_fock_minus': {
            'group': 'GenFockMinus',
            'numer': 'gfockm'
        }
    }

def extract_data(filename, group, estimator, sample=None):
    """Extract data from HDF5 file.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    group : string
        Path to estimator.
    estimator : string
        Estimator to analyse.
    sample : int
        Sample to extract. Optional. Default None (return everything).

    Returns
    -------
    numer : :class:`numpy.ndarray`
        Numerator of estimator.
    denom : :class:`numpy.ndarray`
        Denominator of estimator.
    """
    with h5py.File(filename, 'r') as fh5:
        dsets = list(fh5[group].keys())
        denom_id = [d for d in dsets if 'denominator' in d]
        numer_id = [d for d in dsets if estimator in d]
        if sample is not None:
            assert sample < len(numer_id)
            numer = from_qmcpack_complex(fh5[group][numer_id[sample]][:])
            denom = from_qmcpack_complex(fh5[group][denom_id[sample]][:])
        else:
            numer = numpy.array([from_qmcpack_complex(fh5[group][d][:]) for d in numer_id])
            denom = numpy.array([from_qmcpack_complex(fh5[group][d][:])[0] for d in denom_id])
        return numer, denom

def extract_observable(filename, estimator='back_propagated',
                       name='one_rdm', ix=None, sample=None):
    """Extract observable from HDF5 file.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    name : string
        Name of observable (see estimates.py for list).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).
    sample : int
        Sample to extract. Optional. Default None (return everything).

    Returns
    -------
    obs : :class:`numpy.ndarray`
        Observable for a single sample or the full set of samples. Note if using
        free projection the numerator and denominator are returned separately.
    """
    sym_md = get_metadata(filename)
    free_proj = sym_md['FreeProjection']
    if estimator == 'back_propagated':
        base = 'Observables/BackPropagated/'
        if ix is None:
            bp_md = get_metadata(filename, path=base)
            ix = bp_md['NumAverages'] - 1
        ename = MAP[name]
        base += ename['group'] + '/Average_{}'.format(ix)
    elif estimator == 'mixed':
        base = 'Observables/Mixed/'
        ename = MAP[name]
        base += ename['group']
    else:
        print("Unknown estimator type: {} ".format(estimator))
        return None
    numer, denom = extract_data(filename, base, ename['numer'], sample=sample)
    if free_proj:
        return (numer, denom)
    else:
        # Use array broadcasting to divide by weights.
        return numer / denom[:,None]

def get_metadata(filename, path=''):
    """Extract QMC estimator metadata from h5 file.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    path : string
        Path to metadata in filename.

    Returns
    -------
    md : dict
        Estimator metadata.
    """
    md = {}
    with h5py.File(filename, 'r') as fh5:
        for k, v in fh5[path+'Metadata'].items():
            try:
                # Array valued data.
                md[k] = v[:]
            except ValueError:
                # Scalar data.
                md[k] = v[()]
            except KeyError:
                print("Error: path does not exist.")
                md = None
    return md

def get_estimator_len(filename, name='Mixed'):
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
        dms = fh5[name].keys()
        # Block size of RDM is not necessarily known so just be dumb and count.
        num_dm = len([n for n in dms if 'denominator' in n])
    return num_dm
