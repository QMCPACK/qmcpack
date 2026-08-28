import pytest
from copy import deepcopy
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.PWSCF_ANALYZER)


from . import TEST_DIR
from ..testing import object_eq


def test_empty_init():
    from .. import pwscf_analyzer as pa_module
    from ..pwscf_analyzer import PwscfAnalyzer

    pa = PwscfAnalyzer()
    assert('results_out' in pa)
    assert(pa.results_out is None)
    assert(pa.results_xml is None)
    assert(all(value is False for value in pa.info.data_status.values()))
    free_helpers = ('match_float','read_kpoint_tables')
    for name in free_helpers:
        assert(callable(getattr(pa_module,name)))
        assert(not hasattr(PwscfAnalyzer,name))
    #end for
#end def test_empty_init


@pytest.mark.parametrize('calculation',('scf','nscf','bands','relax','vc-relax','md','vc-md'))
def test_result_initialization(tmp_path,calculation):
    from ..pwscf_analyzer import PwscfAnalyzer

    infile = tmp_path / '{}.in'.format(calculation)
    infile.write_text(
        "&CONTROL\n  calculation = '{}'\n  verbosity = 'low'\n/\n".format(calculation)
        )
    pa = PwscfAnalyzer(infile)

    expected = {
        'Ef','fermi_energies','bands',
        'volume',
        'cputime','walltime','kpoints_cart','kpoints_unit','kweights',
        'K',
        }
    if calculation not in ('nscf','bands'):
        expected.update((
            'E','relax_energies','scf_conv_energy','scf_conv_accuracy',
            'pressure','stress',
            'forces','tot_forces','max_forces',
            ))
    #end if
    if calculation in ('relax','vc-relax','md','vc-md'):
        expected.add('relax_structures')
    #end if
    if calculation in ('md','vc-md'):
        expected.update(('md_data','md_stats'))
    #end if

    assert(set(pa.results_out.keys())==expected)
    assert(all(value is None for value in pa.results_out.values()))
    assert(pa.results_xml is None)
    assert(all(value is False for value in pa.info.data_status.values()))
#end def test_result_initialization


def test_analyze():
    from numpy import array
    from ..developer import obj, to_obj
    from ..pwscf_analyzer import PwscfAnalyzer

    all_result_names = (
        'E','Ef','K','bands','cputime','relax_energies','scf_conv_energy',
        'scf_conv_accuracy','fermi_energies','forces',
        'kpoints_cart','kpoints_unit','kweights','max_forces','md_data',
        'md_stats','pressure','stress','relax_structures','tot_forces','volume',
        'walltime',
        )
    data_status_names = (
        'log','md','fermi_energies','scf_conv_energy','scf_conv_accuracy',
        'relax_energies','bands','relax_structures','pressure','volume',
        'stress','forces','total_forces','timing','kpoints','pw2casino','xml',
        )

    def empty_data_status():
        return obj({name:False for name in data_status_names})
    #end def empty_data_status

    def nest_results(reference,calculation):
        result_names = [
            'Ef','fermi_energies','bands',
            'volume',
            'cputime','walltime','kpoints_cart','kpoints_unit','kweights',
            'K',
            ]
        if calculation not in ('nscf','bands'):
            result_names.extend((
                'E','relax_energies','scf_conv_energy','scf_conv_accuracy',
                'pressure','stress',
                'forces','tot_forces','max_forces',
                ))
        #end if
        if calculation in ('relax','vc-relax','md','vc-md'):
            result_names.append('relax_structures')
        #end if
        if calculation in ('md','vc-md'):
            result_names.extend(('md_data','md_stats'))
        #end if
        results = obj({name:None for name in result_names})
        for name in all_result_names:
            if name in reference:
                results[name] = reference[name]
                del reference[name]
            #end if
        #end for
        reference.results_out = results
        reference.results_xml = None
        return reference
    #end def nest_results

    scf_path = TEST_DIR / "test_pwscf_analyzer_files/scf_output"
    relax_path = TEST_DIR / "test_pwscf_analyzer_files/relax_output"
    nscf_path = TEST_DIR / "test_pwscf_analyzer_files/nscf_output"

    # scf w/o actual analysis
    pa = PwscfAnalyzer(scf_path,'scf.in','scf.out')

    del pa.abspath
    del pa.path

    pa_ref = obj(
        infile_name     = 'scf.in',
        outfile_name    = 'scf.out',
        pw2c_outfile_name = None,
        info = obj(
            md_only         = False,
            warn            = False,
            ),
        input = obj(
            atomic_positions = obj(
                atoms           = ['C','C','C','C','C','C','C','C','C','C','C','C','C','C','C'],
                positions       = array(
                                  [[-0.08993686, -0.08993686, -0.08993686],
                                   [ 3.46309801,  3.46309801, -0.08993686],
                                   [ 5.05974172,  5.05974172,  1.70077347],
                                   [-0.08993686,  3.46309801,  3.46309801],
                                   [ 1.70077347,  5.05974172,  5.05974172],
                                   [ 3.37014984,  6.7493336 ,  3.37014983],
                                   [ 5.05974172,  8.41870996,  5.05974172],
                                   [ 3.46309801, -0.08993686,  3.46309801],
                                   [ 5.05974172,  1.70077347,  5.05974172],
                                   [ 6.7493336 ,  3.37014984,  3.37014984],
                                   [ 8.41870996,  5.05974172,  5.05974172],
                                   [ 3.37014984,  3.37014984,  6.7493336 ],
                                   [ 5.05974172,  5.05974172,  8.41870996],
                                   [ 6.7493336 ,  6.7493336 ,  6.7493336 ],
                                   [ 8.43290286,  8.43290286,  8.43290286]],dtype=float),
                specifier       = 'bohr',
                ),
            atomic_species = obj(
                atoms           = ['C'],
                specifier       = '',
                masses = obj(
                    C               = 12.011,
                    ),
                pseudopotentials = obj(
                    C               = 'C.BFD.upf',
                    ),
                ),
            cell_parameters = obj(
                specifier       = 'bohr',
                vectors         = array(
                                  [[ 6.74632229,  6.74632229,  0.        ],
                                   [-0.        ,  6.74632229,  6.74632229],
                                   [ 6.74632229, -0.        ,  6.74632229]],dtype=float),
                ),
            control = obj(
                calculation     = 'scf',
                outdir          = 'pwscf_output',
                prefix          = 'pwscf',
                pseudo_dir      = './',
                verbosity       = 'high',
                tstress         = True,
                tprnfor         = True,
                ),
            electrons = obj(
                conv_thr        = 1e-07,
                ),
            k_points = obj(
                grid            = array([2.,2.,2.],dtype=float),
                shift           = array([0.,0.,0.],dtype=float),
                specifier       = 'automatic',
                ),
            system = obj(
                ecutwfc         = 75.0,
                ibrav           = 0,
                input_dft       = 'lda',
                nat             = 15,
                nspin           = 1,
                ntyp            = 1,
                tot_charge      = 0.0,
                ),
            ),
        )

    pa_ref.info.data_status = empty_data_status()
    assert(object_eq(to_obj(pa),nest_results(pa_ref,'scf')))

    input_read = deepcopy(pa.input)


    # scf w/ full analysis
    pa = PwscfAnalyzer(scf_path,'scf.in','scf.out',analyze=True)

    assert(object_eq(pa.input,input_read))
    assert('md_data' not in pa.results_out)
    assert(len(pa.results_out.bands.up)==3)
    assert(len(pa.results_out.bands.down)==0)
    for band in pa.results_out.bands.up.values():
        assert(band.eigs.shape==(30,))
        assert(band.occs.shape==(30,))
    #end for
    assert(pa.info.data_status.relax_energies)
    assert(pa.info.data_status.scf_conv_energy)
    assert(pa.info.data_status.scf_conv_accuracy)
    assert(pa.info.data_status.bands)
    assert(pa.info.data_status.forces)
    assert(pa.info.data_status.total_forces)
    assert(pa.info.data_status.stress)
    assert(not pa.info.data_status.xml)
    assert('xml_status_failed' not in pa.info)
    assert(all(isinstance(value,bool) for value in pa.info.data_status.values()))
    del pa.info.data_status

    del pa.input
    del pa.abspath
    del pa.path
    pa.results_out.bands = None

    pa_ref = obj(
        E               = -170.11048381,
        Ef              = None,
        cputime         = 0.001175,
        relax_energies  = array([-170.11048381],dtype=float),
        scf_conv_energy = array(
            [-170.00165599,-170.10149211,-170.10929962,-170.11040360,
             -170.11046415,-170.11048301,-170.11048376],dtype=float),
        scf_conv_accuracy = array(
            [8.8226783e-01,3.0975790e-02,4.5667900e-03,2.2819000e-04,
             8.1210000e-05,1.3500000e-06,2.6000000e-07],dtype=float),
        fermi_energies  = None,
        forces          = array(
                          [[[-0.01852018, -0.01852018, -0.01852018],
                            [ 0.01852018,  0.01852018, -0.01852018],
                            [ 0.        ,  0.        , -0.00189264],
                            [-0.01852018,  0.01852018,  0.01852018],
                            [-0.00189264, -0.        ,  0.        ],
                            [ 0.00046488, -0.00046488,  0.00046488],
                            [-0.        ,  0.00189264,  0.        ],
                            [ 0.01852018, -0.01852018,  0.01852018],
                            [ 0.        , -0.00189264, -0.        ],
                            [-0.00046488,  0.00046488,  0.00046488],
                            [ 0.00189264,  0.        ,  0.        ],
                            [ 0.00046488,  0.00046488, -0.00046488],
                            [ 0.        , -0.        ,  0.00189264],
                            [-0.00046488, -0.00046488, -0.00046488],
                            [-0.        ,  0.        ,  0.        ]]],dtype=float),
        infile_name     = 'scf.in',
        kpoints_cart    = array(
                          [[ 0.       ,  0.       ,  0.       ],
                           [-0.3535534,  0.3535534, -0.3535534],
                           [ 0.       ,  0.       , -0.7071068]],dtype=float),
        kpoints_unit    = array([[ 0.,   0. ,  0. ],
                                 [ 0.,   0. , -0.5],
                                 [ 0.,  -0.5, -0.5]],dtype=float),
        kweights        = array([0.25,1.,  0.75],dtype=float),
        max_forces      = array([0.03207789],dtype=float),
        outfile_name    = 'scf.out',
        pressure        = -170.96,
        pw2c_outfile_name = None,
        stress          = [[-0.00116217, 0.0, 0.0, -170.96, 0.0, 0.0], 
                           [ 0.0, -0.00116217, -0.0, 0.0, -170.96, -0.0], 
                           [ 0.0, 0.0, -0.00116217, 0.0, 0.0, -170.96]],
        tot_forces      = array([0.064343],dtype=float),
        volume          = 614.0889,
        walltime        = 0.00164444444444,
        info = obj(
            md_only         = False,
            warn            = False,
            ),
        )

    assert(object_eq(to_obj(pa),nest_results(pa_ref,'scf')))


    # relax w/ full analysis
    pa = PwscfAnalyzer(relax_path,'relax.in','relax.out',analyze=True)

    assert('md_data' not in pa.results_out)
    assert(len(pa.results_out.bands.up)==3)
    assert(len(pa.results_out.bands.down)==0)
    for band in pa.results_out.bands.up.values():
        assert(band.eigs.shape==(30,))
        assert(band.occs.shape==(0,))
    #end for
    assert(pa.info.data_status.relax_energies)
    assert(pa.info.data_status.scf_conv_energy)
    assert(pa.info.data_status.scf_conv_accuracy)
    assert(pa.info.data_status.relax_structures)
    assert(pa.info.data_status.forces)
    assert(pa.info.data_status.total_forces)
    assert(all(isinstance(value,bool) for value in pa.info.data_status.values()))
    del pa.info.data_status

    del pa.input
    del pa.abspath
    del pa.path
    pa.results_out.bands = None

    pa_ref = obj(
        E               = -168.41267772,
        Ef              = None,
        cputime         = 0.00186111111111,
        relax_energies  = array([-168.38623938,-168.40640935,-168.41263281,-168.41267772],dtype=float),
        scf_conv_energy = array(
            [-168.30366565,-168.37073172,-168.38313345,-168.38614073,
             -168.38622231,-168.33164154,-168.39115701,-168.40610630,
             -168.40637046,-168.40640735,-168.40640922,-168.23494686,
             -168.35091632,-168.41063878,-168.41258264,-168.41263084,
             -168.41263272,-168.41148703,-168.41216839,-168.41266998,
             -168.41267747,-168.41267771],dtype=float),
        scf_conv_accuracy = array(
            [7.1052349e-01,5.1972060e-02,1.3904300e-02,2.6686000e-04,
             7.1900000e-05,1.9229067e-01,6.9241250e-02,2.9424000e-04,
             1.5621000e-04,5.7900000e-06,1.0600000e-06,3.9966070e-01,
             2.7132596e-01,4.7055900e-03,1.3613000e-04,8.9700000e-06,
             9.9000000e-07,2.4786000e-03,2.3529600e-03,9.5100000e-06,
             1.1500000e-06,5.0000000e-08],dtype=float),
        fermi_energies  = None,
        forces          = array(
                          [[[-4.625982e-02, -4.625982e-02, -4.625982e-02],
                            [ 4.625982e-02,  4.625982e-02, -4.625982e-02],
                            [ 0.000000e+00,  0.000000e+00,  2.650527e-02],
                            [-4.625982e-02,  4.625982e-02,  4.625982e-02],
                            [ 2.650527e-02,  0.000000e+00,  0.000000e+00],
                            [-2.041270e-03,  2.041270e-03, -2.041270e-03],
                            [ 0.000000e+00, -2.650527e-02,  0.000000e+00],
                            [ 4.625982e-02, -4.625982e-02,  4.625982e-02],
                            [ 0.000000e+00,  2.650527e-02,  0.000000e+00],
                            [ 2.041270e-03, -2.041270e-03, -2.041270e-03],
                            [-2.650527e-02, -0.000000e+00,  0.000000e+00],
                            [-2.041270e-03, -2.041270e-03,  2.041270e-03],
                            [ 0.000000e+00,  0.000000e+00, -2.650527e-02],
                            [ 2.041270e-03,  2.041270e-03,  2.041270e-03],
                            [ 0.000000e+00,  0.000000e+00,  0.000000e+00]],
                           [[-2.192410e-02, -2.192410e-02, -2.192410e-02],
                            [ 2.192410e-02,  2.192410e-02, -2.192410e-02],
                            [ 0.000000e+00,  0.000000e+00, -1.131570e-02],
                            [-2.192410e-02,  2.192410e-02,  2.192410e-02],
                            [-1.131570e-02,  0.000000e+00,  0.000000e+00],
                            [ 1.693610e-03, -1.693610e-03,  1.693610e-03],
                            [ 0.000000e+00,  1.131570e-02,  0.000000e+00],
                            [ 2.192410e-02, -2.192410e-02,  2.192410e-02],
                            [ 0.000000e+00, -1.131570e-02,  0.000000e+00],
                            [-1.693610e-03,  1.693610e-03,  1.693610e-03],
                            [ 1.131570e-02, -0.000000e+00,  0.000000e+00],
                            [ 1.693610e-03,  1.693610e-03, -1.693610e-03],
                            [-0.000000e+00,  0.000000e+00,  1.131570e-02],
                            [-1.693610e-03, -1.693610e-03, -1.693610e-03],
                            [ 0.000000e+00,  0.000000e+00,  0.000000e+00]],
                           [[ 5.327200e-04,  5.327200e-04,  5.327200e-04],
                            [-5.327200e-04, -5.327200e-04,  5.327200e-04],
                            [ 0.000000e+00,  0.000000e+00, -1.372200e-03],
                            [ 5.327200e-04, -5.327200e-04, -5.327200e-04],
                            [-1.372200e-03, -0.000000e+00,  0.000000e+00],
                            [-3.131240e-03,  3.131240e-03, -3.131240e-03],
                            [ 0.000000e+00,  1.372200e-03, -0.000000e+00],
                            [-5.327200e-04,  5.327200e-04, -5.327200e-04],
                            [-0.000000e+00, -1.372200e-03,  0.000000e+00],
                            [ 3.131240e-03, -3.131240e-03, -3.131240e-03],
                            [ 1.372200e-03,  0.000000e+00, -0.000000e+00],
                            [-3.131240e-03, -3.131240e-03,  3.131240e-03],
                            [ 0.000000e+00,  0.000000e+00,  1.372200e-03],
                            [ 3.131240e-03,  3.131240e-03,  3.131240e-03],
                            [ 0.000000e+00,  0.000000e+00,  0.000000e+00]],
                           [[ 6.653000e-05,  6.653000e-05,  6.653000e-05],
                            [-6.653000e-05, -6.653000e-05,  6.653000e-05],
                            [ 0.000000e+00,  0.000000e+00, -4.582100e-04],
                            [ 6.653000e-05, -6.653000e-05, -6.653000e-05],
                            [-4.582100e-04,  0.000000e+00,  0.000000e+00],
                            [ 8.294900e-04, -8.294900e-04,  8.294900e-04],
                            [ 0.000000e+00,  4.582100e-04,  0.000000e+00],
                            [-6.653000e-05,  6.653000e-05, -6.653000e-05],
                            [ 0.000000e+00, -4.582100e-04,  0.000000e+00],
                            [-8.294900e-04,  8.294900e-04,  8.294900e-04],
                            [ 4.582100e-04,  0.000000e+00, -0.000000e+00],
                            [ 8.294900e-04,  8.294900e-04, -8.294900e-04],
                            [ 0.000000e+00,  0.000000e+00,  4.582100e-04],
                            [-8.294900e-04, -8.294900e-04, -8.294900e-04],
                            [ 0.000000e+00,  0.000000e+00,  0.000000e+00]]],dtype=float),
        infile_name     = 'relax.in',
        max_forces      = array([0.08012436,0.03797366,0.00542347,0.00143672],dtype=float),
        outfile_name    = 'relax.out',
        pressure        = None,
        pw2c_outfile_name = None,
        stress          = None,
        tot_forces      = array([0.173046,0.081060,0.011505,0.003093],dtype=float),
        volume          = 614.0889,
        walltime        = 0.00251388888889,
        info = obj(
            md_only         = False,
            warn            = False,
            ),
        relax_structures = obj({
            0 : obj(
                atoms           = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
                positions       = array(
                                  [[-0.04625982, -0.04625982, -0.04625982],
                                   [ 3.41942097,  3.41942097, -0.04625982],
                                   [ 5.05974172,  5.05974172,  1.71308584],
                                   [-0.04625982,  3.41942097,  3.41942097],
                                   [ 1.71308584,  5.05974172,  5.05974172],
                                   [ 3.37111988,  6.74836356,  3.37111987],
                                   [ 5.05974172,  8.40639759,  5.05974172],
                                   [ 3.41942097, -0.04625982,  3.41942097],
                                   [ 5.05974172,  1.71308584,  5.05974172],
                                   [ 6.74836356,  3.37111988,  3.37111988],
                                   [ 8.40639759,  5.05974172,  5.05974172],
                                   [ 3.37111988,  3.37111988,  6.74836356],
                                   [ 5.05974172,  5.05974172,  8.40639759],
                                   [ 6.74836356,  6.74836356,  6.74836356],
                                   [ 8.43290286,  8.43290286,  8.43290286]],
                                   dtype=float),
                ),
            1 : obj(
                atoms           = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
                positions       = array(
                                  [[-0.09055704, -0.09055704, -0.09055704],
                                   [ 3.46371819,  3.46371819, -0.09055704],
                                   [ 5.05974172,  5.05974172,  1.70201536],
                                   [-0.09055704,  3.46371819,  3.46371819],
                                   [ 1.70201536,  5.05974172,  5.05974172],
                                   [ 3.37322753,  6.74625591,  3.37322752],
                                   [ 5.05974172,  8.41746807,  5.05974172],
                                   [ 3.46371819, -0.09055704,  3.46371819],
                                   [ 5.05974172,  1.70201536,  5.05974172],
                                   [ 6.74625591,  3.37322753,  3.37322753],
                                   [ 8.41746807,  5.05974172,  5.05974172],
                                   [ 3.37322753,  3.37322753,  6.74625591],
                                   [ 5.05974172,  5.05974172,  8.41746807],
                                   [ 6.74625591,  6.74625591,  6.74625591],
                                   [ 8.43290286,  8.43290286,  8.43290286]],
                                   dtype=float),
                ),
            2 : obj(
                atoms           = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
                positions       = array(
                                  [[-0.08993686, -0.08993686, -0.08993686],
                                   [ 3.46309801,  3.46309801, -0.08993686],
                                   [ 5.05974172,  5.05974172,  1.70077347],
                                   [-0.08993686,  3.46309801,  3.46309801],
                                   [ 1.70077347,  5.05974172,  5.05974172],
                                   [ 3.37014984,  6.7493336 ,  3.37014983],
                                   [ 5.05974172,  8.41870996,  5.05974172],
                                   [ 3.46309801, -0.08993686,  3.46309801],
                                   [ 5.05974172,  1.70077347,  5.05974172],
                                   [ 6.7493336 ,  3.37014984,  3.37014984],
                                   [ 8.41870996,  5.05974172,  5.05974172],
                                   [ 3.37014984,  3.37014984,  6.7493336 ],
                                   [ 5.05974172,  5.05974172,  8.41870996],
                                   [ 6.7493336 ,  6.7493336 ,  6.7493336 ],
                                   [ 8.43290286,  8.43290286,  8.43290286]],
                                   dtype=float),
                ),
            3 : obj(
                atoms           = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
                positions       = array(
                                  [[-0.08993686, -0.08993686, -0.08993686],
                                   [ 3.46309801,  3.46309801, -0.08993686],
                                   [ 5.05974172,  5.05974172,  1.70077347],
                                   [-0.08993686,  3.46309801,  3.46309801],
                                   [ 1.70077347,  5.05974172,  5.05974172],
                                   [ 3.37014984,  6.7493336 ,  3.37014983],
                                   [ 5.05974172,  8.41870996,  5.05974172],
                                   [ 3.46309801, -0.08993686,  3.46309801],
                                   [ 5.05974172,  1.70077347,  5.05974172],
                                   [ 6.7493336 ,  3.37014984,  3.37014984],
                                   [ 8.41870996,  5.05974172,  5.05974172],
                                   [ 3.37014984,  3.37014984,  6.7493336 ],
                                   [ 5.05974172,  5.05974172,  8.41870996],
                                   [ 6.7493336 ,  6.7493336 ,  6.7493336 ],
                                   [ 8.43290286,  8.43290286,  8.43290286]],
                                   dtype=float),
                ),
            }),
        )

    assert(object_eq(to_obj(pa),nest_results(pa_ref,'relax')))


    # nscf w/o actual analysis
    pa = PwscfAnalyzer(nscf_path,'nscf.in','nscf.out')

    del pa.abspath
    del pa.path

    pa_ref = obj(
        infile_name     = 'nscf.in',
        outfile_name    = 'nscf.out',
        pw2c_outfile_name = None,
        info = obj(
            md_only         = False,
            warn            = False,
            ),
        input = obj(
            atomic_positions = obj(
                atoms           = ['Sr', 'Co', 'O', 'O', 'O'],
                positions       = array(
                    [[ 0.        ,  0.        ,  0.        ],
                     [ 0.        , -1.06131318,  6.18578654],
                     [ 0.        , -3.62354986,  3.62354986],
                     [-2.56223668,  0.75046175,  4.37401161],
                     [ 2.56223668,  0.75046175,  4.37401161]],dtype=float),
                specifier       = 'alat',
                ),
            atomic_species = obj(
                atoms           = ['Co', 'O', 'Sr'],
                specifier       = '',
                masses = obj(
                    Co              = 58.933,
                    O               = 15.999,
                    Sr              = 87.956,
                    ),
                pseudopotentials = obj(
                    Co              = 'Co.opt.upf',
                    O               = 'O.opt.upf',
                    Sr              = 'Sr.opt.upf',
                    ),
                ),
            cell_parameters = obj(
                specifier       = 'alat',
                vectors         = array(
                    [[-5.12447336, -3.62354986,  3.62354986],
                     [ 5.12447336, -3.62354986,  3.62354986],
                     [ 0.        ,  5.12447336,  5.12447336]],dtype=float),
                ),
            control = obj(
                calculation     = 'nscf',
                outdir          = 'pwscf_output',
                prefix          = 'pwscf',
                pseudo_dir      = './',
                tprnfor         = True,
                tstress         = True,
                verbosity       = 'high',
                wf_collect      = True,
                ),
            electrons = obj(
                conv_thr        = 1e-08,
                electron_maxstep = 1000,
                mixing_beta     = 0.15,
                mixing_mode     = 'local-TF',
                ),
            k_points = obj(
                kpoints         = array(
                    [[ 0. ,  0. ,  0. ],
                     [-0. ,  0.5, -0. ],
                     [ 0.5,  0.5, -0. ],
                     [ 0.5,  0.5,  0.5]],dtype=float),
                nkpoints        = 4,
                specifier       = 'crystal',
                weights         = array([1., 1., 1., 1.],dtype=float),
                ),
            system = obj(
                degauss         = 0.001,
                ecutrho         = 1750.0,
                ecutwfc         = 350.0,
                ibrav           = 0,
                input_dft       = 'lda',
                lda_plus_u      = True,
                nat             = 5,
                nbnd            = 30,
                nosym           = True,
                nspin           = 2,
                ntyp            = 3,
                occupations     = 'smearing',
                smearing        = 'fermi-dirac',
                tot_charge      = 0.0,
                celldm = obj({
                    1               : 1.0,
                    }),
                hubbard_u = obj({
                    1               : 1.0,
                    }),
                starting_magnetization = obj({
                    1               : 1.0,
                    }),
                ),
            ),
        )

    pa_ref.info.data_status = empty_data_status()
    assert(object_eq(to_obj(pa),nest_results(pa_ref,'nscf')))

    input_read = deepcopy(pa.input)

    # nscf w/ analysis
    pa = PwscfAnalyzer(nscf_path,'nscf.in','nscf.out',analyze=True)

    assert(object_eq(pa.input,input_read))
    assert(pa.info.data_status.fermi_energies)
    assert(pa.info.data_status.bands)
    assert(pa.info.data_status.kpoints)
    assert(all(isinstance(value,bool) for value in pa.info.data_status.values()))
    del pa.info.data_status

    del pa.input
    del pa.abspath
    del pa.path

    pa_ref = obj(
        Ef              = 10.1198,
        cputime         = 0.076025,
        fermi_energies  = array([10.1198],dtype=float),
        infile_name     = 'nscf.in',
        kpoints_cart    = array(
            [[ 0.       ,  0.       ,  0.       ],
             [ 0.0487855, -0.0344966,  0.0344966],
             [ 0.       , -0.0689931,  0.0689931],
             [ 0.       , -0.0202076,  0.1177786]],dtype=float),
        kpoints_unit    = array(
            [[0. , 0. , 0. ],
             [0. , 0.5, 0. ],
             [0.5, 0.5, 0. ],
             [0.5, 0.5, 0.5]],dtype=float),
        kweights        = array([0.0041667, 0.0041667, 0.0041667, 0.0041667],dtype=float),
        outfile_name    = 'nscf.out',
        pw2c_outfile_name = None,
        volume          = 380.621,
        walltime        = 0.09731666666666666,
        info = obj(
            md_only         = False,
            warn            = False,
            ),
        bands = obj(
            electronic_structure = 'insulating',
            vbm = obj(
                band_number     = 24,
                energy          = 10.1131,
                index           = 1,
                kpoint_2pi_alat = array([ 0.0487855, -0.0344966,  0.0344966],dtype=float),
                kpoint_rel      = array([0.,  0.5, 0. ],dtype=float),
                pol             = 'up',
                ),
            cbm = obj(
                band_number     = 23,
                energy          = 10.159,
                index           = 0,
                kpoint_2pi_alat = array([0., 0., 0.],dtype=float),
                kpoint_rel      = array([0., 0., 0.],dtype=float),
                pol             = 'down',
                ),
            direct_gap = obj(
                energy          = 0.5689999999999991,
                index           = 2,
                kpoint_2pi_alat = array([ 0.,        -0.0689931,  0.0689931],dtype=float),
                kpoint_rel      = array([0.5, 0.5, 0. ],dtype=float),
                pol             = ['up', 'down'],
                ),
            indirect_gap = obj(
                energy          = 0.046,
                kpoints = obj(
                    cbm = obj(
                        band_number     = 23,
                        energy          = 10.159,
                        index           = 0,
                        kpoint_2pi_alat = array([0., 0., 0.],dtype=float),
                        kpoint_rel      = array([0., 0., 0.],dtype=float),
                        pol             = 'down',
                        ),
                    vbm = obj(
                        band_number     = 24,
                        energy          = 10.1131,
                        index           = 1,
                        kpoint_2pi_alat = array([ 0.0487855, -0.0344966,  0.0344966],dtype=float),
                        kpoint_rel      = array([0.,  0.5, 0. ],dtype=float),
                        pol             = 'up',
                        ),
                    ),
                ),
            up = obj({
                0 : obj(
                    eigs            = array(
                        [-86.6481, -50.4218, -50.4218, -50.4218, -22.2273,  -8.1615,  -6.5642,  -6.5642,
                         -4.3768,  -4.3768,  -4.3768,   5.6456,   5.6456,   5.6456,   6.1126,   6.1126,
                         6.1126,   7.3099,   7.3099,   8.4104,   8.4104,   8.4104,   9.1291,   9.1291,
                         9.1291,  15.3662,  15.3662,  16.3166,  18.0967,  18.0967],dtype=float),
                    index           = 0,
                    kpoint_2pi_alat = array([0., 0., 0.],dtype=float),
                    kpoint_rel      = array([0., 0., 0.],dtype=float),
                    occs            = array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 
                                             1., 1., 1., 1., 1., 1., 1., 1., 1.,
                                             1., 0., 0., 0., 0., 0.],dtype=float),
                    pol             = 'up',
                    ),
                1 : obj(
                    eigs            = array(
                        [-86.6481, -50.4228, -50.4217, -50.4217, -22.1938,  -8.2033,  -6.8047,  -6.5167,
                         -4.2982,  -4.2982,  -3.7851,   3.4922,   5.4009,   5.4009,   5.684 ,   5.684 ,
                         6.1333,   7.0502,   7.3176,   7.5721,   8.602 ,   8.602 ,   9.1755,   9.1755,
                         10.1131,  17.3214,  17.5708,  18.3183,  19.7307,  19.7307],dtype=float),
                    index           = 1,
                    kpoint_2pi_alat = array([ 0.0487855, -0.0344966,  0.0344966],dtype=float),
                    kpoint_rel      = array([0.,  0.5, 0. ],dtype=float),
                    occs            = array(
                        [1.,     1.,     1.,     1.,     1.,     1.,     1.,     1.,     1.,     1.,
                         1.,     1.,     1.,     1.,     1.,     1.,     1.,     1.,     1.,     1.,
                         1.,     1.,     1.,     1.,     0.6219, 0.,     0.,     0.,     0.,     0.],dtype=float),
                    pol             = 'up',
                    ),
                2 : obj(
                    eigs            = array(
                        [-86.6481, -50.4227, -50.4227, -50.4216, -22.1654,  -7.1894,  -7.1894,  -7.1664,
                         -4.2194,  -4.0301,  -4.0301,   3.275 ,   3.6075,   4.0541,   5.1708,   5.1708,
                         6.5471,   8.0257,   8.0257,   8.2978,   8.5616,   8.8243,   8.8243,   9.5042,
                         12.3636,  18.1392,  18.3196,  19.7474,  19.7474,  20.6336],dtype=float),
                    index           = 2,
                    kpoint_2pi_alat = array([ 0.,        -0.0689931,  0.0689931],dtype=float),
                    kpoint_rel      = array([0.5, 0.5, 0. ],dtype=float),
                    occs            = array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 
                                             1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                                             0., 0., 0., 0., 0., 0.],dtype=float),
                    pol             = 'up',
                    ),
                3 : obj(
                    eigs            = array(
                        [-86.6481, -50.4226, -50.4226, -50.4226, -22.1402,  -6.7619,  -6.7619,  -6.7619,
                         -4.7105,  -4.7105,  -4.7105,   3.0227,   3.618 ,   3.618 ,   4.709 ,   4.709 ,
                         4.709 ,   8.8386,   8.8386,   8.8386,   9.5927,   9.5927,   9.5927,  12.3737,
                         12.3737,  17.0409,  17.0409,  17.0409,  20.0775,  20.0775],dtype=float),
                    index           = 3,
                    kpoint_2pi_alat = array([ 0.,        -0.0202076,  0.1177786],dtype=float),
                    kpoint_rel      = array([0.5, 0.5, 0.5],dtype=float),
                    occs            = array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 
                                             1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0.,
                                             0., 0., 0., 0., 0., 0.],dtype=float),
                    pol             = 'up',
                    ),
                }),
            down = obj({
                0 : obj(
                    eigs            = array(
                        [-83.7428, -47.627 , -47.627 , -47.627 , -22.2142,  -7.8414,  -6.1549,  -6.1549, 
                         -4.359 ,  -4.359 ,  -4.359 ,   5.9202,   5.9202,   5.9202,   8.7802,   8.7802,
                         8.7802,   8.8094,   8.8094,   8.8094,   9.4806,   9.4806,   9.4806,  10.159 ,
                         10.159 ,  15.4284,  15.4284,  16.4943,  18.2226,  18.2226],dtype=float),
                    index           = 0,
                    kpoint_2pi_alat = array([0., 0., 0.],dtype=float),
                    kpoint_rel      = array([0., 0., 0.],dtype=float),
                    occs            = array(
                        [1. ,   1. ,   1. ,   1. ,   1. ,   1. ,   1. ,   1. ,   1. ,   1. ,   1. ,   1.,
                         1. ,   1. ,   1. ,   1. ,   1. ,   1. ,   1. ,   1. ,   1. ,   1. ,   1. ,   0.053,
                         0.053, 0. ,   0. ,   0. ,   0. ,   0.   ],dtype=float),
                    pol             = 'down',
                    ),
                1 : obj(
                    eigs            = array(
                        [-83.7428, -47.6283, -47.6269, -47.6269, -22.1812,  -7.915 ,  -6.4635,  -6.1049,
                            -4.2816,  -4.2816,  -3.7107,   4.6303,   5.9542,   5.9542,   7.0588,   7.0588,
                         7.3901,   7.9093,   8.8232,   8.9427,   8.9427,  10.1687,  10.6939,  10.6939,
                         12.1453,  17.4719,  17.6798,  18.3438,  19.8424,  19.8424],dtype=float),
                    index           = 1,
                    kpoint_2pi_alat = array([ 0.0487855, -0.0344966,  0.0344966],dtype=float),
                    kpoint_rel      = array([0.,  0.5, 0. ],dtype=float),
                    occs            = array(
                        [1. ,    1. ,    1. ,    1. ,    1. ,    1. ,    1. ,    1. ,    1. ,    1.,
                         1. ,    1. ,    1. ,    1. ,    1. ,    1. ,    1. ,    1. ,    1. ,    1.,
                         1. ,    0.0268, 0. ,    0. ,    0. ,    0. ,    0. ,    0. ,    0. ,    0. ],dtype=float),
                    pol             = 'down',
                    ),
                2 : obj(
                    eigs            = array(
                        [-83.7428, -47.6282, -47.6282, -47.6268, -22.1532,  -6.905 ,  -6.905 ,  -6.7976,
                         -4.2063,  -3.9506,  -3.9506,   4.0076,   4.9912,   5.0889,   6.6171,   6.6171,
                         6.7985,   8.3962,   8.3962,   9.935 ,  10.504 ,  10.5336,  10.5336,  10.984 ,
                         14.1185,  18.2049,  18.4398,  19.8756,  19.8756,  20.6738],dtype=float),
                    index           = 2,
                    kpoint_2pi_alat = array([ 0. ,       -0.0689931,  0.0689931],dtype=float),
                    kpoint_rel      = array([0.5, 0.5, 0. ],dtype=float),
                    occs            = array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 0.,
                                             0., 0., 0., 0., 0., 0.],dtype=float),
                    pol             = 'down',
                    ),
                3 : obj(
                    eigs            = array(
                        [-83.7428, -47.6281, -47.6281, -47.6281, -22.1282,  -6.418 ,  -6.418 ,  -6.418 ,
                         -4.6851,  -4.6851,  -4.6851,   3.2996,   5.1034,   5.1034,   5.9255,   5.9255,
                         5.9255,  10.0318,  10.0318,  10.0318,  10.8041,  10.8041,  10.8041,  14.1246,
                         14.1246,  17.1249,  17.1249,  17.1249,  20.1029,  20.1029],dtype=float),
                    index           = 3,
                    kpoint_2pi_alat = array([ 0.,        -0.0202076,  0.1177786],dtype=float),
                    kpoint_rel      = array([0.5, 0.5, 0.5],dtype=float),
                    occs            = array(
                        [1. ,    1. ,    1. ,    1. ,    1. ,    1. ,    1. ,    1. ,    1. ,    1. ,
                         1. ,    1. ,    1. ,    1. ,    1. ,    1. ,    1. ,    0.9985, 0.9985, 0.9985,
                         0. ,    0. ,    0. ,    0. ,    0. ,    0. ,    0. ,    0. ,    0. ,    0.    ],dtype=float),
                    pol             = 'down',
                    ),
                }),
            ),
        )

    assert(object_eq(to_obj(pa),nest_results(pa_ref,'nscf')))

#end def test_analyze


def test_modern_output(tmp_path):
    import numpy as np
    from ..pwscf_analyzer import PwscfAnalyzer

    infile = """&CONTROL
  calculation = 'vc-md'
  prefix = 'test'
  outdir = './tmp'
/
&SYSTEM
  ibrav = 1
  celldm(1) = 5.0
  nat = 1
  ntyp = 1
  ecutwfc = 10.0
/
ATOMIC_SPECIES
H 1.0 H.UPF
ATOMIC_POSITIONS crystal
H 0.0 0.0 0.0
K_POINTS gamma
"""
    outfile = """     unit-cell volume          =      125.0000 (a.u.)^3
     number of k points=   1

          k = 0.0000 0.0000 0.0000 (    10 PWs)   bands (ev):

    -1.0000  1.0000   trailing text

!    total energy              =      -1.10000000 Ry
     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =     0.01000000    0.02000000    0.03000000   trailing text

     Total force =     0.037417     Total SCF correction =     0.000000
          total   stress  (Ry/bohr**3)                   (kbar)     P=       10.00
   0.00100000   0.00000000   0.00000000          10.00        0.00        0.00   trailing text
   0.00000000   0.00100000   0.00000000           0.00       10.00        0.00
   0.00000000   0.00000000   0.00100000           0.00        0.00       10.00
     Entering Dynamics;  it =     1   time =  0.00000 pico-seconds
     Ekin =     0.10000000 Ry    T =  100.0 K  Etot =       -1.00000000
     new unit-cell volume =     124.00000 a.u.^3 (    18.00000 Ang^3 )
CELL_PARAMETERS (alat=  5.00000000)
   1.000000000   0.000000000   0.000000000   trailing text
   0.000000000   1.000000000   0.000000000
   0.000000000   0.000000000   1.000000000
ATOMIC_POSITIONS (crystal)
H  0.100000000  0.200000000  0.300000000  1 1 1  trailing text

     PWSCF        :      0.10s CPU      0.20s WALL
"""
    schema = """<?xml version="1.0"?>
<qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0">
  <output>
    <band_structure>
      <lsda>false</lsda>
      <nks>2</nks>
      <ks_energies>
        <k_point weight="0.5">0.0 0.0 0.0</k_point>
        <eigenvalues size="2">-0.5 0.5</eigenvalues>
        <occupations size="2">1.0 0.0</occupations>
      </ks_energies>
      <ks_energies>
        <k_point weight="0.25">0.0 0.0 0.0</k_point>
        <eigenvalues size="2">-0.4 0.6</eigenvalues>
        <occupations size="2">1.0 0.0</occupations>
      </ks_energies>
    </band_structure>
  </output>
</qes:espresso>
"""
    (tmp_path/'pwscf.in').write_text(infile)
    (tmp_path/'pwscf.out').write_text(outfile)
    savedir = tmp_path/'tmp'/'test.save'
    savedir.mkdir(parents=True)
    (savedir/'data-file-schema.xml').write_text(schema)

    pa = PwscfAnalyzer(tmp_path,'pwscf.in','pwscf.out',analyze=True,xml=True)

    assert(np.allclose(pa.results_out.md_data.total_energy,[-1.1]))
    assert(np.allclose(pa.results_out.md_data.time,[0.0]))
    assert(np.allclose(pa.results_out.md_data.kinetic_energy,[0.1]))
    assert(np.allclose(pa.results_out.md_data.temperature,[100.0]))
    assert(np.allclose(pa.results_out.tot_forces,[0.037417]))
    assert(pa.results_out.volume==124.0)
    assert(len(pa.results_out.bands.up)==1)
    assert(pa.results_out.bands.up[0].occs.shape==(0,))
    assert(np.allclose(pa.results_out.relax_structures[0].axes,5*np.eye(3)))
    assert(np.allclose(pa.results_out.relax_structures[0].positions,[[0.5,1.0,1.5]]))
    assert(pa.results_xml is not None)
    assert(not pa.results_xml.failed)
    assert(pa.results_xml.data.root.output.band_structure.nks==2)
    assert(len(pa.results_xml.kpoints)==2)
    assert(np.allclose(pa.results_xml.kpoints[1].up.eigenvalues,[-0.5,0.5]))
    assert(np.allclose(pa.results_xml.kpoints[2].up.eigenvalues,[-0.4,0.6]))
    assert(pa.info.data_status.md)
    assert(not pa.info.data_status.scf_conv_energy)
    assert(not pa.info.data_status.scf_conv_accuracy)
    assert(pa.info.data_status.bands)
    assert(pa.info.data_status.xml)
    assert(not pa.info.xml_status_failed)
    assert(all(isinstance(value,bool) for value in pa.info.data_status.values()))

    count_lines = pa.write_electron_counts().splitlines()
    assert(count_lines[1].split()==['1.50','0.00','0.75','0.75'])
    assert(count_lines[4].split()==['1','0.500000','3.00','0.00','1.50','1.50'])
    assert(count_lines[5].split()==['2','0.250000','1.50','0.00','0.75','0.75'])

    # Recognized but incomplete records are skipped without stopping analysis.
    malformed_tail = """
!    total energy              =      unavailable Ry
     number of k points= unavailable
     Total force = unavailable
          total   stress
   incomplete
     Forces acting on atoms
CELL_PARAMETERS (alat= 5.0)
   1.0 0.0
"""
    (tmp_path/'pwscf.out').write_text(outfile+malformed_tail)
    incomplete_schema = schema.replace(
        '<occupations size="2">1.0 0.0</occupations>',
        '',
        1,
        )
    (savedir/'data-file-schema.xml').write_text(incomplete_schema)

    pa = PwscfAnalyzer(tmp_path,'pwscf.in','pwscf.out',analyze=True,xml=True)

    assert(np.allclose(pa.results_out.md_data.total_energy,[-1.1]))
    assert(np.allclose(pa.results_out.relax_energies,[-1.1]))
    assert(pa.info.data_status.md)
    assert(pa.info.data_status.stress)
    assert(pa.info.data_status.forces)
    assert(pa.info.data_status.total_forces)
    assert(pa.info.data_status.relax_structures)
    assert(pa.results_xml.failed)
    assert(pa.info.data_status.xml)
    assert(pa.info.xml_status_failed)
    assert(len(pa.results_xml.kpoints)==1)

    # XML syntax errors are localized and do not discard parsed log data.
    (savedir/'data-file-schema.xml').write_text('<qes:espresso>')
    pa = PwscfAnalyzer(tmp_path,'pwscf.in','pwscf.out',analyze=True,xml=True)
    assert(np.allclose(pa.results_out.relax_energies,[-1.1]))
    assert(pa.results_xml.failed)
    assert(pa.results_xml.data is None)
    assert(not pa.info.data_status.xml)
    assert(pa.info.xml_status_failed)

    # Missing XML is represented by None rather than an empty XML result.
    (savedir/'data-file-schema.xml').unlink()
    pa = PwscfAnalyzer(tmp_path,'pwscf.in','pwscf.out',analyze=True,xml=True)
    assert(pa.results_xml is None)
    assert(not pa.info.data_status.xml)
    assert('xml_status_failed' not in pa.info)
#end def test_modern_output
