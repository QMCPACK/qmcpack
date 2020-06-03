
import testing
from testing import value_eq,object_eq,text_eq


def test_import():
    from pwscf_analyzer import PwscfAnalyzer
#end def test_import



def test_empty_init():
    from pwscf_analyzer import PwscfAnalyzer

    PwscfAnalyzer()
#end def test_empty_init



def test_analyze():
    import os
    from numpy import array,ndarray
    from generic import obj
    from pwscf_analyzer import PwscfAnalyzer

    tpath = testing.setup_unit_test_output_directory(
        test      = 'pwscf_analyzer',
        subtest   = 'test_analyze',
        file_sets = ['scf_output','relax_output','nscf_output'],
        )

    scf_path = os.path.join(tpath,'scf_output')
    relax_path = os.path.join(tpath,'relax_output')
    nscf_path = os.path.join(tpath,'nscf_output')

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
            xml             = False,
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

    assert(object_eq(pa.to_obj(),pa_ref))

    input_read = pa.input.copy()

    # scf w/ full analysis
    pa = PwscfAnalyzer(scf_path,'scf.in','scf.out',analyze=True)

    assert(object_eq(pa.input,input_read))

    del pa.input
    del pa.abspath
    del pa.path

    # Note: band read is failing for spin unpolarized case.
    # Test needed for spin polarized case w/ bands, then fix unpolarized.
    pa_ref = obj(
        E               = -170.11048381,
        Ef              = 0.0,
        cputime         = 0.001175,
        energies        = array([-170.11048381],dtype=float),
        fermi_energies  = array([],dtype=float),
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
            xml             = False,
            ),
        md_data = obj(
            kinetic_energy   = array([-170.96],dtype=float),
            potential_energy = array([0.84951619],dtype=float),
            pressure         = array([-170.96],dtype=float),
            temperature      = array([0.],dtype=float),
            time             = array([0.],dtype=float),
            total_energy     = array([-170.11048381],dtype=float),
            ),
        md_stats = obj(
            kinetic_energy   = (-170.96, 0.0),
            potential_energy = (0.8495161900000028, 0.0),
            pressure         = (-170.96, 0.0),
            temperature      = (0.0, 0.0),
            time             = (0.0, 0.0),
            total_energy     = (-170.11048381, 0.0),
            ),
        )

    assert(object_eq(pa.to_obj(),pa_ref))


    # relax w/ full analysis
    pa = PwscfAnalyzer(relax_path,'relax.in','relax.out',analyze=True)

    del pa.input
    del pa.abspath
    del pa.path

    pa_ref = obj(
        E               = -168.41267772,
        Ef              = 0.0,
        cputime         = 0.00186111111111,
        energies        = array([-168.38623938,-168.40640935,-168.41263281,-168.41267772],dtype=float),
        fermi_energies  = array([],dtype=float),
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
        pressure        = 0.0,
        pw2c_outfile_name = None,
        stress          = [],
        tot_forces      = array([],dtype=float),
        volume          = 614.0889,
        walltime        = 0.00251388888889,
        info = obj(
            md_only         = False,
            warn            = False,
            xml             = False,
            ),
        structures = obj({
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

    assert(object_eq(pa.to_obj(),pa_ref))


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
            xml             = False,
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
    
    assert(object_eq(pa.to_obj(),pa_ref))

    input_read = pa.input.copy()

    # nscf w/ analysis
    pa = PwscfAnalyzer(nscf_path,'nscf.in','nscf.out',analyze=True)

    assert(object_eq(pa.input,input_read))

    del pa.input
    del pa.abspath
    del pa.path

    pa_ref = obj(
        E               = 0.0,
        Ef              = 10.1198,
        cputime         = 0.076025,
        energies        = array([],dtype=float),
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
        pressure        = 0.0,
        pw2c_outfile_name = None,
        stress          = [],
        volume          = 380.621,
        walltime        = 0.09731666666666666,
        info = obj(
            md_only         = False,
            warn            = False,
            xml             = False,
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

    assert(object_eq(pa.to_obj(),pa_ref))

#end def test_analyze



