
import versions
import testing
from testing import divert_nexus_log,restore_nexus_log
from testing import value_eq,object_eq,text_eq


def test_import():
    from qmcpack_analyzer import QmcpackAnalyzer
#end def test_import



def test_empty_init():
    from generic import obj
    from qmcpack_analyzer import QmcpackAnalyzer

    qa = QmcpackAnalyzer()

    qa_ref = obj(
        info = obj(
            analyzed        = False,
            data_loaded     = False,
            error           = None,
            failed          = False,
            initialized     = False,
            nindent         = 0,
            savefile        = '',
            savefilepath    = './',
            request = obj(
                calculations    = set([]),
                #data_sources    = set(['opt', 'stat', 'dmc', 'storeconfig', 'traces', 'scalar']),
                destination     = '.',
                dm_settings     = None,
                equilibration   = None,
                group_num       = None,
                methods         = set(['opt', 'rmc', 'dmc', 'vmc']),
                ndmc_blocks     = 1000,
                output          = set(['averages', 'samples']),
                quantities      = set(
                    ['mpc', 'localenergy', 'nonlocalecp', 'acceptratio', 
                     'spindensity', 'kinetic', 'blockweight', 'structurefactor',
                     'localecp', 'density', 'kecorr', 'energydensity', 
                     'localenergy_sq', 'blockcpu', 'dm1b', 'localpotential', 
                     'elecelec', 'ionion']),
                savefile        = '',
                source          = './qmcpack.in.xml',
                traces          = False,
                warmup_calculations = set([]),
                ),
            )
        )

    data_sources_ref = set(['opt', 'stat', 'dmc', 'storeconfig', 'traces', 'scalar'])

    req = qa.info.request
    data_sources = req.data_sources
    del req.data_sources

    assert(len(data_sources-data_sources_ref)==0)

    assert(object_eq(qa.to_obj(),qa_ref))
#end def test_empty_init    



def test_vmc_dmc_analysis():
    import os
    from generic import obj
    from qmcpack_analyzer import QmcpackAnalyzer

    tpath = testing.setup_unit_test_output_directory(
        test      = 'qmcpack_analyzer',
        subtest   = 'test_vmc_dmc_analysis',
        file_sets = ['diamond_gamma'],
        )

    # test load of vmc data
    infile = os.path.join(tpath,'diamond_gamma/vmc/vmc.in.xml')

    qa = QmcpackAnalyzer(infile)

    qa_keys = 'info vmc qmc wavefunction'.split()
    assert(set(qa.keys())==set(qa_keys))

    qmc = qa.qmc
    qmc_keys = [0]
    assert(set(qmc.keys())==set(qmc_keys))

    vmc = qa.qmc[0]
    vmc_keys = 'info scalars scalars_hdf'.split()
    assert(len(set(vmc.keys())-set(vmc_keys))==0)
    assert('scalars' in vmc)

    scalars = vmc.scalars
    scalars_keys = 'info method_info data'.split()
    assert(set(scalars.keys())==set(scalars_keys))

    data = scalars.data
    data_keys = '''
        AcceptRatio
        BlockCPU
        BlockWeight
        ElecElec
        IonIon
        Kinetic
        LocalECP
        LocalEnergy
        LocalEnergy_sq
        LocalPotential
        MPC
        NonLocalECP
        '''.split()
    assert(set(data.keys())==set(data_keys))


    # test analysis of vmc data
    qa = QmcpackAnalyzer(infile,analyze=True,equilibration=10)

    scalars = qa.qmc[0].scalars
    skeys = scalars_keys + data_keys + ['LocalEnergyVariance']
    assert(set(scalars.keys())==set(skeys))

    del scalars.data
    del scalars.info
    del scalars.method_info

    scalars_ref = obj(
        AcceptRatio = obj(
            error           = 0.00039521879435689295,
            kappa           = 1.0759731314133185,
            mean            = 0.7692893518516666,
            sample_variance = 1.3065205976562954e-05,
            ),
        BlockCPU = obj(
            error           = 0.00021580706746995528,
            kappa           = 17.761816682736388,
            mean            = 0.035813258793955555,
            sample_variance = 2.3598611606955374e-07,
            ),
        BlockWeight = obj(
            error           = 0.0,
            kappa           = 1.0,
            mean            = 600.0,
            sample_variance = 0.0,
            ),
        ElecElec = obj(
            error           = 0.006371421551195983,
            kappa           = 1.139934660946538,
            mean            = -2.7418612341177777,
            sample_variance = 0.003205053111939163,
            ),
        IonIon = obj(
            error           = 3.744889033148481e-16,
            kappa           = 1.0,
            mean            = -12.775667474000004,
            sample_variance = 1.262177448353619e-29,
            ),
        Kinetic = obj(
            error           = 0.023892472166667816,
            kappa           = 1.2234867804809575,
            mean            = 11.078538966677776,
            sample_variance = 0.04199188841333729,
            ),
        LocalECP = obj(
            error           = 0.028248848774549792,
            kappa           = 1.3238620182010674,
            mean            = -6.548159166556667,
            sample_variance = 0.054250193864959544,
            ),
        LocalEnergy = obj(
            error           = 0.0036722831372396907,
            kappa           = 1.9941956906300935,
            mean            = -10.458057094266668,
            sample_variance = 0.0006086211675753148,
            ),
        LocalEnergyVariance = obj(
            error           = 0.010209447771536196,
            kappa           = 1.0,
            mean            = 0.39878166800218146,
            sample_variance = 0.009380954141975286,
            ),
        LocalEnergy_sq = obj(
            error           = 0.07644640401544833,
            kappa           = 1.9594807148999176,
            mean            = 109.77034847611111,
            sample_variance = 0.26842047376171707,
            ),
        LocalPotential = obj(
            error           = 0.025040457075256945,
            kappa           = 1.2285412284015842,
            mean            = -21.53659606088889,
            sample_variance = 0.045934318559111634,
            ),
        MPC = obj(
            error           = 0.006576052786612818,
            kappa           = 1.094304580426513,
            mean            = -2.4782162843888895,
            sample_variance = 0.0035565987681342825,
            ),
        NonLocalECP = obj(
            error           = 0.009208836158892106,
            kappa           = 1.3807747631543386,
            mean            = 0.5290918141037778,
            sample_variance = 0.005527505216479378,
            ),
        )

    assert(object_eq(scalars.to_obj(),scalars_ref))

    
    # test analysis of dmc data
    infile = os.path.join(tpath,'diamond_gamma/dmc/dmc.in.xml')

    qa = QmcpackAnalyzer(infile,analyze=True,equilibration=5)

    qmc = qa.qmc
    qmc_keys = [0,1,2,3]
    assert(set(qmc.keys())==set(qmc_keys))

    le = obj()
    for n in qmc_keys:
        le[n] = qmc[n].scalars.LocalEnergy.mean
    #end for

    le_ref = obj({
        0 : -10.4891061838,
        1 : -10.531650024088888,
        2 : -10.530189168355555,
        3 : -10.528843362733333,
        })

    assert(object_eq(le,le_ref))

#end def test_vmc_dmc_analysis



def test_optimization_analysis():
    import os
    from numpy import array
    from generic import obj
    from qmcpack_analyzer import QmcpackAnalyzer
    
    tpath = testing.setup_unit_test_output_directory(
        test      = 'qmcpack_analyzer',
        subtest   = 'test_optimization_analysis',
        file_sets = ['diamond_gamma'],
        )

    divert_nexus_log()

    infile = os.path.join(tpath,'diamond_gamma/opt/opt.in.xml')

    qa = QmcpackAnalyzer(infile,analyze=True,equilibration=5)

    qa_keys = 'info wavefunction opt qmc results'.split()
    assert(set(qa.keys())==set(qa_keys))

    qmc = qa.qmc
    qmc_keys = range(6)
    assert(set(qmc.keys())==set(qmc_keys))

    le_opt = obj()
    for n in qmc_keys:
        le_opt[n] = qmc[n].scalars.LocalEnergy.mean
    #end for

    le_opt_ref = obj({
        0 : -10.449225495355556,
        1 : -10.454263886688889,
        2 : -10.459916961266668,
        3 : -10.4583030738,
        4 : -10.462984807955555,
        5 : -10.460860545622223,
        })

    assert(object_eq(le_opt,le_opt_ref))

    results = qa.results
    results_keys = ['optimization']
    assert(set(results.keys())==set(results_keys))

    opt = results.optimization

    opt_keys = '''
        all_complete
        any_complete
        energy
        energy_error
        energy_weight
        failed
        info
        optimal_file
        optimal_series
        optimal_wavefunction
        optimize
        opts
        series
        unstable
        variance
        variance_error
        variance_weight
        '''.split()

    assert(set(opt.keys())==set(opt_keys))

    from qmcpack_input import wavefunction

    assert(isinstance(opt.optimal_wavefunction,wavefunction))

    opt_wf = opt.optimal_wavefunction

    opt_wf_text = opt_wf.write()

    opt_wf_text_ref = '''
<wavefunction name="psi0" target="e">
   <sposet_builder type="bspline" href="../scf/pwscf_output/pwscf.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion0" version="0.1" meshfactor="1.0" precision="float" truncate="no">
      <sposet type="bspline" name="spo_ud" size="4" spindataset="0">                             </sposet>
   </sposet_builder>
   <determinantset>
      <slaterdeterminant>
         <determinant id="updet" group="u" sposet="spo_ud" size="4"/>
         <determinant id="downdet" group="d" sposet="spo_ud" size="4"/>
      </slaterdeterminant>
   </determinantset>
   <jastrow type="Two-Body" name="J2" function="bspline" print="yes">
      <correlation speciesA="u" speciesB="u" size="8" rcut="2.3851851232">
         <coefficients id="uu" type="Array">            
0.2576630369 0.1796686015 0.1326653657 0.09407180823 0.06267013118 0.03899100023 
0.02070235604 0.009229775746
         </coefficients>
      </correlation>
      <correlation speciesA="u" speciesB="d" size="8" rcut="2.3851851232">
         <coefficients id="ud" type="Array">            
0.4385891515 0.3212399072 0.2275448261 0.1558506324 0.1009589176 0.06108433554 
0.03154274436 0.01389485975
         </coefficients>
      </correlation>
   </jastrow>
</wavefunction>
'''

    opt_wf_text_ref = opt_wf_text_ref.replace('"',' " ')
    opt_wf_text     = opt_wf_text.replace('"',' " ')

    assert(text_eq(opt_wf_text,opt_wf_text_ref))

    del opt.info
    del opt.opts
    del opt.optimal_wavefunction

    opt_ref = obj(
        all_complete    = True,
        any_complete    = True,
        energy          = array([-10.45426389,-10.45991696,-10.45830307,
                                 -10.46298481,-10.46086055],
                                dtype=float),
        energy_error    = array([0.00320561,0.00271802,0.00298106,
                                 0.00561322,0.00375811],
                                dtype=float),
        energy_weight   = None,
        failed          = False,
        optimal_file    = 'opt.s003.opt.xml',
        optimal_series  = 3,
        optimize        = (1.0, 0.0),
        series          = array([1,2,3,4,5],dtype=int),
        unstable        = False,
        variance        = array([0.40278743,0.39865602,0.38110459,
                                 0.38927957,0.39354343],
                                dtype=float),
        variance_error  = array([0.01716415,0.00934316,0.00529809,
                                 0.01204068,0.00913372],
                                dtype=float),
        variance_weight = None,
        )

    assert(object_eq(opt.to_obj(),opt_ref,atol=1e-8))

    restore_nexus_log()
#end def test_optimization_analysis



def test_twist_average_analysis():
    import os
    from generic import obj
    from qmcpack_analyzer import QmcpackAnalyzer

    tpath = testing.setup_unit_test_output_directory(
        test      = 'qmcpack_analyzer',
        subtest   = 'test_twist_average_analysis',
        file_sets = ['diamond_twist'],
        )

    # test analysis of twist averaged vmc data
    infile = os.path.join(tpath,'diamond_twist/vmc/vmc.in')

    qa = QmcpackAnalyzer(infile,analyze=True,equilibration=5)

    qa_keys = 'info wavefunction vmc qmc bundled_analyzers'.split()
    assert(set(qa.keys())==set(qa_keys))

    ba = qa.bundled_analyzers
    ba_keys = [0,1,2,3]
    assert(set(ba.keys())==set(ba_keys))
    for n in ba_keys:
        assert(isinstance(ba[n],QmcpackAnalyzer))
    #end for

    qa2 = ba[2]
    qa2_keys = 'info wavefunction vmc qmc'.split()
    assert(set(qa2.keys())==set(qa2_keys))

    le_twists = obj()
    for n in ba_keys:
        le_twists[n] = ba[n].qmc[0].scalars.LocalEnergy.mean
    #end for

    le_twists_ref = obj({
        0 : -10.477644724210526,
        1 : -11.54064069741053,
        2 : -11.547046178357895,
        3 : -11.809361818642103,
        })

    assert(object_eq(le_twists,le_twists_ref))

    le_avg     = qa.qmc[0].scalars.LocalEnergy.mean
    le_avg_ref = -11.343673354655264
    assert(value_eq(le_avg,le_avg_ref))
    
#end def test_twist_average_analysis



if versions.h5py_available:
    def test_density_analysis():
        import os
        from numpy import array
        from generic import obj
        from qmcpack_analyzer import QmcpackAnalyzer

        tpath = testing.setup_unit_test_output_directory(
            test      = 'qmcpack_analyzer',
            subtest   = 'test_density_analysis',
            file_sets = ['diamond_gamma'],
            )

        infile = os.path.join(tpath,'diamond_gamma/dmc/dmc.in.xml')

        qa = QmcpackAnalyzer(infile,analyze=True,equilibration=5)

        qmc = qa.qmc[0]

        qmc_keys = 'info data scalars scalars_hdf SpinDensity'.split()
        assert(set(qmc.keys())==set(qmc_keys))

        sd = qmc.SpinDensity

        sd_keys = 'info method_info data u d'.split()
        assert(set(sd.keys())==set(sd_keys))

        d_keys = 'mean error'.split()
        assert(set(sd.u.keys())==set(d_keys))
        assert(set(sd.d.keys())==set(d_keys))

        d_mean = sd.u.mean + sd.d.mean

        d_mean_ref = array(
                      [[[0.72833333, 0.604     , 0.514     ],
                       [0.63666667, 0.19533333, 0.19683333],
                       [0.51166667, 0.19183333, 0.4175    ]],
                      [[0.6065    , 0.19533333, 0.19      ],
                       [0.19066667, 0.052     , 0.129     ],
                       [0.18783333, 0.13433333, 0.10483333]],
                      [[0.50633333, 0.19633333, 0.42483333],
                       [0.1775    , 0.12966667, 0.09983333],
                       [0.417     , 0.106     , 0.15583333]]],
                       dtype=float)

        assert(d_mean.shape==(3,3,3))
        assert(value_eq(d_mean,d_mean_ref))

    #end def test_density_analysis
#end if

