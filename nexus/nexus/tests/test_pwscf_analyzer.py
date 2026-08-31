import re
import pytest
from copy import deepcopy
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.PWSCF_ANALYZER)


from . import TEST_DIR
from ..testing import object_eq


@pytest.mark.parametrize('text',(
    '0','+7','-12','1.0','3.','.5','-.75','1e3','-2.5E-4',
    '+6D+02','7d-1',
    ))
def test_number_pattern_matches(text):
    from ..pwscf_analyzer import number_pattern

    assert(re.fullmatch(number_pattern,text) is not None)
#end def test_number_pattern_matches


@pytest.mark.parametrize('text',(
    '','.','+','1e','1.2.3','abc123','NaN','Inf','--1','1_000',
    ))
def test_number_pattern_rejects(text):
    from ..pwscf_analyzer import number_pattern

    assert(re.fullmatch(number_pattern,text) is None)
#end def test_number_pattern_rejects


@pytest.mark.parametrize('pattern_name,text,expected',(
    (
        'leading_number_list_pattern',
        '  -5.3120  1.2470  4.9810',
        dict(values='-5.3120  1.2470  4.9810'),
        ),
    (
        'leading_number_list_pattern',
        '0.0488-0.0345 0.0345  convergence achieved',
        dict(values='0.0488-0.0345 0.0345'),
        ),
    (
        'numeric_text_pattern',
        ' 0.0 0.5 -0.5 ',
        dict(values='0.0 0.5 -0.5'),
        ),
    (
        'numeric_text_pattern',
        '1.0D+00\n 0.0D+00',
        dict(values='1.0D+00\n 0.0D+00'),
        ),
    (
        'vector3_pattern',
        '  1.0 0.0 0.0',
        dict(x='1.0',y='0.0',z='0.0'),
        ),
    (
        'vector3_pattern',
        '-2.5D-01  .500000  0.  fourth-column',
        dict(x='-2.5D-01',y='.500000',z='0.'),
        ),
    (
        'kpoint_table_pattern',
        '     k( 1) = ( 0.0000000 0.0000000 0.0000000), wk = 0.2500000',
        dict(kx='0.0000000',ky='0.0000000',kz='0.0000000',weight='0.2500000'),
        ),
    (
        'kpoint_table_pattern',
        'k(12)=(.5 -.5 5.0D-1), wk=1.0D+00',
        dict(kx='.5',ky='-.5',kz='5.0D-1',weight='1.0D+00'),
        ),
    (
        'band_kpoint_pattern',
        '     k = 0.0000 0.0000 -0.7071 ( 2138 PWs) bands (ev):',
        dict(kx='0.0000',ky='0.0000',kz='-0.7071'),
        ),
    (
        'band_kpoint_pattern',
        'k = 0.0488-0.0345 0.0345 ( 42052 PWs) bands (ev):',
        dict(kx='0.0488',ky='-0.0345',kz='0.0345'),
        ),
    (
        'band_kpoint_pattern',
        'k =-0.3536 0.3536-0.3536 ( 2138 PWs) bands (ev):',
        dict(kx='-0.3536',ky='0.3536',kz='-0.3536'),
        ),
    (
        'total_energy_pattern',
        '     total energy              =    -168.12345678 Ry',
        dict(energy='-168.12345678'),
        ),
    (
        'total_energy_pattern',
        '!    total energy = -1.0D+02 Ry',
        dict(energy='-1.0D+02'),
        ),
    (
        'scf_accuracy_pattern',
        '     estimated scf accuracy    <       6.3E-09 Ry',
        dict(accuracy='6.3E-09'),
        ),
    (
        'scf_accuracy_pattern',
        'estimated scf accuracy = 1.2D-10 Ry',
        dict(accuracy='1.2D-10'),
        ),
    (
        'kinetic_energy_pattern',
        '     kinetic energy (Ekin) =           0.00285013 Ry',
        dict(kinetic_energy='0.00285013'),
        ),
    (
        'kinetic_energy_pattern',
        'kinetic energy=2.85D-03 Ry',
        dict(kinetic_energy='2.85D-03'),
        ),
    (
        'temperature_pattern',
        '     temperature           =         300.00000000 K ',
        dict(temperature='300.00000000'),
        ),
    (
        'temperature_pattern',
        'temperature=2.5D+02 K',
        dict(temperature='2.5D+02'),
        ),
    (
        'ekin_temperature_pattern',
        'Ekin = .00285 Ry T = 300.0 K Etot = -15.0 Ry',
        dict(kinetic_energy='.00285',temperature='300.0'),
        ),
    (
        'ekin_temperature_pattern',
        '  Ekin=2.8D-03 Ry  T=0.0 K',
        dict(kinetic_energy='2.8D-03',temperature='0.0'),
        ),
    (
        'md_time_pattern',
        '                           time      =   0.0002 pico-seconds',
        dict(time='0.0002'),
        ),
    (
        'md_time_pattern',
        'Entering Dynamics; it = 1 time = 0.00000 ps',
        dict(time='0.00000'),
        ),
    (
        'pressure_pattern',
        'total stress (Ry/bohr**3) (kbar) P= -170.96',
        dict(pressure='-170.96'),
        ),
    (
        'pressure_pattern',
        'total   stress  (Ry/bohr**3) (kbar) P = 1.2D+03 ',
        dict(pressure='1.2D+03'),
        ),
    (
        'volume_pattern',
        '     unit-cell volume          =     380.6210 (a.u.)^3',
        dict(volume='380.6210'),
        ),
    (
        'volume_pattern',
        'unit-cell volume=3.806210D+02 a.u.^3',
        dict(volume='3.806210D+02'),
        ),
    (
        'alat_pattern',
        'CELL_PARAMETERS (alat= 10.20)',
        dict(alat='10.20'),
        ),
    (
        'alat_pattern',
        'CELL_PARAMETERS (alat = 1.0D+01)',
        dict(alat='1.0D+01'),
        ),
    (
        'atomic_force_pattern',
        'atom 1 type 2 force = -0.001 0.002 0.000',
        dict(atom='1',type='2',fx='-0.001',fy='0.002',fz='0.000'),
        ),
    (
        'atomic_force_pattern',
        'atom 12 type 1 force=1.0D-03 -2.0D-03 .0',
        dict(atom='12',type='1',fx='1.0D-03',fy='-2.0D-03',fz='.0'),
        ),
    (
        'total_force_pattern',
        'Total force = 0.173046 Total SCF correction = 0.001',
        dict(total_force='0.173046'),
        ),
    (
        'total_force_pattern',
        'Total force=1.25D-04',
        dict(total_force='1.25D-04'),
        ),
    (
        'stress_row_pattern',
        ' -.001 0.0 .001 -147.1 0.0 147.1',
        dict(sxx='-.001',sxy='0.0',sxz='.001',kxx='-147.1',kxy='0.0',kxz='147.1'),
        ),
    (
        'stress_row_pattern',
        '1D-3 2D-3 3D-3 1E+2 2E+2 3E+2 trailing',
        dict(sxx='1D-3',sxy='2D-3',sxz='3D-3',kxx='1E+2',kxy='2E+2',kxz='3E+2'),
        ),
    (
        'timing_value_pattern',
        '4m33.69s',
        dict(value='4',unit='m'),
        ),
    (
        'timing_value_pattern',
        '1h 2m 3.5s',
        dict(value='1',unit='h'),
        ),
    (
        'fermi_energies_pattern',
        '     the Fermi energy is    10.1198 ev',
        dict(values='10.1198'),
        ),
    (
        'fermi_energies_pattern',
        'the Fermi energy          =      -3.22772442 eV',
        dict(values='-3.22772442'),
        ),
    (
        'fermi_energies_pattern',
        'the spin up/dw Fermi energies are 5.1 5.2 EV',
        dict(values='5.1 5.2'),
        ),
    ))
def test_tailored_pattern_matches(pattern_name,text,expected):
    from .. import pwscf_analyzer as pa_module

    pattern = getattr(pa_module,pattern_name)
    match   = re.search(pattern,text)
    assert(match is not None)
    for name,value in expected.items():
        assert(match.group(name)==value)
    #end for
#end def test_tailored_pattern_matches


@pytest.mark.parametrize('pattern_name,text',(
    ('leading_number_list_pattern','bands: 1.0 2.0'),
    ('leading_number_list_pattern','occupation numbers'),
    ('numeric_text_pattern','0.5-0.5'),
    ('numeric_text_pattern','1.0 Ry'),
    ('numeric_text_pattern','weight=1.0'),
    ('vector3_pattern','1.0 0.0'),
    ('vector3_pattern','1.0-0.5 0.0'),
    ('kpoint_table_pattern','k(1) = (0.0 0.0), wk = 1.0'),
    ('kpoint_table_pattern','k(1) = (0.0 0.0 0.0) weight = 1.0'),
    ('kpoint_table_pattern','k(1) = (0.0 0.0 0.0), wk = missing'),
    ('band_kpoint_pattern','k = 0.0 0.0 (2138 PWs) bands (ev):'),
    ('band_kpoint_pattern','k-point = 0.0 0.0 0.0 (2138 PWs)'),
    ('band_kpoint_pattern','k = 0.0 0.0 0.0 bands (ev):'),
    ('total_energy_pattern','total energy = -168.1 eV'),
    ('total_energy_pattern','total energy is -168.1 Ry'),
    ('scf_accuracy_pattern','estimated scf accuracy 6.3E-09 Ry'),
    ('scf_accuracy_pattern','estimated accuracy < 6.3E-09 Ry'),
    ('kinetic_energy_pattern','kinetic energy = 0.00285 eV'),
    ('kinetic_energy_pattern','Ekin = 0.00285 Ry'),
    ('temperature_pattern','Starting temperature = 300.0 K'),
    ('temperature_pattern','temperature = 300.0'),
    ('ekin_temperature_pattern','Ekin = .00285 Ry'),
    ('ekin_temperature_pattern','T = 300.0 K Ekin = .00285 Ry'),
    ('ekin_temperature_pattern','Ekin = .00285 Ry T = 300.0'),
    ('md_time_pattern','Entering Dynamics: iteration = 1'),
    ('md_time_pattern','time = 0.0002'),
    ('pressure_pattern','total stress pressure = -170.96'),
    ('pressure_pattern','total stress p = -170.96'),
    ('volume_pattern','new unit-cell volume: 380.6210'),
    ('volume_pattern','cell volume = 380.6210'),
    ('alat_pattern','CELL_PARAMETERS (alat: 10.20)'),
    ('alat_pattern','CELL_PARAMETERS (celldm(1)=10.20)'),
    ('atomic_force_pattern','atom 1 type 2 force = -0.001 0.002'),
    ('atomic_force_pattern','Total force = 0.173046'),
    ('atomic_force_pattern','atom 1 force = -0.001 0.002 0.000'),
    ('total_force_pattern','atom 1 type 2 force = 0.1 0.2 0.3'),
    ('total_force_pattern','total force = 0.173046'),
    ('stress_row_pattern','-.001 0.0 .001 -147.1 0.0'),
    ('stress_row_pattern','stress: -.001 0.0 .001 -147.1 0.0 147.1'),
    ('stress_row_pattern','-.001 0.0.001 -147.1 0.0 147.1'),
    ('timing_value_pattern','CPU time = 12'),
    ('timing_value_pattern','12ms'),
    ('fermi_energies_pattern','the Fermi energy is 10.1198'),
    ('fermi_energies_pattern','the Fermi energies are 5.1 5.2 5.3 eV'),
    ('fermi_energies_pattern','highest occupied level is 10.1198 eV'),
    ('fermi_energies_pattern','Fermi energy convergence was reached'),
    ))
def test_tailored_pattern_rejects(pattern_name,text):
    from .. import pwscf_analyzer as pa_module

    pattern = getattr(pa_module,pattern_name)
    assert(re.search(pattern,text) is None)
#end def test_tailored_pattern_rejects


@pytest.mark.parametrize('text,expected',(
    ('4m33.69s',(('4','m'),('33.69','s'))),
    ('1h 2m 3.5s',(('1','h'),('2','m'),('3.5','s'))),
    ('0.10s',(('0.10','s'),)),
    ))
def test_timing_value_sequence(text,expected):
    from ..pwscf_analyzer import timing_value_pattern

    values = tuple(
        (match.group('value'),match.group('unit'))
        for match in re.finditer(timing_value_pattern,text)
        )
    assert(values==expected)
#end def test_timing_value_sequence


def test_empty_init():
    from .. import pwscf_analyzer as pa_module
    from ..pwscf_analyzer import Pw2CasinoAnalyzer, PwscfAnalyzer, PwscfOutData

    pa = PwscfAnalyzer()
    assert(len(pa)==0)
    free_helpers = ('match_float',)
    for name in free_helpers:
        assert(callable(getattr(pa_module,name)))
        assert(not hasattr(PwscfAnalyzer,name))
    #end for
    assert(not hasattr(pa_module,'read_kpoint_tables'))
    reader_names = (
        'read_calculation','read_md','read_fermi_energies',
        'read_energies','read_scf_convergence','read_bands',
        'read_band_edges','read_structures','read_pressure_volume',
        'read_stress','read_forces','read_timing','read_kpoints',
        )
    assert(all(callable(getattr(PwscfOutData,name)) for name in reader_names))
    assert(not hasattr(PwscfOutData,'read'))
    assert(not any(hasattr(PwscfAnalyzer,'analyze_'+name[5:]) for name in reader_names))
    assert(not hasattr(PwscfOutData,'read_pw2casino'))
    assert(not hasattr(PwscfAnalyzer,'analyze_pw2casino'))
    assert(not hasattr(Pw2CasinoAnalyzer,'read'))
#end def test_empty_init


@pytest.mark.parametrize('calculation,log_text',(
    ('scf',     'Self-consistent Calculation\n'),
    ('nscf',    'Band Structure Calculation\nhighest occupied level (ev): 1.0\n'),
    ('bands',   'Band Structure Calculation\n'),
    ('relax',   'BFGS Geometry Optimization\n'),
    ('vc-relax','BFGS Geometry Optimization\nCELL_PARAMETERS (alat= 1.0)\n'),
    ('md',      'Molecular Dynamics Calculation\nEntering Dynamics: iteration = 1\n'),
    ('vc-md',   'Entering Dynamics; it = 1 time = 0.0 pico-seconds\n'),
    ))
def test_result_initialization(tmp_path,calculation,log_text):
    from ..pwscf_analyzer import PwscfOutData

    outfile = tmp_path / '{}.out'.format(calculation)
    outfile.write_text(log_text)
    out = PwscfOutData(outfile)

    expected = {
        'calculation',
        'Ef','fermi_energies','bands',
        'volume',
        'cputime','walltime','kpoints_cart','kpoints_unit','kweights',
        }
    if calculation in ('scf','relax','vc-relax','md','vc-md'):
        expected.update((
            'E','relax_energies','scf_conv_energy','scf_conv_accuracy',
            'pressure','stress','forces','tot_forces','max_forces',
            ))
    #end if
    if calculation in ('relax','vc-relax','md','vc-md'):
        expected.add('relax_structures')
    #end if
    if calculation in ('md','vc-md'):
        expected.update(('md_data','md_stats'))
    #end if

    assert(set(out.keys())==expected)
    assert(out.calculation==calculation)
    assert(all(value is None for name,value in out.items() if name!='calculation'))
#end def test_result_initialization


def test_pw2casino_analyzer_read(tmp_path):
    from ..pwscf_analyzer import Pw2CasinoAnalyzer, PwscfAnalyzer, PwscfOutData

    (tmp_path/'pwscf.in').write_text("&CONTROL\n  calculation = 'scf'\n/\n")
    (tmp_path/'pwscf.out').write_text('Self-consistent Calculation\n')
    (tmp_path/'pw2casino.out').write_text(
        'Kinetic energy from orbitals = 1.25D+01\n'
        )
    pw2casino = Pw2CasinoAnalyzer(tmp_path/'pw2casino.out')
    assert(pw2casino.K==12.5)

    pa = PwscfAnalyzer(
        tmp_path,
        'pwscf.in',
        'pwscf.out',
        'pw2casino.out',
        analyze = True,
        )

    assert(isinstance(pa.results_out,PwscfOutData))
    assert('K' not in pa.results_out)
    assert(isinstance(pa.pw2casino,Pw2CasinoAnalyzer))
    assert(pa.pw2casino.K==12.5)
    assert('data_status' not in pa.info)

    (tmp_path/'pw2casino.out').write_text('Kinetic energy is unavailable\n')
    pw2casino = Pw2CasinoAnalyzer(tmp_path/'pw2casino.out')
    assert(pw2casino.K is None)

    pw2casino = Pw2CasinoAnalyzer(tmp_path/'missing.out')
    assert(pw2casino.K is None)
#end def test_pw2casino_analyzer_read


def test_analyze():
    from numpy import array
    from ..developer import obj, to_obj
    from ..pwscf_analyzer import PwscfAnalyzer, PwscfOutData

    all_result_names = (
        'E','Ef','bands','cputime','relax_energies','scf_conv_energy',
        'scf_conv_accuracy','fermi_energies','forces',
        'kpoints_cart','kpoints_unit','kweights','max_forces','md_data',
        'md_stats','pressure','stress','relax_structures','tot_forces','volume',
        'walltime',
        )
    def nest_results(reference,calculation):
        reference.pw2casino = None
        if calculation is None:
            reference.results_out = None
            reference.results_xml = None
            return reference
        #end if
        result_names = [
            'Ef','fermi_energies','bands','volume','cputime','walltime',
            'kpoints_cart','kpoints_unit','kweights',
            ]
        if calculation in ('scf','relax','vc-relax','md','vc-md'):
            result_names.extend((
                'E','relax_energies','scf_conv_energy','scf_conv_accuracy',
                'pressure','stress','forces','tot_forces','max_forces',
                ))
        #end if
        if calculation in ('relax','vc-relax','md','vc-md'):
            result_names.append('relax_structures')
        #end if
        if calculation in ('md','vc-md'):
            result_names.extend(('md_data','md_stats'))
        #end if
        results = obj({name:None for name in result_names})
        results.calculation = calculation
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
    assert(pa.make_movie('unused.xyz') is None)

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

    assert(object_eq(to_obj(pa),nest_results(pa_ref,None)))

    input_read = deepcopy(pa.input)


    # scf w/ full analysis
    pa = PwscfAnalyzer(scf_path,'scf.in','scf.out',analyze=True)

    assert(object_eq(pa.input,input_read))
    assert(isinstance(pa.results_out,PwscfOutData))
    assert('md_data' not in pa.results_out)
    assert(len(pa.results_out.bands.up)==3)
    assert(len(pa.results_out.bands.down)==0)
    for band in pa.results_out.bands.up.values():
        assert(band.eigs.shape==(30,))
        assert(band.occs.shape==(30,))
    #end for
    assert(pa.results_out.relax_energies is not None)
    assert(pa.results_out.scf_conv_energy is not None)
    assert(pa.results_out.scf_conv_accuracy is not None)
    assert(pa.results_out.bands is not None)
    assert(pa.results_out.forces is not None)
    assert(pa.results_out.tot_forces is not None)
    assert(pa.results_out.stress is not None)
    assert('xml_status_failed' not in pa.info)

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
    assert(pa.results_out.relax_energies is not None)
    assert(pa.results_out.scf_conv_energy is not None)
    assert(pa.results_out.scf_conv_accuracy is not None)
    assert(pa.results_out.relax_structures is not None)
    assert(pa.results_out.forces is not None)
    assert(pa.results_out.tot_forces is not None)

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

    assert(object_eq(to_obj(pa),nest_results(pa_ref,None)))

    input_read = deepcopy(pa.input)

    # nscf w/ analysis
    pa = PwscfAnalyzer(nscf_path,'nscf.in','nscf.out',analyze=True)

    assert(object_eq(pa.input,input_read))
    assert(pa.results_out.fermi_energies is not None)
    assert(pa.results_out.bands is not None)
    assert(pa.results_out.kpoints_cart is not None)

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
    assert(pa.results_out.md_data is not None)
    assert(pa.results_out.scf_conv_energy is None)
    assert(pa.results_out.scf_conv_accuracy is None)
    assert(pa.results_out.bands is not None)
    assert(pa.results_xml is not None)
    assert(not pa.info.xml_status_failed)
    assert('data_status' not in pa.info)

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
    assert(pa.results_out.md_data is not None)
    assert(pa.results_out.stress is not None)
    assert(pa.results_out.forces is not None)
    assert(pa.results_out.tot_forces is not None)
    assert(pa.results_out.relax_structures is not None)
    assert(pa.results_xml.failed)
    assert(pa.info.xml_status_failed)
    assert(len(pa.results_xml.kpoints)==1)

    # XML syntax errors are localized and do not discard parsed log data.
    (savedir/'data-file-schema.xml').write_text('<qes:espresso>')
    pa = PwscfAnalyzer(tmp_path,'pwscf.in','pwscf.out',analyze=True,xml=True)
    assert(np.allclose(pa.results_out.relax_energies,[-1.1]))
    assert(pa.results_xml.failed)
    assert(pa.results_xml.data is None)
    assert(pa.info.xml_status_failed)

    # Missing XML is represented by None rather than an empty XML result.
    (savedir/'data-file-schema.xml').unlink()
    pa = PwscfAnalyzer(tmp_path,'pwscf.in','pwscf.out',analyze=True,xml=True)
    assert(pa.results_xml is None)
    assert('xml_status_failed' not in pa.info)

    # A total force is retained even when no atomic-force block is present.
    (tmp_path/'pwscf.out').write_text('     Total force = 0.125 Ry/au\n')
    pa = PwscfAnalyzer(tmp_path,'pwscf.in','pwscf.out',analyze=True)
    assert(pa.results_out.forces is None)
    assert(np.allclose(pa.results_out.tot_forces,[0.125]))
    assert('data_status' not in pa.info)
#end def test_modern_output
