import pytest

from . import TEST_DIR, NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.PWSCF_ANALYZER)


@pytest.mark.parametrize(
    argnames='text,expected',
    argvalues=(
        ('0',0.0),('+7',7.0),('-12',-12.0),('3.',3.0),('.5',0.5),
        ('-2.5E-4',-2.5e-4),('+6D+02',600.0),(' 7d-1 ',0.7),
        ),
    )
def test_parse_float(text,expected):
    from ..pwscf_analyzer import parse_float

    assert(parse_float(text)==expected)
#end def test_parse_float


@pytest.mark.parametrize(
    argnames='text',
    argvalues=(
        '','.','+','1e','1.2.3','abc123','NaN','Inf','--1','1_000',
        ),
    )
def test_parse_float_rejects(text):
    from ..pwscf_analyzer import parse_float

    assert(parse_float(text) is None)
#end def test_parse_float_rejects


@pytest.mark.parametrize(
    argnames='text,expected',
    argvalues=(
        ('the Fermi energy is 10.1198 eV',[10.1198]),
        ('the Fermi energy = -3.22772442 eV',[-3.22772442]),
        ('the spin up/dw Fermi energies are 5.1 5.2 EV',[5.1,5.2]),
        ),
    )
def test_fermi_energy_parsing(tmp_path,text,expected):
    import numpy as np

    from ..pwscf_analyzer import PwscfOutData

    outfile = tmp_path/'pwscf.out'
    outfile.write_text(f'Self-consistent Calculation\n{text}\n')
    assert(np.allclose(PwscfOutData(outfile).fermi_energies,expected))
#end def test_fermi_energy_parsing


@pytest.mark.parametrize(
    argnames='text',
    argvalues=(
        'the Fermi energy is 10.1198',
        'the Fermi energies are 5.1 5.2 5.3 eV',
        'highest occupied level is 10.1198 eV',
        'Fermi energy convergence was reached',
        ),
    )
def test_fermi_energy_rejects(tmp_path,text):
    from ..pwscf_analyzer import PwscfOutData

    outfile = tmp_path/'pwscf.out'
    outfile.write_text(f'Self-consistent Calculation\n{text}\n')
    assert(PwscfOutData(outfile).fermi_energies is None)
#end def test_fermi_energy_rejects


def test_empty_init():
    from .. import pwscf_analyzer as pa_module
    from ..pwscf_analyzer import Pw2CasinoAnalyzer, PwscfAnalyzer, PwscfOutData

    pa = PwscfAnalyzer()
    assert(len(pa)==0)
    with pytest.raises(
        RuntimeError,
        match=r'PWSCF output file name is not available',
        ):
        pa.analyze()
    free_helpers = ('parse_float',)
    for name in free_helpers:
        assert(callable(getattr(pa_module,name)))
        assert(not hasattr(PwscfAnalyzer,name))
    #end for
    assert(not hasattr(pa_module,'read_kpoint_tables'))
    reader_names = (
        'read_calculation','read_fermi_energies',
        'read_energies','read_scf_convergence','read_bands','read_band_edges',
        'read_structures','read_pressure','read_volume','read_stress',
        'read_forces','read_timing','read_kpoints',
        )
    assert(all(callable(getattr(PwscfOutData,name)) for name in reader_names))
    assert(not hasattr(PwscfOutData,'read'))
    assert(not any(hasattr(PwscfAnalyzer,'analyze_'+name[5:]) for name in reader_names))
    assert(not hasattr(PwscfAnalyzer,'analyze_schema_xml'))
    assert(callable(Pw2CasinoAnalyzer))
    assert(not hasattr(Pw2CasinoAnalyzer,'read'))
    assert(not hasattr(pa_module,'PwscfXmlData'))
#end def test_empty_init


@pytest.mark.parametrize(
    argnames='calculation,log_text',
    argvalues=(
        ('scf',     'Self-consistent Calculation\n'),
        ('nscf',    'Band Structure Calculation\nhighest occupied level (ev): 1.0\n'),
        ('relax',   'BFGS Geometry Optimization\n'),
        ('vc-relax','BFGS Geometry Optimization\nCELL_PARAMETERS (alat= 1.0)\n'),
        ),
    )
def test_result_initialization(tmp_path,calculation,log_text):
    from ..pwscf_analyzer import PwscfOutData

    outfile = tmp_path / f'{calculation}.out'
    outfile.write_text(log_text)
    out = PwscfOutData(outfile)

    expected = {
        'calculation',
        'Ef','fermi_energies','bands',
        'volume','cputime','walltime',
        'kpoints_cart','kpoints_unit','kweights',
        }
    if calculation in {'scf','relax','vc-relax'}:
        expected.update((
            'E','relax_energies','scf_conv_energy','scf_conv_accuracy',
            'pressure','stress','forces','tot_forces','max_forces',
            ))
    #end if
    if calculation in {'relax','vc-relax'}:
        expected.add('relax_structures')
    #end if

    assert(set(out.keys())==expected)
    assert(out.calculation==calculation)
    assert(all(value is None for name,value in out.items() if name!='calculation'))
#end def test_result_initialization


def test_tokenized_log_parsing(tmp_path):
    import numpy as np

    from ..pwscf_analyzer import PwscfOutData

    scf_file = tmp_path/'scf.out'
    scf_file.write_text('''\
Self-consistent Calculation
number of atoms/cell   = 2 trailing tokens
number of k points = 1 trailing tokens
cart. coord.
k(1) = (0.0 0.0 0.0), wk = 1.0 trailing
cryst. coord.
k(1) = (0.0 0.0 0.0), wk = 1.0 trailing
!  total   energy = -1.0D+02 Ry trailing tokens
     total energy              =    -9.9D+01 Ry
     estimated scf accuracy < 6.3E-09 Ry
unit-cell volume = 3.806210D+02 (a.u.)^3
total stress (Ry/bohr**3) (kbar) P = 1.2D+03 trailing tokens
 1D-3 2D-3 3D-3 1E+2 2E+2 3E+2 trailing
 4D-3 5D-3 6D-3 4E+2 5E+2 6E+2 trailing
 7D-3 8D-3 9D-3 7E+2 8E+2 9E+2 trailing
Forces acting on atoms
atom 1 type 1 force = 0.1 0.2 0.3
atom 2 type 1 force = 0.4 0.5 0.6
Total force = 1.25D-04 Total SCF correction = 0.0
PWSCF        : 1h 2m 3.5s CPU 4m33.69s WALL
''')
    scf = PwscfOutData(scf_file)

    assert(scf.E==-100.0)
    assert(np.allclose(scf.relax_energies,[-100.0]))
    assert(np.allclose(scf.scf_conv_energy,[-99.0]))
    assert(np.allclose(scf.scf_conv_accuracy,[6.3e-9]))
    assert(scf.volume==380.621)
    assert(np.isclose(scf.cputime,1+2/60+3.5/3600))
    assert(np.isclose(scf.walltime,4/60+33.69/3600))
    assert(scf.pressure==1200.0)
    assert(np.allclose(scf.stress,[[[100.,200.,300.],
                                    [400.,500.,600.],
                                    [700.,800.,900.]]]))
    assert(scf.forces.shape==(1,2,3))
    assert(np.allclose(scf.tot_forces,[1.25e-4]))
    assert(np.allclose(scf.max_forces,[np.linalg.norm([.4,.5,.6])]))
    assert(scf.kpoints_cart.shape==(1,3))
    assert(scf.kpoints_unit.shape==(1,3))
    assert(scf.kweights.shape==(1,))

    relax_file = tmp_path/'relax.out'
    relax_file.write_text('''\
BFGS Geometry Optimization
CELL_PARAMETERS (alat = 2.0D+00) trailing tokens
1.0 0.0 0.0 trailing
0.0 1.0 0.0 trailing
0.0 0.0 1.0 trailing
ATOMIC_POSITIONS (crystal)
H .25 .25 .25 0 0 0
H .75 .75 .75 1 1 1
End final coordinates
''')
    relax = PwscfOutData(relax_file)
    assert(isinstance(relax.relax_structures,list))
    structure = relax.relax_structures[0]

    assert(np.allclose(structure.axes,2*np.eye(3)))
    assert(np.allclose(structure.positions,[[.5,.5,.5],[1.5,1.5,1.5]]))

    malformed_file = tmp_path/'malformed.out'
    malformed_file.write_text('''\
Self-consistent Calculation
! total energy = -168.1 eV
total stress (Ry/bohr**3) (kbar) p = -170.96
stress: -.001 0.0 .001 -147.1 0.0 147.1
-.001 0.0 .001 -147.1 0.0
-.001 0.0.001 -147.1 0.0 147.1
''')
    malformed = PwscfOutData(malformed_file)

    assert(malformed.E is None)
    assert(malformed.pressure is None)
    assert(malformed.stress is None)
#end def test_tokenized_log_parsing


def test_pw2casino_analyzer_read(tmp_path):
    from ..pwscf_analyzer import Pw2CasinoAnalyzer, PwscfAnalyzer, PwscfOutData

    (tmp_path/'pwscf.in').write_text("&CONTROL\n  calculation = 'scf'\n/\n")
    (tmp_path/'pwscf.out').write_text('Self-consistent Calculation\n')
    (tmp_path/'pw2casino.out').write_text('Kinetic energy from orbitals = 1.25D+01\n')

    pw2casino = Pw2CasinoAnalyzer(tmp_path/'pw2casino.out')
    assert(pw2casino.K==12.5)
    analyzer = PwscfAnalyzer(
        tmp_path,
        'pwscf.in',
        'pwscf.out',
        'pw2casino.out',
        analyze = True,
        )
    assert(isinstance(analyzer.results_out,PwscfOutData))
    assert(isinstance(analyzer.pw2casino,Pw2CasinoAnalyzer))
    assert(analyzer.pw2casino.K==12.5)

    (tmp_path/'pw2casino.out').write_text('Kinetic energy is unavailable\n')
    assert(Pw2CasinoAnalyzer(tmp_path/'pw2casino.out').K is None)
    with pytest.raises(FileNotFoundError):
        Pw2CasinoAnalyzer(tmp_path/'missing.out')
    (tmp_path/'pw2casino.out').unlink()
    with pytest.raises(FileNotFoundError):
        PwscfAnalyzer(tmp_path,'pwscf.in','pwscf.out','pw2casino.out',analyze=True)
#end def test_pw2casino_analyzer_read


@pytest.mark.parametrize(
    argnames='verbosity,calculation',
    argvalues=tuple(
        (verbosity,calculation)
        for verbosity in ('default','high','low')
        for calculation in ('scf','nscf','relax','vc-relax')
        ),
    )
def test_qe_7_0_calculation_modes(verbosity,calculation):
    from ..pwscf_analyzer import PwscfAnalyzer, PwscfOutData

    fixture_path = (
        TEST_DIR/'test_pwscf_analyzer_files'/'qe_7_0'/verbosity/calculation
        )
    analyzer = PwscfAnalyzer(
        fixture_path,
        'pwscf.in',
        'pwscf.out',
        analyze = True,
        )
    out = analyzer.results_out

    kpoint_counts = {
        'scf'      : 3,
        'nscf'     : 4,
        'relax'    : 6,
        'vc-relax' : 6,
        }
    band_counts = {
        'scf'      : 4,
        'nscf'     : 8,
        'relax'    : 4,
        'vc-relax' : 4,
        }
    scf_cycle_counts = {
        'scf'      : 6,
        'relax'    : 13,
        'vc-relax' : 23,
        }
    npoints = kpoint_counts[calculation]
    nbands  = band_counts[calculation]

    assert(analyzer.input.control.calculation==calculation)
    assert(isinstance(out,PwscfOutData))
    assert(out.calculation==calculation)
    assert(out.bands is not None)
    assert(out.volume is not None)
    assert(out.cputime is not None)
    assert(out.walltime is not None)
    assert(analyzer.results_xml is None)
    if verbosity=='high':
        assert(out.kpoints_cart.shape==(npoints,3))
        assert(out.kpoints_unit.shape==(npoints,3))
        assert(out.kweights.shape==(npoints,))
        assert(analyzer.kpoints('B').shape==(npoints,3))
        assert(analyzer.kweights().shape==(npoints,))
    else:
        assert(out.kpoints_cart is None)
        assert(out.kpoints_unit is None)
        assert(out.kweights is None)
    #end if
    assert(analyzer.eigenvalues('Ha').shape==(npoints,nbands))
    band = out.bands.up[0]
    assert({'index','kpoint_2pi_alat','kpoint_rel','eigs','occs','pol'}<=set(band))
    if verbosity=='high':
        assert(analyzer.occupations().shape==(npoints,nbands))
    else:
        assert(analyzer.occupations() is None)
    #end if
    if calculation=='nscf':
        assert('E' not in out)
        assert('forces' not in out)
        assert('stress' not in out)
        if verbosity=='high':
            assert(out.bands.electronic_structure=='insulating')
            assert(out.bands.direct_gap.energy>0)
            assert(out.bands.indirect_gap.energy>0)
            assert(out.bands.vbm.index!=out.bands.cbm.index)
        with pytest.raises(
            RuntimeError,
            match=r'not supported for calculation "nscf"',
            ):
            analyzer.forces()
    else:
        nsteps = 1 if calculation=='scf' else 3
        assert(out.E is not None)
        assert(out.relax_energies.shape==(nsteps,))
        assert(out.scf_conv_energy.shape==(scf_cycle_counts[calculation],))
        assert(out.scf_conv_accuracy.shape==(scf_cycle_counts[calculation],))
        assert(out.forces.shape==(nsteps,2,3))
        assert(out.tot_forces.shape==(nsteps,))
        assert(out.max_forces.shape==(nsteps,))
        assert(out.pressure is not None)
        assert(out.stress is not None)
        assert(analyzer.energy('Ry')==out.E)
        assert(analyzer.forces('Ry/B').shape==(nsteps,2,3))
        assert(analyzer.stress('kbar').shape==(nsteps,3,3))
        assert(analyzer.pressure('kbar') is not None)
    #end if
    if calculation in {'relax','vc-relax'}:
        assert(len(out.relax_structures)==3)
        assert(out.relax_structures[2].positions.shape==(2,3))
        if calculation=='vc-relax':
            assert(analyzer.relaxed_structure('B').pos.shape==(2,3))
    else:
        assert('relax_structures' not in out)
    #end if
#end def test_qe_7_0_calculation_modes


@pytest.mark.parametrize(
    argnames='calculation',
    argvalues=('bands','md','vc-md'),
    )
def test_unsupported_calculation_modes(calculation):
    from ..pwscf_analyzer import PwscfAnalyzer

    fixture_path = (
        TEST_DIR/'test_pwscf_analyzer_files'/'qe_7_0'/'default'/calculation
        )
    with pytest.raises(
        RuntimeError,
        match=r'PWSCF .* calculations are not supported',
        ):
        PwscfAnalyzer(
            fixture_path,
            'pwscf.in',
            'pwscf.out',
            analyze = True,
            )
#end def test_unsupported_calculation_modes


def test_quantity_accessors():
    import numpy as np

    from ..pwscf_analyzer import PwscfAnalyzer

    fixture_root = TEST_DIR/'test_pwscf_analyzer_files'/'qe_7_0'/'high'
    scf = PwscfAnalyzer(
        fixture_root/'scf',
        'pwscf.in',
        'pwscf.out',
        analyze = True,
        )
    assert(np.isclose(scf.energy('Ha'),scf.energy('Ry')/2))
    assert(scf.kpoints('B').shape==(3,3))
    assert(scf.kweights().shape==(3,))
    assert(scf.eigenvalues('eV').shape==(3,4))
    assert(scf.occupations().shape==(3,4))
    assert(scf.forces('Ry/B').shape==(1,2,3))
    assert(scf.stress('kbar').shape==(1,3,3))
    assert(scf.pressure('kbar') is not None)

    input_fixture = TEST_DIR/'test_pwscf_analyzer_files'/'scf_output'
    input_scf = PwscfAnalyzer(
        input_fixture,
        'scf.in',
        'scf.out',
        analyze = True,
        )
    assert(input_scf.initial_structure('A').units=='A')

    nscf = PwscfAnalyzer(
        fixture_root/'nscf',
        'pwscf.in',
        'pwscf.out',
        analyze = True,
        )
    assert(nscf.energy() is None)
    assert(nscf.Ef() is None)
    assert(nscf.Evbm() is not None)
    assert(nscf.Ecbm() is not None)
    assert(nscf.band_gap() is not None)
    assert(not nscf.fractional_occs())
    with pytest.raises(
        RuntimeError,
        match=r'not supported for calculation "nscf"',
        ):
        nscf.forces()
    with pytest.raises(
        RuntimeError,
        match=r'not supported for calculation "nscf"',
        ):
        nscf.relaxed_structure()
    with pytest.raises(
        ValueError,
        match=r'energy units must be one of',
        ):
        nscf.energy('invalid')
    with pytest.raises(
        ValueError,
        match=r'kpoints units must be one of',
        ):
        nscf.kpoints('invalid')
    with pytest.raises(
        ValueError,
        match=r'eigenvalues units must be one of',
        ):
        nscf.eigenvalues('invalid')

    fermi_fixture = TEST_DIR/'test_pwscf_analyzer_files'/'nscf_output'
    fermi_nscf = PwscfAnalyzer(
        fermi_fixture,
        'nscf.in',
        'nscf.out',
        analyze = True,
        )
    assert(fermi_nscf.results_out.fermi_energies.shape==(1,))
#end def test_quantity_accessors


def test_legacy_xml(tmp_path,monkeypatch):
    from .. import pwscf_analyzer as pa_module
    from ..developer import obj
    from ..pwscf_analyzer import PwscfAnalyzer

    savedir = tmp_path/'pwscf.save'
    savedir.mkdir()
    datafile = savedir/'data-file.xml'
    eigfile  = savedir/'eigenval.xml'
    datafile.write_text('<legacy/>')
    eigfile.write_text('<eigenvalues/>')

    data = obj(root=obj(eigenvalues=obj(k_point=obj({
        1:obj(
            k_point_coords = [0.0,0.0,0.0],
            weight         = 1.0,
            datafile       = obj({1:obj(iotk_link='eigenval.xml')}),
            ),
        }))))
    eigenvalues = obj(root=obj(
        units_for_energies = obj(units='Hartree'),
        eigenvalues        = [-0.5,0.5],
        occupations        = [1.0,0.0],
        ))

    def read_qexml(filepath):
        return data if str(filepath).endswith('data-file.xml') else eigenvalues
    #end def read_qexml

    monkeypatch.setattr(pa_module,'read_qexml',read_qexml)
    analyzer = PwscfAnalyzer()
    analyzer.path        = str(tmp_path)
    analyzer.input       = obj(control=obj(outdir='.',prefix='pwscf'))
    analyzer.results_xml = None
    analyzer.analyze_xml()

    assert(analyzer.results_xml is not None)
    assert(analyzer.results_xml.data is data)
    assert(analyzer.results_xml.kpoints[1].weight==1.0)
    assert(analyzer.results_xml.kpoints[1].up.units=='Ha')
    assert(analyzer.results_xml.kpoints[1].up.eigenvalues==[-0.5,0.5])
#end def test_legacy_xml
