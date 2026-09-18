import pytest

from . import TEST_DIR, NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.PWSCF_ANALYZER)

PWSCF_ANALYZER_FIXTURES = TEST_DIR/'test_pwscf_analyzer_files'
SCHEMA_FIXTURES = tuple(sorted(
    PWSCF_ANALYZER_FIXTURES.rglob('data-file-schema.xml'),
    ))
TRIPLE_SOURCE_CASES = (
    (
        'scf',
        'qe_7_0/high/scf',
        (
            'initial_structure','energy','kpoints','kweights','eigenvalues',
            'occupations','fractional_occs','forces','stress','pressure',
            ),
        ),
    (
        'nscf',
        'qe_7_0/high/nscf',
        (
            'initial_structure','kpoints','kweights','eigenvalues',
            'occupations','Evbm','Ecbm','band_gap','fractional_occs',
            ),
        ),
    (
        'relax',
        'qe_7_0/high/relax',
        (
            'initial_structure','energy','kpoints','kweights','eigenvalues',
            'occupations','fractional_occs','relaxed_structure','forces',
            'stress','pressure',
            ),
        ),
    (
        'vc-relax',
        'qe_7_0/high/vc-relax',
        (
            'initial_structure','energy','kweights','eigenvalues',
            'occupations','fractional_occs','relaxed_structure','forces',
            'stress','pressure',
            ),
        ),
    (
        'spin',
        'qe_7_6/supplemental/cbn_spin',
        (
            'initial_structure','energy','kpoints','kweights','eigenvalues',
            'occupations','Ef','Evbm','Ecbm','band_gap','fractional_occs',
            'forces','stress','pressure',
            ),
        ),
    )


def write_input(directory,calculation='scf',prefix='pwscf',outdir='.'):
    """Write a minimal PWSCF input used by analyzer policy tests."""
    infile = directory/'pwscf.in'
    infile.write_text(
        "&CONTROL\n"
        f"  calculation = '{calculation}'\n"
        f"  prefix = '{prefix}'\n"
        f"  outdir = '{outdir}'\n"
        "/\n"
        )
    return infile
#end def write_input


def write_schema(directory,calculation='scf',energy=-1.0,prefix='pwscf'):
    """Write a minimal modern QE schema file and return its path."""
    savedir = directory/f'{prefix}.save'
    savedir.mkdir(parents=True,exist_ok=True)
    schema_file = savedir/'data-file-schema.xml'
    schema_file.write_text(
        '<espresso>\n'
        '  <input><control_variables>'
        f'<calculation>{calculation}</calculation>'
        '</control_variables></input>\n'
        '  <output><total_energy>'
        f'<etot>{energy}</etot>'
        '</total_energy></output>\n'
        '</espresso>\n'
        )
    return schema_file
#end def write_schema


def test_analyzer_input_object():
    from ..pwscf_analyzer import PwscfAnalyzer
    from ..pwscf_input import PwscfInput

    fixture = PWSCF_ANALYZER_FIXTURES/'qe_7_0/default/scf'
    input_data = PwscfInput(fixture/'pwscf.in')
    analyzer = PwscfAnalyzer(
        input   = input_data,
        outfile = 'pwscf.out',
        path    = fixture,
        analyze = True,
        )

    assert analyzer.input is input_data
    assert analyzer.calculation=='scf'
#end def test_analyzer_input_object


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
    from ..pwscf_analyzer import (
        Pw2CasinoAnalyzer,
        PwscfAnalyzer,
        PwscfOutData,
        PwscfXmlData,
    )

    pa = PwscfAnalyzer()
    assert(pa.read_all)
    assert(pa.strict)
    assert(pa.required==set())
    with pytest.raises(
        FileNotFoundError,
        match=r'PWSCF schema XML file and output file are not available',
        ):
        pa.analyze()
    with pytest.raises(
        RuntimeError,
        match=r'PWSCF output has not been analyzed',
        ):
        pa.make_movie('movie.xyz')
    free_helpers = ('parse_float',)
    for name in free_helpers:
        assert(callable(getattr(pa_module,name)))
        assert(not hasattr(PwscfAnalyzer,name))
    #end for
    assert(not hasattr(pa_module,'read_kpoint_tables'))
    reader_names = (
        'read_calculation','read_fermi_energies',
        'read_energies','read_scf_convergence','read_bands','read_initial_structure',
        'read_structures','read_pressure','read_volume','read_stress',
        'read_forces','read_timing','read_kpoints',
        )
    assert(all(callable(getattr(PwscfOutData,name)) for name in reader_names))
    assert(not hasattr(PwscfOutData,'read'))
    assert(not hasattr(PwscfOutData,'read_band_edges'))
    assert(not any(hasattr(PwscfAnalyzer,'analyze_'+name[5:]) for name in reader_names))
    assert(not hasattr(PwscfAnalyzer,'analyze_schema_xml'))
    assert(callable(Pw2CasinoAnalyzer))
    assert(not hasattr(Pw2CasinoAnalyzer,'read'))
    assert(pa_module.PwscfXmlData is PwscfXmlData)
#end def test_empty_init


@pytest.mark.parametrize(
    argnames='schema_file',
    argvalues=SCHEMA_FIXTURES,
    ids=lambda path:str(path.relative_to(PWSCF_ANALYZER_FIXTURES)),
    )
def test_all_schema_fixtures_are_parsed(schema_file):
    from ..pwscf_analyzer import PwscfAnalyzer, PwscfXmlData

    fixture_path = schema_file.parent
    while (
        fixture_path!=PWSCF_ANALYZER_FIXTURES
        and not (fixture_path/'pwscf.in').is_file()
        ):
        fixture_path = fixture_path.parent
    assert((fixture_path/'pwscf.in').is_file())

    analyzer = PwscfAnalyzer(
        input = 'pwscf.in',
        path  = fixture_path,
        analyze = True,
        read_all = False,
        )
    selected,status = analyzer._schema_file()
    assert(status=='found')
    assert(schema_file.samefile(selected))
    assert(isinstance(analyzer.results_xml,PwscfXmlData))
    assert(analyzer.source_status.xml=='parsed')
#end def test_all_schema_fixtures_are_parsed


@pytest.mark.parametrize(
    argnames='case,relative_path,quantities',
    argvalues=TRIPLE_SOURCE_CASES,
    ids=[case[0] for case in TRIPLE_SOURCE_CASES],
    )
def test_query_agreement_across_sources(case,relative_path,quantities):
    import numpy as np

    from ..pwscf_analyzer import PwscfAnalyzer, PwscfXmlData
    from ..structure import Structure

    fixture_path = PWSCF_ANALYZER_FIXTURES/relative_path
    analyzers = {
        mode:PwscfAnalyzer(
            input   = 'pwscf.in',
            outfile = 'pwscf.out',
            path    = fixture_path,
            analyze = True,
            read_all = mode=='all',
            )
        for mode in ('all','xml')
        }
    both = analyzers['all']
    xml  = analyzers['xml']

    assert(isinstance(both.results_xml,PwscfXmlData))
    assert(both.results_out is not None)
    assert(isinstance(xml.results_xml,PwscfXmlData))
    assert(xml.results_out is None)
    assert(both.calculation==xml.calculation)

    tolerances = {
        'initial_structure' : 1e-7,
        'relaxed_structure' : 1e-7,
        'energy'            : 1e-7,
        'kpoints'           : 1e-7,
        'kweights'          : 1e-7,
        'eigenvalues'       : 1e-4,
        'occupations'       : 1e-10,
        'Ef'                : 1e-4,
        'Evbm'              : 1e-4,
        'Ecbm'              : 1e-4,
        'band_gap'          : 1e-4,
        'forces'            : 1e-5,
        'stress'            : 1e-3,
        'pressure'          : 1e-3,
        }

    for quantity in quantities:
        assert(all(
            analyzer.available(quantity)
            for analyzer in analyzers.values()
            )), f'{case}: {quantity} is not available from every source'
        values = {
            source:getattr(analyzer,quantity)()
            for source,analyzer in analyzers.items()
            }
        both_value = values['all']
        xml_value  = values['xml']
        if isinstance(xml_value,Structure):
            assert(isinstance(both_value,Structure))
            assert(both_value.units==xml_value.units)
            assert(np.array_equal(both_value.elem,xml_value.elem))
            assert(np.array_equal(both_value.axes,xml_value.axes))
            assert(np.array_equal(both_value.pos,xml_value.pos))
        elif isinstance(xml_value,(bool,np.bool_)):
            assert(both_value==xml_value)
        else:
            both_array = np.asarray(both_value)
            xml_array  = np.asarray(xml_value)
            assert(both_array.shape==xml_array.shape)
            assert(both_array.dtype.kind==xml_array.dtype.kind)
            assert(np.array_equal(both_array,xml_array))
        #end if
    #end for
#end def test_query_agreement_across_sources


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
        'calculation','run_type_detected',
            'Ef','fermi_energies','bands',
            'volume','cputime','walltime',
            'kpoints_cart','kpoints_unit','kweights','initial_structure_data',
            'E','relax_energies','scf_conv_energy','scf_conv_accuracy',
            'pressure','stress','forces','tot_forces','max_forces',
            'md_data','md_stats','relax_structures',
        }

    assert(set(out.keys())==expected)
    assert(out.calculation==calculation)
    assert(out.run_type_detected)
    assert(all(
        value is None
        for name,value in out.items()
        if name not in {'calculation','run_type_detected'}
        ))
#end def test_result_initialization


def test_tokenized_log_parsing(tmp_path):
    import numpy as np

    from ..pwscf_analyzer import PwscfOutData

    scf_file = tmp_path/'scf.out'
    scf_file.write_text('''\
Self-consistent Calculation
number of atoms/cell   = 2 trailing tokens
number of k points = 1 trailing tokens

unrelated informational line
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

    angstrom_file = tmp_path/'angstrom.out'
    angstrom_file.write_text('''\
BFGS Geometry Optimization
CELL_PARAMETERS (angstrom)
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
ATOMIC_POSITIONS (angstrom)
H 1.0 2.0 3.0
End final coordinates
''')
    angstrom = PwscfOutData(angstrom_file)
    structure = angstrom.relax_structures[0]
    bohr_per_angstrom = 1.0/0.529177210903
    assert(np.allclose(structure.axes,bohr_per_angstrom*np.eye(3)))
    assert(np.allclose(
        structure.positions,
        np.array([1.0,2.0,3.0])*bohr_per_angstrom,
        ))
#end def test_tokenized_log_parsing


def test_malformed_log_records_are_skipped(tmp_path):
    import numpy as np

    from ..pwscf_analyzer import PwscfOutData

    md_file = tmp_path/'md.out'
    md_file.write_text('''\
! total energy =
total stress (Ry/bohr**3) (kbar) P=
time =
Ekin T
temperature =
! total energy = -1.0 Ry
total stress (Ry/bohr**3) (kbar) P= 2.0
time = 0.5
kinetic energy = 0.25 Ry
temperature = 300 K
''')
    md = PwscfOutData(md_file,'md',md_only=True)
    assert(md.md_data is not None)
    assert(np.allclose(md.md_data.total_energy,[-1.0]))
    assert(np.allclose(md.md_data.pressure,[2.0]))
    assert(np.allclose(md.md_data.time,[0.5]))

    forces_file = tmp_path/'forces.out'
    forces_file.write_text('''\
Self-consistent Calculation
Forces acting on atoms
atom 1 type 1 force = 0.1 0.2 0.3
Forces acting on atoms
atom 1 type 1 force = 0.4 0.5 0.6
    atom 2 type 1 force = 0.7 0.8 0.9
''')
    forces = PwscfOutData(forces_file)
    assert(forces.forces.shape==(1,2,3))
    assert(np.allclose(forces.forces[0,0],[0.4,0.5,0.6]))

    timing_file = tmp_path/'timing.out'
    timing_file.write_text(
        'Self-consistent Calculation\n'
        'PWSCF : 1e999s CPU 2s WALL\n'
        )
    timing = PwscfOutData(timing_file)
    assert(timing.cputime is None)
    assert(np.isclose(timing.walltime,2/3600))

    bands_file = tmp_path/'bands.out'
    bands_file.write_text('''\
Band Structure Calculation
highest occupied level (ev): 1.0
bands (ev):  k = 0.0 0.0 0.0 ( 1 PWs) bands (ev):
1e999 2.0
''')
    assert(PwscfOutData(bands_file).bands is None)

    structure_file = tmp_path/'structure.out'
    structure_file.write_text('''\
BFGS Geometry Optimization
number of atoms/cell = 2
ATOMIC_POSITIONS (bohr)
H 0.0 0.0 0.0
End final coordinates
''')
    assert(PwscfOutData(structure_file).relax_structures is None)
#end def test_malformed_log_records_are_skipped


def test_calculation_sources_are_independent(tmp_path):
    from ..pwscf_analyzer import PwscfAnalyzer

    (tmp_path/'pwscf.in').write_text('''\
&CONTROL
  calculation = 'relax'
/
''')
    (tmp_path/'pwscf.out').write_text('''\
BFGS Geometry Optimization
CELL_PARAMETERS (alat= 5.0)
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
''')
    analyzer = PwscfAnalyzer(
        input   = 'pwscf.in',
        outfile = 'pwscf.out',
        path    = tmp_path,
        analyze = True,
        strict = False,
        )
    assert(analyzer.results_out.calculation=='vc-relax')
    assert(analyzer.calculation==('relax',None,'vc-relax'))

    (tmp_path/'pwscf.out').write_text('Self-consistent Calculation\n')
    analyzer = PwscfAnalyzer(
        input   = 'pwscf.in',
        outfile = 'pwscf.out',
        path    = tmp_path,
        analyze = True,
        )
    assert(analyzer.results_out.calculation=='scf')
    assert(analyzer.calculation==('relax',None,'scf'))
#end def test_calculation_sources_are_independent


def test_pw2casino_analyzer_read(tmp_path):
    from ..pwscf_analyzer import Pw2CasinoAnalyzer, PwscfAnalyzer, PwscfOutData

    (tmp_path/'pwscf.in').write_text("&CONTROL\n  calculation = 'scf'\n/\n")
    (tmp_path/'pwscf.out').write_text('Self-consistent Calculation\n')
    (tmp_path/'pw2casino.out').write_text('Kinetic energy from orbitals = 1.25D+01\n')

    pw2casino = Pw2CasinoAnalyzer(tmp_path/'pw2casino.out')
    assert(pw2casino.K==12.5)
    analyzer = PwscfAnalyzer(
        input   = 'pwscf.in',
        outfile = 'pwscf.out',
        path    = tmp_path,
        pw2c_outfile_name = 'pw2casino.out',
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
        PwscfAnalyzer(
            input   = 'pwscf.in',
            outfile = 'pwscf.out',
            path    = tmp_path,
            pw2c_outfile_name = 'pw2casino.out',
            analyze = True,
            )
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
    import numpy as np

    from ..pwscf_analyzer import PwscfAnalyzer, PwscfOutData

    fixture_path = (
        TEST_DIR/'test_pwscf_analyzer_files'/'qe_7_0'/verbosity/calculation
        )
    analyzer = PwscfAnalyzer(
        input   = 'pwscf.in',
        outfile = 'pwscf.out',
        path    = fixture_path,
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
    assert(analyzer.results_xml is not None)
    if verbosity=='high':
        assert(out.kpoints_cart.shape==(npoints,3))
        assert(out.kpoints_unit.shape==(npoints,3))
        assert(out.kweights.shape==(npoints,))
    else:
        assert(out.kpoints_cart is None)
        assert(out.kpoints_unit is None)
        assert(out.kweights is None)
    #end if
    assert(analyzer.kpoints('B').shape==(npoints,3))
    assert(analyzer.kweights().shape==(npoints,))
    assert(analyzer.eigenvalues('Ha').shape==(npoints,nbands))
    band = out.bands.up[0]
    assert({'index','kpoint_2pi_alat','kpoint_rel','eigs','occs','pol'}<=set(band))
    assert(analyzer.occupations().shape==(npoints,nbands))
    if calculation=='nscf':
        assert('E' in out)
        assert('forces' in out)
        assert('stress' in out)
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
        assert(np.isclose(analyzer.energy('Ry'),out.E))
        assert(type(analyzer.energy()) is float)
        assert(analyzer.forces('Ry/B').shape==(2,3))
        assert(analyzer.stress('kbar').shape==(3,3))
        assert(analyzer.pressure('kbar') is not None)
        assert(type(analyzer.pressure()) is float)
    #end if
    if calculation in {'relax','vc-relax'}:
        assert(len(out.relax_structures)==3)
        assert(out.relax_structures[2].positions.shape==(2,3))
        if calculation=='vc-relax':
            assert(analyzer.relaxed_structure('B').pos.shape==(2,3))
    else:
        assert(out.relax_structures is None)
    #end if
#end def test_qe_7_0_calculation_modes


@pytest.mark.parametrize(
    argnames='verbosity,calculation',
    argvalues=tuple(
        (verbosity,calculation)
        for verbosity in ('default','high','low')
        for calculation in ('bands','md','vc-md')
        ),
    )
def test_supported_calculation_modes(verbosity,calculation):
    from ..pwscf_analyzer import PwscfAnalyzer, PwscfXmlData

    fixture_path = (
        TEST_DIR/'test_pwscf_analyzer_files'/'qe_7_0'/verbosity/calculation
        )
    analyzer = PwscfAnalyzer(
        input   = 'pwscf.in',
        outfile = 'pwscf.out',
        path    = fixture_path,
        analyze = True,
        )
    assert(analyzer.results_out.calculation==calculation)
    assert(isinstance(analyzer.results_xml,PwscfXmlData))
#end def test_supported_calculation_modes


@pytest.mark.parametrize(
    argnames='version,case,calculation,eigenvalue_shape,has_schema',
    argvalues=tuple(
        (version,case,calculation,eigenvalue_shape,has_schema)
        for version in ('qe_7_0','qe_7_6')
        for case,calculation,eigenvalue_shape,has_schema in (
            ('cbn_crystal_kpoints','scf',(4,4),True),
            ('cbn_relax','relax',(6,4),True),
            ('cbn_scf','scf',(3,4),True),
            ('cbn_smearing','scf',(4,8),True),
            ('cbn_spin','scf',(4,2,8),True),
            ('cbn_vc_relax','vc-relax',(6,4),True),
            ('scf_crystal_kpoints','scf',(4,4),True),
            ('scf_no_symmetry','scf',(27,4),True),
            ('scf_no_xml','scf',(3,4),False),
            ('scf_smearing','scf',(4,8),True),
            ('scf_spin','scf',(4,2,8),True),
            ('scf_symmetry','scf',(4,4),True),
            )
        ),
    )
def test_supplemental_qe_runs(
    version,
    case,
    calculation,
    eigenvalue_shape,
    has_schema,
    ):
    import numpy as np

    from ..pwscf_analyzer import PwscfAnalyzer, PwscfXmlData

    fixture_path = (
        TEST_DIR/'test_pwscf_analyzer_files'/version/'supplemental'/case
        )
    schema_file = fixture_path/'pwscf.save'/'data-file-schema.xml'
    assert(schema_file.is_file()==has_schema)
    analyzer = PwscfAnalyzer(
        input   = 'pwscf.in',
        outfile = 'pwscf.out',
        path    = fixture_path,
        analyze = True,
        )
    out = analyzer.results_out

    assert(out.calculation==calculation)
    assert(isinstance(analyzer.results_xml,PwscfXmlData)==has_schema)
    assert(analyzer.source_status.xml==('parsed' if has_schema else 'missing'))
    assert(analyzer.eigenvalues().shape==eigenvalue_shape)
    assert(analyzer.occupations().shape==eigenvalue_shape)
    assert(out.forces.shape[1:]==(2,3))
    assert(out.stress.shape[1:]==(3,3))
    if has_schema or case in {'cbn_smearing','cbn_spin','scf_smearing','scf_spin'}:
        assert(analyzer.Ef() is not None)
    else:
        assert(analyzer.Ef() is None)
    if 'spin' in case:
        assert(out.bands.up[0].pol=='up')
        assert(out.bands.down[0].pol=='down')
    else:
        assert(out.bands.up[0].pol=='none')
        assert(len(out.bands.down)==0)
    if 'crystal_kpoints' in case:
        assert(np.isclose(analyzer.kweights().sum(),2.0))
    if case in {'cbn_relax','cbn_vc_relax'}:
        assert(out.relax_structures[-1].atoms==['B','N'])
        assert(analyzer.relaxed_structure('B').pos.shape==(2,3))
    if case.startswith('cbn_'):
        initial = analyzer.initial_structure('B')
        assert(initial.elem.tolist()==['B','N'])
        assert(initial.axes.shape==(3,3))
    #end if
#end def test_supplemental_qe_runs


@pytest.mark.parametrize(
    argnames='version,case',
    argvalues=tuple(
        (version,case)
        for version in ('qe_7_0','qe_7_6')
        for case in ('md_iprint','vc_md_iprint')
        ),
    )
def test_supplemental_md_runs(version,case):
    from ..pwscf_analyzer import PwscfAnalyzer, PwscfXmlData

    fixture_path = (
        TEST_DIR/'test_pwscf_analyzer_files'/version/'supplemental'/case
        )
    analyzer = PwscfAnalyzer(
        input   = 'pwscf.in',
        outfile = 'pwscf.out',
        path    = fixture_path,
        analyze = True,
        )
    assert(isinstance(analyzer.results_xml,PwscfXmlData))
    assert(analyzer.results_out.calculation in {'md','vc-md'})
    assert(analyzer.results_out.md_data is not None)
    assert(len(analyzer.results_out.md_data.time)>0)
    stats = analyzer.md_statistics()
    assert(stats is not None)
    assert(set(stats)==set(analyzer.results_out.md_data))
    md_only = PwscfAnalyzer(
        input   = 'pwscf.in',
        outfile = 'pwscf.out',
        path    = fixture_path,
        analyze = True,
        md_only = True,
        )
    assert(isinstance(md_only.results_xml,PwscfXmlData))
    assert(md_only.results_out.md_data is not None)
    assert(md_only.results_out.E is None)
#end def test_supplemental_md_runs


def test_quantity_accessors():
    import numpy as np

    from ..pwscf_analyzer import PwscfAnalyzer
    from ..unit_converter import UnitConverter

    fixture_root = TEST_DIR/'test_pwscf_analyzer_files'/'qe_7_0'/'high'
    scf = PwscfAnalyzer(
        input   = 'pwscf.in',
        outfile = 'pwscf.out',
        path    = fixture_root/'scf',
        analyze = True,
        )
    assert(np.isclose(scf.energy('Ha'),scf.energy('Ry')/2))
    assert(scf.kpoints('B').shape==(3,3))
    assert(scf.kweights().shape==(3,))
    assert(scf.eigenvalues('eV').shape==(3,4))
    assert(scf.occupations().shape==(3,4))
    assert('nspin' not in scf.input.system)
    assert(np.all(scf.occupations()==1.0))
    assert(not scf.fractional_occs())
    scf.results_xml.occupations[0,0] = .9995
    assert(not scf.fractional_occs())
    assert(scf.fractional_occs(tol=1e-4))
    assert(scf.forces('Ry/B').shape==(2,3))
    assert(scf.stress('kbar').shape==(3,3))
    assert(scf.pressure('kbar') is not None)
    unit_scales = {
        'eV/A^3'    : UnitConverter.A**3/UnitConverter.eV,
        'Ha/Bohr^3' : UnitConverter.B**3/UnitConverter.Ha,
        'Ry/Bohr^3' : UnitConverter.B**3/UnitConverter.Ry,
        }
    for units,scale in unit_scales.items():
        assert(np.allclose(scf.stress(units),scf.stress('kbar')*1e8*scale))
        assert(np.isclose(scf.pressure(units),scf.pressure('kbar')*1e8*scale))
    #end for

    input_fixture = TEST_DIR/'test_pwscf_analyzer_files'/'scf_output'
    input_scf = PwscfAnalyzer(
        input   = 'scf.in',
        outfile = 'scf.out',
        path    = input_fixture,
        analyze = True,
        )
    assert(input_scf.initial_structure('A').units=='A')

    nscf = PwscfAnalyzer(
        input   = 'pwscf.in',
        outfile = 'pwscf.out',
        path    = fixture_root/'nscf',
        analyze = True,
        )
    assert(nscf.energy() is not None)
    assert(nscf.Ef() is not None)
    assert(nscf.Evbm() is not None)
    assert(nscf.Ecbm() is not None)
    assert(nscf.band_gap() is not None)
    for quantity in ('energy','Ef','Evbm','Ecbm','band_gap'):
        assert(type(getattr(nscf,quantity)()) is float)
    assert(not nscf.fractional_occs())
    with pytest.raises(
        RuntimeError,
        match=r'not supported for calculation "nscf"',
        ):
        nscf.forces()
    assert(nscf.relaxed_structure() is not None)
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
        input   = 'nscf.in',
        outfile = 'nscf.out',
        path    = fermi_fixture,
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

    browse_only = PwscfAnalyzer(
        path = tmp_path,
        analyze = True,
        read_all = False,
        strict = False,
        )
    assert(browse_only.results_xml.data is data)
    assert(browse_only.energy() is None)
    assert(browse_only.calculation is None)

    def broken_read_qexml(filepath):
        raise RuntimeError('malformed legacy XML')
    #end def broken_read_qexml

    analyzer.results_out = obj(marker='preserved')
    monkeypatch.setattr(pa_module,'read_qexml',broken_read_qexml)
    analyzer.analyze_xml()
    assert(analyzer.results_xml is None)
    assert(analyzer.results_out.marker=='preserved')
#end def test_legacy_xml


def test_schema_queries_are_preferred(tmp_path):
    """Modern XML supplies query data even when the text output is sparse."""
    import numpy as np

    from ..pwscf_analyzer import PwscfAnalyzer, PwscfXmlData

    (tmp_path/'pwscf.in').write_text("""&CONTROL
 calculation = 'relax'
 prefix = 'pwscf'
 outdir = '.'
/
""")
    (tmp_path/'pwscf.out').write_text('Self-consistent Calculation\n')
    savedir = tmp_path/'pwscf.save'
    savedir.mkdir()
    (savedir/'data-file-schema.xml').write_text("""<espresso>
  <input><control_variables><calculation>relax</calculation></control_variables>
    <atomic_structure alat="2"><atomic_positions><atom name="H">0 0 0</atom></atomic_positions><cell><a1>2 0 0</a1><a2>0 2 0</a2><a3>0 0 2</a3></cell></atomic_structure>
  </input>
  <output><atomic_structure alat="2"><atomic_positions><atom name="H">0.5 0 0</atom></atomic_positions><cell><a1>2 0 0</a1><a2>0 2 0</a2><a3>0 0 2</a3></cell></atomic_structure>
    <total_energy><etot>-1</etot></total_energy><forces>0.1 0.2 0.3</forces><stress>1 0 0 0 2 0 0 0 3</stress>
    <band_structure><lsda>false</lsda><fermi_energy>0</fermi_energy><ks_energies><k_point weight="1">0 0 0</k_point><eigenvalues>-0.5 0.25</eigenvalues><occupations>1 0</occupations></ks_energies></band_structure>
  </output>
</espresso>""")

    analyzer = PwscfAnalyzer(input='pwscf.in',outfile='pwscf.out',path=tmp_path,analyze=True)
    assert(isinstance(analyzer.results_xml,PwscfXmlData))
    assert(np.isclose(analyzer.energy('Ha'),-1.0))
    assert(analyzer.initial_structure('B').pos.shape==(1,3))
    assert(analyzer.relaxed_structure('B').pos[0,0]==0.5)
    assert(analyzer.kpoints('B').shape==(1,3))
    assert(analyzer.kweights().shape==(1,))
    assert(analyzer.eigenvalues('Ha').shape==(1,2))
    assert(analyzer.occupations().shape==(1,2))
    assert(np.isclose(analyzer.Ef('Ha'),0.0))
    assert(np.isclose(analyzer.Evbm('Ha'),-0.5))
    assert(np.isclose(analyzer.Ecbm('Ha'),0.25))
    assert(np.isclose(analyzer.band_gap('Ha'),0.75))
    assert(not analyzer.fractional_occs())
    assert(analyzer.forces('Ha/B').shape==(1,3))
    stress = analyzer.stress('GPa')
    assert(stress.shape==(3,3))
    assert(np.isclose(analyzer.pressure('GPa'),np.trace(stress)/3))

    xml_only = PwscfAnalyzer(
        input = 'pwscf.in',
        path  = tmp_path,
        analyze = True,
        read_all = False,
        )
    assert(xml_only.results_out is None)
    assert(xml_only.results_xml is not None)
    assert(np.isclose(xml_only.energy('Ha'),-1.0))
    with pytest.raises(TypeError,match='read_all must be a bool'):
        PwscfAnalyzer(input='pwscf.in',path=tmp_path,read_all='invalid')

    malformed = tmp_path/'bad-schema.xml'
    malformed.write_text('<espresso>')
    assert(PwscfXmlData(malformed).parse_failed)

    bad_encoding = tmp_path/'bad-encoding-schema.xml'
    bad_encoding.write_text(
        '<?xml version="1.0" encoding="DTF-8"?>\n<espresso/>\n'
        )
    assert(PwscfXmlData(bad_encoding).parse_failed)
#end def test_schema_queries_are_preferred


def test_schema_semantic_errors_are_field_local(tmp_path):
    import numpy as np

    from ..pwscf_analyzer import PwscfXmlData

    schema_file = tmp_path/'data-file-schema.xml'
    schema_file.write_text('''\
<espresso>
  <input>
    <control_variables><calculation>scf</calculation></control_variables>
    <atomic_structure alat="bad">
      <atomic_positions><atom name="H">bad coordinates</atom></atomic_positions>
      <cell><a1>1 0 0</a1><a2>0 1 0</a2><a3>0 0 1</a3></cell>
    </atomic_structure>
  </input>
  <output>
    <atomic_structure alat="2">
      <atomic_positions><atom name="H">0 0 0</atom></atomic_positions>
      <cell><a1>2 0 0</a1><a2>0 2 0</a2><a3>0 0 2</a3></cell>
    </atomic_structure>
    <total_energy><etot>not-a-number</etot></total_energy>
    <forces>not-a-force</forces>
    <stress>1 0 0 0 2 0 0 0 3</stress>
    <band_structure>
      <lsda>false</lsda><fermi_energy>0.1</fermi_energy>
      <ks_energies><k_point weight="0.5">0 0 0</k_point><eigenvalues>-1 1</eigenvalues><occupations>1 0</occupations></ks_energies>
      <ks_energies><k_point weight="0.5">0.5 0 0</k_point><eigenvalues>-0.5 0.5 1</eigenvalues><occupations>1 0 0</occupations></ks_energies>
    </band_structure>
  </output>
</espresso>
''')
    xml = PwscfXmlData(schema_file)
    assert(not xml.parse_failed)
    assert(xml.calculation=='scf')
    assert(xml.initial_atoms is None)
    assert(xml.initial_alat is None)
    assert(xml.positions.shape==(1,3))
    assert(xml.total_energy is None)
    assert(xml.forces is None)
    assert(np.allclose(xml.stress,np.diag([1.,2.,3.])))
    assert(xml.eigenvalues is None)
    assert(xml.occupations is None)
    assert(xml.fermi_energy==0.1)
#end def test_schema_semantic_errors_are_field_local


def test_default_source_requires_one_file(tmp_path):
    from ..pwscf_analyzer import PwscfAnalyzer, PwscfXmlData

    (tmp_path/'pwscf.in').write_text("""&CONTROL
 calculation = 'scf'
 prefix = 'pwscf'
 outdir = '.'
/
""")
    savedir = tmp_path/'pwscf.save'
    savedir.mkdir()
    (savedir/'data-file-schema.xml').write_text('''\
<espresso>
  <input><control_variables><calculation>scf</calculation></control_variables></input>
  <output><total_energy><etot>-1</etot></total_energy></output>
</espresso>
''')

    analyzer = PwscfAnalyzer(
        input = 'pwscf.in',
        path  = tmp_path,
        analyze = True,
        )
    assert(analyzer.results_out is None)
    assert(isinstance(analyzer.results_xml,PwscfXmlData))
    assert(analyzer.energy('Ha')==-1.0)

    with pytest.raises(FileNotFoundError,match='PWSCF output file is not available'):
        PwscfAnalyzer(
            input   = 'pwscf.in',
            outfile = 'missing.out',
            path    = tmp_path,
            analyze = True,
            )
#end def test_default_source_requires_one_file


def test_source_required_and_strict_validation(tmp_path):
    from ..pwscf_analyzer import PwscfAnalyzer

    quantities = {
        'initial_structure','energy','kpoints','kweights','eigenvalues',
        'occupations','Ef','Evbm','Ecbm','band_gap','fractional_occs',
        'relaxed_structure','forces','stress','pressure',
        }
    analyzer = PwscfAnalyzer(
        path = tmp_path,
        analyze = False,
        required = 'energy',
        )
    assert(analyzer.read_all)
    assert(analyzer.strict is True)
    assert(analyzer.required=={'energy'})
    assert(set(analyzer.quantity_names)==quantities)

    analyzer = PwscfAnalyzer(
        path = tmp_path,
        analyze = False,
        required = ['energy','forces','energy'],
        strict = False,
        )
    assert(analyzer.required=={'energy','forces'})
    with pytest.raises(TypeError,match='strict must be a bool'):
        PwscfAnalyzer(path=tmp_path,strict=1)
    with pytest.raises(TypeError,match='read_all must be a bool'):
        PwscfAnalyzer(path=tmp_path,read_all=1)
    with pytest.raises(ValueError,match='unknown PWSCF quantity'):
        PwscfAnalyzer(path=tmp_path,required='timing')
    with pytest.raises(ValueError,match='unknown PWSCF quantity'):
        PwscfAnalyzer(path=tmp_path,required=[None])
#end def test_source_required_and_strict_validation


def test_file_requirements_are_checked_at_analysis_time(tmp_path):
    from ..pwscf_analyzer import PwscfAnalyzer

    analyzer = PwscfAnalyzer(
        input   = 'missing.in',
        outfile = 'missing.out',
        path    = tmp_path,
        analyze = False,
        )
    with pytest.raises(FileNotFoundError,match='PWSCF input file is not available'):
        analyzer.analyze()

    direct = PwscfAnalyzer(input=tmp_path/'direct.in',analyze=False)
    with pytest.raises(FileNotFoundError,match='PWSCF input file is not available'):
        direct.analyze()
    permissive = PwscfAnalyzer(
        input = 'missing.in',
        path  = tmp_path,
        analyze = True,
        strict = False,
        )
    assert(permissive.input is None)
    assert(permissive.energy() is None)

    for read_all in (False,True):
        analyzer = PwscfAnalyzer(
            path = tmp_path,
            analyze = False,
            read_all = read_all,
            strict = False,
            )
        with pytest.raises(RuntimeError,match='has not been analyzed'):
            analyzer.energy()
        with pytest.raises(RuntimeError,match='has not been analyzed'):
            analyzer.available('energy')
        analyzer.analyze()
        assert(analyzer.energy() is None)
        assert(not analyzer.available('energy'))
    #end for

    with pytest.raises(FileNotFoundError,match='schema XML'):
        PwscfAnalyzer(path=tmp_path,xmlfile='missing.xml',analyze=True)
    with pytest.raises(FileNotFoundError,match='output file'):
        PwscfAnalyzer(path=tmp_path,outfile='missing.out',analyze=True)
    with pytest.raises(FileNotFoundError) as error:
        PwscfAnalyzer(path=tmp_path,analyze=True)
    message = str(error.value)
    assert('schema XML' in message)
    assert('output file' in message)
#end def test_file_requirements_are_checked_at_analysis_time


def test_strict_source_file_matrix(tmp_path):
    from ..pwscf_analyzer import PwscfAnalyzer

    write_input(tmp_path)
    write_schema(tmp_path)
    (tmp_path/'pwscf.out').write_text(
        'Self-consistent Calculation\n'
        '! total energy = -2.0 Ry\n'
        )

    xml = PwscfAnalyzer(input='pwscf.in',path=tmp_path,analyze=True,read_all = False)
    assert(xml.results_xml is not None)
    assert(xml.results_out is None)
    both = PwscfAnalyzer(input='pwscf.in',path=tmp_path,analyze=True)
    assert(both.results_xml is not None)
    assert(both.results_out is not None)
    assert(both.calculation=='scf')

    (tmp_path/'pwscf.save'/'data-file-schema.xml').unlink()
    out_only = PwscfAnalyzer(input='pwscf.in',path=tmp_path,analyze=True)
    assert(out_only.results_xml is None)
    assert(out_only.results_out is not None)
    assert(out_only.calculation=='scf')
#end def test_strict_source_file_matrix


def test_ambiguous_source_discovery(tmp_path):
    from ..pwscf_analyzer import PwscfAnalyzer

    write_schema(tmp_path,energy=-1.0,prefix='one')
    write_schema(tmp_path,energy=-2.0,prefix='two')
    with pytest.raises(RuntimeError,match='multiple.*schema XML'):
        PwscfAnalyzer(path=tmp_path,analyze=True,read_all=False)
    permissive_xml = PwscfAnalyzer(
        path = tmp_path,
        analyze = True,
        read_all = False,
        strict = False,
        )
    assert(permissive_xml.results_xml is None)
    assert(permissive_xml.energy() is None)

    (tmp_path/'one.out').write_text('Self-consistent Calculation\n')
    (tmp_path/'two.out').write_text('Self-consistent Calculation\n')
    with pytest.raises(RuntimeError,match='multiple.*output'):
        PwscfAnalyzer(path=tmp_path,analyze=True)
    permissive_out = PwscfAnalyzer(
        path = tmp_path,
        analyze = True,
        strict = False,
        )
    assert(permissive_out.results_out is None)

    (tmp_path/'one.out').unlink()
    (tmp_path/'two.out').unlink()
    write_input(tmp_path,prefix='expected')
    with pytest.raises(FileNotFoundError,match='schema XML'):
        PwscfAnalyzer(
            input   = 'pwscf.in',
            xmlfile = 'missing.xml',
            path    = tmp_path,
            analyze = True,
            read_all = False,
            )
    with pytest.raises(FileNotFoundError,match='output file'):
        PwscfAnalyzer(
            input   = 'pwscf.in',
            outfile = 'missing.out',
            path    = tmp_path,
            analyze = True,
            )
#end def test_ambiguous_source_discovery


def test_unique_output_discovery_excludes_auxiliary(tmp_path):
    from ..pwscf_analyzer import PwscfAnalyzer

    (tmp_path/'actual.out').write_text(
        'Self-consistent Calculation\n'
        '! total energy = -4.0 Ry\n'
        )
    (tmp_path/'pw2casino.out').write_text(
        'Kinetic energy from orbitals = 2.0\n'
        )
    analyzer = PwscfAnalyzer(
        path = tmp_path,
        pw2c_outfile_name = 'pw2casino.out',
        analyze = True,
        )
    assert(analyzer.results_out is not None)
    assert(analyzer.energy('Ry')==-4.0)
    assert(analyzer.pw2casino.K==2.0)
#end def test_unique_output_discovery_excludes_auxiliary


def test_required_controls_both_source_fallback(tmp_path):
    import numpy as np

    from ..pwscf_analyzer import PwscfAnalyzer

    write_input(tmp_path)
    write_schema(tmp_path,energy=-1.0)
    (tmp_path/'pwscf.out').write_text('''\
Self-consistent Calculation
number of atoms/cell = 1
! total energy = -3.0 Ry
Forces acting on atoms
atom 1 type 1 force = 0.1 0.2 0.3
''')

    xml_suffices = PwscfAnalyzer(
        input = 'pwscf.in',
        path  = tmp_path,
        analyze = True,
        read_all = False,
        required = 'energy',
        )
    assert(xml_suffices.results_xml is not None)
    assert(xml_suffices.results_out is None)
    assert(xml_suffices.energy('Ha')==-1.0)
    assert(xml_suffices.source_status.xml=='parsed')
    assert(xml_suffices.source_status.out=='skipped')

    (tmp_path/'pwscf.out').unlink()
    xml_without_output = PwscfAnalyzer(
        input = 'pwscf.in',
        path  = tmp_path,
        analyze = True,
        read_all = False,
        required = 'energy',
        )
    assert(xml_without_output.results_out is None)
    assert(xml_without_output.energy('Ha')==-1.0)
    (tmp_path/'pwscf.out').write_text('''\
Self-consistent Calculation
number of atoms/cell = 1
! total energy = -3.0 Ry
Forces acting on atoms
atom 1 type 1 force = 0.1 0.2 0.3
''')

    fallback = PwscfAnalyzer(
        input = 'pwscf.in',
        path  = tmp_path,
        analyze = True,
        read_all = False,
        required = ('energy','forces'),
        )
    assert(fallback.results_out is not None)
    assert(np.allclose(fallback.forces('Ry/B'),[[.1,.2,.3]]))
    assert(fallback.source_status.out=='parsed')

    no_requirements = PwscfAnalyzer(input='pwscf.in',path=tmp_path,analyze=True)
    assert(no_requirements.results_out is not None)

    xml_only = PwscfAnalyzer(
        input = 'pwscf.in',
        path  = tmp_path,
        analyze = True,
        read_all = False,
        )
    xml_only.require('forces')
    assert(xml_only.results_out is None)
    with pytest.raises(RuntimeError,match='required PWSCF quantity "forces"'):
        xml_only.forces()
#end def test_required_controls_both_source_fallback


def test_require_only_updates_policy(tmp_path):
    from ..pwscf_analyzer import PwscfAnalyzer

    write_input(tmp_path)
    write_schema(tmp_path)
    (tmp_path/'pwscf.out').write_text('Self-consistent Calculation\n')
    analyzer = PwscfAnalyzer(
        input = 'pwscf.in',
        path  = tmp_path,
        analyze = True,
        read_all = False,
        required = 'energy',
        )
    assert(analyzer.results_out is None)
    analyzer.require('forces')
    assert(analyzer.required=={'energy','forces'})
    assert(analyzer.results_out is None)
    with pytest.raises(RuntimeError,match='required PWSCF quantity "forces"'):
        analyzer.forces()

    before = set(analyzer.required)
    with pytest.raises(ValueError,match='unknown PWSCF quantity'):
        analyzer.require('stress','timing')
    assert(analyzer.required==before)
    assert(analyzer.require() is None)
    analyzer.analyze()
    assert(analyzer.results_out is not None)
#end def test_require_only_updates_policy


def test_available_is_policy_free(tmp_path):
    from ..pwscf_analyzer import PwscfAnalyzer

    write_schema(tmp_path,energy=0.0)
    analyzer = PwscfAnalyzer(
        path = tmp_path,
        analyze = True,
        read_all = False,
        required = ('energy','forces'),
        )
    assert(analyzer.available())
    assert(analyzer.available('energy'))
    assert(not analyzer.available('energy','forces'))
    assert(not analyzer.available('relaxed_structure'))
    with pytest.raises(ValueError,match='unknown PWSCF quantity'):
        analyzer.available('energy','timing')
    with pytest.raises(RuntimeError,match='required PWSCF quantity "forces"'):
        analyzer.forces()

    fixture = TEST_DIR/'test_pwscf_analyzer_files'/'qe_7_0'/'high'/'scf'
    scf = PwscfAnalyzer(
        input   = 'pwscf.in',
        outfile = 'pwscf.out',
        path    = fixture,
        analyze = True,
        )
    assert(scf.fractional_occs() is False)
    assert(scf.available('fractional_occs'))
#end def test_available_is_policy_free


def test_internal_query_dependencies_ignore_requirements(tmp_path):
    from ..pwscf_analyzer import PwscfAnalyzer

    outfile = tmp_path/'pwscf.out'
    outfile.write_text('''\
Self-consistent Calculation
total stress (Ry/bohr**3) (kbar) P = 12.0
malformed stress row
''')
    analyzer = PwscfAnalyzer(
        outfile = outfile,
        analyze = True,
        required = ('stress','occupations'),
        )
    assert(analyzer.pressure('kbar')==12.0)
    assert(analyzer.fractional_occs() is None)
    with pytest.raises(RuntimeError,match='required PWSCF quantity "stress"'):
        analyzer.stress()
#end def test_internal_query_dependencies_ignore_requirements


def test_malformed_xml_field_falls_back_to_output(tmp_path):
    from ..pwscf_analyzer import PwscfAnalyzer

    write_input(tmp_path)
    schema_file = write_schema(tmp_path)
    schema_file.write_text('''\
<espresso>
  <input><control_variables><calculation>scf</calculation></control_variables></input>
  <output><total_energy><etot>malformed</etot></total_energy></output>
</espresso>
''')
    (tmp_path/'pwscf.out').write_text(
        'Self-consistent Calculation\n'
        '! total energy = -5.0 Ry\n'
        )
    analyzer = PwscfAnalyzer(input='pwscf.in',path=tmp_path,analyze=True)
    assert(analyzer.results_xml.total_energy is None)
    assert(analyzer.energy('Ry')==-5.0)
#end def test_malformed_xml_field_falls_back_to_output


def test_legacy_xml_does_not_satisfy_modern_xml_source(tmp_path):
    from ..pwscf_analyzer import PwscfAnalyzer

    savedir = tmp_path/'pwscf.save'
    savedir.mkdir()
    (savedir/'data-file.xml').write_text('<legacy/>')
    with pytest.raises(FileNotFoundError,match='schema XML'):
        PwscfAnalyzer(
            path    = tmp_path,
            xmlfile = 'missing.xml',
            analyze = True,
            read_all = False,
            )

    analyzer = PwscfAnalyzer(
        path = tmp_path,
        analyze = True,
        read_all = False,
        strict = False,
        )
    assert(analyzer.results_xml is None)
    assert(analyzer.energy() is None)
#end def test_legacy_xml_does_not_satisfy_modern_xml_source


def test_calculation_disagreement_tuple(tmp_path):
    from ..pwscf_analyzer import PwscfAnalyzer

    write_input(tmp_path,calculation='relax')
    write_schema(tmp_path,calculation='nscf')
    (tmp_path/'pwscf.out').write_text('Self-consistent Calculation\n')

    with pytest.raises(RuntimeError,match='top-level calculation types disagree'):
        PwscfAnalyzer(
            input = 'pwscf.in',
            path  = tmp_path,
            analyze = True,
            read_all = False,
            required = 'energy',
            )

    skipped_out = PwscfAnalyzer(
        input = 'pwscf.in',
        path  = tmp_path,
        analyze = True,
        read_all = False,
        required = 'energy',
        strict = False,
        )
    assert(skipped_out.calculation==('relax','nscf',None))

    all_sources = PwscfAnalyzer(
        input = 'pwscf.in',
        path  = tmp_path,
        analyze = True,
        strict = False,
        )
    assert(all_sources.calculation==('relax','nscf','scf'))
#end def test_calculation_disagreement_tuple


def test_output_reader_selection_uses_detected_run_type(tmp_path):
    import numpy as np

    from ..pwscf_analyzer import PwscfAnalyzer, PwscfOutData

    outfile = tmp_path/'pwscf.out'
    outfile.write_text('''\
Band Structure Calculation
highest occupied level (ev): 1.0
number of atoms/cell = 1
! total energy = -2.0 Ry
Forces acting on atoms
atom 1 type 1 force = 0.1 0.2 0.3
ATOMIC_POSITIONS (bohr)
H 0.0 0.0 0.0
End final coordinates
    ''')
    out = PwscfOutData(outfile)
    assert(out.calculation=='nscf')
    assert(out.run_type_detected)
    assert(out.E is None)
    assert(out.forces is None)
    assert(out.relax_structures is None)

    outfile.write_text('''\
highest occupied level (ev): 1.0
number of atoms/cell = 1
! total energy = -2.0 Ry
Forces acting on atoms
atom 1 type 1 force = 0.1 0.2 0.3
ATOMIC_POSITIONS (bohr)
H 0.0 0.0 0.0
End final coordinates
''')
    out = PwscfOutData(outfile)
    assert(out.calculation is None)
    assert(not out.run_type_detected)
    assert(out.E==-2.0)
    assert(np.allclose(out.forces,[[[.1,.2,.3]]]))
    assert(len(out.relax_structures)==1)

    analyzer = PwscfAnalyzer(
        outfile = outfile,
        analyze = True,
        )
    assert(analyzer.calculation is None)
    assert(np.allclose(analyzer.forces('Ry/B'),[[.1,.2,.3]]))
#end def test_output_reader_selection_uses_detected_run_type


def test_undetected_run_type_permissively_reads_md(tmp_path):
    import numpy as np

    from ..pwscf_analyzer import PwscfOutData

    outfile = tmp_path/'unknown.out'
    outfile.write_text('''\
! total energy = -1.0 Ry
total stress (Ry/bohr**3) (kbar) P= 2.0
time = 0.5
kinetic energy = 0.25 Ry
temperature = 300 K
''')
    out = PwscfOutData(outfile)
    assert(out.calculation is None)
    assert(not out.run_type_detected)
    assert(out.md_data is not None)
    assert(np.allclose(out.md_data.total_energy,[-1.0]))
    assert(out.E==-1.0)

    md_only = PwscfOutData(outfile,md_only=True)
    assert(md_only.md_data is not None)
    assert(md_only.E is None)
#end def test_undetected_run_type_permissively_reads_md
