import re

import pytest

from . import NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.PWSCF_ANALYZER)


from . import TEST_DIR


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
        {'values': '-5.3120  1.2470  4.9810'},
        ),
    (
        'leading_number_list_pattern',
        '0.0488-0.0345 0.0345  convergence achieved',
        {'values': '0.0488-0.0345 0.0345'},
        ),
    (
        'vector3_pattern',
        '  1.0 0.0 0.0',
        {'x': '1.0','y': '0.0','z': '0.0'},
        ),
    (
        'vector3_pattern',
        '-2.5D-01  .500000  0.  fourth-column',
        {'x': '-2.5D-01','y': '.500000','z': '0.'},
        ),
    (
        'kpoint_table_pattern',
        '     k( 1) = ( 0.0000000 0.0000000 0.0000000), wk = 0.2500000',
        {'kx': '0.0000000','ky': '0.0000000','kz': '0.0000000','weight': '0.2500000'},
        ),
    (
        'kpoint_table_pattern',
        'k(12)=(.5 -.5 5.0D-1), wk=1.0D+00',
        {'kx': '.5','ky': '-.5','kz': '5.0D-1','weight': '1.0D+00'},
        ),
    (
        'band_kpoint_pattern',
        '     k = 0.0000 0.0000 -0.7071 ( 2138 PWs) bands (ev):',
        {'kx': '0.0000','ky': '0.0000','kz': '-0.7071'},
        ),
    (
        'band_kpoint_pattern',
        'k = 0.0488-0.0345 0.0345 ( 42052 PWs) bands (ev):',
        {'kx': '0.0488','ky': '-0.0345','kz': '0.0345'},
        ),
    (
        'band_kpoint_pattern',
        'k =-0.3536 0.3536-0.3536 ( 2138 PWs) bands (ev):',
        {'kx': '-0.3536','ky': '0.3536','kz': '-0.3536'},
        ),
    (
        'total_energy_pattern',
        '     total energy              =    -168.12345678 Ry',
        {'energy': '-168.12345678'},
        ),
    (
        'total_energy_pattern',
        '!    total energy = -1.0D+02 Ry',
        {'energy': '-1.0D+02'},
        ),
    (
        'scf_accuracy_pattern',
        '     estimated scf accuracy    <       6.3E-09 Ry',
        {'accuracy': '6.3E-09'},
        ),
    (
        'scf_accuracy_pattern',
        'estimated scf accuracy = 1.2D-10 Ry',
        {'accuracy': '1.2D-10'},
        ),
    (
        'pressure_pattern',
        'total stress (Ry/bohr**3) (kbar) P= -170.96',
        {'pressure': '-170.96'},
        ),
    (
        'pressure_pattern',
        'total   stress  (Ry/bohr**3) (kbar) P = 1.2D+03 ',
        {'pressure': '1.2D+03'},
        ),
    (
        'volume_pattern',
        '     unit-cell volume          =     380.6210 (a.u.)^3',
        {'volume': '380.6210'},
        ),
    (
        'volume_pattern',
        'unit-cell volume=3.806210D+02 a.u.^3',
        {'volume': '3.806210D+02'},
        ),
    (
        'alat_pattern',
        'CELL_PARAMETERS (alat= 10.20)',
        {'alat': '10.20'},
        ),
    (
        'alat_pattern',
        'CELL_PARAMETERS (alat = 1.0D+01)',
        {'alat': '1.0D+01'},
        ),
    (
        'atomic_force_pattern',
        'atom 1 type 2 force = -0.001 0.002 0.000',
        {'atom': '1','type': '2','fx': '-0.001','fy': '0.002','fz': '0.000'},
        ),
    (
        'atomic_force_pattern',
        'atom 12 type 1 force=1.0D-03 -2.0D-03 .0',
        {'atom': '12','type': '1','fx': '1.0D-03','fy': '-2.0D-03','fz': '.0'},
        ),
    (
        'total_force_pattern',
        'Total force = 0.173046 Total SCF correction = 0.001',
        {'total_force': '0.173046'},
        ),
    (
        'total_force_pattern',
        'Total force=1.25D-04',
        {'total_force': '1.25D-04'},
        ),
    (
        'stress_row_pattern',
        ' -.001 0.0 .001 -147.1 0.0 147.1',
        {'sxx': '-.001','sxy': '0.0','sxz': '.001','kxx': '-147.1','kxy': '0.0','kxz': '147.1'},
        ),
    (
        'stress_row_pattern',
        '1D-3 2D-3 3D-3 1E+2 2E+2 3E+2 trailing',
        {'sxx': '1D-3','sxy': '2D-3','sxz': '3D-3','kxx': '1E+2','kxy': '2E+2','kxz': '3E+2'},
        ),
    (
        'timing_value_pattern',
        '4m33.69s',
        {'value': '4','unit': 'm'},
        ),
    (
        'timing_value_pattern',
        '1h 2m 3.5s',
        {'value': '1','unit': 'h'},
        ),
    (
        'fermi_energies_pattern',
        '     the Fermi energy is    10.1198 ev',
        {'values': '10.1198'},
        ),
    (
        'fermi_energies_pattern',
        'the Fermi energy          =      -3.22772442 eV',
        {'values': '-3.22772442'},
        ),
    (
        'fermi_energies_pattern',
        'the spin up/dw Fermi energies are 5.1 5.2 EV',
        {'values': '5.1 5.2'},
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
    from ..pwscf_analyzer import (
        Pw2CasinoAnalyzer,
        PwscfAnalyzer,
        PwscfOutData,
        )

    pa = PwscfAnalyzer()
    assert(len(pa)==0)
    with pytest.raises(RuntimeError):
        pa.analyze()
    free_helpers = ('match_float',)
    for name in free_helpers:
        assert(callable(getattr(pa_module,name)))
        assert(not hasattr(PwscfAnalyzer,name))
    #end for
    assert(not hasattr(pa_module,'read_kpoint_tables'))
    reader_names = (
        'read_calculation','read_fermi_energies',
        'read_energies','read_scf_convergence','read_bands',
        'read_band_edges','read_structures','read_pressure_volume',
        'read_stress','read_forces','read_timing','read_kpoints',
        )
    assert(all(callable(getattr(PwscfOutData,name)) for name in reader_names))
    assert(not hasattr(PwscfOutData,'read'))
    assert(not any(hasattr(PwscfAnalyzer,'analyze_'+name[5:]) for name in reader_names))
    assert(not hasattr(PwscfOutData,'read_pw2casino'))
    assert(not hasattr(PwscfAnalyzer,'analyze_pw2casino'))
    assert(not hasattr(PwscfAnalyzer,'analyze_schema_xml'))
    assert(not hasattr(Pw2CasinoAnalyzer,'read'))
    assert(not hasattr(pa_module,'PwscfXmlData'))
#end def test_empty_init


@pytest.mark.parametrize('calculation,log_text',(
    ('scf',     'Self-consistent Calculation\n'),
    ('nscf',    'Band Structure Calculation\nhighest occupied level (ev): 1.0\n'),
    ('relax',   'BFGS Geometry Optimization\n'),
    ('vc-relax','BFGS Geometry Optimization\nCELL_PARAMETERS (alat= 1.0)\n'),
    ))
def test_result_initialization(tmp_path,calculation,log_text):
    from ..pwscf_analyzer import PwscfOutData

    outfile = tmp_path / f'{calculation}.out'
    outfile.write_text(log_text)
    out = PwscfOutData(outfile)

    expected = {
        'calculation',
        'Ef','fermi_energies','bands',
        'volume',
        'cputime','walltime','kpoints_cart','kpoints_unit','kweights',
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


@pytest.mark.parametrize(
    'verbosity,calculation',
    tuple(
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
    if verbosity=='high':
        assert(analyzer.occupations().shape==(npoints,nbands))
    else:
        assert(analyzer.occupations() is None)
    #end if
    if calculation=='nscf':
        assert('E' not in out)
        assert('forces' not in out)
        assert('stress' not in out)
        with pytest.raises(RuntimeError):
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


def test_pw2casino_analyzer_read(tmp_path):
    from ..pwscf_analyzer import Pw2CasinoAnalyzer, PwscfAnalyzer, PwscfOutData

    (tmp_path/'pwscf.in').write_text("&CONTROL\n  calculation = 'scf'\n/\n")
    (tmp_path/'pwscf.out').write_text('Self-consistent Calculation\n')
    (tmp_path/'pw2casino.out').write_text(
        'Kinetic energy from orbitals = 1.25D+01\n'
        )
    with pytest.raises(FileNotFoundError):
        PwscfAnalyzer(tmp_path,'missing.in','pwscf.out')
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

    (tmp_path/'pw2casino.out').write_text('Kinetic energy is unavailable\n')
    pw2casino = Pw2CasinoAnalyzer(tmp_path/'pw2casino.out')
    assert(pw2casino.K is None)

    with pytest.raises(FileNotFoundError):
        Pw2CasinoAnalyzer(tmp_path/'missing.out')

    (tmp_path/'pw2casino.out').unlink()
    with pytest.raises(FileNotFoundError):
        PwscfAnalyzer(
            tmp_path,
            'pwscf.in',
            'pwscf.out',
            'pw2casino.out',
            analyze = True,
            )
#end def test_pw2casino_analyzer_read


@pytest.mark.parametrize('calculation',('bands','md','vc-md'))
def test_unsupported_calculation_modes(calculation):
    from ..pwscf_analyzer import PwscfAnalyzer

    fixture_path = (
        TEST_DIR/'test_pwscf_analyzer_files'/'qe_7_0'/'default'/calculation
        )
    with pytest.raises(RuntimeError):
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
    scf = PwscfAnalyzer(fixture_root/'scf','pwscf.in','pwscf.out',analyze=True)
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
        input_fixture,'scf.in','scf.out',analyze=True
        )
    assert(input_scf.initial_structure('A').units=='A')

    nscf = PwscfAnalyzer(fixture_root/'nscf','pwscf.in','pwscf.out',analyze=True)
    assert(nscf.energy() is None)
    assert(nscf.Ef() is None)
    assert(nscf.Evbm() is not None)
    assert(nscf.Ecbm() is not None)
    assert(nscf.band_gap() is not None)
    assert(not nscf.fractional_occs())
    with pytest.raises(RuntimeError):
        nscf.forces()
    with pytest.raises(RuntimeError):
        nscf.relaxed_structure()
    with pytest.raises(ValueError):
        nscf.energy('invalid')
    with pytest.raises(ValueError):
        nscf.kpoints('invalid')
    with pytest.raises(ValueError):
        nscf.eigenvalues('invalid')

    fermi_fixture = TEST_DIR/'test_pwscf_analyzer_files'/'nscf_output'
    fermi_nscf = PwscfAnalyzer(
        fermi_fixture,'nscf.in','nscf.out',analyze=True
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
