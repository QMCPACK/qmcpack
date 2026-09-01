import pytest

from . import TEST_DIR, NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.RMG_ANALYZER)


representative_root = TEST_DIR/'test_rmg_analyzer_files'
representative_runs = (
    pytest.param(
        'electronic/input.scf.02.log',
        'scf',
        -10.93257703,
        id = 'scf',
        ),
    pytest.param(
        'electronic/input.nscf.05.log',
        'nscf',
        0.0,
        id = 'nscf',
        ),
    pytest.param(
        'ionic/relax/input.01.log',
        'relax',
        -10.93250973,
        id = 'relax',
        ),
    )



def test_empty_init():
    from ..developer import obj
    from ..rmg_analyzer import RmgAnalyzer

    analyzer = RmgAnalyzer()

    expected_members = {
        'abspath','info','input','outfile_name','path','results','run_mode',
        }
    assert set(analyzer.keys())==expected_members
    assert analyzer.path is None
    assert analyzer.abspath is None
    assert analyzer.outfile_name is None
    assert analyzer.input is None
    assert analyzer.run_mode is None
    assert analyzer.results is None
    assert isinstance(analyzer.info,obj)
    assert len(analyzer.info)==0
#end def test_empty_init


def rmg_log(calculation_type,body=''):
    return f'''
Files
   Control input file:        input
   Data output file:          Waves/wave.out

Run Setup
    Calculation type:         {calculation_type}
    Description:              analyzer test

Grid Points (Anisotropy: 1.000)
    X: Total: 8 Per PE: 4 Spacing: 0.5 a0
    Y: Total: 8 Per PE: 4 Spacing: 0.5 a0
    Z: Total: 8 Per PE: 8 Spacing: 0.5 a0
    Equivalent energy cutoffs: 20.0 80.0 Ha

Lattice Setup (a0)
    X Basis Vector: 4.0 0.0 0.0
    Y Basis Vector: 0.0 4.0 0.0
    Z Basis Vector: 0.0 0.0 4.0

Initial Ionic Positions And Displacements (Bohr)
Species      X           Y           Z           dX          dY          dZ
  H        1.0000      1.0000      1.0000      0.0000      0.0000      0.0000

Initial Ionic Positions And Displacements (Angstrom)
{body}

--------TIMING INFORMATION FOR Processor owned the most atoms----------------
                                        Total time               Per SCF/step
1-TOTAL                                             3.00                0.50
'''
#end def rmg_log


@pytest.mark.parametrize(
    argnames='calculation_type,short_mode',
    argvalues=[
        ('Quench electrons','scf'),
        ('NSCF calculate','nscf'),
        ('Structure Optimization.','relax'),
        ],
    )
def test_run_modes(tmp_path,calculation_type,short_mode):
    from ..rmg_analyzer import RmgAnalyzer, RmgOutData
    from ..structure import Structure

    logfile = tmp_path/'rmg.log'
    logfile.write_text(rmg_log(calculation_type))
    outdata = RmgOutData(str(logfile))
    assert outdata.run_mode==short_mode
    assert outdata.setup_info.run_mode==short_mode
    assert isinstance(outdata.setup_info.structure,Structure)

    analyzer = RmgAnalyzer(str(logfile),analyze=True)

    assert analyzer.results.run_mode==short_mode
    assert analyzer.initial_structure() is not None
#end def test_run_modes


@pytest.mark.parametrize(
    argnames='relative_path,run_mode,energy',
    argvalues=representative_runs,
    )
def test_representative_outputs(relative_path,run_mode,energy):
    import numpy as np

    from ..rmg_analyzer import RmgAnalyzer, RmgOutData
    from ..structure import Structure

    logfile  = representative_root/relative_path
    outdata  = RmgOutData(str(logfile))
    analyzer = RmgAnalyzer(str(logfile),analyze=True)

    assert outdata.run_mode==run_mode
    assert analyzer.run_mode==run_mode
    assert isinstance(analyzer.initial_structure(),Structure)
    assert np.isclose(analyzer.energy(),energy)
#end def test_representative_outputs


def test_representative_physical_results():
    from ..rmg_analyzer import RmgAnalyzer
    from ..structure import Structure

    electronic = representative_root/'electronic'
    scf         = RmgAnalyzer(str(electronic/'input.scf.02.log'),analyze=True)
    nscf        = RmgAnalyzer(str(electronic/'input.nscf.05.log'),analyze=True)

    assert scf.kpoints().shape==(1,3)
    assert scf.kweights() is None
    assert scf.eigenvalues().shape==(1,14)
    assert scf.occupations().shape==(1,14)
    assert scf.forces().shape==(1,2,3)
    assert nscf.eigenvalues().shape==(1,14)
    assert nscf.occupations().shape==(1,14)

    relax = RmgAnalyzer(
        str(representative_root/'ionic/relax/input.01.log'),
        analyze = True,
        )

    assert isinstance(relax.relaxed_structure(),Structure)
    assert relax.forces().shape==(3,2,3)
#end def test_representative_physical_results


def test_physical_results(tmp_path):
    import numpy as np

    from ..rmg_analyzer import RmgAnalyzer

    body = '''
K-points
    Kx Ky Kz Weight in crystal unit
    0.0 0.0 0.0 0.25
    0.5 0.0 0.0 0.75

KOHN SHAM EIGENVALUES [eV] AT K-POINT [  0]: 0.0 0.0 0.0
[kpt 0 0 0] -2.0 [2.000] 1.0 [0.000]
KOHN SHAM EIGENVALUES [eV] AT K-POINT [  1]: 0.5 0.0 0.0
[kpt 1 0 0] -1.0 [1.500] 2.0 [0.500]
FERMI ENERGY = 5.25 eV
spinup: valence band maximum = 4.0 eV, conduction band minumm = 6.0 eV
spinup: Band gap = 2.0 eV
final total energy from eig sum = -1.23450000 Ha

@ION  Ion  Species       X           Y           Z       Charge  Mag       FX          FY         FZ      Movable
@ION    1     H      1.1000000   1.2000000   1.3000000    0.050  0.100   0.0100000   0.0200000   0.0300000  1 1 1

 stress total in unit of kbar
 1.0 0.1 0.2
 0.1 2.0 0.3
 0.2 0.3 3.0
'''
    logfile = tmp_path/'scf.log'
    logfile.write_text(rmg_log('Quench electrons',body))
    analyzer = RmgAnalyzer(str(logfile),analyze=True)

    from ..structure import Structure
    from ..unit_converter import convert

    assert isinstance(analyzer.initial_structure(),Structure)
    assert analyzer.initial_structure().units=='A'
    assert analyzer.initial_structure(units='B').units=='B'
    assert analyzer.energy()==-1.2345
    assert np.isclose(analyzer.energy(units='Ry'),-2.469)
    assert analyzer.kpoints().shape==(2,3)
    assert np.allclose(
        analyzer.kpoints(units='A'),analyzer.kpoints()*convert(1.0,'A','B'))
    assert np.allclose(analyzer.kweights(),[0.25,0.75])
    assert analyzer.eigenvalues().shape==(2,2)
    assert np.allclose(
        analyzer.eigenvalues(units='Ha'),convert(analyzer.eigenvalues(),'eV','Ha'))
    assert analyzer.occupations().shape==(2,2)
    assert analyzer.Ef()==5.25
    assert analyzer.Evbm()==4.0
    assert analyzer.Ecbm()==6.0
    assert analyzer.band_gap()==2.0
    assert analyzer.fractional_occs()
    force_factor = convert(1.0,'Ha','eV')/convert(1.0,'B','A')
    assert np.allclose(analyzer.forces(),analyzer.results.forces*force_factor)
    assert np.allclose(analyzer.forces(units='Ha/B'),analyzer.results.forces)
    assert np.allclose(analyzer.forces(units='Ry/B'),2*analyzer.results.forces)
    assert np.allclose(analyzer.stress(),analyzer.results.stress*0.1)
    assert np.allclose(analyzer.stress(units='kbar'),analyzer.results.stress)
    assert np.isclose(analyzer.pressure(),-0.2)
    assert np.isclose(analyzer.pressure(units='kbar'),-2.0)
    with pytest.raises(RuntimeError,match='relaxed_structure'):
        analyzer.relaxed_structure()
    with pytest.raises(ValueError,match='energy units'):
        analyzer.energy(units='J')

    relax_log = tmp_path/'relax.log'
    relax_log.write_text(rmg_log('Structure Optimization.',body))
    relax = RmgAnalyzer(str(relax_log),analyze=True)
    assert isinstance(relax.relaxed_structure(),Structure)
    assert relax.relaxed_structure().units=='A'
    assert np.allclose(
        relax.relaxed_structure(units='B').pos[0],[1.1,1.2,1.3])
#end def test_physical_results


def test_missing_property_data(tmp_path):
    from ..rmg_analyzer import RmgAnalyzer

    logfile = tmp_path/'missing.log'
    logfile.write_text(rmg_log('Quench electrons'))
    analyzer = RmgAnalyzer(str(logfile),analyze=True)

    for name in (
        'energy','kpoints','kweights','eigenvalues','occupations','Ef','Evbm','Ecbm',
        'band_gap','fractional_occs','forces','stress','pressure',
        ):
        assert getattr(analyzer,name)() is None
    with pytest.raises(RuntimeError,match='relaxed_structure'):
        analyzer.relaxed_structure()
    with pytest.raises(RuntimeError,match='has not been analyzed'):
        RmgAnalyzer().energy()
#end def test_missing_property_data


def test_whitespace_and_trailing_fields(tmp_path):
    import numpy as np

    from ..rmg_analyzer import RmgAnalyzer

    log = rmg_log('Quench electrons').replace(
        'Calculation type:', 'Calculation\t type   :').replace(
        '1-TOTAL                                             3.00                0.50',
        '1 - TOTAL\t3.00\t0.50\tnew timing annotation')
    body = '''
FERMI   ENERGY : 5.25 eV trailing diagnostic 77
spin0: conduction band minimum = 6.0 eV, valence band maximum = 4.0 eV extra 88
spin0: Band gap : 2.0 eV extra 99
final   total energy from eigenvalue sum : -1.2345 Ry trailing 42

@ION\tIon\tSpecies X Y Z Charge Mag FX FY FZ Movable
@ION 1 H 1.1 1.2 1.3 0.05 0.10 0.01 0.02 0.03 1 1 1 trailing 77

stress total in unit of kbar

1 1.0 0.1 0.2 trailing
2 0.1 2.0 0.3 trailing
3 0.2 0.3 3.0 trailing
'''
    log = log.replace(
        '\n\n--------TIMING INFORMATION',body+'\n--------TIMING INFORMATION')
    logfile = tmp_path/'spacing.log'
    logfile.write_text(log)

    analyzer = RmgAnalyzer(str(logfile),analyze=True)

    assert analyzer.results.run_mode=='scf'
    assert analyzer.energy()==-0.61725
    assert analyzer.Ef()==5.25
    assert analyzer.Evbm()==4.0
    assert analyzer.Ecbm()==6.0
    assert analyzer.band_gap()==2.0
    assert np.allclose(analyzer.forces(units='Ha/B')[0,0],[0.01,0.02,0.03])
    assert np.allclose(
        analyzer.stress(units='kbar')[0],
        [[1.0,0.1,0.2],[0.1,2.0,0.3],[0.2,0.3,3.0]])
#end def test_whitespace_and_trailing_fields


def test_malformed_sections_do_not_stop_analysis(tmp_path):
    from ..rmg_analyzer import RmgAnalyzer

    body = '''
@@ TOTAL ENERGY = malformed
final total energy from eig sum = malformed Ha
@ION Ion Species X Y Z Charge Mag FX FY FZ Movable
@ION 1 H malformed row
@ION 2 H 2.1 2.2 2.3 0.0 0.0 0.1 0.2 0.3 1 1 1 valid trailing data
'''
    logfile = tmp_path/'malformed.log'
    logfile.write_text(rmg_log('Quench electrons',body))

    analyzer = RmgAnalyzer(str(logfile),analyze=True)

    assert analyzer.energy() is None
    assert analyzer.forces(units='Ha/B').shape==(1,1,3)
    assert len(analyzer.info)==0
#end def test_malformed_sections_do_not_stop_analysis
