import pytest
import numpy as np
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.RMG_ANALYZER)




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
    return '''
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
'''.format(calculation_type=calculation_type,body=body)
#end def rmg_log


@pytest.mark.parametrize('calculation_type,short_mode',[
    ('Quench electrons','scf'),
    ('NSCF calculate','nscf'),
    ('Band structure calculation.','band'),
    ("calculate Exx integral's from saved wave functions",'exx'),
    ('Structure Optimization.','relax'),
    ('Molecular dynamics - CVE','md_VE'),
    ('Molecular dynamics - CVT','md_TE'),
    ('Time dependent DFT (TDDFT) calculation','tddft'),
    ('calculate STM charge density','stm'),
    ('Molecular dynamics using Nudged Elastic Band.','neb'),
    ])
def test_run_modes(tmp_path,calculation_type,short_mode):
    from ..rmg_analyzer import RmgAnalyzer, RmgOutData

    logfile = tmp_path/'rmg.log'
    logfile.write_text(rmg_log(calculation_type))
    outdata = RmgOutData(str(logfile))
    metadata_fields = {
        'abspath','input','outfile_name','path','run_mode','setup_info',
        }
    result_fields = set(outdata.keys())

    assert outdata.run_mode==short_mode
    assert outdata.setup_info.run_mode==short_mode
    assert metadata_fields<result_fields
    assert {'geometry','convergence','timing'}<result_fields

    electronic_run = short_mode in {
        'scf','nscf','relax','md_VE','md_TE','tddft','neb',
        }
    field_applicability = dict(
        energy         = electronic_run,
        ionic_steps    = electronic_run,
        pressure       = electronic_run,
        electronic     = electronic_run or short_mode=='band',
        bands          = short_mode=='band',
        md             = short_mode in {'md_VE','md_TE'},
        tddft          = short_mode=='tddft',
        produced_files = short_mode in {'scf','exx','stm'},
        )
    for name,applies in field_applicability.items():
        assert (name in result_fields)==applies
    #end for

    analyzer = RmgAnalyzer(str(logfile),analyze=True)

    assert analyzer.results.run_mode==short_mode
    assert set(analyzer.results.keys())==result_fields
    assert analyzer.results.timing is not None
    assert analyzer.results.timing.total==3.0
    assert analyzer.results.setup_info.grid_points.grid.tolist()==[8,8,8]
    assert len(analyzer.results.setup_info.structure.elem)==1
#end def test_run_modes


def test_physical_results(tmp_path):
    from ..rmg_analyzer import RmgAnalyzer

    body = '''
KOHN SHAM EIGENVALUES [eV] AT K-POINT [  0]: 0.0 0.0 0.0
[kpt 0 0 0] -2.0 [2.000] 1.0 [0.000]
KOHN SHAM EIGENVALUES [eV] AT K-POINT [  1]: 0.5 0.0 0.0
[kpt 1 0 0] -1.0 [1.500] 2.0 [0.500]
FERMI ENERGY = 5.25 eV
spinup: valence band maximum = 4.0 eV, conduction band minumm = 6.0 eV
spinup: Band gap = 2.0 eV
@@ EIGENVALUE SUM     =       -0.500000 Ha
@@ ION_ION            =        0.100000 Ha
@@ ELECTROSTATIC      =       -0.200000 Ha
@@ VXC                =       -0.300000 Ha
@@ EXC                =       -0.250000 Ha
@@ TOTAL ENERGY       =       -1.250000 Ha
@@ estimated error    =        1.00e-08 Ha
 quench: [md: 0/2 scf: 3/20 step time: 0.20 scf time: 0.80 secs RMS[dV]: 2.00e-05 ]
final total energy from eig sum = -1.23450000 Ha
final total energy from direct = -1.23440000 Ha
Total charge in supercell = 1.0
SUM FORCE = 0.1 0.2 0.3
 volume and energy per atom = 64.0 -33.0 eV

@ION  Ion  Species       X           Y           Z       Charge  Mag       FX          FY         FZ      Movable
@ION    1     H      1.1000000   1.2000000   1.3000000    0.050  0.100   0.0100000   0.0200000   0.0300000  1 1 1

 stress total in unit of kbar
 1.0 0.1 0.2
 0.1 2.0 0.3
 0.2 0.3 3.0
potential convergence has been achieved. stopping ...
'''
    logfile = tmp_path/'scf.log'
    logfile.write_text(rmg_log('Quench electrons',body))
    analyzer = RmgAnalyzer(str(logfile),analyze=True)

    from ..structure import Structure
    from ..unit_converter import convert

    assert analyzer.results.timing is not None
    assert isinstance(analyzer.input_structure(),Structure)
    assert analyzer.input_structure().units=='A'
    assert analyzer.input_structure('B').units=='B'
    assert analyzer.energy()==-1.2345
    assert np.isclose(analyzer.energy('Ry'),-2.469)
    assert analyzer.kpoints().shape==(2,3)
    assert np.allclose(
        analyzer.kpoints('A'),analyzer.kpoints()*convert(1.0,'A','B'))
    analyzer.results.geometry.kweights = np.array([0.25,0.75])
    assert np.allclose(analyzer.kweights(),[0.25,0.75])
    assert analyzer.eigenvalues().shape==(2,2)
    assert np.allclose(
        analyzer.eigenvalues('Ha'),convert(analyzer.eigenvalues(),'eV','Ha'))
    assert analyzer.occupations().shape==(2,2)
    assert analyzer.Ef()==5.25
    assert analyzer.Evbm()==4.0
    assert analyzer.Ecbm()==6.0
    assert analyzer.band_gap()==2.0
    assert analyzer.fractional_occs()
    assert analyzer.results.energy==-1.2345
    assert analyzer.results.electronic.fermi_energies[-1]==5.25
    assert analyzer.results.energy_units=='Ha'
    assert analyzer.results.scf.scf_steps.tolist()==[3]
    assert np.allclose(analyzer.results.forces[0,0],[0.01,0.02,0.03])
    assert np.allclose(analyzer.results.positions[0,0],[1.1,1.2,1.3])
    assert np.allclose(
        analyzer.results.stress[0],[[1.0,0.1,0.2],[0.1,2.0,0.3],[0.2,0.3,3.0]])
    assert analyzer.results.pressure==-2.0
    assert analyzer.results.convergence.electronic_converged
    force_factor = convert(1.0,'Ha','eV')/convert(1.0,'B','A')
    assert np.allclose(analyzer.forces(),analyzer.results.forces*force_factor)
    assert np.allclose(analyzer.forces('Ha/B'),analyzer.results.forces)
    assert np.allclose(analyzer.forces('Ry/B'),2*analyzer.results.forces)
    assert np.allclose(analyzer.stress(),analyzer.results.stress*0.1)
    assert np.allclose(analyzer.stress('kbar'),analyzer.results.stress)
    assert np.isclose(analyzer.pressure(),-0.2)
    assert np.isclose(analyzer.pressure('kbar'),-2.0)
    with pytest.raises(RuntimeError,match='relaxed_structure'):
        analyzer.relaxed_structure()
    with pytest.raises(ValueError,match='energy units'):
        analyzer.energy('J')

    relax_log = tmp_path/'relax.log'
    relax_log.write_text(rmg_log('Structure Optimization.',body))
    relax = RmgAnalyzer(str(relax_log),analyze=True)
    assert isinstance(relax.relaxed_structure(),Structure)
    assert relax.relaxed_structure().units=='A'
    assert np.allclose(relax.relaxed_structure('B').pos[0],[1.1,1.2,1.3])
#end def test_physical_results


def test_missing_property_data(tmp_path):
    from ..rmg_analyzer import RmgAnalyzer

    logfile = tmp_path/'missing.log'
    logfile.write_text(rmg_log('Quench electrons'))
    analyzer = RmgAnalyzer(str(logfile),analyze=True)

    for name in {
        'energy','kpoints','kweights','eigenvalues','occupations','Ef','Evbm','Ecbm',
        'band_gap','fractional_occs','forces','stress','pressure',
        }:
        assert getattr(analyzer,name)() is None
    #end for
    with pytest.raises(RuntimeError,match='relaxed_structure'):
        analyzer.relaxed_structure()
    with pytest.raises(RuntimeError,match='has not been analyzed'):
        RmgAnalyzer().energy()
#end def test_missing_property_data


def test_mode_specific_results(tmp_path):
    from ..rmg_analyzer import RmgAnalyzer

    # Band structure companion file
    band_log = tmp_path/'band.log'
    band_log.write_text(rmg_log('Band structure calculation.'))
    (tmp_path/'band_spin0.bandstructure.dat').write_text(
        '0.0 -1.0\n1.0 -0.5\n&&\n0.0 1.0\n1.0 1.5\n&&\n')
    band = RmgAnalyzer(str(band_log),analyze=True)
    assert band.results.bands[0].energies.shape==(2,2)
    assert band.eigenvalues().shape==(2,2)
    with pytest.raises(RuntimeError,match='energy'):
        band.energy()
    with pytest.raises(RuntimeError,match='occupations'):
        band.occupations()

    # Molecular dynamics thermodynamic record
    md_log = tmp_path/'md.log'
    md_log.write_text(rmg_log('Molecular dynamics - CVE',
        '@CVE 1 -1.0 0.1 -0.9 300.0 2.5e-4'))
    md = RmgAnalyzer(str(md_log),analyze=True)
    assert md.results.md.step.tolist()==[1]
    assert md.results.md.temperature.tolist()==[300.0]
    assert md.results.md_stats.temperature.mean==300.0

    # TDDFT energy and dipole time series
    td_log = tmp_path/'td.log'
    td_log.write_text(rmg_log('Time dependent DFT (TDDFT) calculation'))
    (tmp_path/'td_totalE').write_text(
        '&& time kinetic hartree xc total\n0.0 1.0 2.0 3.0 4.0\n')
    (tmp_path/'td_spin0_dipole.dat').write_text(
        '&&electric field in cartesian unit: 0.0 0.0 0.1\n'
        '&&dipole at ground state: 1.0 2.0 3.0\n'
        '0.0 1.1 2.1 3.1\n')
    td = RmgAnalyzer(str(td_log),analyze=True)
    assert td.results.tddft.energy.total_energy_change.tolist()==[4.0]
    assert np.allclose(td.results.tddft.dipoles[0].dipole,[[1.1,2.1,3.1]])

    # Exact-exchange and STM products
    exx_log = tmp_path/'exx.log'
    exx_log.write_text(rmg_log("calculate Exx integral's from saved wave functions"))
    (tmp_path/'exx_integrals.h5').touch()
    exx = RmgAnalyzer(str(exx_log),analyze=True)
    assert len(exx.results.produced_files.exx_integrals)==1

    stm_dir = tmp_path/'STM'
    stm_dir.mkdir()
    (stm_dir/'bias.stm').touch()
    (stm_dir/'bias.cube').touch()
    stm_log = tmp_path/'stm.log'
    stm_log.write_text(rmg_log('calculate STM charge density'))
    stm = RmgAnalyzer(str(stm_log),analyze=True)
    assert len(stm.results.produced_files.stm)==1
    assert len(stm.results.produced_files.stm_cube)==1
#end def test_mode_specific_results


def test_whitespace_and_trailing_fields(tmp_path):
    from ..rmg_analyzer import RmgAnalyzer

    log = rmg_log('Quench electrons').replace(
        'Calculation type:', 'Calculation\t type   :').replace(
        '1-TOTAL                                             3.00                0.50',
        '1 - TOTAL\t3.00\t0.50\tnew timing annotation')
    body = '''
FERMI   ENERGY : 5.25 eV trailing diagnostic 77
spin0: valence band maximum = 4.0 eV, conduction band minimum = 6.0 eV extra 88
spin0: Band gap : 2.0 eV extra 99
@@   TOTAL ENERGY : -1.250000 Ha extra 123
@@ estimated error = 1.0D-8 Ha extra 456
quench: [ RMS [ dV ] : 2.0D-5 scf time: 0.80 md: 0/2 extra 44 step time: 0.20 scf: 3/20 ]
final   total energy from eigenvalue sum : -1.2345 Ry trailing 42
SUM FORCE = 0.1 0.2 0.3 trailing 55
volume and energy per atom = 64.0 -33.0 eV trailing 66

@ION\tIon\tSpecies X Y Z Charge Mag FX FY FZ Movable
@ION 1 H 1.1 1.2 1.3 0.05 0.10 0.01 0.02 0.03 1 1 1 trailing 77

stress total in unit of kbar

1 1.0 0.1 0.2 trailing
2 0.1 2.0 0.3 trailing
3 0.2 0.3 3.0 trailing
potential   convergence has been achieved trailing status
'''
    log = log.replace(
        '\n\n--------TIMING INFORMATION',body+'\n--------TIMING INFORMATION')
    logfile = tmp_path/'spacing.log'
    logfile.write_text(log)

    analyzer = RmgAnalyzer(str(logfile),analyze=True)

    assert analyzer.results.run_mode=='scf'
    assert analyzer.results.timing is not None
    assert analyzer.results.energy_units=='Ry'
    assert analyzer.results.energy==-1.2345
    assert analyzer.results.electronic.fermi_energies[-1]==5.25
    assert analyzer.results.electronic.valence_band_maxima.tolist()==[4.0]
    assert analyzer.results.electronic.conduction_band_minima.tolist()==[6.0]
    assert analyzer.results.electronic.sum_forces.tolist()==[[0.1,0.2,0.3]]
    assert analyzer.results.scf.scf_steps.tolist()==[3]
    assert analyzer.results.scf.step_times.tolist()==[0.2]
    assert np.allclose(analyzer.results.forces[0,0],[0.01,0.02,0.03])
    assert np.allclose(
        analyzer.results.stress[0],[[1.0,0.1,0.2],[0.1,2.0,0.3],[0.2,0.3,3.0]])
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

    assert analyzer.results.timing is not None
    assert len(analyzer.results.ionic_steps)==1
    assert analyzer.results.positions.shape==(1,1,3)
    assert analyzer.results.scf is None
    assert len(analyzer.info)==0
#end def test_malformed_sections_do_not_stop_analysis
