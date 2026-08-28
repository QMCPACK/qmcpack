import pytest
import numpy as np
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.RMG_ANALYZER)




def test_empty_init():
    from ..rmg_analyzer import RmgAnalyzer

    RmgAnalyzer()
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
    from ..rmg_analyzer import RmgAnalyzer

    logfile = tmp_path/'rmg.log'
    logfile.write_text(rmg_log(calculation_type))
    analyzer = RmgAnalyzer(str(logfile),analyze=True)

    assert analyzer.calculation_shortmode==short_mode
    assert analyzer.run_completed
    assert analyzer.results.timing.total==3.0
    assert analyzer.setup_info.grid_points.grid.tolist()==[8,8,8]
    assert len(analyzer.setup_info.structure.elem)==1
#end def test_run_modes


def test_physical_results(tmp_path):
    from ..rmg_analyzer import RmgAnalyzer

    body = '''
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

    assert analyzer.analysis_succeeded
    assert analyzer.run_completed
    assert analyzer.E==-1.2345
    assert analyzer.Ef==5.25
    assert analyzer.results.energy_units=='Ha'
    assert analyzer.scf_data.scf_steps.tolist()==[3]
    assert np.allclose(analyzer.forces[0,0],[0.01,0.02,0.03])
    assert np.allclose(analyzer.positions[0,0],[1.1,1.2,1.3])
    assert np.allclose(analyzer.stress[0],[[1.0,0.1,0.2],[0.1,2.0,0.3],[0.2,0.3,3.0]])
    assert analyzer.pressure==-2.0
    assert analyzer.convergence.electronic_converged
#end def test_physical_results


def test_mode_specific_results(tmp_path):
    from ..rmg_analyzer import RmgAnalyzer

    # Band structure companion file
    band_log = tmp_path/'band.log'
    band_log.write_text(rmg_log('Band structure calculation.'))
    (tmp_path/'band_spin0.bandstructure.dat').write_text(
        '0.0 -1.0\n1.0 -0.5\n&&\n0.0 1.0\n1.0 1.5\n&&\n')
    band = RmgAnalyzer(str(band_log),analyze=True)
    assert band.bands[0].energies.shape==(2,2)

    # Molecular dynamics thermodynamic record
    md_log = tmp_path/'md.log'
    md_log.write_text(rmg_log('Molecular dynamics - CVE',
        '@CVE 1 -1.0 0.1 -0.9 300.0 2.5e-4'))
    md = RmgAnalyzer(str(md_log),analyze=True)
    assert md.md_data.step.tolist()==[1]
    assert md.md_data.temperature.tolist()==[300.0]
    assert md.md_stats.temperature.mean==300.0

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
    assert td.tddft.energy.total_energy_change.tolist()==[4.0]
    assert np.allclose(td.tddft.dipoles[0].dipole,[[1.1,2.1,3.1]])

    # Exact-exchange and STM products
    exx_log = tmp_path/'exx.log'
    exx_log.write_text(rmg_log("calculate Exx integral's from saved wave functions"))
    (tmp_path/'exx_integrals.h5').touch()
    exx = RmgAnalyzer(str(exx_log),analyze=True)
    assert len(exx.artifacts.exx_integrals)==1

    stm_dir = tmp_path/'STM'
    stm_dir.mkdir()
    (stm_dir/'bias.stm').touch()
    (stm_dir/'bias.cube').touch()
    stm_log = tmp_path/'stm.log'
    stm_log.write_text(rmg_log('calculate STM charge density'))
    stm = RmgAnalyzer(str(stm_log),analyze=True)
    assert len(stm.artifacts.stm)==1
    assert len(stm.artifacts.stm_cube)==1
#end def test_mode_specific_results
