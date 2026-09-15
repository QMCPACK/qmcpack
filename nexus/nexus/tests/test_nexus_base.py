import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.NEXUS_BASE)

from . import isolate_nexus_core, TEST_DIR
from ..testing import object_eq


TEST_FILES = {
    "old_nxs_pwscf_input.p": (TEST_DIR / "test_generic_files/old_nxs_pwscf_input.p").resolve(),
    "old_nxs_sim.p":         (TEST_DIR / "test_generic_files/old_nxs_sim.p").resolve(),
    "old_nxs_pwscf_input_numpy_1.p": (TEST_DIR / "test_generic_files/old_nxs_pwscf_input_numpy_1.p").resolve(),
    "old_nxs_sim_numpy_1.p":         (TEST_DIR / "test_generic_files/old_nxs_sim_numpy_1.p").resolve(),
    }

def test_namespaces():
    from ..nexus_base import nexus_core,nexus_core_defaults
    from ..nexus_base import nexus_noncore,nexus_noncore_defaults
    from ..nexus_base import nexus_core_noncore,nexus_core_noncore_defaults
    
    assert('runs' in nexus_core_defaults)
    assert('basis_dir' in nexus_noncore_defaults)
    assert('pseudo_dir' in nexus_core_noncore_defaults)

    assert(object_eq(nexus_core,nexus_core_defaults))
    assert(object_eq(nexus_noncore,nexus_noncore_defaults))
    assert(object_eq(nexus_core_noncore,nexus_core_noncore_defaults))
#end def test_namespaces



def test_empty_init():
    from ..nexus_base import NexusCore
    nc = NexusCore()
#end def test_empty_init


def test_write_splash(capsys):
    from ..nexus_base import write_splash

    assert(not hasattr(write_splash, "wrote_splash"))
    write_splash()
    captured = capsys.readouterr()
    assert('Nexus' in captured.out)
    assert('Please cite:' in captured.out)
    assert(hasattr(write_splash, "wrote_splash"))
    assert(write_splash.wrote_splash)
#end def test_write_splash
    

def test_enter_leave(tmp_path, capsys):
    import os
    from ..nexus_base import NexusCore
    cwd = os.getcwd()

    nc = NexusCore()

    nc.enter(tmp_path)
    tcwd = os.getcwd()
    assert(tcwd==str(tmp_path))
    captured = capsys.readouterr()
    assert('Entering' in captured.out)
    assert(str(tmp_path) in captured.out)

    nc.leave()
    assert(os.getcwd()==cwd)
#end def test_enter_leave


def test_old_nexus_unpickle():
    import numpy as np
    from ..nexus_base import NexusCore

    sim_obj = NexusCore()
    if np.lib.NumpyVersion(np.__version__) >= '2.0.0b1':
        sim_obj.load(TEST_FILES["old_nxs_sim.p"])
    else:
        sim_obj.load(TEST_FILES["old_nxs_sim_numpy_1.p"])

    assert(sim_obj.analyzed            is True)
    assert(sim_obj.analyzer_image      == "analyzer.p")
    assert(sim_obj.app_name            == "pw.x")
    assert(set(sim_obj.app_props)      == set(["serial", "mpi"]))
    assert(sim_obj.block               is False)
    assert(sim_obj.block_subcascade    is False)
    assert(sim_obj.errfile             == "relax.err")
    assert(sim_obj.failed              is False)
    assert(sim_obj.files               == {"Ge.pbe-kjpaw.UPF", "relax.in"})
    assert(sim_obj.finished            is True)
    assert(sim_obj.got_output          is True)
    assert(sim_obj.identifier          == "relax")
    assert(sim_obj.image_dir           == "sim_relax")
    assert(sim_obj.imlocdir            == "./runs/relax/kgrid_111/sim_relax")
    assert(sim_obj.imremdir            == "./runs/relax/kgrid_111/sim_relax")
    assert(sim_obj.imresdir            == "./runs/relax/kgrid_111/sim_relax")
    assert(sim_obj.infile              == "relax.in")
    assert(sim_obj.input_image         == "input.p")
    assert(sim_obj.locdir              == "./runs/relax/kgrid_111")
    assert(sim_obj.outfile             == "relax.out")
    assert(sim_obj.outputs             is None)
    assert(sim_obj.path                == "relax/kgrid_111")
    assert(sim_obj.process_id          == 1)
    assert(sim_obj.remdir              == "./runs/relax/kgrid_111")
    assert(sim_obj.resdir              == "./runs/relax/kgrid_111")
    assert(sim_obj.sent_files          is True)
    assert(sim_obj.setup               is True)
    assert(sim_obj.sim_image           == "sim.p")
    assert(sim_obj.subcascade_finished is False)
    assert(sim_obj.submitted           is True)


    ref_atomic_positions_atoms = [
        "Ge", "Ge", "Ge", "Ge", "Ge", "Ge", "Ge", "Ge", "Ge",
        "Ge", "Ge", "Ge", "Ge", "Ge", "Ge", "Ge", "Ge",
        ]

    ref_atomic_positions_positions = np.array([
        [ 5.34792496,  5.34792496,  5.34792496],
        [ 0.        ,  0.        ,  0.        ],
        [ 2.67396248,  2.67396248,  2.67396248],
        [ 5.34792496,  5.34792496,  0.        ],
        [ 8.02188743,  8.02188743,  2.67396248],
        [ 0.        ,  5.34792496,  5.34792496],
        [ 2.67396248,  8.02188743,  8.02188743],
        [ 5.34792496, 10.69584991,  5.34792496],
        [ 8.02188743, 13.36981239,  8.02188743],
        [ 5.34792496,  0.        ,  5.34792496],
        [ 8.02188743,  2.67396248,  8.02188743],
        [10.69584991,  5.34792496,  5.34792496],
        [13.36981239,  8.02188743,  8.02188743],
        [ 5.34792496,  5.34792496, 10.69584991],
        [ 8.02188743,  8.02188743, 13.36981239],
        [10.69584991, 10.69584991, 10.69584991],
        [13.36981239, 13.36981239, 13.36981239],
        ], dtype=np.float64)

    ref_cell_parameters_vectors = np.array([
        [10.69584991, 10.69584991,  0.        ],
        [ 0.        , 10.69584991, 10.69584991],
        [10.69584991,  0.        , 10.69584991],
        ], dtype=np.float64)

    inp_obj = NexusCore()
    if np.lib.NumpyVersion(np.__version__) >= '2.0.0b1':
        inp_obj.load(TEST_FILES["old_nxs_pwscf_input.p"])
    else:
        inp_obj.load(TEST_FILES["old_nxs_pwscf_input_numpy_1.p"])

    assert(inp_obj.atomic_positions.atoms == ref_atomic_positions_atoms)
    np.testing.assert_allclose(inp_obj.atomic_positions.positions, ref_atomic_positions_positions, rtol=1e-7)
    assert(inp_obj.atomic_positions.specifier == "bohr")

    assert(inp_obj.atomic_species.atoms               == ["Ge"])
    assert(inp_obj.atomic_species.specifier           == "")
    np.testing.assert_allclose(inp_obj.atomic_species.masses.Ge, 72.61, rtol=1e-7)
    assert(inp_obj.atomic_species.pseudopotentials.Ge == "Ge.pbe-kjpaw.UPF")

    assert(inp_obj.cell_parameters.specifier == "bohr")
    np.testing.assert_allclose(inp_obj.cell_parameters.vectors, ref_cell_parameters_vectors, rtol=1e-7)

    assert(inp_obj.control.calculation  == "relax")
    assert(inp_obj.control.disk_io      == "low")
    assert(inp_obj.control.outdir       == "pwscf_output")
    assert(inp_obj.control.prefix       == "pwscf")
    assert(inp_obj.control.pseudo_dir   == "./")
    assert(inp_obj.control.restart_mode == "from_scratch")
    assert(inp_obj.control.verbosity    == "high")
    assert(inp_obj.control.wf_collect   == False)

    np.testing.assert_allclose(inp_obj.electrons.conv_thr, 1e-06, rtol=1e-7)
    assert(inp_obj.electrons.diagonalization  == "david")
    assert(inp_obj.electrons.electron_maxstep == 1000)
    np.testing.assert_allclose(inp_obj.electrons.mixing_beta, 0.7, rtol=1e-7)
    assert(inp_obj.electrons.mixing_mode      == "plain")

    assert(inp_obj.ions.ion_dynamics      == "bfgs")
    assert(inp_obj.ions.pot_extrapolation == "second_order")
    assert(inp_obj.ions.upscale           == 100)
    assert(inp_obj.ions.wfc_extrapolation == "second_order")

    assert(inp_obj.k_points.grid      == (1, 1, 1))
    assert(inp_obj.k_points.shift     == (1, 1, 1))
    assert(inp_obj.k_points.specifier == "automatic")

    np.testing.assert_allclose(inp_obj.system.degauss, 0.0001, rtol=1e-7)
    assert(inp_obj.system.ecutrho     == 200)
    assert(inp_obj.system.ecutwfc     == 50)
    assert(inp_obj.system.ibrav       == 0)
    assert(inp_obj.system.input_dft   == "pbe")
    assert(inp_obj.system.nat         == 17)
    assert(inp_obj.system.nosym       == True)
    assert(inp_obj.system.ntyp        == 1)
    assert(inp_obj.system.occupations == "smearing")
    assert(inp_obj.system.smearing    == "fermi-dirac")
    assert(inp_obj.system.tot_charge  == 0)
#end def test_old_nexus_unpickle
