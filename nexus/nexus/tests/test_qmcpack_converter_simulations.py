import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.QMCPACK_CONVERTER_SIMULATIONS)

from ..generic import generic_settings
generic_settings.raise_error = True

from pathlib import Path

from . import isolate_nexus_core, create_pseudo_files
from nexus.nexus_base import nexus_core
from ..testing import clear_all_sims
from ..testing import failed,FailedTest
from ..testing import object_eq


#====================================================================#
#                         pw2qmcpack tests                           #
#====================================================================#

def get_pw2qmcpack_sim(**kwargs):
    from ..machines import job
    from ..qmcpack_converters import Pw2qmcpack,generate_pw2qmcpack
    
    sim = generate_pw2qmcpack(
        job = job(machine='ws1',cores=1),
        **kwargs
        )

    assert(isinstance(sim,Pw2qmcpack))

    return sim
#end def get_pw2qmcpack_sim



def test_pw2qmcpack_minimal_init():
    from ..machines import job
    from ..qmcpack_converters import Pw2qmcpack,generate_pw2qmcpack
    
    sim = generate_pw2qmcpack(
        job = job(machine='ws1',cores=1),
        )

    assert(isinstance(sim,Pw2qmcpack))

    clear_all_sims()
#end def test_pw2qmcpack_minimal_init



def test_pw2qmcpack_check_result():
    sim = get_pw2qmcpack_sim()
    
    assert(not sim.check_result('unknown',None))
    assert(sim.check_result('orbitals',None))

    clear_all_sims()
#end def test_pw2qmcpack_check_result



def test_pw2qmcpack_get_result():
    from ..developer import NexusError, obj

    sim = get_pw2qmcpack_sim()
    
    try:
        sim.get_result('unknown',None)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    result = sim.get_result('orbitals',None)

    result_ref = obj(
        h5file   = './runs/pwscf_output/pwscf.pwscf.h5',
        ptcl_xml = './runs/pwscf_output/pwscf.ptcl.xml',
        wfs_xml  = './runs/pwscf_output/pwscf.wfs.xml',
        )

    assert(object_eq(result,result_ref))

    clear_all_sims()
#end def test_pw2qmcpack_get_result


@isolate_nexus_core(needs_tmp_path=True)
def test_pw2qmcpack_incorporate_result(tmp_path):
    from ..developer import NexusError
    from ..simulation import Simulation
    from .test_pwscf_simulation import get_pwscf_sim

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = tmp_path
    nexus_core.file_locations = nexus_core.file_locations + [tmp_path]

    create_pseudo_files(tmp_path, ["C.BFD.upf"])

    other = Simulation()

    scf = get_pwscf_sim('scf')

    sim = get_pw2qmcpack_sim(path='scf')

    try:
        sim.incorporate_result('unknown',None,other)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    sim.incorporate_result('orbitals',None,scf)

    clear_all_sims()
#end def test_pw2qmcpack_incorporate_result


@isolate_nexus_core(needs_tmp_path=True)
def test_pw2qmcpack_check_sim_status(tmp_path):
    from ..nexus_base import nexus_core

    nexus_core.runs = ''
    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = tmp_path
    nexus_core.file_locations = nexus_core.file_locations + [tmp_path]

    sim = get_pw2qmcpack_sim()

    assert(Path(sim.locdir).resolve() == tmp_path)

    assert(not sim.finished)
    assert(not sim.failed)

    try:
        sim.check_sim_status()
        raise FailedTest
    except IOError:
        None
    except Exception as e:
        failed(str(e))
    #end try

    sim.create_directories()
    outfile = Path(sim.locdir) / sim.outfile
    outfile_text = 'JOB DONE'
    outfile.write_text(outfile_text)
    assert(outfile_text in outfile.read_text())
    prefix = sim.input.inputpp.prefix
    outdir = sim.input.inputpp.outdir
    outdir_path = tmp_path / outdir
    outdir_path.mkdir(parents=True)

    filenames = [prefix+'.pwscf.h5',prefix+'.ptcl.xml',prefix+'.wfs.xml']
    for filename in filenames:
        filepath = outdir_path / filename
        filepath.touch()
        assert(filepath.exists())
    #end for
    sim.job.finished = True
    
    sim.check_sim_status()
    
    assert(sim.finished)
    assert(not sim.failed)

    clear_all_sims()
#end def test_pw2qmcpack_check_sim_status


#====================================================================#
#                        convert4qmc tests                           #
#====================================================================#

def get_convert4qmc_sim(**kwargs):
    from ..machines import job
    from ..qmcpack_converters import Convert4qmc,generate_convert4qmc
    
    sim = generate_convert4qmc(
        job = job(machine='ws1',cores=1),
        **kwargs
        )

    assert(isinstance(sim,Convert4qmc))

    return sim
#end def get_convert4qmc_sim



def test_convert4qmc_minimal_init():
    from ..machines import job
    from ..qmcpack_converters import Convert4qmc,generate_convert4qmc
    
    sim = generate_convert4qmc(
        job = job(machine='ws1',cores=1),
        )

    assert(isinstance(sim,Convert4qmc))

    clear_all_sims()
#end def test_convert4qmc_minimal_init



def test_convert4qmc_check_result():
    sim = get_convert4qmc_sim()
    
    assert(not sim.check_result('unknown',None))
    assert(sim.check_result('orbitals',None))
    assert(sim.check_result('particles',None))

    clear_all_sims()
#end def test_convert4qmc_check_result



def test_convert4qmc_get_result():
    from ..developer import NexusError, obj

    sim = get_convert4qmc_sim()
    
    try:
        sim.get_result('unknown',None)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    result = sim.get_result('orbitals',None)

    result_ref = obj(
        location = './runs/sample.wfj.xml',
        orbfile  = './runs/sample.orbs.h5',
        )

    assert(object_eq(result,result_ref))

    result = sim.get_result('particles',None)

    result_ref = obj(
        location = './runs/sample.structure.xml',
        )

    assert(object_eq(result,result_ref))

    clear_all_sims()
#end def test_convert4qmc_get_result



def test_convert4qmc_incorporate_result():
    from ..developer import NexusError, obj
    from ..simulation import Simulation
    from ..gamess import Gamess
    from ..quantum_package import QuantumPackage
    from .test_gamess_simulation import get_gamess_sim
    from .test_pyscf_simulation import get_pyscf_sim
    from .test_quantum_package_simulation import get_quantum_package_sim

    other = Simulation()

    gms = get_gamess_sim('rhf')
    Gamess.ericfmt = None
    gms_result = obj(
        location        = './rhf/rhf.out',
        mos             = 0,
        norbitals       = 0,
        outfile         = './rhf/rhf.out',
        scftyp          = 'rohf',
        vec             = None,
        )

    pscf = get_pyscf_sim()
    pscf_result = obj(
        h5_file = './scf.h5',
        )

    qp = get_quantum_package_sim()
    QuantumPackage.qprc = None
    qp_result = obj(
        outfile = './qp_savewf.out',
        )

    sim_start = get_convert4qmc_sim()

    sim = sim_start.copy()
    try:
        sim.incorporate_result('unknown',None,other)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # incorporate orbitals from gamess
    sim = sim_start.copy()

    assert(sim.input_code is None)
    assert(sim.input.gamess_ascii is None)
    assert(sim.job.app_command=='convert4qmc')

    sim.incorporate_result('orbitals',gms_result,gms)

    assert(sim.input_code=='gamess')
    assert(sim.input.gamess_ascii=='../rhf/rhf.out')
    assert(sim.job.app_command=='convert4qmc -gamess ../rhf/rhf.out')

    # incorporate orbitals from pyscf
    sim = sim_start.copy()

    assert(sim.input_code is None)
    assert(sim.input.pyscf is None)

    sim.incorporate_result('orbitals',pscf_result,pscf)

    assert(sim.input_code=='pyscf')
    assert(sim.input.orbitals=='../scf.h5')
    
    # incorporate orbitals from quantum package
    sim = sim_start.copy()

    assert(sim.input_code is None)
    assert(sim.input.qp is None)

    sim.incorporate_result('orbitals',qp_result,qp)

    assert(sim.input_code=='qp')
    #assert(sim.input.qp=='../qp_savewf.out')
    assert(sim.input.orbitals=='../qp_savewf.out')

    clear_all_sims()
#end def test_convert4qmc_incorporate_result


@isolate_nexus_core(needs_tmp_path=True)
def test_convert4qmc_check_sim_status(tmp_path):
    from ..nexus_base import nexus_core

    nexus_core.runs = ''
    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = tmp_path
    nexus_core.file_locations = nexus_core.file_locations + [tmp_path]

    sim = get_convert4qmc_sim()

    assert(Path(sim.locdir).resolve() == tmp_path)

    assert(not sim.finished)
    assert(not sim.failed)

    try:
        sim.check_sim_status()
        raise FailedTest
    except IOError:
        None
    except Exception as e:
        failed(str(e))
    #end try

    sim.create_directories()

    outfile = Path(sim.locdir).resolve() / sim.outfile
    outfile_text = 'QMCGaussianParserBase::dump'
    outfile.write_text(outfile_text)

    assert(outfile_text in open(outfile,'r').read())
    for filename in sim.list_output_files():
        filepath = Path(sim.locdir).resolve() / filename
        filepath.touch()
        assert(filepath.exists())
    #end for
    sim.job.finished = True

    sim.check_sim_status()

    assert(sim.finished)
    assert(not sim.failed)

    clear_all_sims()
#end def test_convert4qmc_check_sim_status



#====================================================================#
#                       pyscf_to_afqmc tests                         #
#====================================================================#

def get_pyscf_to_afqmc_sim(**kwargs):
    from ..machines import job
    from ..qmcpack_converters import PyscfToAfqmc,generate_pyscf_to_afqmc
    
    sim = generate_pyscf_to_afqmc(
        job = job(machine='ws1',cores=1),
        **kwargs
        )

    assert(isinstance(sim,PyscfToAfqmc))

    return sim
#end def get_pyscf_to_afqmc_sim



def test_pyscf_to_afqmc_minimal_init():
    from ..machines import job
    from ..qmcpack_converters import PyscfToAfqmc,generate_pyscf_to_afqmc
    
    sim = generate_pyscf_to_afqmc(
        job = job(machine='ws1',cores=1),
        )

    assert(isinstance(sim,PyscfToAfqmc))

    clear_all_sims()
#end def test_pyscf_to_afqmc_minimal_init



def test_pyscf_to_afqmc_check_result():
    sim = get_pyscf_to_afqmc_sim()

    assert(not sim.check_result('unknown',None))
    assert(not sim.check_result('wavefunction',None))
    assert(not sim.check_result('hamiltonian',None))

    sim.input.output = 'afqmc.h5'

    assert(sim.check_result('wavefunction',None))
    assert(sim.check_result('hamiltonian',None))

    clear_all_sims()
#end def test_pyscf_to_afqmc_check_result



def test_pyscf_to_afqmc_get_result():
    from ..developer import NexusError, obj

    sim = get_pyscf_to_afqmc_sim()
    
    sim.input.output = 'afqmc.h5'

    try:
        sim.get_result('unknown',None)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    result_ref = obj(
        h5_file = './runs/afqmc.h5',
        )

    result = sim.get_result('wavefunction',None)

    assert(object_eq(result,result_ref))

    result = sim.get_result('hamiltonian',None)

    assert(object_eq(result,result_ref))

    clear_all_sims()
#end def test_pyscf_to_afqmc_get_result



def test_pyscf_to_afqmc_incorporate_result():
    import os
    from ..developer import NexusError, obj
    from ..simulation import Simulation
    from .test_pyscf_simulation import get_pyscf_sim

    other = Simulation()

    scf = get_pyscf_sim()

    sim = get_pyscf_to_afqmc_sim()

    try:
        sim.incorporate_result('unknown',None,other)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    result = obj(
        chkfile = os.path.join(scf.locdir,'scf.chk'),
        )

    assert(sim.input.input==None)

    sim.incorporate_result('wavefunction',result,scf)

    assert(sim.input.input=='scf.chk')

    clear_all_sims()
#end def test_pyscf_to_afqmc_incorporate_result


@isolate_nexus_core(needs_tmp_path=True)
def test_pyscf_to_afqmc_check_sim_status(tmp_path):
    from ..nexus_base import nexus_core

    nexus_core.runs = ''
    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = tmp_path
    nexus_core.file_locations = nexus_core.file_locations + [tmp_path]

    sim = get_pyscf_to_afqmc_sim()

    assert(Path(sim.locdir).resolve()==tmp_path)

    assert(not sim.finished)
    assert(not sim.failed)

    try:
        sim.check_sim_status()
        raise FailedTest
    except IOError:
        None
    except Exception as e:
        failed(str(e))
    #end try

    sim.input.output = 'afqmc.h5'

    sim.create_directories()
    outfile = Path(sim.locdir).resolve() / sim.outfile
    outfile_text = '# Finished.'
    outfile.write_text(outfile_text)

    assert(outfile_text in outfile.read_text())
    output_file = Path(sim.locdir).resolve() / sim.input.output
    output_file.touch()
    assert(output_file.exists())
    sim.job.finished = True

    sim.check_sim_status()

    assert(sim.finished)
    assert(not sim.failed)

    clear_all_sims()
#end def test_pyscf_to_afqmc_check_sim_status
