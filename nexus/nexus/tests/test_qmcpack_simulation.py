import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.QMCPACK_SIMULATION)

from ..generic import generic_settings
generic_settings.raise_error = True

from pathlib import Path

from . import isolate_nexus_core, create_pseudo_files, TEST_DIR

from ..testing import clear_all_sims
from ..testing import failed,FailedTest
from ..testing import value_eq,text_eq



def get_system(tiling=(1,1,1)):
    from ..physical_system import generate_physical_system

    system = generate_physical_system(
        units  = 'A',
        axes   = [[ 1.785,  1.785,  0.   ],
                  [ 0.   ,  1.785,  1.785],
                  [ 1.785,  0.   ,  1.785]],
        elem   = ['C','C'],
        pos    = [[ 0.    ,  0.    ,  0.    ],
                  [ 0.8925,  0.8925,  0.8925]],
        tiling = tiling,
        kgrid  = (1,1,1),
        kshift = (0,0,0),
        #C      = 4
        )

    return system
#end def get_system



def get_qmcpack_sim(type='rsqmc',**kwargs):
    from ..machines import job
    from ..qmcpack import Qmcpack,generate_qmcpack

    if type=='rsqmc':
        tiling = kwargs.pop('tiling',(1,1,1))

        sim = generate_qmcpack(
            job    = job(machine='ws1',cores=1),
            system = get_system(tiling=tiling),
            **kwargs
            )
    elif type=='afqmc':
        sim = generate_qmcpack(
            job        = job(machine='ws1',cores=1),
            input_type = 'basic_afqmc',
            **kwargs
            )
    else:
        failed()
    #end if

    assert(isinstance(sim,Qmcpack))

    return sim
#end def get_qmcpack_sim



def test_minimal_init():
    from ..machines import job
    from ..qmcpack import Qmcpack,generate_qmcpack

    sim = generate_qmcpack(
        job    = job(machine='ws1',cores=1),
        system = get_system(),
        )

    assert(isinstance(sim,Qmcpack))

    sim = generate_qmcpack(
        identifier = 'afqmc',
        job        = job(machine='ws1',cores=1),
        input_type = 'basic_afqmc',
        )

    assert(isinstance(sim,Qmcpack))

    clear_all_sims()
#end def test_minimal_init



def test_check_result():
    sim = get_qmcpack_sim()

    assert(not sim.check_result('unknown',None))
    assert(not sim.check_result('jastrow',None))
    assert(not sim.check_result('wavefunction',None))
    assert(not sim.check_result('cuspcorr',None))

    ds = sim.input.get('determinantset')
    ds.cuspcorrection = True

    assert(not sim.check_result('jastrow',None))
    assert(not sim.check_result('wavefunction',None))
    assert(sim.check_result('cuspcorr',None))

    opt = get_qmcpack_sim(identifier='opt',qmc='opt')

    assert(opt.check_result('jastrow',None))
    assert(opt.check_result('wavefunction',None))
    assert(not opt.check_result('cuspcorr',None))

    clear_all_sims()
#end def test_check_result


@isolate_nexus_core
def test_get_result(tmp_path):
    from ..developer import NexusError, obj
    from ..nexus_base import nexus_core
    from ..qmcpack_analyzer import QmcpackAnalyzer

    nexus_core.runs    = ''
    nexus_core.results = ''
    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    sim = get_qmcpack_sim()

    assert(Path(sim.locdir).resolve() == tmp_path)

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

    result = sim.get_result('cuspcorr',None)

    result_ref = obj(
        updet_cusps  = 'updet.cuspInfo.xml',
        spo_dn_cusps = 'qmcpack.spo-dn.cuspInfo.xml',
        spo_up_cusps = 'qmcpack.spo-up.cuspInfo.xml',
        dndet_cusps  = 'downdet.cuspInfo.xml',
        )

    assert(set(result.keys())==set(result_ref.keys()))
    for k,v in result.items():
        assert(Path(v).relative_to(tmp_path)==Path(result_ref[k]))
    #end for

    opt_infile = TEST_DIR / "test_qmcpack_analyzer_files/diamond_gamma/opt/opt.in.xml"
    assert(opt_infile.exists())

    qa = QmcpackAnalyzer(opt_infile,analyze=True,equilibration=5)

    assert(qa.results.optimization.optimal_file=='opt.s003.opt.xml')

    sim.create_directories()

    qa.save(Path(sim.imresdir).resolve() / sim.analyzer_image)

    result = sim.get_result('jastrow',None)

    result_ref = obj(
        opt_file = 'opt.s003.opt.xml',
        )

    assert(set(result.keys())==set(result_ref.keys()))
    for k,v in result.items():
        assert(Path(v).relative_to(tmp_path)==Path(result_ref[k]))
    #end for

    clear_all_sims()
#end def test_get_result


@isolate_nexus_core
def test_incorporate_result(tmp_path):
    import shutil
    from numpy import array
    from ..developer import obj
    from ..nexus_base import nexus_core
    from .test_vasp_simulation import setup_vasp_sim as get_vasp_sim
    from .test_qmcpack_converter_simulations import get_pw2qmcpack_sim
    from .test_qmcpack_converter_simulations import get_convert4qmc_sim
    from .test_qmcpack_converter_simulations import get_pyscf_to_afqmc_sim

    nexus_core.runs    = ''
    nexus_core.results = ''
    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    create_pseudo_files(
        tmp_dir=tmp_path,
        pseudos=["C.POTCAR"],
        pseudo_strs=[(
            "This is not a real POTCAR file.\n"
            "\n"
            "End of Dataset\n"
        )]
    )

    # incorporate vasp structure
    sim = get_qmcpack_sim(identifier='qmc_vasp_structure',tiling=(2,2,2))

    vasp_struct = get_vasp_sim(tmp_path,identifier='vasp_structure',copy_files=True)

    assert(Path(sim.locdir).resolve()==tmp_path)

    pc_file = 'diamond_POSCAR'
    cc_file = vasp_struct.identifier+'.CONTCAR'
    shutil.copy2(
        src = tmp_path / pc_file,
        dst = tmp_path / cc_file,
    )

    result = vasp_struct.get_result('structure',None)

    ion0 = sim.input.get('ion0')
    c = ion0.groups.C
    cp0 = c.position[0].copy()
    s = result.structure.copy()
    s.change_units('B')
    rp0 = s.pos[0]
    zero = array([0,0,0],dtype=float)
    assert(value_eq(c.position[0],cp0,atol=1e-8))
    c.position[0] += 0.1
    assert(not value_eq(c.position[0],cp0,atol=1e-8))
    assert(not value_eq(c.position[0],rp0,atol=1e-8))

    sim.incorporate_result('structure',result,vasp_struct)

    ion0 = sim.input.get('ion0')
    c = ion0.groups.C
    assert(value_eq(c.position[0],rp0,atol=1e-8))


    # incorporate pw2qmcpack orbitals
    sim = get_qmcpack_sim(identifier='qmc_p2q_orbitals')
    
    p2q_orb = get_pw2qmcpack_sim(identifier='p2q_orbitals')

    result = p2q_orb.get_result('orbitals',None)

    p2q_output_path = tmp_path / 'pwscf_output'
    p2q_output_path.mkdir()

    p2q_h5file = p2q_output_path / 'pwscf.pwscf.h5'
    p2q_h5file.touch()
    assert(p2q_h5file.exists())

    spo = sim.input.get('bspline')
    assert(spo.href=='MISSING.h5')

    sim.incorporate_result('orbitals',result,p2q_orb)

    assert(spo.href=='pwscf_output/pwscf.pwscf.h5')


    # incorporate convert4qmc orbitals
    sim = get_qmcpack_sim(identifier='qmc_c4q_orbitals')
    
    c4q_orb = get_convert4qmc_sim(identifier='c4q_orbitals')

    result = c4q_orb.get_result('orbitals',None)

    wfn_file  = tmp_path / 'c4q_orbitals.wfj.xml'
    wfn_file2 = tmp_path / 'c4q_orbitals.orbs.h5'
    input = sim.input.copy()
    dset = input.get('determinantset')
    dset.href = 'orbs.h5'
    qs = input.simulation.qmcsystem
    del input.simulation
    input.qmcsystem = qs
    input.write(wfn_file)
    assert(wfn_file.exists())
    wfn_file2.write_text('fake')
    assert(wfn_file2.exists())

    from ..qmcpack_input import QmcpackInput
    inp = QmcpackInput(wfn_file)

    dset = sim.input.get('determinantset')
    assert('href' not in dset)

    sim.incorporate_result('orbitals',result,c4q_orb)

    dset = sim.input.get('determinantset')
    assert(dset.href=='c4q_orbitals.orbs.h5')


    # incorporate qmcpack jastrow
    sim = get_qmcpack_sim(identifier='qmc_jastrow')

    opt_path = TEST_DIR / "test_qmcpack_analyzer_files/diamond_gamma/opt"
    opt_infile = opt_path / 'opt.in.xml'
    assert(opt_infile.exists())
    opt_file = opt_path / 'opt.s004.opt.xml'
    assert(opt_file.exists())

    result = obj(opt_file=opt_file)

    sim.incorporate_result('jastrow',result,sim)

    j = sim.input.get('jastrow')

    j_text = j.write().replace('"',' " ')

    j_text_ref = '''
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
        '''.replace('"',' " ')
    
    assert(text_eq(j_text,j_text_ref))


    # incorporate qmcpack wavefunction
    sim = get_qmcpack_sim(identifier='qmc_wavefunction')

    sim.incorporate_result('wavefunction',result,sim)

    j = sim.input.get('jastrow')

    j_text = j.write().replace('"',' " ')

    assert(text_eq(j_text,j_text_ref))


    # incorporate qmcpack wavefunction afqmc
    sim = get_qmcpack_sim(type='afqmc',identifier='qmc_wf_afqmc')

    p2a_wf = get_pyscf_to_afqmc_sim(identifier='p2a_wavefunction')

    result = p2a_wf.get_result('wavefunction',None)

    del result.xml

    wfn = sim.input.get('wavefunction')
    ham = sim.input.get('hamiltonian')

    assert(wfn.filename=='MISSING.h5')
    assert(ham.filename=='MISSING.h5')

    sim.incorporate_result('wavefunction',result,p2a_wf)

    assert(wfn.filename=='p2a_wavefunction.afqmc.h5')
    assert(ham.filename=='p2a_wavefunction.afqmc.h5')
    
    clear_all_sims()
#end def test_incorporate_result()


@isolate_nexus_core
def test_check_sim_status(tmp_path):
    from ..nexus_base import nexus_core

    nexus_core.runs = ''
    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    sim = get_qmcpack_sim(identifier='qmc')

    assert(Path(sim.locdir).resolve()==tmp_path)

    assert(not sim.finished)
    assert(not sim.failed)

    try:
        sim.check_sim_status()
    except IOError:
        None
    #end try

    assert(not sim.finished)
    assert(not sim.failed)

    outfile = tmp_path / sim.outfile
    out_text = 'Total Execution'
    outfile.write_text(out_text)
    assert(outfile.exists())
    assert(out_text in outfile.read_text())
    errfile = tmp_path / sim.errfile
    errfile.touch()
    assert(errfile.exists())

    sim.check_sim_status()

    assert(sim.finished)
    assert(not sim.failed)

    clear_all_sims()
#end def test_check_sim_status()

