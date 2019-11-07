
import testing
from testing import divert_nexus,restore_nexus
from testing import failed,FailedTest
from testing import value_eq,object_eq,text_eq


def clear_all_sims():
    from gamess import Gamess

    Gamess.ericfmt = None

    testing.clear_all_sims()
#end def clear_all_sims


def get_gamess_sim(type='rhf'):
    from machines import job
    from gamess import Gamess,generate_gamess,GamessInput
    from test_gamess_input import get_files

    sim = None

    Gamess.ericfmt = ''

    files = get_files()

    if type=='rhf':
        gi_input = GamessInput(files['rhf.inp'])

        sim = generate_gamess(
            identifier = 'rhf',
            path       = 'rhf',
            job        = job(machine='ws1',cores=1),
            input      = gi_input,
            )
    else:
        failed()
    #end if

    assert(sim is not None)
    assert(isinstance(sim,Gamess))

    return sim
#end def get_gamess_sim



def test_import():
    from gamess import Gamess,generate_gamess
#end def test_import



def test_minimal_init():
    from machines import job
    from gamess import Gamess,generate_gamess

    Gamess.ericfmt = ''

    sim = generate_gamess(
        job = job(machine='ws1',cores=1),
        )

    assert(isinstance(sim,Gamess))

    clear_all_sims()
#end def test_minimal_init



def test_check_result():
    tpath = testing.setup_unit_test_output_directory('gamess_simulation','test_check_result')

    sim = get_gamess_sim('rhf')
    
    assert(not sim.check_result('unknown',None))
    assert(sim.check_result('orbitals',None))

    clear_all_sims()
#end def test_check_result



def test_get_result():
    import os
    from generic import NexusError,obj
    from nexus_base import nexus_core

    tpath = testing.setup_unit_test_output_directory('gamess_simulation','test_get_result',divert=True)

    nexus_core.runs = ''

    sim = get_gamess_sim('rhf')

    assert(sim.locdir.rstrip('/')==os.path.join(tpath,'rhf').rstrip('/'))

    sim.create_directories()
    sim.write_inputs()
    if not os.path.exists(sim.imresdir):
        os.makedirs(sim.imresdir)
    #end if
    analyzer = sim.analyzer_type(sim)
    analyzer.save(os.path.join(sim.imresdir,sim.analyzer_image))

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
        location        = 'rhf/rhf.out',
        mos             = 0,
        norbitals       = 0,
        outfile         = 'rhf/rhf.out',
        scftyp          = 'rohf',
        vec             = None,
        )

    result.location = result.location.replace(tpath,'').lstrip('/')
    result.outfile  = result.outfile.replace(tpath,'').lstrip('/')

    assert(object_eq(result,result_ref))

    clear_all_sims()
    restore_nexus()
#end def test_get_result



def test_incorporate_result():
    import os
    from generic import NexusError,obj

    tpath = testing.setup_unit_test_output_directory('gamess_simulation','test_incorporate_result')

    sim = get_gamess_sim('rhf')

    try:
        sim.incorporate_result('unknown',None,None)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    result = obj(
        vec       = 'vec text',
        norbitals = 10,
        )

    input = sim.input

    assert('vec' not in input)
    assert('norb' not in input.guess)
    assert('prtmo' not in input.guess)

    sim.incorporate_result('orbitals',result,None)

    assert(input.vec.text=='vec text')
    assert(input.guess.guess=='moread')
    assert(input.guess.norb==10)
    assert(input.guess.prtmo==True)

    clear_all_sims()
#end def test_incorporate_result



def test_check_sim_status():
    import os
    from generic import NexusError,obj
    from nexus_base import nexus_core

    tpath = testing.setup_unit_test_output_directory('gamess_simulation','test_check_sim_status',divert=True)

    nexus_core.runs = ''

    sim = get_gamess_sim('rhf')

    assert(sim.locdir.rstrip('/')==os.path.join(tpath,'rhf').rstrip('/'))

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
    outfile = os.path.join(sim.locdir,sim.outfile)
    outfile_text = 'EXECUTION OF GAMESS TERMINATED NORMALLY'
    out = open(outfile,'w')
    out.write(outfile_text)
    out.close()
    assert(outfile_text in open(outfile,'r').read())

    sim.check_sim_status()

    assert(sim.finished)
    assert(not sim.failed)

    clear_all_sims()
    restore_nexus()
#end def test_check_sim_status
