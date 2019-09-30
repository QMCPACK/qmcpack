
import testing
from testing import value_eq,object_eq
from testing import TestFailed,failed


def get_test_simulation_class():

    from simulation import Simulation,SimulationInput,SimulationAnalyzer

    class TestSimulationInput(SimulationInput):
        def is_valid(self):
            return True
        #end def is_valid

        def read(self,filepath):
            None
        #end def read

        def write(self,filepath=None):
            None
        #end def write

        def read_text(self,text,filepath=None):
            None
        #end def read_text

        def write_text(self,filepath=None):
            None
        #end def write_text

        def incorporate_system(self,system):
            None
        #end def incorporate_system

        def return_system(self):
            self.not_implemented()
        #end def return_system
    #end class TestSimulationInput


    class TestSimulationAnalyzer(SimulationAnalyzer):
        def __init__(self,sim):
            None
        #end def __init__

        def analyze(self):
            None
        #end def analyze
    #end class TestSimulationAnalyzer


    class TestSimulation(Simulation):
        input_type    = TestSimulationInput
        analyzer_type = TestSimulationAnalyzer

        def check_sim_status(self):
            self.finished = True
        #end def check_sim_status

        def get_output_files(self):
            return []
        #end def get_output_files
    #end class TestSimulation

    return TestSimulation
#end def get_test_simulation_class





def test_import():
    import simulation
    from simulation import Simulation,SimulationInput,SimulationAnalyzer
    from simulation import SimulationImage
    from simulation import NullSimulationInput,NullSimulationAnalyzer
    from simulation import GenericSimulation
    from simulation import SimulationInputTemplate
    from simulation import SimulationInputMultiTemplate
    from simulation import input_template,multi_input_template
    from simulation import generate_simulation

    print('\nend test')
#end def test_import



def test_simulation_input():
    import os
    from simulation import SimulationInput

    tpath = testing.setup_unit_test_output_directory('simulation','test_simulation_input')

    # empty init
    si = SimulationInput()

    # write
    infile = os.path.join(tpath,'sim_input.in')
    wtext = 'simulation input'
    si.write_file_text(infile,wtext)
    assert(os.path.exists(infile))

    # read
    rtext = si.read_file_text(infile)
    assert(rtext==wtext)

    # virtuals
    cls = SimulationInput
    virts = [
        cls.is_valid,
        cls.return_structure,
        cls.read_text,
        cls.write_text,
        (cls.incorporate_system,[None]),
        cls.return_system,
        ]
    for v in virts:
        args = []
        if isinstance(v,tuple):
            v,args = v
        #end if
        try:
            v(*args)
            raise TestFailed
        except TestFailed:
            failed(str(v))
        except:
            None
        #end try
    #end for
#end def test_simulation_input



def test_simulation_analyzer():
    import os
    from simulation import SimulationAnalyzer

    # empty init
    try:
        SimulationAnalyzer()
        raise TestFailed
    except TestFailed:
        failed()
    except:
        None
    #end try

    # virtuals
    try:
        SimulationAnalyzer(None)
        raise TestFailed
    except TestFailed:
        failed()
    except:
        None
    #end try
#end def test_simulation_analyzer



def test_simulation_code_name():
    from simulation import Simulation

    cn = Simulation.code_name()
    assert(isinstance(cn,str))
    assert(' ' not in cn)
#end def test_simulation_code_name



def test_simulation_init():
    from generic import obj
    from machines import job,Job
    from simulation import Simulation,SimulationInput

    # empty init, tests set(), set_directories(), set_files()
    se = Simulation()

    se_ref = obj(
        analyzed             = False,
        analyzer_image       = 'analyzer.p',
        app_name             = 'simapp',
        app_props            = ['serial'],
        block                = False,
        block_subcascade     = False,
        bundleable           = True,
        bundled              = False,
        bundler              = None,
        created_directories  = False,
        dependency_ids       = set([]),
        errfile              = 'sim.err',
        failed               = False,
        fake_sim             = False,
        files                = set([]),
        finished             = False,
        force_restart        = False,
        force_write          = False,
        got_dependencies     = False,
        got_output           = False,
        identifier           = 'sim',
        image_dir            = 'sim_sim',
        imlocdir             = './runs/sim_sim',
        imremdir             = './runs/sim_sim',
        imresdir             = './results/runs/sim_sim',
        infile               = 'sim.in',
        input_image          = 'input.p',
        locdir               = './runs/',
        job                  = None,
        loaded               = False,
        ordered_dependencies = [],
        outfile              = 'sim.out',
        outputs              = None,
        path                 = '',
        process_id           = None,
        restartable          = False,
        remdir               = './runs/',
        resdir               = './results/runs/',
        sent_files           = False,
        setup                = False,
        sim_image            = 'sim.p',
        simlabel             = None,
        skip_submit          = False,
        subcascade_finished  = False,
        submitted            = False,
        system               = None,
        wait_ids             = set([]),
        dependencies         = obj(),
        dependents           = obj(),
        input                = SimulationInput(),
        )

    assert(object_eq(se.obj(se_ref.keys()),se_ref))
    assert(isinstance(se.simid,int))
    assert(se.simid>=0)
    assert(se.simid<Simulation.sim_count)

    Simulation.clear_all_sims()


    # make a test job
    test_job = job(machine='ws1',app_command='test.x')

    
    # minimal non-empty init, tests init_job()
    sm = Simulation(job=test_job)

    sm_ref = se_ref.copy()
    del sm_ref.job
    assert(object_eq(sm.obj(sm_ref.keys()),sm_ref))
    assert(isinstance(se.simid,int))
    assert(se.simid>=0)
    assert(se.simid<Simulation.sim_count)
    assert(isinstance(sm.job,Job))
    assert(id(sm.job)!=id(test_job))


    # initialization tests for set_directories()
    # two sims in same directory w/ same identifier should fail
    try:
        s1 = Simulation(
            identifier = 'same_identifier',
            path       = 'same_directory',
            job        = test_job,
            )
        s2 = Simulation(
            identifier = 'same_identifier',
            path       = 'same_directory',
            job        = test_job,
            )
        raise TestFailed
    except TestFailed:
        failed()
    except:
        None
    #end try

    # two sims in same directory w/ different identifiers should be ok
    s1 = Simulation(
        identifier = 'identifier1',
        path       = 'same_directory',
        job        = test_job,
        )
    s2 = Simulation(
        identifier = 'identifier2',
        path       = 'same_directory',
        job        = test_job,
        )

    # two sims in different directories w/ same identifier should be ok
    s1 = Simulation(
        identifier = 'same_identifier',
        path       = 'directory1',
        job        = test_job,
        )
    s2 = Simulation(
        identifier = 'same_identifier',
        path       = 'directory2',
        job        = test_job,
        )


#end def test_simulation_init
