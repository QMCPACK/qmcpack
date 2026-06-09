import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.SIMULATION)

from ..generic import generic_settings
generic_settings.raise_error = True


from pathlib import Path
from . import isolate_nexus_core
from nexus.nexus_base import nexus_core

from ..testing import value_eq,object_eq
from ..testing import FailedTest,failed

from ..developer import obj
from ..simulation import Simulation,SimulationInput,SimulationAnalyzer


class SimulationInputForTests(SimulationInput):
    def __init__(self,*args,**kwargs):
        SimulationInput.__init__(self,*args,**kwargs)

        self.result_data = obj()
    #end def __init__

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
#end class SimulationInputForTests


class SimulationAnalyzerForTests(SimulationAnalyzer):
    def __init__(self,sim):
        self.analysis_performed = False
    #end def __init__

    def analyze(self):
        self.analysis_performed = True
    #end def analyze
#end class SimulationAnalyzerForTests


class SimulationForTests(Simulation):
    input_type    = SimulationInputForTests
    analyzer_type = SimulationAnalyzerForTests

    application_results = set(['quant1','quant2','quant3'])

    def check_sim_status(self):
        self.finished = True
    #end def check_sim_status

    def get_output_files(self):
        return []
    #end def get_output_files

    def check_result(self,result_name,sim):
        return result_name in self.application_results
    #end def check_result

    def get_result(self,result_name,sim):
        result = obj()
        result.name  = result_name
        result.simid = sim.simid
        return result
    #end def get_result

    def incorporate_result(self,result_name,result,sim):
        self.input.result_data[result.simid] = result.name
    #end def incorporate_result
#end class SimulationForTests




get_sim_simulations = []

def get_sim(**kwargs):
    from ..machines import job
    from ..simulation import Simulation
    test_job = job(machine='ws1',app_command='test.x')

    n = len(get_sim_simulations)

    sim = Simulation(identifier='sim'+str(n),job=test_job,**kwargs)

    get_sim_simulations.append(sim)

    return sim
#end def get_sim


get_test_sim_simulations = []

def get_test_sim(**kwargs):
    from ..machines import job

    test_job = kwargs.pop('job',None)
    if test_job is None:
        test_job = job(machine='ws128',serial=True,app_command='echo run')
    #end if

    n = len(get_test_sim_simulations)

    test_sim = SimulationForTests(
        identifier = 'test_sim'+str(n),
        job        = test_job,
        **kwargs
        )

    get_test_sim_simulations.append(test_sim)

    return test_sim
#end def get_test_sim



n_test_workflows = 9

    
def generate_network():
    from numpy.random import randint

    nsims        = 10
    nconnections = 3
    nheads       = 3

    sims = []

    for n in range(randint(nheads)+1):
        sims.append([])
    #end for

    for isim in range(nsims):
        deps = []
        for idep in range(randint(nconnections)+1):
            if len(sims)>0:
                i = randint(len(sims))
                deps.append(i)
            #end if
        #end for
        sims.append(list(sorted(set(deps))))
    #end for

    sims_dict = {}
    for i,d in enumerate(sims):
        sims_dict[i] = d
    #end for

    return sims_dict
#end def generate_network


def get_test_workflow(index,**kwargs):
    from ..developer import obj

    def make_network(network,**kwargs):
        sims = obj()
        for i in range(len(network)):
            s = get_test_sim(**kwargs)
            for j in network[i]:
                s.depends((sims['s'+str(j)],'other'))
            #end for
            sims['s'+str(i)] = s
        #end for
        return sims
    #end def make_network

    nsims_bef = len(get_test_sim_simulations)

    sims = obj()

    if index==0:
        # single simulation
        sims.s = get_test_sim(**kwargs)
    elif index==1:
        # linear simulation chain
        sims.s1 = get_test_sim(**kwargs)
        sims.s2 = get_test_sim(dependencies=(sims.s1,'other'),**kwargs)
        sims.s3 = get_test_sim(dependencies=(sims.s2,'other'),**kwargs)
    elif index==2:
        # linear chain with split
        sims.s1  = get_test_sim(**kwargs)
        sims.s2  = get_test_sim(dependencies=(sims.s1,'other'),**kwargs)
        sims.s3  = get_test_sim(dependencies=(sims.s2,'other'),**kwargs)
        sims.s41 = get_test_sim(dependencies=(sims.s3,'other'),**kwargs)
        sims.s51 = get_test_sim(dependencies=(sims.s41,'other'),**kwargs)
        sims.s42 = get_test_sim(dependencies=(sims.s3,'other'),**kwargs)
        sims.s52 = get_test_sim(dependencies=(sims.s42,'other'),**kwargs)
    elif index==3:
        # linear chains with join
        sims.s11 = get_test_sim(**kwargs)
        sims.s21 = get_test_sim(dependencies=(sims.s11,'other'),**kwargs)
        sims.s12 = get_test_sim(**kwargs)
        sims.s22 = get_test_sim(dependencies=(sims.s12,'other'),**kwargs)
        sims.s3  = get_test_sim(dependencies=[(sims.s21,'other'),(sims.s22,'other')])
        sims.s4  = get_test_sim(dependencies=(sims.s3,'other'),**kwargs)
        sims.s5  = get_test_sim(dependencies=(sims.s4,'other'),**kwargs)
    elif index==4:
        # net-like workflow
        sims.s11 = get_test_sim(**kwargs)
        sims.s12 = get_test_sim(**kwargs)
        sims.s13 = get_test_sim(**kwargs)

        sims.s21 = get_test_sim(
            dependencies = [
                (sims.s11,'other'),
                (sims.s12,'other'),
                ],
            **kwargs
            )
        sims.s22 = get_test_sim(
            dependencies = [
                (sims.s12,'other'),
                (sims.s13,'other'),
                ],
            **kwargs
            )

        sims.s31 = get_test_sim(
            dependencies = [
                (sims.s21,'other'),
                (sims.s22,'other'),
                ],
            **kwargs
            )
        sims.s32 = get_test_sim(
            dependencies = [
                (sims.s21,'other'),
                (sims.s22,'other'),
                ],
            **kwargs
            )
        sims.s33 = get_test_sim(
            dependencies = [
                (sims.s21,'other'),
                (sims.s22,'other'),
                ],
            **kwargs
            )

        sims.s41 = get_test_sim(
            dependencies = [
                (sims.s11,'other'),
                (sims.s22,'other'),
                (sims.s32,'other'),
                ],
            **kwargs
            )
    elif index==5:
        # random network 1
        network = {
            0  : [],
            1  : [],
            2  : [1],
            3  : [2],
            4  : [1, 2],
            5  : [0, 4],
            6  : [4],
            7  : [1],
            8  : [7],
            9  : [5],
            10 : [0, 3, 9],
            11 : [4, 5, 7],
            }
        sims = make_network(network,**kwargs)
    elif index==6:
        # random network 2
        network = {
            0  : [],
            1  : [0],
            2  : [0],
            3  : [1, 2],
            4  : [1, 2, 3],
            5  : [0, 1, 2],
            6  : [1, 3],
            7  : [2],
            8  : [1, 3, 7],
            9  : [1, 8],
            10 : [2, 4],
            }
        sims = make_network(network,**kwargs)
    elif index==7:
        # random network 3
        network = {
            0  : [],
            1  : [],
            2  : [],
            3  : [0, 1, 2],
            4  : [2, 3],
            5  : [1],
            6  : [3, 5],
            7  : [0],
            8  : [3],
            9  : [6],
            10 : [0],
            11 : [10],
            12 : [8, 10, 11],
            }
        sims = make_network(network,**kwargs)
    elif index==8:
        # larger random network
        network = {
            0  : [], 
            1  : [], 
            2  : [], 
            3  : [0], 
            4  : [2], 
            5  : [0, 2], 
            6  : [2, 5], 
            7  : [1], 
            8  : [1, 5, 7], 
            9  : [0, 8], 
            10 : [5], 
            11 : [0, 10], 
            12 : [1], 
            13 : [6, 9, 11], 
            14 : [2, 4], 
            15 : [8, 9, 10, 12], 
            16 : [3], 
            17 : [7], 
            18 : [6, 13], 
            19 : [7, 11, 14], 
            20 : [8, 16, 17], 
            21 : [0, 8], 
            22 : [14], 
            23 : [2], 
            24 : [16], 
            25 : [13, 16, 22, 24], 
            26 : [6, 10], 
            27 : [12, 17, 24], 
            28 : [8], 
            29 : [10, 12, 23, 26], 
            30 : [1], 
            31 : [28], 
            32 : [20],
            }
        sims = make_network(network,**kwargs)
    else:
        failed('index exceeds available workflows')
    #end if

    nsims_aft = len(get_test_sim_simulations)

    assert(len(sims)==nsims_aft-nsims_bef)

    return sims
#end def get_test_workflow



def test_simulation_input(tmp_path):
    from ..developer import NexusError
    from ..simulation import SimulationInput

    # empty init
    si = SimulationInput()

    # write
    infile = tmp_path / 'sim_input.in'
    wtext = 'simulation input'
    si.write_file_text(infile,wtext)
    assert(infile.exists())

    # read
    rtext = si.read_file_text(infile)
    assert(rtext==wtext)

    # virtuals
    virts = [
        si.is_valid,
        si.return_structure,
        (si.read_text,[None]),
        si.write_text,
        (si.incorporate_system,[None]),
        si.return_system,
        ]
    for v in virts:
        args = []
        if isinstance(v,tuple):
            v,args = v
        #end if
        try:
            v(*args)
            raise FailedTest
        except NexusError:
            None
        except FailedTest:
            failed(str(v))
        except Exception as e:
            failed(str(e))
        #end try
    #end for
#end def test_simulation_input



def test_simulation_analyzer():
    from ..developer import NexusError
    from ..simulation import SimulationAnalyzer

    # empty init
    try:
        SimulationAnalyzer()
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try

    # virtuals
    try:
        SimulationAnalyzer(None)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try
#end def test_simulation_analyzer



def test_simulation_input_template(tmp_path):
    from string import Template
    from ..developer import obj, NexusError
    from ..simulation import SimulationInput
    from ..simulation import GenericSimulationInput
    from ..simulation import SimulationInputTemplate
    from ..simulation import input_template

    # empty init
    si_empty = input_template()
    assert(isinstance(si_empty,SimulationInput))
    assert(isinstance(si_empty,GenericSimulationInput))
    assert(isinstance(si_empty,SimulationInputTemplate))

    si_empty_ref = obj(
        template      = None,
        keywords      = set(),
        values        = obj(),
        allow_not_set = set(),
        )

    assert(len(si_empty)==4)
    assert(object_eq(si_empty.to_obj(),si_empty_ref))


    # template reference data
    template_text = '''
a     = "$a"
b     = $b
file1 = "$file.$ext1"
file2 = "$file.$ext2"
'''

    template_filepath = tmp_path / 'template_file.txt'

    template_filepath.write_text(template_text)

    # read
    si_read = input_template(template_filepath)

    assert(isinstance(si_read.template,Template))
    assert(si_read.keywords==set(['a','b','ext1','ext2','file']))


    # assign
    si_assign = input_template()
    try:
        si_assign.assign(b=1)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    si_assign = input_template(template_filepath)
    try:
        si_assign.assign(c=1)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    values = obj(
        a    = 'name',
        b    = 1,
        file = 'my_file',
        ext1 = 'txt',
        ext2 = 'dat',
        )

    si_assign.assign(**values)

    assert(object_eq(si_assign.values,values))


    # write
    def try_write(si):
        try:
            si.write()
            raise FailedTest
        except NexusError:
            None
        except FailedTest:
            failed()
        except Exception as e:
            failed(str(e))
        #end try
    #end def try_write

    si_write = input_template()
    try_write(si_write)

    si_write = input_template(template_filepath)
    try_write(si_write)

    si_write.assign(b=1)
    try_write(si_write)


    text_ref = '''
a     = "name"
b     = 1
file1 = "my_file.txt"
file2 = "my_file.dat"
'''

    si_write.assign(**values)
    text = si_write.write()
    assert(text==text_ref)
    
    input_filepath = tmp_path / 'input_file.txt'
    si_write.write(input_filepath)
    assert(input_filepath.read_text()==text_ref)
#end def test_simulation_input_template



def test_simulation_input_multi_template(tmp_path):
    from string import Template
    from ..developer import obj
    from ..simulation import SimulationInput
    from ..simulation import GenericSimulationInput
    from ..simulation import SimulationInputMultiTemplate
    from ..simulation import multi_input_template

    # make template files
    template1_filepath = tmp_path / 'template1.txt'
    template2_filepath = tmp_path / 'template2.txt'
    template3_filepath = tmp_path / 'template3.txt'

    template1_filepath.write_text('''
name = "$name"
a    = $a
''')
    template2_filepath.write_text('''
name = "$name"
b    = $b
''')
    template3_filepath.write_text('''
name = "$name"
c    = $c
''')

    input1_filepath = tmp_path / 'input_file1.txt'
    input2_filepath = tmp_path / 'input_file2.txt'
    input3_filepath = tmp_path / 'input_file3.txt'


    # empty init
    si_empty = multi_input_template()

    assert(isinstance(si_empty,SimulationInput))
    assert(isinstance(si_empty,GenericSimulationInput))
    assert(isinstance(si_empty,SimulationInputMultiTemplate))

    si_empty_ref = obj(
        filenames = obj(),
        )

    assert(len(si_empty)==1)
    assert(object_eq(si_empty.to_obj(),si_empty_ref))


    # filename init
    filenames = obj(
        input1 = 'input_file1.txt',
        input2 = 'input_file2.txt',
        input3 = 'input_file3.txt',
        )

    si = multi_input_template(**filenames)

    assert(len(si)==1)
    assert(len(si.filenames)==3)
    assert(object_eq(si.filenames,filenames))

    
    # init read
    si_init = multi_input_template(
        input1 = ('input_file1.txt',template1_filepath),
        input2 = ('input_file2.txt',template2_filepath),
        input3 = ('input_file3.txt',template3_filepath),
        )
    si = si_init
    assert(len(si)==4)
    assert(len(si.filenames)==3)
    assert(object_eq(si.filenames,filenames))
    keywords_ref = dict(
        input1 = set(['a', 'name']),
        input2 = set(['b', 'name']),
        input3 = set(['c', 'name']),
        )
    for name,keyword_set in keywords_ref.items():
        assert(name in si)
        sit = si[name]
        assert(sit.keywords==keyword_set)
        assert(isinstance(sit.template,Template))
        assert(object_eq(sit.values,obj()))
        assert(sit.allow_not_set==set())
    #end for


    # write
    write_ref = obj(
        input1 = '''
name = "name1"
a    = 1
''',
        input2 = '''
name = "name2"
b    = 2
''',
        input3 = '''
name = "name3"
c    = 3
''',
        )

    si_write = multi_input_template(
        input1 = ('input_file1.txt',template1_filepath),
        input2 = ('input_file2.txt',template2_filepath),
        input3 = ('input_file3.txt',template3_filepath),
        )
    si_write.input1.assign(
        name = 'name1',
        a    = 1,
        )
    si_write.input2.assign(
        name = 'name2',
        b    = 2,
        )
    si_write.input3.assign(
        name = 'name3',
        c    = 3,
        )
    assert(object_eq(si_write.write(),write_ref))

    si_write.write(input1_filepath)
    assert(input1_filepath.exists())
    assert(input2_filepath.exists())
    assert(input3_filepath.exists())


    # read
    si_read = multi_input_template(**filenames)

    si_read.read(input1_filepath)

    si = si_read
    assert(len(si)==4)
    assert(len(si.filenames)==3)
    assert(object_eq(si.filenames,filenames))
    for name,keyword_set in keywords_ref.items():
        assert(name in si)
        sit = si[name]
        assert(sit.keywords==set())
        assert(isinstance(sit.template,Template))
        assert(object_eq(sit.values,obj()))
        assert(sit.allow_not_set==set())
    #end for
    assert(object_eq(si_read.write(),write_ref))
#end def test_simulation_input_multi_template



def test_code_name():
    from ..simulation import Simulation

    cn = Simulation.code_name()
    assert(isinstance(cn,str))
    assert(' ' not in cn)
#end def test_code_name



def test_init():
    from ..developer import obj
    from ..machines import job,Job
    from ..simulation import Simulation,SimulationInput

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
        imresdir             = './runs/sim_sim',
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
        resdir               = './runs/',
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

    assert(object_eq(se.obj(list(se_ref.keys())),se_ref))
    assert(isinstance(se.simid,int))
    assert(se.simid>=0)
    assert(se.simid<Simulation.sim_count)

    Simulation.clear_all_sims()
    assert(len(Simulation.all_sims)==0)
    assert(len(Simulation.sim_directories)==0)
    assert(Simulation.sim_count==0)


    # make a test job
    test_job = job(machine='ws1',app_command='test.x')

    
    # minimal non-empty init, tests init_job()
    sm = Simulation(job=test_job)

    sm_ref = se_ref.copy()
    del sm_ref.job
    assert(object_eq(sm.obj(list(sm_ref.keys())),sm_ref))
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
        raise FailedTest
    except FailedTest:
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

    Simulation.clear_all_sims()

#end def test_init



def test_virtuals():
    from ..developer import NexusError
    from ..simulation import Simulation

    s = Simulation()

    virts = [
        (s.check_result,[None,None]),
        (s.get_result,[None,None]),
        (s.incorporate_result,[None,None,None]),
        s.app_command,
        s.check_sim_status,
        s.get_output_files,
        ]
    for v in virts:
        args = []
        if isinstance(v,tuple):
            v,args = v
        #end if
        try:
            v(*args)
            raise FailedTest
        except NexusError:
            None
        except FailedTest:
            failed(str(v))
        except Exception as e:
            failed(str(e))
        #end try
    #end for

    vacuous_virts = [
        s.propagate_identifier,
        s.pre_init,
        s.post_init,
        s.pre_create_directories,
        s.write_prep,
        (s.pre_write_inputs,[None]),
        (s.pre_send_files,[None]),
        s.post_submit,
        s.pre_check_status,
        (s.post_analyze,[None]),
        ]
    for v in vacuous_virts:
        args = []
        if isinstance(v,tuple):
            v,args = v
        #end if
        v(*args)
    #end for

    Simulation.clear_all_sims()

#end def test_virtuals



def test_reset_indicators():
    from ..simulation import Simulation

    indicators = '''
        got_dependencies
        setup     
        sent_files
        submitted 
        finished  
        failed    
        got_output
        analyzed  
        '''.split()

    s = Simulation()

    for i in indicators:
        s[i] = True
    #end for

    s.reset_indicators()

    for i in indicators:
        ind = s[i]
        assert(isinstance(ind,bool))
        assert(not ind)
    #end for

    Simulation.clear_all_sims()
#end def test_reset_indicators



def test_indicator_checks():
    from ..machines import job
    from ..simulation import Simulation

    def complete(sim):
        sim.setup      = True
        sim.sent_files = True
        sim.submitted  = True
        sim.finished   = True
        sim.got_output = True
        sim.analyzed   = True
        sim.failed     = False
    #end def complete

    # test completed()
    s = Simulation()
    assert(not s.completed())
    complete(s)
    assert(s.completed())
    s.reset_indicators()
    Simulation.clear_all_sims()

    # test ready() and active()
    test_job = job(machine='ws1',app_command='test.x')

    simdeps = []
    for i in range(5):
        s = Simulation(identifier='id'+str(i),job=test_job)
        complete(s)
        simdeps.append((s,'other'))
    #end for

    s = Simulation(job=test_job,dependencies=simdeps)
    assert(s.ready())
    assert(s.active())
    s.submitted = True
    assert(not s.ready())
    assert(s.active())

    Simulation.clear_all_sims()

#end def test_indicator_checks


@isolate_nexus_core
def test_create_directories(tmp_path):
    import os
    from ..simulation import Simulation

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    s = Simulation()

    assert(not os.path.exists(s.locdir))
    assert(not os.path.exists(s.imlocdir))
    assert(not s.created_directories)

    s.create_directories()

    assert(os.path.exists(s.locdir))
    assert(os.path.exists(s.imlocdir))
    assert(s.created_directories)

    Simulation.clear_all_sims()
#end def test_create_directories


@isolate_nexus_core
def test_file_text(tmp_path):
    from ..simulation import Simulation

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    s = Simulation()
    s.create_directories()

    outfile = Path(s.locdir).resolve() / s.outfile
    errfile = Path(s.locdir).resolve() / s.errfile

    out_text = 'output'
    err_text = 'error'

    outfile.write_text(out_text)
    errfile.write_text(err_text)

    assert(s.outfile_text()==out_text)
    assert(s.errfile_text()==err_text)

    Simulation.clear_all_sims()
#end def test_file_text



def check_dependency_objects(*sims,**kwargs):
    from ..developer import obj
    from ..simulation import Simulation
    empty    = kwargs.get('empty',False)
    wait_ids = kwargs.get('wait_ids',True)
    if len(sims)==1 and isinstance(sims[0],list):
        sims = sims[0]
    #end if
    for sim in sims:
        if empty:
            assert(value_eq(sim.ordered_dependencies,[]))
            assert(isinstance(sim.dependencies,obj))
            assert(len(sim.dependencies)==0)
            assert(sim.dependency_ids==set())
            assert(sim.wait_ids==set())
        else:
            # check dependencies object
            for simid,dep in sim.dependencies.items():
                assert(isinstance(simid,int))
                assert(isinstance(dep,obj))
                assert('result_names' in dep)
                assert('results' in dep)
                assert('sim' in dep)
                assert(len(dep)==3)
                assert(isinstance(dep.sim,Simulation))
                assert(simid==dep.sim.simid)
                assert(isinstance(dep.result_names,list))
                for name in dep.result_names:
                    assert(isinstance(name,str))
                #end for
                assert(isinstance(dep.results,obj))
                assert(len(dep.results)==0)
            #end for
            # check ordered_dependencies object
            for dep in sim.ordered_dependencies:
                dep2 = sim.dependencies[dep.sim.simid]
                assert(id(dep2)==id(dep))
            #end for
            # check dependents object
            for dsimid,dsim in sim.dependents.items():
                assert(isinstance(dsimid,int))
                assert(isinstance(dsim,Simulation))
                assert(dsimid==dsim.simid)
                assert(sim.simid in dsim.dependency_ids)
                assert(sim.simid in dsim.dependencies)
                assert(id(sim)==id(dsim.dependencies[sim.simid].sim))
                found = False
                for s in dsim.ordered_dependencies:
                    found |= id(s.sim)==id(sim)
                #end for
                assert(found)
            #end for
            # check dependency_ids
            for simid in sim.dependency_ids:
                assert(isinstance(simid,int))
                assert(simid in sim.dependencies)
            #end for
            # check wait_ids
            if wait_ids:
                assert(sim.wait_ids==sim.dependency_ids)
            #end if
        #end if
    #end if
#end def check_dependency_objects



def check_dependency(sim2,sim1,quants=['other'],only=False,objects=False):
    # sim2 depends on sim1 for all quantities
    if objects:
        check_dependency_objects(sim1)
        check_dependency_objects(sim2)
    #end if
    assert(sim2.simid in sim1.dependents)
    assert(id(sim1.dependents[sim2.simid])==id(sim2))
    assert(sim1.simid in sim2.dependency_ids)
    assert(sim1.simid in sim2.dependencies)
    assert(id(sim2.dependencies[sim1.simid].sim)==id(sim1))
    assert(set(sim2.dependencies[sim1.simid].result_names)==set(quants))
    if only:
        assert(len(sim1.dependents)==1)
        assert(len(sim2.dependency_ids)==1)
        assert(len(sim2.dependencies)==1)
        assert(len(sim2.ordered_dependencies)==1)
    #end if
#end def check_dependency


def test_depends():
    from ..developer import NexusError
    from ..simulation import Simulation

    # single dependency, single quantity
    s1 = get_sim()
    s2 = get_sim()

    check_dependency_objects(s1,empty=True)
    check_dependency_objects(s2,empty=True)

    s2.depends(s1,'other')

    check_dependency(s2,s1,objects=True,only=True)
    del s1,s2

    s1 = get_sim()
    s2 = get_sim()
    s2.depends((s1,'other'))
    check_dependency(s2,s1,objects=True,only=True)
    del s1,s2

    s1 = get_sim()
    s2 = get_sim(
        dependencies = (s1,'other'),
        )
    check_dependency(s2,s1,objects=True,only=True)
    del s1,s2

    s1 = get_sim()
    s2 = get_sim(
        dependencies = [(s1,'other')],
        )
    check_dependency(s2,s1,objects=True,only=True)
    del s1,s2


    # single dependency, multiple quantities
    s1 = get_test_sim()
    s2 = get_test_sim()

    quants = ['quant1','quant2','quant3']

    check_dependency_objects(s1,empty=True)
    check_dependency_objects(s2,empty=True)

    s2.depends(s1,'quant1','quant2','quant3')

    check_dependency(s2,s1,quants,objects=True,only=True)
    del s1,s2

    s1 = get_test_sim()
    s2 = get_test_sim()
    s2.depends(s1,'quant1')
    s2.depends(s1,'quant2')
    s2.depends(s1,'quant3')
    check_dependency(s2,s1,quants,objects=True,only=True)
    del s1,s2

    s1 = get_test_sim()
    s2 = get_test_sim()
    s2.depends(s1,'quant1','quant2')
    s2.depends(s1,'quant3')
    check_dependency(s2,s1,quants,objects=True,only=True)
    del s1,s2

    s1 = get_test_sim()
    s2 = get_test_sim()
    s2.depends((s1,'quant1','quant2','quant3'))
    check_dependency(s2,s1,quants,objects=True,only=True)
    del s1,s2

    s1 = get_test_sim()
    s2 = get_test_sim()
    s2.depends((s1,'quant1'))
    s2.depends((s1,'quant2'))
    s2.depends((s1,'quant3'))
    check_dependency(s2,s1,quants,objects=True,only=True)
    del s1,s2

    s1 = get_test_sim()
    s2 = get_test_sim(
        dependencies = (s1,'quant1','quant2','quant3'),
        )
    check_dependency(s2,s1,quants,objects=True,only=True)
    del s1,s2

    s1 = get_test_sim()
    s2 = get_test_sim(
        dependencies = [(s1,'quant1','quant2','quant3')],
        )
    check_dependency(s2,s1,quants,objects=True,only=True)
    del s1,s2

    s1 = get_test_sim()
    s2 = get_test_sim(
        dependencies = [
            (s1,'quant1'),
            (s1,'quant2'),
            (s1,'quant3'),
            ],
        )
    check_dependency(s2,s1,quants,objects=True,only=True)
    del s1,s2


    # multiple dependencies
    s11 = get_test_sim()
    s12 = get_test_sim()
    s13 = get_test_sim()

    s21 = get_test_sim(
        dependencies = [
            (s11,'quant1'),
            (s12,'quant2'),
            ]
        )
    s22 = get_test_sim(
        dependencies = [
            (s12,'quant2'),
            (s13,'quant3'),
            ]
        )

    s31 = get_test_sim(
        dependencies = [
            (s21,'quant1'),
            (s22,'quant2'),
            ]
        )
    s32 = get_test_sim(
        dependencies = [
            (s21,'quant1'),
            (s22,'quant2'),
            ]
        )
    s33 = get_test_sim(
        dependencies = [
            (s21,'quant1'),
            (s22,'quant2'),
            ]
        )

    s41 = get_test_sim(
        dependencies = [
            (s11,'quant1'),
            (s22,'quant2'),
            (s32,'quant3'),
            ]
        )

    check_dependency_objects(s11,s12,s13,s21,s22,s31,s32,s33,s41)

    check_dependency(s21,s11,['quant1'])
    check_dependency(s21,s12,['quant2'])

    check_dependency(s22,s12,['quant2'])
    check_dependency(s22,s13,['quant3'])

    check_dependency(s31,s21,['quant1'])
    check_dependency(s31,s22,['quant2'])

    check_dependency(s32,s21,['quant1'])
    check_dependency(s32,s22,['quant2'])

    check_dependency(s33,s21,['quant1'])
    check_dependency(s33,s22,['quant2'])

    check_dependency(s41,s11,['quant1'])
    check_dependency(s41,s22,['quant2'])
    check_dependency(s41,s32,['quant3'])

    del s11,s12,s13,s21,s22,s31,s32,s33,s41


    # fail when dependency does not exist
    try:
        s1 = get_sim()
        s2 = get_sim(
            dependencies = [(s1,'quant1')],
            )
        raise FailedTest
    except NexusError:
        None
    except:
        failed()
    #end try

    try:
        s1 = get_sim()
        s2 = get_sim(
            dependencies = [(s1,'other','quant2')],
            )
        raise FailedTest
    except NexusError:
        None
    except:
        failed()
    #end try

    try:
        s1 = get_test_sim()
        s2 = get_test_sim(
            dependencies = [(s1,'quant1','apple')],
            )
        raise FailedTest
    except NexusError:
        None
    except:
        failed()
    #end try

    Simulation.clear_all_sims()
#end def test_depends



def test_undo_depends():
    from ..simulation import Simulation

    # single dependency, single quantity
    s1 = get_sim()
    s2 = get_sim()

    check_dependency_objects(s1,empty=True)
    check_dependency_objects(s2,empty=True)

    s2.depends(s1,'other')

    check_dependency(s2,s1,objects=True,only=True)

    s2.undo_depends(s1)

    check_dependency_objects(s1,empty=True)
    check_dependency_objects(s2,empty=True)

    Simulation.clear_all_sims()
#end def test_undo_depends



def test_has_generic_input():
    from ..simulation import Simulation
    from ..simulation import SimulationInput,GenericSimulationInput

    s = get_sim()
    assert(not s.has_generic_input())
    del s

    class GenInput(SimulationInput,GenericSimulationInput):
        None
    #end class GenInput

    s = get_sim(
        input = GenInput(),
        )
    assert(s.has_generic_input())
    del s

    Simulation.clear_all_sims()
#end def test_has_generic_input


@isolate_nexus_core
def test_check_dependencies():
    from ..developer import obj, NexusError
    from ..simulation import Simulation
    from ..simulation import SimulationInput,GenericSimulationInput

    result = obj()
    result.dependencies_satisfied = True

    s11 = get_test_sim()
    s12 = get_test_sim()
    s13 = get_test_sim()

    s21 = get_test_sim(
        dependencies = [
            (s11,'quant1'),
            (s12,'quant2'),
            ]
        )
    s22 = get_test_sim(
        dependencies = [
            (s12,'quant2'),
            (s13,'quant3'),
            ]
        )

    s31 = get_test_sim(
        dependencies = [
            (s21,'quant1'),
            (s22,'quant2'),
            ]
        )
    s32 = get_test_sim(
        dependencies = [
            (s21,'quant1'),
            (s22,'quant2'),
            ]
        )
    s33 = get_test_sim(
        dependencies = [
            (s21,'quant1'),
            (s22,'quant2'),
            ]
        )

    s41 = get_test_sim(
        dependencies = [
            (s11,'quant1'),
            (s22,'quant2'),
            (s32,'quant3'),
            ]
        )

    sims = [s11,s12,s13,s21,s22,s31,s32,s33,s41]
    for s in sims:
        s.check_dependencies(result)
    #end for
    assert(result.dependencies_satisfied)


    # non-existent dependency
    try:
        s  = get_test_sim()
        s2 = get_test_sim(dependencies=((s,'nonexistent')))
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try


    # existent dependency but generic input
    class GenInput(SimulationInput,GenericSimulationInput):
        None
    #end class GenInput

    s = get_test_sim(
        input = GenInput(),
        )

    s2 = get_test_sim(
        input = GenInput(),
        dependencies = (s,'quant1')
        )

    result = obj(dependencies_satisfied=True)
    s.check_dependencies(result)
    assert(result.dependencies_satisfied)

    try:
        s2.check_dependencies(result)
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    Simulation.clear_all_sims()
#end def test_check_dependencies



def test_get_dependencies():
    from ..developer import obj
    from ..simulation import Simulation

    simdeps = obj()

    deps = []
    s11 = get_test_sim()
    simdeps[s11.simid] = deps

    s12 = get_test_sim()
    simdeps[s12.simid] = deps

    s13 = get_test_sim()
    simdeps[s13.simid] = deps

    deps = [
        (s11,'quant1'),
        (s12,'quant2'),
        ]
    s21 = get_test_sim(dependencies=deps)
    simdeps[s21.simid] = deps

    deps = [
        (s12,'quant2'),
        (s13,'quant3'),
        ]
    s22 = get_test_sim(dependencies=deps)
    simdeps[s22.simid] = deps

    deps = [
        (s21,'quant1'),
        (s22,'quant2'),
        ]
    s31 = get_test_sim(dependencies=deps)
    simdeps[s31.simid] = deps

    deps = [
        (s21,'quant1'),
        (s22,'quant2'),
        ]
    s32 = get_test_sim(dependencies=deps)
    simdeps[s32.simid] = deps

    deps = [
        (s21,'quant1'),
        (s22,'quant2'),
        ]
    s33 = get_test_sim(dependencies=deps)
    simdeps[s33.simid] = deps

    deps = [
        (s11,'quant1'),
        (s22,'quant2'),
        (s32,'quant3'),
        ]
    s41 = get_test_sim(dependencies=deps)
    simdeps[s41.simid] = deps

    sims = [s11,s12,s13,s21,s22,s31,s32,s33,s41]
    assert(len(simdeps)==len(sims))
    for s in sims:
        assert(not s.got_dependencies)
        assert(len(s.input.result_data)==0)
        s.get_dependencies()
        resdata = s.input.result_data
        deps = simdeps[s.simid]
        for sim,resname in deps:
            assert(sim.simid in resdata)
            assert(resdata[sim.simid]==resname)
        #end for
    #end for

    Simulation.clear_all_sims()
#end def test_get_dependencies



def test_downstream_simids():
    from ..developer import obj
    from ..simulation import Simulation

    s11 = get_test_sim()
    s12 = get_test_sim()
    s13 = get_test_sim()

    s21 = get_test_sim(
        dependencies = [
            (s11,'quant1'),
            (s12,'quant2'),
            ]
        )
    s22 = get_test_sim(
        dependencies = [
            (s12,'quant2'),
            (s13,'quant3'),
            ]
        )

    s31 = get_test_sim(
        dependencies = [
            (s21,'quant1'),
            (s22,'quant2'),
            ]
        )
    s32 = get_test_sim(
        dependencies = [
            (s21,'quant1'),
            (s22,'quant2'),
            ]
        )
    s33 = get_test_sim(
        dependencies = [
            (s21,'quant1'),
            (s22,'quant2'),
            ]
        )

    s41 = get_test_sim(
        dependencies = [
            (s11,'quant1'),
            (s22,'quant2'),
            (s32,'quant3'),
            ]
        )

    sims = obj(
        s11 = s11,
        s12 = s12,
        s13 = s13,
        s21 = s21,
        s22 = s22,
        s31 = s31,
        s32 = s32,
        s33 = s33,
        s41 = s41,
        )

    downstream_sims = obj(
        s11 = [s21,s31,s32,s33,s41],
        s12 = [s21,s22,s31,s32,s33,s41],
        s13 = [s22,s31,s32,s33,s41],
        s21 = [s31,s32,s33,s41],
        s22 = [s31,s32,s33,s41],
        s31 = [],
        s32 = [s41],
        s33 = [],
        s41 = [],
        )

    n = 0
    for sname in sorted(sims.keys()):
        s = sims[sname]
        ds_ids = s.downstream_simids()
        ds_ids_ref = set([sd.simid for sd in downstream_sims[sname]])
        assert(ds_ids==ds_ids_ref)
        n+=1
    #end for
    assert(n==9)

    Simulation.clear_all_sims()
#end def test_downstream_simids



def test_copy_file(tmp_path):
    import os
    from ..simulation import Simulation
    
    opath = tmp_path / 'other'
    if not os.path.exists(opath):
        os.makedirs(opath)
    #end if

    file1 = tmp_path / 'file.txt'
    file2 = opath / 'file.txt'

    file1.write_text('text')
    assert(file1.exists())

    s = get_sim()

    s.copy_file(file1,opath)

    assert(file2.exists())
    assert(file2.read_text().strip()=='text')
    
    Simulation.clear_all_sims()
#end def test_copy_file


@isolate_nexus_core
def test_save_load_image(tmp_path):
    from ..developer import obj
    from ..simulation import Simulation,SimulationImage

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    nsave = 30
    nload = 22

    assert(len(SimulationImage.save_fields)==nsave)
    assert(len(SimulationImage.load_fields)==nload)
    assert(len(SimulationImage.save_only_fields&SimulationImage.load_fields)==0)

    sim = get_sim()

    sim.create_directories()

    sim.save_image()

    imagefile = Path(sim.imlocdir) / sim.sim_image
    assert(imagefile.exists())

    image = obj()
    image.load(imagefile)
    assert(len(image)==nsave)
    for field in SimulationImage.save_fields:
        assert(field in image)
        assert(field in sim)
        assert(value_eq(image[field],sim[field]))
    #end for

    orig = obj()
    for field in SimulationImage.load_fields:
        orig[field] = sim[field]
        del sim[field]
    #end for
    sim.sim_image = orig.sim_image
    sim.load_image()
    for field in SimulationImage.load_fields:
        assert(field in sim)
        assert(value_eq(sim[field],orig[field]))
    #end for
    Simulation.clear_all_sims()
#end def test_save_load_image


@isolate_nexus_core
def test_load_analyzer_image(tmp_path):
    from ..simulation import Simulation

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    sim = get_test_sim()

    if not Path(sim.imresdir).exists():
        Path(sim.imresdir).mkdir(parents=True)
    #end if

    analyzer_file = Path(sim.imresdir).resolve() / sim.analyzer_image

    a = sim.analyzer_type(None)
    assert(not a.analysis_performed)
    a.analyze()
    assert(a.analysis_performed)
    a.save(analyzer_file)
    assert(Path(analyzer_file).exists())

    a2 = sim.load_analyzer_image()
    assert(isinstance(a2,sim.analyzer_type))
    assert(a2.analysis_performed)
    assert(object_eq(a2,a))

    Simulation.clear_all_sims()
#end def test_load_analyzer_image


@isolate_nexus_core
def test_save_attempt(tmp_path):
    from ..simulation import Simulation

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    sim = get_test_sim()

    sim.create_directories()

    files = (sim.infile,sim.outfile,sim.errfile)

    assert(sim.attempt_files()==files)
    for file in files:
        (Path(sim.locdir).resolve() / file).write_text('made an attempt')
    #end for

    attempt_dir = Path(sim.locdir).resolve() / f'{sim.identifier}_attempt1'
    assert(not attempt_dir.exists())
    sim.save_attempt()
    assert(attempt_dir.exists())
    for file in files:
        assert(not (Path(sim.locdir).resolve() / file).exists())
        assert((Path(attempt_dir).resolve() / file).exists())
    #end for

    Simulation.clear_all_sims()
#end def test_save_attempt


@isolate_nexus_core
def test_write_inputs(tmp_path):
    from ..simulation import Simulation,input_template

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    template = '''
name = "$name"
a    = $a
'''

    input_ref = '''
name = "input_name"
a    = 1
'''
    
    si = input_template(text=template)
    si.assign(name='input_name',a=1)

    s = get_test_sim(
        input = si,
        )
    s.create_directories()

    input_file       = Path(s.locdir).resolve() / s.infile
    image_file       = Path(s.imlocdir).resolve() / s.sim_image
    input_image_file = Path(s.imlocdir).resolve() / s.input_image

    assert(not s.setup)
    assert(not input_file.exists())
    assert(not image_file.exists())
    assert(not input_image_file.exists())

    s.write_inputs()

    assert(s.setup)
    assert(input_file.exists())
    assert(image_file.exists())
    assert(input_image_file.exists())

    assert(input_file.read_text()==input_ref)

    s.setup = False
    s.load_image()
    assert(s.setup)

    Simulation.clear_all_sims()
#end def test_write_inputs


@isolate_nexus_core
def test_send_files(tmp_path):
    from ..nexus_base import nexus_core
    from ..simulation import Simulation

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    # make fake data files
    data_file1 = 'data_file1.txt'
    data_file2 = 'data_file2.txt'

    (tmp_path / data_file1).write_text('data1')
    (tmp_path / data_file2).write_text('data2')

    data_files = [data_file1,data_file2] 

    s = get_test_sim(
        files = data_files,
        )
    s.infile = None

    assert(s.locdir==s.remdir)
    assert(s.files==set(data_files))

    s.create_directories()

    loc_data_file1 = Path(s.locdir).resolve() / data_file1
    loc_data_file2 = Path(s.locdir).resolve() / data_file2

    assert(not s.sent_files)
    assert(not loc_data_file1.exists())
    assert(not loc_data_file2.exists())

    s.send_files()

    assert(s.sent_files)
    assert(loc_data_file1.exists())
    assert(loc_data_file2.exists())

    assert(loc_data_file1.read_text()=='data1')
    assert(loc_data_file2.read_text()=='data2')

    s.sent_files = False
    s.load_image()
    assert(s.sent_files)

    Simulation.clear_all_sims()
#end def test_send_files


@isolate_nexus_core
def test_submit(tmp_path):
    from ..machines import job
    from ..simulation import Simulation

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    s = get_test_sim(
        job = job(machine='ws1',app_command='echo run'),
        )

    j = s.job
    m = j.get_machine()

    s.create_directories()

    assert(not s.submitted)
    assert(not s.job.submitted)
    assert(j.internal_id not in m.jobs)
    assert(j.internal_id not in m.waiting)

    s.submit()

    assert(s.submitted)
    assert(s.job.submitted)
    assert(j.internal_id in m.jobs)
    assert(j.internal_id in m.waiting)

    Simulation.clear_all_sims()
#end def test_submit


@isolate_nexus_core
def test_update_process_id(tmp_path):
    from ..simulation import Simulation

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    s = get_test_sim()
    j = s.job

    s.create_directories()

    ref_pid = 0

    assert(s.process_id is None)
    assert(j.system_id is None)

    j.system_id = ref_pid

    s.update_process_id()

    assert(s.process_id==ref_pid)

    s.process_id = None
    s.load_image()
    assert(s.process_id==ref_pid)

    Simulation.clear_all_sims()
#end def test_update_process_id


@isolate_nexus_core
def test_check_status(tmp_path):
    from ..simulation import Simulation

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    s = get_test_sim()
    j = s.job

    assert(not s.finished)
    assert(not j.finished)

    s.check_status()
    assert(not s.finished)

    s.create_directories()

    (Path(s.locdir).resolve() / s.outfile).write_text('out')
    (Path(s.locdir).resolve() / s.errfile).write_text('err')
    j.finished = True

    s.check_status()

    assert(s.finished)

    s.finished = False
    s.load_image()
    assert(s.finished)

    Simulation.clear_all_sims()
#end def test_check_status


@isolate_nexus_core
def test_get_output(tmp_path):
    from ..simulation import Simulation

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    s = get_test_sim()

    s.create_directories()

    remote_image  = Path(s.imremdir).resolve() / s.sim_image
    results_image = Path(s.imresdir).resolve() / s.sim_image

    assert(not remote_image.exists())
    assert(not results_image.exists())
    assert(not s.finished)

    assert(value_eq(s.get_output_files(),[]))

    files = [s.infile,s.outfile,s.errfile]

    loc_files = []
    res_files = []

    for file in files:
        loc_files.append(Path(s.locdir).resolve() / file)
        res_files.append(Path(s.resdir).resolve() / file)
    #end for

    for loc_file in loc_files:
        loc_file.write_text('contents')
    #end for
    s.finished = True

    assert(not s.got_output)
    for loc_file,res_file in zip(loc_files,res_files):
        assert(loc_file.exists())
        if s.resdir!=s.locdir:
            assert(not res_file.exists())
        else:
            assert(res_file.exists())
        #end if
    #end for

    s.get_output()

    assert(s.got_output)
    for loc_file,res_file in zip(loc_files,res_files):
        assert(loc_file.exists())
        assert(res_file.exists())
    #end for

    s.got_output = False
    s.load_image()
    assert(s.got_output)

    Simulation.clear_all_sims()
#end def test_get_output


@isolate_nexus_core
def test_analyze(tmp_path):
    from ..simulation import Simulation

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    s = get_test_sim()

    s.create_directories()

    assert(not s.finished)
    s.finished = True

    analyzer_image = Path(s.imresdir).resolve() / s.analyzer_image

    assert(not s.analyzed)
    assert(not analyzer_image.exists())

    s.analyze()

    assert(s.analyzed)
    assert(analyzer_image.exists())

    s.analyzed = False
    s.load_image()
    assert(s.analyzed)

    Simulation.clear_all_sims()
#end def test_analyze


@isolate_nexus_core
def test_progress(tmp_path):
    from ..nexus_base import nexus_core
    from ..simulation import Simulation,input_template

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    assert(nexus_core.mode==nexus_core.modes.stages)
    assert(len(nexus_core.stages)==0)

    nexus_core.stages     = list(nexus_core.primary_modes)
    nexus_core.stages_set = set(nexus_core.stages)

    primary_modes = ['setup','send_files','submit','get_output','analyze']
    assert(value_eq(nexus_core.stages,primary_modes))
    assert(value_eq(nexus_core.stages_set,set(primary_modes)))


    template = '''
name = "$name"
a    = $a
'''
    si = input_template(text=template)
    si.assign(name='input_name',a=1)

    s = get_test_sim(input=si)

    indicators = [
        'created_directories',
        'got_dependencies',
        'setup',
        'sent_files',
        'submitted',
        'finished',
        'failed',
        'got_output',
        'analyzed',
        ]

    inds = obj()

    # first progression
    #   directory creation
    #   dependency processing
    #   input file write
    #   auxilliary file transfer
    #   job submission
    assert(not s.created_directories)
    assert(not s.got_dependencies)
    assert(not s.setup)
    assert(not s.sent_files)
    assert(not s.submitted)
    assert(not s.finished)
    assert(not s.got_output)
    assert(not s.analyzed)
    assert(s.files==set())
    assert(s.job.status==0)
    assert(not Path(s.locdir).exists())
    assert(not Path(s.remdir).exists())
    assert(not Path(s.resdir).exists())
    assert(not Path(s.imlocdir).exists())
    assert(not Path(s.imremdir).exists())
    assert(not Path(s.imresdir).exists())

    s.progress()

    assert(s.created_directories)
    assert(s.got_dependencies)
    assert(s.setup)
    assert(s.sent_files)
    assert(s.submitted)
    assert(not s.finished)
    assert(not s.got_output)
    assert(not s.analyzed)
    assert(s.files==set([s.infile]))
    assert(s.job.status==1)
    assert(Path(s.locdir).exists())
    assert(Path(s.remdir).exists())
    assert(Path(s.imlocdir).exists())
    assert(Path(s.imremdir).exists())
    assert((Path(s.locdir).resolve() / s.infile).exists())
    assert((Path(s.imlocdir).resolve() / s.sim_image).exists())
    assert((Path(s.imlocdir).resolve() / s.input_image).exists())
    assert(not (Path(s.locdir).resolve() / s.outfile).exists())
    assert(not (Path(s.locdir).resolve() / s.errfile).exists())
    assert(not (Path(s.imlocdir).resolve() / s.analyzer_image).exists())
    if s.resdir!=s.locdir:
        assert(not Path(s.resdir).exists())
        assert(not Path(s.imresdir).exists())
    else:
        assert(Path(s.resdir).exists())
        assert(Path(s.imresdir).exists())
    #end if
    
    # check image
    inds.transfer_from(s,indicators)
    s.reset_indicators()
    s.load_image()
    assert(s.setup)
    assert(s.sent_files)
    assert(not s.submitted) # submitted is not stored yet
    assert(not s.finished)
    assert(not s.got_output)
    assert(not s.analyzed)
    s.transfer_from(inds,indicators)


    # simulate job completion
    #   create output and error files
    #   set job status to finished
    (Path(s.locdir).resolve() / s.outfile).write_text('out')
    (Path(s.locdir).resolve() / s.errfile).write_text('err')
    s.job.finished = True

    assert((Path(s.locdir) / s.outfile).exists())
    assert((Path(s.locdir) / s.errfile).exists())


    # second progression
    #   status check
    #   output file transfer
    #   output analysis
    assert(s.created_directories)
    assert(s.got_dependencies)
    assert(s.setup)
    assert(s.sent_files)
    assert(s.submitted)
    assert(not s.finished)
    assert(not s.got_output)
    assert(not s.analyzed)

    s.progress()

    assert(s.created_directories)
    assert(s.got_dependencies)
    assert(s.setup)
    assert(s.sent_files)
    assert(s.submitted)
    assert(s.finished)
    assert(s.got_output)
    assert(s.analyzed)
    assert(Path(s.resdir).exists())
    assert(Path(s.imresdir).exists())
    assert((Path(s.resdir) / s.infile).exists())
    assert((Path(s.resdir) / s.errfile).exists())
    assert((Path(s.resdir) / s.outfile).exists())
    assert((Path(s.imresdir) / s.sim_image).exists())
    assert((Path(s.imresdir) / s.input_image).exists())
    assert((Path(s.imresdir) / s.analyzer_image).exists())
    if s.resdir!=s.locdir:
        assert(not (Path(s.imlocdir) / s.analyzer_image).exists())
    else:
        assert((Path(s.imlocdir) / s.analyzer_image).exists())
    #end if

    # check image
    inds.transfer_from(s,indicators)
    s.reset_indicators()
    s.load_image()
    assert(s.setup)
    assert(s.sent_files)
    assert(s.submitted)
    assert(s.finished)
    assert(s.got_output)
    assert(s.analyzed)
    s.transfer_from(inds,indicators)

    
    # attempt third progression
    #   nothing should happen
    sbef = s.copy()
    sbef.input.template = s.input.template

    s.progress()

    assert(object_eq(s,sbef))

    Simulation.clear_all_sims()
#end def test_progress


@isolate_nexus_core
def test_execute(tmp_path):
    from ..machines import job
    from ..simulation import Simulation

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    import shutil
    serial = shutil.which('mpirun') is None

    s = get_test_sim(
        job = job(machine='ws1',app_command='echo run',serial=serial),
        )

    s.create_directories()

    outfile = Path(s.locdir) / s.outfile
    errfile = Path(s.locdir) / s.errfile

    assert(not s.submitted)
    assert(not s.job.finished)
    assert(s.job.status==0)
    assert(not outfile.exists())
    assert(not errfile.exists())

    s.execute()

    assert(s.submitted)
    assert(s.job.finished)
    assert(s.job.status==4)
    assert(outfile.exists())
    assert(errfile.exists())
    assert(outfile.read_text().strip()=='run')
    err_contents = errfile.read_text().strip()
    # Handle spurious error message from OpenMPI
    #   see also: https://github.com/QMCPACK/qmcpack/pull/4339#discussion_r1033813856
    err_contents = err_contents.replace('Invalid MIT-MAGIC-COOKIE-1 key','').strip()
    assert(err_contents=='')

    Simulation.clear_all_sims()
#end def test_execute



def test_reset_wait_ids():
    from ..simulation import Simulation

    for i in range(n_test_workflows):
        sims = get_test_workflow(i)
        for s in sims:
            s.wait_ids = None
        #end for
        for s in sims:
            if len(s.dependencies)==0:
                s.reset_wait_ids()
            #end if
        #end for
        for s in sims:
            assert(isinstance(s.wait_ids,set))
            assert(s.wait_ids==s.dependency_ids)
        #end for
    #end for

    Simulation.clear_all_sims()
#end def test_reset_wait_ids



def test_check_subcascade():
    from ..simulation import Simulation

    def finish(sim):
        sim.finished = True
    #end def finish

    for i in range(n_test_workflows):
        sims = get_test_workflow(i)

        # no cascades are finished
        for s in sims:
            assert(not s.finished)
        #end for
        for s in sims:
            if len(s.dependencies)==0:
                finished = s.check_subcascade()
                assert(isinstance(finished,bool))
                assert(not finished)
            #end if
        #end for

        # all cascades are finished
        for s in sims:
            s.finished = True
        #end for
        for s in sims:
            if len(s.dependencies)==0:
                finished = s.check_subcascade()
            #end if
        #end for

        # only a single cascade is finished
        for s in sims:
            s.finished = False
        #end for
        single = None
        for s in sims:
            if len(s.dependencies)==0:
                if single is None:
                    single = s
                #end if
            #end if
        #end for
        single.traverse_full_cascade(finish)
        for s in sims:
            if len(s.dependencies)==0:
                finished = s.check_subcascade()
                if id(s)==id(single):
                    if not finished:
                        from ..simulation import graph_sims
                        for sim in sims:
                            if sim.finished:
                                sim.block = True
                            #end if
                        #end for
                        graph_sims(sims.list())
                    #end if
                    assert(finished)
                else:
                    assert(not finished)
                #end if
            #end if
        #end for

        # all simulations are finished except one
        # not all cascades are finished
        for s in sims:
            s.finished = True
        #end for
        n = 0
        for key in sorted(sims.keys()):
            n+=1
            if 2*n>len(sims):
                sims[key].finished = False
            #end if
        #end for
        finished = True
        for s in sims:
            if len(s.dependencies)==0:
                finished &= s.check_subcascade()
            #end if
        #end for
        assert(not finished)
    #end for

    Simulation.clear_all_sims()
#end def test_check_subcascade



def test_block_dependents():
    from ..simulation import Simulation

    def assert_blocked(sim):
        assert(sim.block)
        assert(sim.block_subcascade)
    #end def assert_blocked

    for i in range(n_test_workflows):
        sims = get_test_workflow(i)
        for s in sims:
            if len(s.dependencies)==0:
                s.block_dependents()
                s.traverse_full_cascade(assert_blocked)
            #end if
        #end for
    #end for

    Simulation.clear_all_sims()
#end def test_block_dependents


@isolate_nexus_core
def test_reconstruct_cascade(tmp_path):
    from ..simulation import Simulation

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    sims = get_test_workflow(2)
    assert(len(sims)==7)

    for s in sims:
        imagefile = Path(s.imlocdir) / s.sim_image
        assert(not imagefile.exists())
        assert(not s.loaded)
        assert(not s.submitted)
        assert(not s.finished)
        assert(s.process_id is None)
        assert(s.job.system_id is None)
    #end for

    for s in sims:
        s.create_directories()
        s.save_image()
    #end for

    for s in sims:
        imagefile = Path(s.imlocdir) / s.sim_image
        assert(imagefile.exists())
        assert(not s.loaded)
        assert(not s.submitted)
        assert(not s.finished)
        assert(s.process_id is None)
        assert(s.job.system_id is None)
    #end for

    sims.s1.reconstruct_cascade()

    for s in sims:
        imagefile = Path(s.imlocdir) / s.sim_image
        assert(imagefile.exists())
        assert(s.loaded)
        assert(not s.submitted)
        assert(not s.finished)
        assert(s.process_id is None)
        assert(s.job.system_id is None)
    #end for

    Simulation.clear_all_sims()


    sims = get_test_workflow(2)

    def get_process_id():
        get_process_id.current += 1
        return get_process_id.current
    #end def get_process_id
    get_process_id.current = 0

    def finish(s):
        s.setup            = True
        s.sent_files       = True
        s.submitted        = True
        s.finished         = True
        s.failed           = True
        s.got_output       = True
        s.analyzed         = True
        s.process_id       = get_process_id()
    #end def finish

    def clear(s):
        s.loaded           = False
        s.setup            = False
        s.sent_files       = False
        s.submitted        = False
        s.finished         = False
        s.failed           = False
        s.got_output       = False
        s.analyzed         = False
        s.process_id       = None
        s.job.system_id    = None
    #end def clear

    def finished(s):
        f = True
        f &= s.setup            
        f &= s.sent_files       
        f &= s.submitted        
        f &= s.finished         
        f &= s.failed           
        f &= s.got_output       
        f &= s.analyzed         
        f &= isinstance(s.process_id,int)
        return f
    #end def finished

    def empty(s):
        e = True
        e &= not s.setup            
        e &= not s.sent_files       
        e &= not s.submitted        
        e &= not s.finished         
        e &= not s.failed           
        e &= not s.got_output       
        e &= not s.analyzed         
        e &= s.process_id is None
        e &= s.job.system_id is None
        return e
    #end def empty

    def cleared(s):
        return empty(s) and not s.loaded
    #end def cleared

    for s in sims:
        assert(cleared(s))
    #end for

    finish(sims.s1)
    finish(sims.s2)
    finish(sims.s3)

    s = sims.s41
    s.got_dependencies = True
    s.setup            = True
    s.sent_files       = True
    s.submitted        = True
    s.finished         = True
    s.process_id       = get_process_id()

    s = sims.s42
    s.got_dependencies = True
    s.setup            = True
    s.sent_files       = True
    s.submitted        = True
    s.process_id       = get_process_id()
    
    for s in sims:
        s.create_directories()
        s.save_image()
        clear(s)
    #end for

    for s in sims:
        assert(cleared(s))
    #end for

    sims.s1.reconstruct_cascade()

    for s in sims:
        assert(s.loaded)
    #end for

    assert(finished(sims.s1))
    assert(finished(sims.s2))
    assert(finished(sims.s3))

    s = sims.s41
    assert(not finished(s) and not empty(s))
    assert(s.got_dependencies)
    assert(s.setup           )
    assert(s.sent_files      )
    assert(s.submitted       )
    assert(s.finished        )
    assert(not s.failed      )     
    assert(not s.got_output  )     
    assert(not s.analyzed    )
    assert(s.process_id==4   )
    assert(s.job.system_id is None )

    s = sims.s42
    assert(not finished(s) and not empty(s))
    assert(s.got_dependencies)
    assert(s.setup           )
    assert(s.sent_files      )
    assert(s.submitted       )
    assert(not s.finished    )
    assert(not s.failed      )     
    assert(not s.got_output  )     
    assert(not s.analyzed    )
    assert(s.process_id==5   )
    assert(s.job.system_id==5)

    assert(empty(sims.s51))
    assert(empty(sims.s52))

    Simulation.clear_all_sims()
#end def test_reconstruct_cascade



def test_traverse_cascade():
    from ..simulation import Simulation

    def count_visits(sim,visit_counts):
        i = sim.simid
        if i not in visit_counts:
            visit_counts[i] = 1
        else:
            visit_counts[i] += 1
        #end if
    #end def count_visits

    for i in range(n_test_workflows):
        sims = get_test_workflow(i)
        counts = dict()
        for s in sims:
            if len(s.dependencies)==0:
                s.traverse_cascade(count_visits,counts)
            #end if
        #end for
        assert(len(counts)==len(sims))
        for s in sims:
            assert(s.simid in counts)
            assert(counts[s.simid]==1)
        #end for
    #end for

    Simulation.clear_all_sims()
#end def test_traverse_cascade



def test_traverse_full_cascade():
    from ..simulation import Simulation

    def finish(sim):
        sim.finished = True
    #end def finish

    for i in range(n_test_workflows):
        sims = get_test_workflow(i)
        for s in sims:
            assert(not s.finished)
        #end for
        for s in sims:
            if len(s.dependencies)==0:
                s.traverse_full_cascade(finish)
            #end if
        #end for
        for s in sims:
            assert(s.finished)
        #end for
    #end for
    
    Simulation.clear_all_sims()
#end def test_traverse_full_cascade


@isolate_nexus_core
def test_write_dependents():
    from ..simulation import Simulation

    for i in range(n_test_workflows):
        sims = get_test_workflow(i)
        for s in sims:
            if len(s.dependencies)==0:
                s.write_dependents()
            #end if
        #end for
    #end for

    Simulation.clear_all_sims()
#end def test_write_dependents



def test_generate_simulation():
    from ..developer import NexusError
    from ..simulation import Simulation,GenericSimulation
    from ..simulation import generate_simulation

    try:
        sim = generate_simulation(sim_type='unknown')
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    sim = generate_simulation()
    assert(isinstance(sim,Simulation))
    assert(isinstance(sim,GenericSimulation))

    Simulation.clear_all_sims()
#end def test_generate_simulation


def test_generic_simulation(tmp_path):
    from ..simulation import Simulation,GenericSimulation
    from ..simulation import generate_simulation,SimulationInputTemplate
    from ..machines import job

    # Test 1: Create GenericSimulation with string input
    script_text = 'print("Hello from GenericSimulation")'
    sim1 = generate_simulation(
        identifier = 'test_generic_string',
        path       = str(tmp_path / 'test1'),
        job        = job(machine='ws1', serial=True, app='python3'),
        input      = str(script_text),
        outfiles   = ['output.txt'],
        )
    assert(isinstance(sim1,GenericSimulation))
    assert(isinstance(sim1,Simulation))
    assert(sim1.outfiles == ['output.txt'])
    assert(isinstance(sim1.input,SimulationInputTemplate))
    
    sim1_input_text = sim1.input.write_text()
    assert(sim1_input_text == script_text)

    assert(sim1.get_output_files() == ['output.txt'])

    assert(sim1.app_command() is not None)
    assert('python3' in sim1.app_command() or sim1.job.app_name == 'python3')


    script_file = tmp_path / 'test_script.py'
    script_file_content = 'print("Hello from file")\n'
    script_file.write_text(script_file_content)
    #end with
    
    # Test 2: GenericSimulation with file path input
    sim2 = generate_simulation(
        identifier = 'test_generic_file',
        path       = str(tmp_path / 'test2'),
        job        = job(machine='ws1', serial=True, app='python3'),
        input      = str(script_file),
        outfiles   = ['result.txt'],
        )
    assert(isinstance(sim2,GenericSimulation))
    assert(sim2.outfiles == ['result.txt'])
    assert(isinstance(sim2.input,SimulationInputTemplate))
    
    sim2_input_text = sim2.input.write_text()
    assert(sim2_input_text == script_file_content)

    assert(sim2.get_output_files() == ['result.txt'])

    assert(sim2.app_command() is not None)
    assert('python3' in sim2.app_command() or sim2.job.app_name == 'python3')

    Simulation.clear_all_sims()
#end def test_generic_simulation



def test_graph_sims():
    _ = pytest.importorskip("matplotlib")
    _ = pytest.importorskip("pydot")
    from ..simulation import Simulation,graph_sims

    sims = get_test_workflow(3)

    graph_sims(sims.list(),display=False,exit=False)

    Simulation.clear_all_sims()
#end def test_graph_sims
