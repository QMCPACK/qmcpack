
import testing
from testing import value_eq,object_eq,failed,FailedTest
from testing import divert_nexus_log,restore_nexus_log


all_machines = []
machines_data = dict()


def get_machine_data():
    from generic import obj
    from machines import Machine,Workstation,Supercomputer
    if len(machines_data)==0:
        workstations   = obj()
        supercomputers = obj()
        for machine in Machine.machines:
            if isinstance(machine,Workstation):
                workstations.append(machine)
            elif isinstance(machine,Supercomputer):
                supercomputers[machine.name] = machine
            else:
                failed()
            #end if
        #end for
        machines_data['ws'] = workstations
        machines_data['sc'] = supercomputers
    #end if
    ws = machines_data['ws']
    sc = machines_data['sc']
    return ws,sc
#end def get_machine_data


def get_all_machines():
    from machines import Machine
    if len(all_machines)==0:
        for m in Machine.machines:
            all_machines.append(m)
        #end for
    #end if
    return all_machines
#end def get_all_machines


def get_supercomputers():
    ws,sc = get_machine_data()
    return sc
#end def get_supercomputers




def test_import():
    from machines import Machine,Workstation,InteractiveCluster,Supercomputer
    from machines import Job,job
    from machines import get_machine,get_machine_name
#end def test_import



def test_cpu_count():
    from machines import cpu_count
    assert(isinstance(cpu_count(),int))
#end def test_cpu_count



def test_options():
    from generic import obj
    from machines import Options

    # empty init
    o = Options()
    assert(len(o)==0)

    # filled init
    inputs = dict(
        n = '-n 1',
        d = '-d 2',
        exe = '--exe',
        )
    oi = Options(**inputs)
    assert(oi.to_dict()==inputs)

    # add
    oa = Options()
    oa.add(**inputs)
    assert(object_eq(oa,oi))

    # read
    opts = '-a -b 1 -c=2 --dname="0 1 2" --evar -fval -gval other --hval'
    ref = obj(
        a     = '-a',
        b     = '-b 1',
        c     = '-c=2',
        dname = '--dname="0 1 2"',
        evar  = '--evar',
        fval  = '-fval',
        gval  = '-gval other',
        hval  = '--hval',
        )
    o = Options()
    o.read(opts)
    assert(object_eq(o.to_obj(),ref))

    # write
    opts_write = o.write()
    o2 = Options()
    o2.read(opts_write)
    assert(object_eq(o2.to_obj(),ref))
#end def test_options



def test_job_init():
    from machines import Job,job
    from machines import job_defaults
    from machines import job_defaults_assign,job_defaults_nonassign

    jda  = set(job_defaults_assign.keys())
    jdna = set(job_defaults_nonassign.keys())
    jd   = set(job_defaults.keys())
    assert('nodes' in jda)
    assert('app' in jdna)
    assert(len(jda&jdna)==0)
    assert(jd==(jda|jdna))

    assert(id(job)==id(Job))

    # empty init should fail w/o implicit or explicit machine
    try:
        job()
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try

    # empty init should succeed if machine is bypassed
    j = job(skip_machine=True)

    assert(len(jda-set(j.keys()))==0)
    assert(len(j.app_props)==0)
    j.app_props = None # set back to default
    assert(object_eq(j.obj(*jda),job_defaults_assign))
#end def test_job_init



def test_job_time():
    from generic import obj
    from machines import job,Job

    def check_time(j,d,h,m,s,td,th,tm,ts):
        t = obj(days=d,hours=h,minutes=m,seconds=s)
        assert(object_eq(j.get_time(),t))
        assert(j.total_seconds()==ts)
        assert(j.total_minutes()==tm)
        assert(j.total_hours()==th)
        assert(j.total_days()==td)
    #end def check_time

    sm = dict(skip_machine=True)

    j1 = job(**sm)
    j2 = job(hours=4,minutes=30,**sm)
    j3 = job(days=2,hours=8,**sm)
    j4 = job(hours=53,minutes=127,**sm)

    check_time(j1, 0, 0,  0, 0, 0,  0,    0,      0)
    check_time(j2, 0, 4, 30, 0, 0,  4,  270,  16200)
    check_time(j3, 2, 8,  0, 0, 2, 56, 3360, 201600)
    check_time(j4, 2, 7,  7, 0, 2, 55, 3307, 198420)

    # pbs walltime
    assert(j1.pbs_walltime()=='00:00:00')
    assert(j2.pbs_walltime()=='04:30:00')
    assert(j3.pbs_walltime()=='2:08:00:00')
    assert(j4.pbs_walltime()=='2:07:07:00')

    # sbatch walltime
    assert(j1.sbatch_walltime()=='00:00:00')
    assert(j2.sbatch_walltime()=='04:30:00')
    assert(j3.sbatch_walltime()=='56:00:00')
    assert(j4.sbatch_walltime()=='55:07:00')

    # ll walltime
    assert(j1.ll_walltime()=='00:00:00')
    assert(j2.ll_walltime()=='04:30:00')
    assert(j3.ll_walltime()=='56:00:00')
    assert(j4.ll_walltime()=='55:07:00')

    # lsf walsftime
    assert(j1.lsf_walltime()=='00:00')
    assert(j2.lsf_walltime()=='04:30')
    assert(j3.lsf_walltime()=='56:00')
    assert(j4.lsf_walltime()=='55:07')

    t = Job.zero_time()
    t = j1.max_time(t)
    t = j2.max_time(t)
    t = j3.max_time(t)
    t = j4.max_time(t)
    assert(object_eq(t,j3.get_time()))

#end def test_job_time



def test_job_set_id():
    from machines import job

    j1 = job(skip_machine=True)
    j2 = job(skip_machine=True)

    j1.set_id()
    j2.set_id()

    assert(isinstance(j1.internal_id,int))
    assert(isinstance(j2.internal_id,int))
    assert(j2.internal_id-j1.internal_id==1)
#end def test_job_set_id



def test_job_get_machine():
    from machines import job

    machines = get_all_machines()
    for m in machines:
        j = job(machine=m.name,skip_machine=True)
        mj = j.get_machine()
        assert(id(mj)==id(m))
    #end for

#end def test_job_get_machine



def test_job_set_environment():
    from machines import job

    workstations,supercomputers = get_machine_data()

    machines = []
    machines.append(workstations.first())
    machines.append(supercomputers.first())

    for m in machines:
        j = job(machine=m.name,skip_machine=True)
        j.set_environment(OMP_NUM_THREADS=1)
        assert(j.env['OMP_NUM_THREADS']=='1')
    #end for

#end def test_job_set_environment



def test_job_clone():
    from machines import job

    j1 = job(skip_machine=True)
    j1.set_id()

    j2 = j1.clone()
    assert(id(j2)!=id(j1))
    assert(j2.internal_id-j1.internal_id==1)
    del j1.internal_id
    del j2.internal_id
    assert(object_eq(j2,j1))
#end def_test_job_clone



def test_job_serial_clone():
    from machines import job

    j1 = job(skip_machine=True)
    assert(not j1.serial)

    j2 = j1.serial_clone()
    assert(j2.serial)
    assert(j2.cores==1)
    assert(id(j2)!=id(j1))
    keys = 'serial cores init_info'.split()
    j1.delete(keys)
    j2.delete(keys)
    assert(object_eq(j2,j1))
#end def test_job_serial_clone



def test_machine_virtuals():
    from machines import Machine
    arg0 = None
    arg1 = None
    try:
        Machine.query_queue(arg0)
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try
    try:
        Machine.submit_jobs(arg0)
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try
    try:
        Machine.process_job(arg0,arg1)
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try
    try:
        Machine.process_job_options(arg0,arg1)
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try
    try:
        Machine.write_job(arg0,arg1,file=False)
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try
    try:
        Machine.submit_job(arg0,arg1)
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try
#end def test_machine_virtuals



def test_machine_list():
    from machines import Machine

    assert(len(Machine.machines)>0)
    for m in Machine.machines:
        assert(isinstance(m,Machine))
        exists = m.name in Machine.machines
        assert(exists)
        assert(Machine.exists(m.name))
        assert(Machine.is_unique(m))
        m.validate()
    #end for
#end def test_machine_list



def test_machine_add():
    from machines import Machine
    mtest = Machine.machines.first()
    assert(isinstance(mtest,Machine))
    try:
        Machine.add(mtest)
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try
    try:
        Machine.add('my_machine')
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try
#end def test_machine_add



def test_machine_get():
    from machines import Machine
    mtest = Machine.machines.first()
    assert(isinstance(mtest,Machine))

    m = Machine.get(mtest.name)
    assert(isinstance(m,Machine))
    assert(id(m)==id(mtest))
    try:
        Machine.get(m)
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try
    try:
        Machine.get('some_nonexistant_machine')
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try
#end def test_machine_get



def test_machine_instantiation():
    from machines import Machine
    # test guards against empty/invalid instantiation
    try:
        Machine()
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try
    try:
        Machine(123)
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try
    # test creation of a new machine
    test_name = 'test_machine'
    assert(not Machine.exists(test_name))
    m = Machine(name=test_name)
    assert(isinstance(m,Machine))
    assert(Machine.exists(m.name))
    assert(Machine.is_unique(m))
    m.validate()

    # test guards against multiple instantiation
    try:
        Machine(name=test_name)
        raise FailedTest
    except FailedTest:
        failed()
    except:
        None
    #end try

    # remove test machine
    del Machine.machines.test_machine
    assert(not Machine.exists(test_name))
#end def test_machine_instantiation



def test_workstation_init():
    from generic import obj
    from machines import Workstation

    ws = Workstation('wsi',16,'mpirun')

    refws = obj(
        account         = None,
        app_directories = None,
        app_directory   = None,
        app_launcher    = 'mpirun',
        cores           = 16,
        finished        = set([]),
        local_directory = None,
        name            = 'wsi',
        process_granularity = 1,
        queue_size      = 16,
        running         = set([]),
        user            = None,
        waiting         = set([]),
        jobs            = obj(),
        processes       = obj(),
        )

    assert(object_eq(ws.to_obj(),refws))

#end def test_workstation_init



# imitate Job.initialize w/o involving simulation object
def init_job(j,
             id  = 'job_ident',
             dir = './',
             ):
    import os
    identifier = id
    directory  = dir
    j.set_id()
    j.identifier  = identifier
    j.directory   = directory
    j.abs_dir     = os.path.abspath(directory)
    j.subdir      = j.directory
    j.abs_subdir  = j.abs_dir
    j.name        = identifier
    j.simid       = j.internal_id+1000
    j.outfile     = identifier+'.out'
    j.errfile     = identifier+'.err'
    j.app_command = 'echo run'
    j.process()
#end def init_job



def test_workstation_scheduling():
    import time
    from machines import Workstation
    from machines import job,Job

    tpath = testing.setup_unit_test_output_directory('machines','test_workstation_scheduling')

    # create workstation for testing
    ws = Workstation('wss',16,'mpirun')


    # test process_job(), process_job_options(), write_job()
    j = job(machine=ws.name,cores=4,threads=2,local=True)
    assert(j.machine==ws.name)
    assert(j.cores==4)
    assert(j.threads==2)
    assert(j.local)
    assert(j.processes==2)
    assert(j.run_options.np=='-np 2')
    assert(j.batch_mode==False)
    init_job(j,dir=tpath) # imitate interaction w/ simulation object
    assert(ws.write_job(j)=='export OMP_NUM_THREADS=2\nmpirun -np 2 echo run')

    j = job(machine=ws.name,serial=True)
    assert(j.machine==ws.name)
    assert(j.cores==1)
    assert(j.threads==1)
    assert(j.serial)
    assert(j.processes==1)
    assert(j.run_options.np=='-np 1')
    assert(j.batch_mode==False)
    init_job(j,dir=tpath) # imitate interaction w/ simulation object
    assert(ws.write_job(j)=='export OMP_NUM_THREADS=1\necho run')


    # test add_job()
    assert(j.status==Job.states.none)
    assert(not j.submitted)
    assert(len(ws.waiting)==0)

    j.submit() # calls ws.add_job(j)

    assert(j.status==Job.states.waiting)
    assert(j.submitted)
    assert(ws.waiting==set([j.internal_id]))
    assert(set(ws.jobs.keys())==set([j.internal_id]))
    assert(id(ws.jobs[j.internal_id])==id(j))


    # test submit_jobs() and submit_job()
    divert_nexus_log()
    assert(j.system_id is None)
    assert(len(ws.running)==0)
    assert(len(ws.processes)==0)

    ws.submit_jobs() # will call system echo

    assert(j.status==Job.states.running)
    assert(isinstance(j.system_id,int))
    assert(len(ws.waiting)==0)
    assert(ws.running==set([j.internal_id]))
    assert(set(ws.processes.keys())==set([j.system_id]))
    p = ws.processes[j.system_id]
    assert(p.popen.pid==j.system_id)
    assert(id(p.job)==id(j))
    assert(set(ws.jobs.keys())==set([j.internal_id]))
    restore_nexus_log()

    # allow a moment for all system calls to resolve
    time.sleep(0.1)

    # test query_queue()
    assert(not j.finished)

    ws.query_queue() # echo process complete

    assert(j.finished)
    assert(j.status==Job.states.finished)
    assert(len(ws.running)==0)
    assert(len(ws.processes)==0)
    assert(ws.finished==set([j.internal_id]))
    assert(set(ws.jobs.keys())==set([j.internal_id]))

#end def test_workstation_scheduling



def test_supercomputer_init():
    from generic import obj
    from machines import Theta

    class ThetaInit(Theta):
        name = 'theta_init'
    #end class ThetaInit

    sc = ThetaInit(4392,1,64,192,1000,'aprun','qsub','qstata','qdel')

    refsc = obj(
        account         = None,
        app_directories = None,
        app_directory   = None,
        app_launcher    = 'aprun',
        cores           = 281088,
        cores_per_node  = 64,
        cores_per_proc  = 64,
        finished        = set([]),
        job_remover     = 'qdel',
        local_directory = None,
        name            = 'theta_init',
        nodes           = 4392,
        procs           = 4392,
        procs_per_node  = 1,
        queue_querier   = 'qstata',
        queue_size      = 1000,
        ram             = 843264,
        ram_per_node    = 192,
        running         = set([]),
        sub_launcher    = 'qsub',
        user            = None,
        waiting         = set([]),
        jobs            = obj(),
        processes       = obj(),
        system_queue    = obj(),
        )

    assert(object_eq(sc.to_obj(),refsc))

#end def test_supercomputer_init



def test_supercomputer_scheduling():
    import os
    import time
    from generic import obj
    from machines import Theta
    from machines import job,Job

    tpath = testing.setup_unit_test_output_directory('machines','test_supercomputer_scheduling')

    # create supercomputer for testing
    class ThetaSched(Theta):
        name = 'theta_sched'
    #end class ThetaSched

    sc = ThetaSched(4392,1,64,192,1000,'aprun','echo','test_query','qdel')


    # test process_job() and process_job_options()
    j = job(machine=sc.name,nodes=2,threads=8,hours=3,minutes=30,account='ABC123')
    assert(j.machine==sc.name)
    assert(j.nodes==2)
    assert(j.threads==8)
    assert(j.processes==16)
    assert(j.processes_per_node==8)
    assert(j.cores==128)
    assert(j.hours==3)
    assert(j.minutes==30)
    assert(j.account=='ABC123')
    refro = obj(
        N               = '-N 8',
        cc              = '-cc depth',
        d               = '-d 8',
        e               = '-e OMP_NUM_THREADS=8',
        j               = '-j 1',
        n               = '-n 16',
        )
    assert(object_eq(j.run_options.to_obj(),refro))
    assert(j.batch_mode==True)


    # test write_job()
    init_job(j,id='123',dir=tpath) # imitate interaction w/ simulation object
    ref_wj = '''#!/bin/bash
#COBALT -q default
#COBALT -A ABC123
#COBALT -n 2
#COBALT -t 210
#COBALT -O 123
#COBALT --attrs mcdram=cache:numa=quad

export OMP_NUM_THREADS=8
aprun -e OMP_NUM_THREADS=8 -d 8 -cc depth -j 1 -n 16 -N 8 echo run'''
    wj = sc.write_job(j)
    for flag in refro:
        assert(flag in wj)
    #end for
    assert('aprun ' in wj)
    assert(' echo run' in wj)
    def scomp(s):
        return s.strip().rsplit('\n',1)[0].strip()
    #end if
    assert(scomp(wj)==scomp(ref_wj))


    # test add_job()
    assert(j.status==Job.states.none)
    assert(not j.submitted)
    assert(len(sc.waiting)==0)

    j.submit() # calls sc.add_job(j)

    assert(j.status==Job.states.waiting)
    assert(j.submitted)
    assert(sc.waiting==set([j.internal_id]))
    assert(set(sc.jobs.keys())==set([j.internal_id]))
    assert(id(sc.jobs[j.internal_id])==id(j))

    
    # test write_job() to file
    sc.write_job(j,file=True)

    subfile_path = os.path.join(tpath,j.subfile)
    assert(os.path.exists(subfile_path))
    wj = open(subfile_path,'r').read().strip()
    assert('aprun ' in wj)
    assert(' echo run' in wj)
    def scomp(s):
        return s.strip().rsplit('\n',1)[0].strip()
    #end if
    assert(scomp(wj)==scomp(ref_wj))


    # test sub_command()
    assert(sc.sub_command(j)=='echo 123.echo.in')


    # test submit_jobs() and submit_job()
    divert_nexus_log()
    assert(j.system_id is None)

    sc.submit_jobs() # will call system echo

    assert(j.status==Job.states.running)
    assert(j.system_id==123)
    assert(len(sc.waiting)==0)
    assert(sc.running==set([j.internal_id]))
    assert(set(sc.processes.keys())==set([123]))
    assert(set(sc.jobs.keys())==set([j.internal_id]))
    restore_nexus_log()


    # allow a moment for all system calls to resolve
    time.sleep(0.1)


    # test query_queue()
    assert(not j.finished)

    sc.query_queue() # echo process complete

    assert(j.finished)
    assert(j.status==Job.states.finished)
    assert(len(sc.running)==0)
    assert(len(sc.processes)==0)
    assert(sc.finished==set([j.internal_id]))
    assert(set(sc.jobs.keys())==set([j.internal_id]))

#end def test_supercomputer_scheduling



def test_process_job():
    from random import randint
    from generic import obj
    from machines import Machine,Job

    nw  = 5
    nwj = 5
    nsj = 5
    nij = 5

    workstations,supercomputers = get_machine_data()

    allow_warn = Machine.allow_warnings
    Machine.allow_warnings = False

    not_idempotent = obj()

    # check workstations
    nworkstations = nw
    if nworkstations is None:
        nworkstations=len(workstations)
    #end if
    nworkstations = min(nworkstations,len(workstations))
    njobs = nwj
    for nm in range(nworkstations):
        if nworkstations<len(workstations):
            machine = workstations.select_random() # select machine at random
        else:
            machine = workstations[nm]
        #end if
        cores_min     = 1
        cores_max     = machine.cores
        processes_min = 1
        processes_max = machine.cores
        threads_min   = 1
        threads_max   = machine.cores
        job_inputs = []
        job_inputs_base = []
        for nj in range(njobs): # vary cores
            cores   = randint(cores_min,cores_max)
            threads = randint(threads_min,threads_max)
            job_inputs_base.append(obj(cores=cores,threads=threads))
        #end for
        for nj in range(njobs): # vary processes
            processes   = randint(processes_min,processes_max)
            threads = randint(threads_min,threads_max)
            job_inputs_base.append(obj(processes=processes,threads=threads))
        #end for
        job_inputs.extend(job_inputs_base)
        for job_input in job_inputs_base: # run in serial
            ji = job_input.copy()
            ji.serial = True
            job_inputs.append(ji)
        #end for
        # perform idempotency test
        machine_idempotent = True
        for job_input in job_inputs:
            job = Job(machine=machine.name,**job_input)
            job2 = obj.copy(job)
            machine.process_job(job2)
            machine_idempotent &= job==job2
        #end for
        if not machine_idempotent:
            not_idempotent[machine.name] = machine
        #end if
    #end for

    # check supercomputers
    njobs = nsj
    small_node_ceiling = 20
    nodes_min   = 1
    cores_min   = 1
    threads_min = 1
    shared_job_inputs = obj(name='some_job',account='some_account')
    for machine in supercomputers:
        job_inputs = []
        job_inputs_base = []
        threads_max = 2*machine.cores_per_node
        # sample small number of nodes more heavily
        nodes_max   = min(small_node_ceiling,machine.nodes)
        cores_max   = min(small_node_ceiling*machine.cores_per_node,machine.cores)
        for nj in range(njobs): # nodes alone
            nodes   = randint(nodes_min,nodes_max)
            threads = randint(threads_min,threads_max)
            job_input = obj(nodes=nodes,threads=threads,**shared_job_inputs)
            job_inputs_base.append(job_input)
        #end for
        for nj in range(njobs): # cores alone
            cores   = randint(cores_min,cores_max)
            threads = randint(threads_min,threads_max)
            job_input = obj(cores=cores,threads=threads,**shared_job_inputs)
            job_inputs_base.append(job_input)
        #end for
        for nj in range(njobs): # nodes and cores
            nodes   = randint(nodes_min,nodes_max)
            cores   = randint(cores_min,cores_max)
            threads = randint(threads_min,threads_max)
            job_input = obj(nodes=nodes,cores=cores,threads=threads,**shared_job_inputs)
            job_inputs_base.append(job_input)
        #end for
        # sample full node set
        nodes_max = machine.nodes
        cores_max = machine.cores
        for nj in range(njobs): # nodes alone
            nodes   = randint(nodes_min,nodes_max)
            threads = randint(threads_min,threads_max)
            job_input = obj(nodes=nodes,threads=threads,**shared_job_inputs)
            job_inputs_base.append(job_input)
        #end for
        for nj in range(njobs): # cores alone
            cores   = randint(cores_min,cores_max)
            threads = randint(threads_min,threads_max)
            job_input = obj(cores=cores,threads=threads,**shared_job_inputs)
            job_inputs_base.append(job_input)
        #end for
        for nj in range(njobs): # nodes and cores
            nodes   = randint(nodes_min,nodes_max)
            cores   = randint(cores_min,cores_max)
            threads = randint(threads_min,threads_max)
            job_input = obj(nodes=nodes,cores=cores,threads=threads,**shared_job_inputs)
            job_inputs_base.append(job_input)
        #end for
        job_inputs.extend(job_inputs_base)
        # now add serial jobs
        for job_input in job_inputs_base:
            ji = job_input.copy()
            ji.serial = True
            job_inputs.append(ji)
        #end for
        # now add local, serial jobs
        for job_input in job_inputs_base:
            ji = job_input.copy()
            ji.serial = True
            ji.local  = True
            job_inputs.append(ji)
        #end for
        # perform idempotency test
        machine_idempotent = True
        for job_input in job_inputs:
            job = Job(machine=machine.name,**job_input)
            assert(isinstance(job.processes,int))
            assert(isinstance(job.nodes,int))
            if job.processes_per_node is not None:
                assert(job.processes==job.nodes*job.processes_per_node)
            #end if
            if not job.serial:
                nodes_input = 'nodes' in job_input
                cores_input = 'cores' in job_input
                if job.processes%job.nodes==0:
                    assert(job.processes_per_node is not None)
                    assert(job.processes==job.nodes*job.processes_per_node)
                #end if
                if nodes_input:
                    assert(job.nodes==job_input.nodes)
                #end if
                if nodes_input and cores_input:
                    assert(job.cores<=job_input.cores)
                elif nodes_input:
                    assert(job.cores%machine.cores_per_node==0)
                    assert(job.cores==job.nodes*machine.cores_per_node)
                elif cores_input:
                    assert(job.cores==job_input.cores)
                #end if
            #end if
            job2 = obj.copy(job)
            machine.process_job(job2)
            job_idempotent = object_eq(job,job2)
            if not job_idempotent:
                d,d1,d2 = object_diff(job,job2,full=True)
                change = obj(job_before=obj(d1),job_after=obj(d2))
                msg = machine.name+'\n'+str(change)
                failed(msg)
            #end if
            machine_idempotent &= job_idempotent
        #end for
        if not machine_idempotent:
            not_idempotent[machine.name] = machine
        #end if
    #end for

    if len(not_idempotent)>0:
        mlist = ''
        for name in sorted(not_idempotent.keys()):
            mlist+= '\n  '+name
        #end for
        msg='\n\nsome machines failed process_job idempotency test:{0}'.format(mlist)
        failed(msg)
    #end if
    Machine.allow_warnings = allow_warn

#end def test_process_job



def test_job_run_command():
    from generic import obj
    from machines import Machine,Job

    workstations,supercomputers = get_machine_data()

    allow_warn = Machine.allow_warnings
    Machine.allow_warnings = False

    def parse_job_command(command):
        tokens = command.split()
        #tokens = command.replace(':',' ').split()
        launcher = tokens[0]
        exe = tokens[-1]
        args = []
        options = obj()
        last_option = None
        for t in tokens[1:-1]:
            if t.startswith('-'):
                options[t] = None
                last_option = t
            elif last_option is not None:
                if options[last_option] is None:
                    options[last_option] = {t}
                else:
                    options[last_option].add(t)
                #end if
                if last_option!='--envs':
                    last_option = None
                #end if
            else:
                args.append(t)
            #end if
        #end for
        jc = obj(
            launcher   = launcher,
            executable = exe,
            args       = args,
            options    = options,
            )
        return jc
    #end def parse_job_command

    def job_commands_equal(c1,c2):
        jc1 = parse_job_command(c1)
        jc2 = parse_job_command(c2)
        return object_eq(jc1,jc2)
    #end def job_command_equal

    job_run_ref = obj({
        ('amos'           , 'n1'            ) : 'srun test.x',
        ('amos'           , 'n1_p1'         ) : 'srun test.x',
        ('amos'           , 'n2'            ) : 'srun test.x',
        ('amos'           , 'n2_t2'         ) : 'srun test.x',
        ('amos'           , 'n2_t2_e'       ) : 'srun test.x',
        ('amos'           , 'n2_t2_p2'      ) : 'srun test.x',
        ('andes'          , 'n1'            ) : 'srun -N 1 -n 32 test.x',
        ('andes'          , 'n1_p1'         ) : 'srun -N 1 -n 1 test.x',
        ('andes'          , 'n2'            ) : 'srun -N 2 -n 64 test.x',
        ('andes'          , 'n2_t2'         ) : 'srun -N 2 -c 2 --cpu-bind=cores -n 32 test.x',
        ('andes'          , 'n2_t2_e'       ) : 'srun -N 2 -c 2 --cpu-bind=cores -n 32 test.x',
        ('andes'          , 'n2_t2_p2'      ) : 'srun -N 2 -c 2 --cpu-bind=cores -n 4 test.x',
        ('archer2'        , 'n1'            ) : 'srun --distribution=block:block --hint=nomultithread -N 1 -n 128 test.x',
        ('archer2'        , 'n1_p1'         ) : 'srun --distribution=block:block --hint=nomultithread -N 1 -n 1 test.x',
        ('archer2'        , 'n2'            ) : 'srun --distribution=block:block --hint=nomultithread -N 2 -n 256 test.x',
        ('archer2'        , 'n2_t2'         ) : 'srun --distribution=block:block --hint=nomultithread -N 2 -c 2 -n 128 test.x',
        ('archer2'        , 'n2_t2_e'       ) : 'srun --distribution=block:block --hint=nomultithread -N 2 -c 2 -n 128 test.x',
        ('archer2'        , 'n2_t2_p2'      ) : 'srun --distribution=block:block --hint=nomultithread -N 2 -c 2 -n 4 test.x',
        ('attaway'        , 'n1'            ) : 'srun test.x',
        ('attaway'        , 'n1_p1'         ) : 'srun test.x',
        ('attaway'        , 'n2'            ) : 'srun test.x',
        ('attaway'        , 'n2_t2'         ) : 'srun test.x',
        ('attaway'        , 'n2_t2_e'       ) : 'srun test.x',
        ('attaway'        , 'n2_t2_p2'      ) : 'srun test.x',
        ('bluewaters_xe'  , 'n1'            ) : 'aprun -n 32 test.x',
        ('bluewaters_xe'  , 'n1_p1'         ) : 'aprun -n 1 test.x',
        ('bluewaters_xe'  , 'n2'            ) : 'aprun -n 64 test.x',
        ('bluewaters_xe'  , 'n2_t2'         ) : 'aprun -d 2 -n 32 test.x',
        ('bluewaters_xe'  , 'n2_t2_e'       ) : 'aprun -d 2 -n 32 test.x',
        ('bluewaters_xe'  , 'n2_t2_p2'      ) : 'aprun -d 2 -n 4 test.x',
        ('bluewaters_xk'  , 'n1'            ) : 'aprun -n 16 test.x',
        ('bluewaters_xk'  , 'n1_p1'         ) : 'aprun -n 1 test.x',
        ('bluewaters_xk'  , 'n2'            ) : 'aprun -n 32 test.x',
        ('bluewaters_xk'  , 'n2_t2'         ) : 'aprun -d 2 -n 16 test.x',
        ('bluewaters_xk'  , 'n2_t2_e'       ) : 'aprun -d 2 -n 16 test.x',
        ('bluewaters_xk'  , 'n2_t2_p2'      ) : 'aprun -d 2 -n 4 test.x',
        ('cades'          , 'n1'            ) : 'mpirun -np 36 test.x',
        ('cades'          , 'n1_p1'         ) : 'mpirun -np 1 test.x',
        ('cades'          , 'n2'            ) : 'mpirun -np 72 test.x',
        ('cades'          , 'n2_t2'         ) : 'mpirun -np 36 test.x',
        ('cades'          , 'n2_t2_e'       ) : 'mpirun -np 36 test.x',
        ('cades'          , 'n2_t2_p2'      ) : 'mpirun -np 4 test.x',
        ('cades_moab'     , 'n1'            ) : 'mpirun -np 36 test.x',
        ('cades_moab'     , 'n1_p1'         ) : 'mpirun -np 1 test.x',
        ('cades_moab'     , 'n2'            ) : 'mpirun -np 72 test.x',
        ('cades_moab'     , 'n2_t2'         ) : 'mpirun --npersocket 9 -np 36 test.x',
        ('cades_moab'     , 'n2_t2_e'       ) : 'mpirun --npersocket 9 -np 36 test.x',
        ('cades_moab'     , 'n2_t2_p2'      ) : 'mpirun --npersocket 1 -np 4 test.x',
        ('cetus'          , 'n1'            ) : 'runjob --envs OMP_NUM_THREADS=1 --np 16 -p 16 --verbose=INFO $LOCARGS : test.x',
        ('cetus'          , 'n1_p1'         ) : 'runjob --envs OMP_NUM_THREADS=1 --np 1 -p 1 --verbose=INFO $LOCARGS : test.x',
        ('cetus'          , 'n2'            ) : 'runjob --envs OMP_NUM_THREADS=1 --np 32 -p 16 --verbose=INFO $LOCARGS : test.x',
        ('cetus'          , 'n2_t2'         ) : 'runjob --envs OMP_NUM_THREADS=2 --np 16 -p 8 --verbose=INFO $LOCARGS : test.x',
        ('cetus'          , 'n2_t2_e'       ) : 'runjob --envs OMP_NUM_THREADS=2 ENV_VAR=1 --np 16 -p 8 --verbose=INFO $LOCARGS : test.x',
        ('cetus'          , 'n2_t2_p2'      ) : 'runjob --envs OMP_NUM_THREADS=2 --np 4 -p 2 --verbose=INFO $LOCARGS : test.x',
        ('chama'          , 'n1'            ) : 'srun test.x',
        ('chama'          , 'n1_p1'         ) : 'srun test.x',
        ('chama'          , 'n2'            ) : 'srun test.x',
        ('chama'          , 'n2_t2'         ) : 'srun test.x',
        ('chama'          , 'n2_t2_e'       ) : 'srun test.x',
        ('chama'          , 'n2_t2_p2'      ) : 'srun test.x',
        ('cooley'         , 'n1'            ) : 'mpirun -np 12 test.x',
        ('cooley'         , 'n1_p1'         ) : 'mpirun -np 1 test.x',
        ('cooley'         , 'n2'            ) : 'mpirun -np 24 test.x',
        ('cooley'         , 'n2_t2'         ) : 'mpirun -np 12 test.x',
        ('cooley'         , 'n2_t2_e'       ) : 'mpirun -np 12 test.x',
        ('cooley'         , 'n2_t2_p2'      ) : 'mpirun -np 4 test.x',
        ('cori'           , 'n1'            ) : 'srun test.x',
        ('cori'           , 'n1_p1'         ) : 'srun test.x',
        ('cori'           , 'n2'            ) : 'srun test.x',
        ('cori'           , 'n2_t2'         ) : 'srun test.x',
        ('cori'           , 'n2_t2_e'       ) : 'srun test.x',
        ('cori'           , 'n2_t2_p2'      ) : 'srun test.x',
        ('eclipse'        , 'n1'            ) : 'srun test.x',
        ('eclipse'        , 'n1_p1'         ) : 'srun test.x',
        ('eclipse'        , 'n2'            ) : 'srun test.x',
        ('eclipse'        , 'n2_t2'         ) : 'srun test.x',
        ('eclipse'        , 'n2_t2_e'       ) : 'srun test.x',
        ('eclipse'        , 'n2_t2_p2'      ) : 'srun test.x',
        ('edison'         , 'n1'            ) : 'srun test.x',
        ('edison'         , 'n1_p1'         ) : 'srun test.x',
        ('edison'         , 'n2'            ) : 'srun test.x',
        ('edison'         , 'n2_t2'         ) : 'srun test.x',
        ('edison'         , 'n2_t2_e'       ) : 'srun test.x',
        ('edison'         , 'n2_t2_p2'      ) : 'srun test.x',
        ('eos'            , 'n1'            ) : 'aprun -n 16 test.x',
        ('eos'            , 'n1_p1'         ) : 'aprun -n 1 test.x',
        ('eos'            , 'n2'            ) : 'aprun -n 32 test.x',
        ('eos'            , 'n2_t2'         ) : 'aprun -ss -cc numa_node -d 2 -n 16 test.x',
        ('eos'            , 'n2_t2_e'       ) : 'aprun -ss -cc numa_node -d 2 -n 16 test.x',
        ('eos'            , 'n2_t2_p2'      ) : 'aprun -ss -cc numa_node -d 2 -n 4 test.x',
        ('jaguar'         , 'n1'            ) : 'aprun -n 16 test.x',
        ('jaguar'         , 'n1_p1'         ) : 'aprun -n 1 test.x',
        ('jaguar'         , 'n2'            ) : 'aprun -n 32 test.x',
        ('jaguar'         , 'n2_t2'         ) : 'aprun -d 2 -n 16 test.x',
        ('jaguar'         , 'n2_t2_e'       ) : 'aprun -d 2 -n 16 test.x',
        ('jaguar'         , 'n2_t2_p2'      ) : 'aprun -d 2 -n 4 test.x',
        ('komodo'         , 'n1'            ) : 'mpirun -np 12 test.x',
        ('komodo'         , 'n1_p1'         ) : 'mpirun -np 1 test.x',
        ('komodo'         , 'n2'            ) : 'mpirun -np 24 test.x',
        ('komodo'         , 'n2_t2'         ) : 'mpirun -np 12 test.x',
        ('komodo'         , 'n2_t2_e'       ) : 'mpirun -np 12 test.x',
        ('komodo'         , 'n2_t2_p2'      ) : 'mpirun -np 4 test.x',
        ('kraken'         , 'n1'            ) : 'aprun -n 12 test.x',
        ('kraken'         , 'n1_p1'         ) : 'aprun -n 1 test.x',
        ('kraken'         , 'n2'            ) : 'aprun -n 24 test.x',
        ('kraken'         , 'n2_t2'         ) : 'aprun -d 2 -n 12 test.x',
        ('kraken'         , 'n2_t2_e'       ) : 'aprun -d 2 -n 12 test.x',
        ('kraken'         , 'n2_t2_p2'      ) : 'aprun -d 2 -n 4 test.x',
        ('lonestar'       , 'n1'            ) : 'ibrun -n 12 -o 0 test.x',
        ('lonestar'       , 'n1_p1'         ) : 'ibrun -n 1 -o 0 test.x',
        ('lonestar'       , 'n2'            ) : 'ibrun -n 24 -o 0 test.x',
        ('lonestar'       , 'n2_t2'         ) : 'ibrun -n 12 -o 0 test.x',
        ('lonestar'       , 'n2_t2_e'       ) : 'ibrun -n 12 -o 0 test.x',
        ('lonestar'       , 'n2_t2_p2'      ) : 'ibrun -n 4 -o 0 test.x',
        ('matisse'        , 'n1'            ) : 'mpirun -np 16 test.x',
        ('matisse'        , 'n1_p1'         ) : 'mpirun -np 1 test.x',
        ('matisse'        , 'n2'            ) : 'mpirun -np 32 test.x',
        ('matisse'        , 'n2_t2'         ) : 'mpirun -np 16 test.x',
        ('matisse'        , 'n2_t2_e'       ) : 'mpirun -np 16 test.x',
        ('matisse'        , 'n2_t2_p2'      ) : 'mpirun -np 4 test.x',
        ('mira'           , 'n1'            ) : 'runjob --envs OMP_NUM_THREADS=1 --np 16 -p 16 --verbose=INFO $LOCARGS : test.x',
        ('mira'           , 'n1_p1'         ) : 'runjob --envs OMP_NUM_THREADS=1 --np 1 -p 1 --verbose=INFO $LOCARGS : test.x',
        ('mira'           , 'n2'            ) : 'runjob --envs OMP_NUM_THREADS=1 --np 32 -p 16 --verbose=INFO $LOCARGS : test.x',
        ('mira'           , 'n2_t2'         ) : 'runjob --envs OMP_NUM_THREADS=2 --np 16 -p 8 --verbose=INFO $LOCARGS : test.x',
        ('mira'           , 'n2_t2_e'       ) : 'runjob --envs OMP_NUM_THREADS=2 ENV_VAR=1 --np 16 -p 8 --verbose=INFO $LOCARGS : test.x',
        ('mira'           , 'n2_t2_p2'      ) : 'runjob --envs OMP_NUM_THREADS=2 --np 4 -p 2 --verbose=INFO $LOCARGS : test.x',
        ('oic5'           , 'n1'            ) : 'mpirun -np 32 test.x',
        ('oic5'           , 'n1_p1'         ) : 'mpirun -np 1 test.x',
        ('oic5'           , 'n2'            ) : 'mpirun -np 64 test.x',
        ('oic5'           , 'n2_t2'         ) : 'mpirun -np 32 test.x',
        ('oic5'           , 'n2_t2_e'       ) : 'mpirun -np 32 test.x',
        ('oic5'           , 'n2_t2_p2'      ) : 'mpirun -np 4 test.x',
        ('rhea'           , 'n1'            ) : 'srun -N 1 -n 16 test.x',
        ('rhea'           , 'n1_p1'         ) : 'srun -N 1 -n 1 test.x',
        ('rhea'           , 'n2'            ) : 'srun -N 2 -n 32 test.x',
        ('rhea'           , 'n2_t2'         ) : 'srun -N 2 -n 16 -c 2 --cpu-bind=cores test.x',
        ('rhea'           , 'n2_t2_e'       ) : 'srun -N 2 -n 16 -c 2 --cpu-bind=cores test.x',
        ('rhea'           , 'n2_t2_p2'      ) : 'srun -N 2 -n 4 -c 2 --cpu-bind=cores test.x',
        ('skybridge'      , 'n1'            ) : 'srun test.x',
        ('skybridge'      , 'n1_p1'         ) : 'srun test.x',
        ('skybridge'      , 'n2'            ) : 'srun test.x',
        ('skybridge'      , 'n2_t2'         ) : 'srun test.x',
        ('skybridge'      , 'n2_t2_e'       ) : 'srun test.x',
        ('skybridge'      , 'n2_t2_p2'      ) : 'srun test.x',
        ('solo'           , 'n1'            ) : 'srun test.x',
        ('solo'           , 'n1_p1'         ) : 'srun test.x',
        ('solo'           , 'n2'            ) : 'srun test.x',
        ('solo'           , 'n2_t2'         ) : 'srun test.x',
        ('solo'           , 'n2_t2_e'       ) : 'srun test.x',
        ('solo'           , 'n2_t2_p2'      ) : 'srun test.x',
        ('stampede2'      , 'n1'            ) : 'ibrun -n 68 -o 0 test.x',
        ('stampede2'      , 'n1_p1'         ) : 'ibrun -n 1 -o 0 test.x',
        ('stampede2'      , 'n2'            ) : 'ibrun -n 136 -o 0 test.x',
        ('stampede2'      , 'n2_t2'         ) : 'ibrun -n 68 -o 0 test.x',
        ('stampede2'      , 'n2_t2_e'       ) : 'ibrun -n 68 -o 0 test.x',
        ('stampede2'      , 'n2_t2_p2'      ) : 'ibrun -n 4 -o 0 test.x',
        ('summit'         , 'n1'            ) : 'jsrun -a 21 -r 2 -b rs -c 21 -d packed -n 2 -g 0 test.x',
        ('summit'         , 'n1_g6'         ) : 'jsrun -a 7 -r 6 -b rs -c 7 -d packed -n 6 -g 1 test.x',
        ('summit'         , 'n2'            ) : 'jsrun -a 21 -r 2 -b rs -c 21 -d packed -n 4 -g 0 test.x',
        ('summit'         , 'n2_g6'         ) : 'jsrun -a 7 -r 6 -b rs -c 7 -d packed -n 12 -g 1 test.x',
        ('summit'         , 'n2_t2'         ) : 'jsrun -a 10 -r 2 -b rs -c 20 -d packed -n 4 -g 0 test.x',
        ('summit'         , 'n2_t2_e'       ) : 'jsrun -a 10 -r 2 -b rs -c 20 -d packed -n 4 -g 0 test.x',
        ('summit'         , 'n2_t2_e_g6'    ) : 'jsrun -a 3 -r 6 -b rs -c 6 -d packed -n 12 -g 1 test.x',
        ('summit'         , 'n2_t2_g6'      ) : 'jsrun -a 3 -r 6 -b rs -c 6 -d packed -n 12 -g 1 test.x',
        ('supermuc'       , 'n1'            ) : 'mpiexec -n 28 test.x',
        ('supermuc'       , 'n1_p1'         ) : 'mpiexec -n 1 test.x',
        ('supermuc'       , 'n2'            ) : 'mpiexec -n 56 test.x',
        ('supermuc'       , 'n2_t2'         ) : 'mpiexec -n 28 test.x',
        ('supermuc'       , 'n2_t2_e'       ) : 'mpiexec -n 28 test.x',
        ('supermuc'       , 'n2_t2_p2'      ) : 'mpiexec -n 4 test.x',
        ('supermucng'     , 'n1'            ) : 'mpiexec -n 48 test.x',
        ('supermucng'     , 'n1_p1'         ) : 'mpiexec -n 1 test.x',
        ('supermucng'     , 'n2'            ) : 'mpiexec -n 96 test.x',
        ('supermucng'     , 'n2_t2'         ) : 'mpiexec -n 48 test.x',
        ('supermucng'     , 'n2_t2_e'       ) : 'mpiexec -n 48 test.x',
        ('supermucng'     , 'n2_t2_p2'      ) : 'mpiexec -n 4 test.x',
        ('golub'           , 'n1'            ) : 'mpirun -np 12 test.x',
        ('golub'           , 'n1_p1'         ) : 'mpirun -np 1 test.x',
        ('golub'           , 'n2'            ) : 'mpirun -np 24 test.x',
        ('golub'           , 'n2_t2'         ) : 'mpirun -np 12 test.x',
        ('golub'           , 'n2_t2_e'       ) : 'mpirun -np 12 test.x',
        ('golub'           , 'n2_t2_p2'      ) : 'mpirun -np 4 test.x',
        ('theta'          , 'n1'            ) : 'aprun -e OMP_NUM_THREADS=1 -d 1 -cc depth -j 1 -n 64 -N 64 test.x',
        ('theta'          , 'n1_p1'         ) : 'aprun -e OMP_NUM_THREADS=1 -d 1 -cc depth -j 1 -n 1 -N 1 test.x',
        ('theta'          , 'n2'            ) : 'aprun -e OMP_NUM_THREADS=1 -d 1 -cc depth -j 1 -n 128 -N 64 test.x',
        ('theta'          , 'n2_t2'         ) : 'aprun -e OMP_NUM_THREADS=2 -d 2 -cc depth -j 1 -n 64 -N 32 test.x',
        ('theta'          , 'n2_t2_e'       ) : 'aprun -e OMP_NUM_THREADS=2 -d 2 -cc depth -j 1 -n 64 -N 32 test.x',
        ('theta'          , 'n2_t2_p2'      ) : 'aprun -e OMP_NUM_THREADS=2 -d 2 -cc depth -j 1 -n 4 -N 2 test.x',
        ('titan'          , 'n1'            ) : 'aprun -n 16 test.x',
        ('titan'          , 'n1_p1'         ) : 'aprun -n 1 test.x',
        ('titan'          , 'n2'            ) : 'aprun -n 32 test.x',
        ('titan'          , 'n2_t2'         ) : 'aprun -d 2 -n 16 test.x',
        ('titan'          , 'n2_t2_e'       ) : 'aprun -d 2 -n 16 test.x',
        ('titan'          , 'n2_t2_p2'      ) : 'aprun -d 2 -n 4 test.x',
        ('tomcat3'        , 'n1'            ) : 'mpirun -np 64 test.x',
        ('tomcat3'        , 'n1_p1'         ) : 'mpirun -np 1 test.x',
        ('tomcat3'        , 'n2'            ) : 'mpirun -np 128 test.x',
        ('tomcat3'        , 'n2_t2'         ) : 'mpirun -np 64 test.x',
        ('tomcat3'        , 'n2_t2_e'       ) : 'mpirun -np 64 test.x',
        ('tomcat3'        , 'n2_t2_p2'      ) : 'mpirun -np 4 test.x',
        ('uno'            , 'n1'            ) : 'srun test.x',
        ('uno'            , 'n1_p1'         ) : 'srun test.x',
        ('uno'            , 'n2'            ) : 'srun test.x',
        ('uno'            , 'n2_t2'         ) : 'srun test.x',
        ('uno'            , 'n2_t2_e'       ) : 'srun test.x',
        ('uno'            , 'n2_t2_p2'      ) : 'srun test.x',
        ('vesta'          , 'n1'            ) : 'runjob --envs OMP_NUM_THREADS=1 --np 16 -p 16 --verbose=INFO $LOCARGS : test.x',
        ('vesta'          , 'n1_p1'         ) : 'runjob --envs OMP_NUM_THREADS=1 --np 1 -p 1 --verbose=INFO $LOCARGS : test.x',
        ('vesta'          , 'n2'            ) : 'runjob --envs OMP_NUM_THREADS=1 --np 32 -p 16 --verbose=INFO $LOCARGS : test.x',
        ('vesta'          , 'n2_t2'         ) : 'runjob --envs OMP_NUM_THREADS=2 --np 16 -p 8 --verbose=INFO $LOCARGS : test.x',
        ('vesta'          , 'n2_t2_e'       ) : 'runjob --envs OMP_NUM_THREADS=2 ENV_VAR=1 --np 16 -p 8 --verbose=INFO $LOCARGS : test.x',
        ('vesta'          , 'n2_t2_p2'      ) : 'runjob --envs OMP_NUM_THREADS=2 --np 4 -p 2 --verbose=INFO $LOCARGS : test.x',
        })

    if testing.global_data['job_ref_table']:
        print('\n\n')
    #end if
    job_inputs_orig = obj(
        n1        = obj(nodes=1),
        n1_p1     = obj(nodes=1,processes_per_node=1),
        n2        = obj(nodes=2),
        n2_t2     = obj(nodes=2,threads=2),
        n2_t2_p2  = obj(nodes=2,threads=2,processes_per_node=2),
        n2_t2_e   = obj(nodes=2,threads=2,env=obj(ENV_VAR=1)),
        )
    for name in sorted(supercomputers.keys()):
        m = supercomputers[name]
        if m.requires_account:
            acc = 'ABC123'
        else:
            acc = None
        #end if
        job_inputs = job_inputs_orig
        if name=='summit': # exceptional treatment for summit nodes
            job_inputs = job_inputs_orig.copy()
            jtypes = list(job_inputs.keys())
            for jtype in jtypes:
                if 'p' in jtype:
                    del job_inputs[jtype]
                else:
                    jcpu = job_inputs[jtype]
                    jcpu.gpus = 0
                    jgpu = jcpu.copy()
                    jgpu.gpus = 6
                    job_inputs[jtype+'_g6'] = jgpu
                #end if
            #end for
        #end if
        for jtype in sorted(job_inputs.keys()):
            job = Job(app_command = 'test.x',
                      machine     = name,
                      account     = acc,
                      **job_inputs[jtype]
                      )
            command = job.run_command()
            if testing.global_data['job_ref_table']:
                sname = "'{0}'".format(name)
                stype = "'{0}'".format(jtype)
                print("        ({0:<16} , {1:<16}) : '{2}',".format(sname,stype,command))
                continue
            #end if
            ref_command = job_run_ref[name,jtype]
            if not job_commands_equal(command,ref_command):
                failed('Job.run_command for machine "{0}" does not match the reference\njob inputs:\n{1}\nreference command: {2}\nincorrect command: {3}'.format(name,job_inputs[jtype],ref_command,command))
            #end for
        #end for
    #end for


    # test split_nodes
    for name in sorted(supercomputers.keys()):
        m = supercomputers[name]
        if m.app_launcher=='srun': # no slurm support yet
            continue
        #end if
        if name=='summit': # no summit support
            continue
        #end if
        if m.requires_account:
            acc = 'ABC123'
        else:
            acc = None
        #end if
        # make a 4 node job
        job = Job(app_command = 'test.x',
                  machine     = name,
                  account     = acc,
                  nodes       = 4,
                  threads     = m.cores_per_node,
                  )
        # split the job into 1 and 3 nodes
        job1,job2 = job.split_nodes(1)
        # get the split run commands
        rc  = job.run_command()
        rc1 = job1.run_command()
        rc2 = job2.run_command()
        ns  = ' {0} '.format(job.nodes)
        ns1 = ' {0} '.format(job1.nodes)
        ns2 = ' {0} '.format(job2.nodes)
        # verify that node count is in each command
        assert(ns  in rc )
        assert(ns1 in rc1)
        assert(ns2 in rc2)
        # verify that text on either side of node count 
        # agrees for original and split commands
        assert(len(rc1)==len(rc))
        assert(len(rc2)==len(rc))

        loc = rc.find(ns)+1
        assert(rc[loc]=='4')
        assert(rc1[loc]=='1')
        assert(rc2[loc]=='3')
        rcf  = rc[:loc]+' '+rc[loc+1:]
        rc1f = rc1[:loc]+' '+rc1[loc+1:]
        rc2f = rc2[:loc]+' '+rc2[loc+1:]

        assert(job_commands_equal(rcf,rc1f))
        assert(job_commands_equal(rcf,rc2f))
    #end for

    Machine.allow_warnings = allow_warn
#end def test_job_run_command



def test_write_job():
    from generic import obj
    from machines import Machine,job

    supercomputers = get_supercomputers()

    allow_warn = Machine.allow_warnings
    Machine.allow_warnings = False

    job_write_ref = dict(
        amos = '''#!/bin/bash -x
#SBATCH --export=ALL
#SBATCH -J None
#SBATCH -p debug
#SBATCH -o test.out
#SBATCH -e test.err
#SBATCH --nodes 2
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 01:00:00

export OMP_NUM_THREADS=1
export ENV_VAR=1
srun test.x''',
        andes = '''#!/bin/bash
#SBATCH --job-name jobname
#SBATCH --account=ABC123
#SBATCH -N 2
#SBATCH -t 06:30:00
#SBATCH -o test.out
#SBATCH -e test.err

cd $SLURM_SUBMIT_DIR

echo JobID : $SLURM_JOBID 
echo Number of nodes requested: $SLURM_JOB_NUM_NODES 
echo List of nodes assigned to the job: $SLURM_NODELIST 


export ENV_VAR=1
export OMP_NUM_THREADS=1
srun -N 2 -n 64 test.x''',
        archer2 = '''#!/bin/bash
#SBATCH --job-name jobname
#SBATCH --account=ABC123
#SBATCH -N 2
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH -t 06:30:00
#SBATCH -o test.out
#SBATCH -e test.err
#SBATCH --partition=standard
#SBATCH --qos=standard

echo JobID : $SLURM_JOBID
echo Number of nodes requested: $SLURM_JOB_NUM_NODES
echo List of nodes assigned to the job: $SLURM_NODELIST

export ENV_VAR=1
export OMP_NUM_THREADS=1

srun --distribution=block:block --hint=nomultithread -N 2 -n 256 test.x''',
        attaway = '''#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name jobname
#SBATCH --account=ABC123
#SBATCH -N 2
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH -t 06:30:00
#SBATCH -o test.out
#SBATCH -e test.err

export OMP_NUM_THREADS=1
export ENV_VAR=1
srun test.x''',
        bluewaters_xe = '''#!/bin/bash
#PBS -N jobname
#PBS -l walltime=06:30:00
#PBS -l nodes=2:ppn=32:xe
#PBS -o test.out
#PBS -e test.err
#PBS -V

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=1
export ENV_VAR=1
aprun -n 64 test.x''',
        bluewaters_xk = '''#!/bin/bash
#PBS -N jobname
#PBS -l walltime=06:30:00
#PBS -l nodes=2:ppn=16:xk
#PBS -o test.out
#PBS -e test.err
#PBS -V

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=1
export ENV_VAR=1
aprun -n 32 test.x''',
        cades = '''#!/bin/bash
#SBATCH -A ABC123
#SBATCH -p skylake
#SBATCH -J jobname
#SBATCH -t 06:30:00
#SBATCH -N 2
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH -o test.out
#SBATCH -e test.err
#SBATCH --exclusive
#SBATCH --export=ALL

export ENV_VAR=1
export OMP_NUM_THREADS=1
mpirun -np 72 test.x''',
        cades_moab = '''#!/bin/bash
#PBS -A ABC123
#PBS -W group_list=cades-ABC123
#PBS -q skylake
#PBS -N jobname
#PBS -o test.out
#PBS -e test.err
#PBS -l qos=std
#PBS -l walltime=06:30:00
#PBS -l nodes=2:ppn=36

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

export ENV_VAR=1
export OMP_NUM_THREADS=1
mpirun -np 72 test.x''',
        cetus = '''#!/bin/bash
#COBALT -q default
#COBALT -A ABC123
#COBALT -n 2
#COBALT -t 390
#COBALT -O None

LOCARGS="--block $COBALT_PARTNAME ${COBALT_CORNER:+--corner} $COBALT_CORNER ${COBALT_SHAPE:+--shape} $COBALT_SHAPE"
echo "Cobalt location args: $LOCARGS" >&2


runjob --np 32 -p 16 $LOCARGS --verbose=INFO --envs OMP_NUM_THREADS=1 ENV_VAR=1 : test.x''',
        chama = '''#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name jobname
#SBATCH --account=ABC123
#SBATCH -N 2
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 06:30:00
#SBATCH -o test.out
#SBATCH -e test.err

export OMP_NUM_THREADS=1
export ENV_VAR=1
srun test.x''',
        cooley = '''#!/bin/bash
#COBALT -q default
#COBALT -A ABC123
#COBALT -n 2
#COBALT -t 390
#COBALT -O None

export OMP_NUM_THREADS=1
export ENV_VAR=1
mpirun -np 24 test.x''',
        cori = '''#!/bin/bash
#SBATCH -p regular
#SBATCH -C knl
#SBATCH -J jobname
#SBATCH -t 06:30:00
#SBATCH -N 2
#SBATCH --tasks-per-node=68
#SBATCH --cpus-per-task=4
#SBATCH -o test.out
#SBATCH -e test.err
#SBATCH --export=ALL

echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=1
export ENV_VAR=1
srun test.x''',
        eclipse = '''#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name jobname
#SBATCH --account=ABC123
#SBATCH -N 2
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH -t 06:30:00
#SBATCH -o test.out
#SBATCH -e test.err

export OMP_NUM_THREADS=1
export ENV_VAR=1
srun test.x''',
        edison = '''#!/bin/bash
#SBATCH -p regular
#SBATCH -J jobname
#SBATCH -t 06:30:00
#SBATCH -N 2
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH -o test.out
#SBATCH -e test.err
#SBATCH --export=ALL

echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=1
export ENV_VAR=1
srun test.x''',
        eos = '''#!/bin/bash
#PBS -A ABC123
#PBS -q batch
#PBS -N jobname
#PBS -o test.out
#PBS -e test.err
#PBS -l walltime=06:30:00
#PBS -l nodes=2
#PBS -l gres=atlas1
#PBS -V

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=1
export ENV_VAR=1
aprun -n 32 test.x''',
        jaguar = '''#!/bin/bash
#PBS -A ABC123
#PBS -q batch
#PBS -N jobname
#PBS -l walltime=06:30:00
#PBS -l size=32
#PBS -l gres=widow2%widow3
#PBS -o test.out
#PBS -e test.err
#PBS -V

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
export MPICH_PTL_SEND_CREDITS=-1
export MPICH_MAX_SHORT_MSG_SIZE=1024
export MPICH_PTL_UNEX_EVENTS=800000
export MPICH_UNEX_BUFFER_SIZE=16M
export MPI_MSGS_PER_PROC=32768

export OMP_NUM_THREADS=1
export ENV_VAR=1
aprun -n 32 test.x''',
        komodo = '''#!/bin/bash -x
#SBATCH --export=ALL
#SBATCH -J None
#SBATCH -p defq
#SBATCH -o test.out
#SBATCH -e test.err
#SBATCH --nodes 2
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH -t 06:30:00

export OMP_NUM_THREADS=1
export ENV_VAR=1
mpirun -np 24 test.x''',
        kraken = '''#!/bin/bash
#PBS -A ABC123
#PBS -N jobname
#PBS -l walltime=06:30:00
#PBS -l size=24
#PBS -o test.out
#PBS -e test.err
#PBS -V

cd $PBS_O_WORKDIR
export MPICH_PTL_SEND_CREDITS=-1
export MPICH_MAX_SHORT_MSG_SIZE=1024
export MPICH_PTL_UNEX_EVENTS=800000
export MPICH_UNEX_BUFFER_SIZE=16M
export MPI_MSGS_PER_PROC=32768

export OMP_NUM_THREADS=1
export ENV_VAR=1
aprun -n 24 test.x''',
        lonestar = '''#!/bin/bash
#$ -q batch
#$ -N jobname
#$ -o test.out
#$ -e test.err
#$ -l h_rt=06:30:00
#$ -pe 12way 24
#$ -cwd
#$ -V

export OMP_NUM_THREADS=1
export ENV_VAR=1
ibrun -n 24 -o 0 test.x''',
        matisse = '''#!/bin/bash -x
#SBATCH --export=ALL
#SBATCH -J None
#SBATCH -p defq
#SBATCH -o test.out
#SBATCH -e test.err
#SBATCH --nodes 2
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 06:30:00

export OMP_NUM_THREADS=1
export ENV_VAR=1
mpirun -np 32 test.x''',
        mira = '''#!/bin/bash
#COBALT -q default
#COBALT -A ABC123
#COBALT -n 2
#COBALT -t 390
#COBALT -O None

LOCARGS="--block $COBALT_PARTNAME ${COBALT_CORNER:+--corner} $COBALT_CORNER ${COBALT_SHAPE:+--shape} $COBALT_SHAPE"
echo "Cobalt location args: $LOCARGS" >&2


runjob --np 32 -p 16 $LOCARGS --verbose=INFO --envs OMP_NUM_THREADS=1 ENV_VAR=1 : test.x''',
        oic5 = '''#!/bin/bash
#PBS -q mstqmc13q
#PBS -N jobname
#PBS -l walltime=06:30:00
#PBS -l nodes=2:ppn=32
#PBS -W x="NACCESSPOLICY:SINGLEJOB"
#PBS -o test.out
#PBS -e test.err
#PBS -V

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=1
export ENV_VAR=1
mpirun -np 64 test.x''',
        rhea = '''#!/bin/bash
#SBATCH --job-name jobname
#SBATCH --account=ABC123
#SBATCH -N 2
#SBATCH -t 06:30:00
#SBATCH -o test.out
#SBATCH -e test.err

cd $SLURM_SUBMIT_DIR

echo JobID : $SLURM_JOBID 
echo Number of nodes requested: $SLURM_JOB_NUM_NODES 
echo List of nodes assigned to the job: $SLURM_NODELIST 


export ENV_VAR=1
export OMP_NUM_THREADS=1
srun -N 2 -n 32 test.x''',
        skybridge = '''#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name jobname
#SBATCH --account=ABC123
#SBATCH -N 2
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 06:30:00
#SBATCH -o test.out
#SBATCH -e test.err

export OMP_NUM_THREADS=1
export ENV_VAR=1
srun test.x''',
        solo = '''#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name jobname
#SBATCH --account=ABC123
#SBATCH -N 2
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH -t 06:30:00
#SBATCH -o test.out
#SBATCH -e test.err

export OMP_NUM_THREADS=1
export ENV_VAR=1
srun test.x''',
        stampede2 = '''#!/bin/bash
#SBATCH --job-name jobname
#SBATCH --account=ABC123
#SBATCH -N 2
#SBATCH --ntasks-per-node=68
#SBATCH --cpus-per-task=1
#SBATCH -t 06:30:00
#SBATCH -o test.out
#SBATCH -e test.err
#SBATCH -p normal


export OMP_NUM_THREADS=1
export ENV_VAR=1
ibrun -n 136 -o 0 test.x''',
        summit = '''#!/bin/bash
#BSUB -P ABC123
#BSUB -J jobname
#BSUB -o test.out
#BSUB -e test.err
#BSUB -W 06:30
#BSUB -nnodes 2
#BSUB -alloc_flags "smt1"

export OMP_NUM_THREADS=1
export ENV_VAR=1
jsrun -a 7 -r 6 -b rs -c 7 -d packed -n 12 -g 1 test.x''',
        supermuc = '''#!/bin/bash
#@ job_name         = jobname
#@ job_type         = MPICH
#@ class            = general
#@ node             = 2
#@ island_count     = 1,1
#@ total_tasks      = 56
#@ wall_clock_limit = 06:30:00
#@ network.MPI      = sn_all,not_shared,us
#@ initialdir       = /path/on/supermuc
#@ output           = test.out
#@ error            = test.err
#@ energy_policy_tag = my_energy_tag
#@ minimize_time_to_solution = yes
#@ notification     = never
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
module unload mpi.ibm
module load mpi.intel

export OMP_NUM_THREADS=1
export ENV_VAR=1
mpiexec -n 56 test.x''',
        supermucng = '''#!/bin/bash
#SBATCH --account=ABC123
#SBATCH --partition=general
#SBATCH -J jobname
#SBATCH --time=06:30:00
#SBATCH -o ./test.out
#SBATCH -e ./test.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH -D ./
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env

export OMP_NUM_THREADS=1
export ENV_VAR=1
mpiexec -n 96 test.x''',
        golub = '''#PBS -q secondary
#PBS -N jobname
#PBS -l nodes=2:ppn=12
#PBS -l walltime=06:30:00
#PBS -e test.err
#PBS -o test.out
#PBS -V

cd ${PBS_O_WORKDIR}


export ENV_VAR=1
export OMP_NUM_THREADS=1
mpirun -np 24 test.x''',
        theta = '''#!/bin/bash
#COBALT -q default
#COBALT -A ABC123
#COBALT -n 2
#COBALT -t 390
#COBALT -O None
#COBALT --attrs mcdram=cache:numa=quad

export OMP_NUM_THREADS=1
export ENV_VAR=1
aprun -e OMP_NUM_THREADS=1 -d 1 -cc depth -j 1 -n 128 -N 64 test.x''',
        titan = '''#!/bin/bash
#PBS -A ABC123
#PBS -q batch
#PBS -N jobname
#PBS -o test.out
#PBS -e test.err
#PBS -l walltime=06:30:00
#PBS -l nodes=2
#PBS -l gres=atlas1
#PBS -V

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=1
export ENV_VAR=1
aprun -n 32 test.x''',
        tomcat3 = '''#!/bin/bash -l
#SBATCH -J jobname
#SBATCH -N 2
#SBATCH -t 06:30:00
#SBATCH -p tomcat
#. /home/rcohen/.bashrc
unalias cd; source /mnt/beegfs/intel/parallel_studio_xe_2019.3.062/bin/psxevars.sh
ulimit -a

export OMP_NUM_THREADS=1
export ENV_VAR=1
mpirun -np 128 test.x >test.out 2>test.err''',
        uno = '''#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name jobname
#SBATCH --account=ABC123
#SBATCH -N 2
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 06:30:00
#SBATCH -o test.out
#SBATCH -e test.err

export OMP_NUM_THREADS=1
export ENV_VAR=1
srun test.x''',
        vesta = '''#!/bin/bash
#COBALT -q default
#COBALT -A ABC123
#COBALT -n 2
#COBALT -t 390
#COBALT -O None

LOCARGS="--block $COBALT_PARTNAME ${COBALT_CORNER:+--corner} $COBALT_CORNER ${COBALT_SHAPE:+--shape} $COBALT_SHAPE"
echo "Cobalt location args: $LOCARGS" >&2


runjob --np 32 -p 16 $LOCARGS --verbose=INFO --envs OMP_NUM_THREADS=1 ENV_VAR=1 : test.x''',
        )

    def process_job_file(jf):
        jf = jf.strip()
        lines = []
        tokens = []
        for line in jf.splitlines():
            if len(line.strip())==0:
                continue
            #end if
            if line.startswith('export'):
                tokens.append(line)
            else:
                lines.append(line)
            #end if
        #end for
        tokens.extend(lines.pop().split())
        jo = obj(
            lines  = lines,
            tokens = set(tokens),
            )
        return jo
    #end def process_job_file

    def job_files_same(jf1,jf2):
        jf1 = process_job_file(jf1)
        jf2 = process_job_file(jf2)
        if not object_eq(jf1,jf2): print(f"compare --------------------\n * wj *\n{jf1}\n * ref_wj *\n{jf2}\n")
        return object_eq(jf1,jf2)
    #end def job_files_same

    if testing.global_data['job_ref_table']:
        print('\n\n')
    #end if
    job_inputs = obj(
        nodes       = 2,
        hours       = 6,
        minutes     = 30,
        env         = obj(ENV_VAR=1),
        identifier  = 'test',
        outfile     = 'test.out',
        errfile     = 'test.err',
        app_command = 'test.x',
        )
    for name in sorted(supercomputers.keys()):
        m = supercomputers[name]
        if m.requires_account:
            acc = 'ABC123'
        else:
            acc = None
        #end if
        ji = job_inputs.copy()
        if name=='summit': # exceptional treatment for summit nodes
            ji.gpus = 6
        #end if
        j = job(
            machine = name,
            account = acc,
            **ji
            )
        j.abs_dir = '/path/on/'+name
        wj = m.write_job(j)
        if testing.global_data['job_ref_table']:
            print("        {} = '''{}''',".format(name,wj.strip()))
            continue
        #end if
        ref_wj = job_write_ref[name]
        assert(job_files_same(wj,ref_wj))
    #end for

    Machine.allow_warnings = allow_warn
#end def test_write_job
