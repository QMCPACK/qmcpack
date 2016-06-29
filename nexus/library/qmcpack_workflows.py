
import os
from developer import obj,ci,error as dev_error,devlog
from pwscf import generate_pwscf
from qmcpack_converters import generate_pw2qmcpack
from qmcpack_input import generate_jastrow,loop,linear,cslinear,vmc,dmc
from qmcpack import generate_qmcpack


def error(msg,loc=None,exit=True,trace=True,indent='    ',logfile=devlog):
    header = 'qmcpack_workflows'
    if loc!=None:
        msg+='\nfunction location: {0}'.format(loc)
    #end if
    dev_error(msg,header,exit,trace,indent,logfile)
#end def error



defaults_version = 'v1'





jastrow_factor_defaults = obj(
    v1 = obj(
        J1       = True,
        J2       = True,
        J3       = False,
        system   = None,
        J1_size  = 10,
        J1_rcut  = None,
        J2_size  = 10,
        J2_rcut  = None,
        J2_init  = 'zero',
        J3_isize = 3,
        J3_esize = 3,
        J3_rcut  = 5.0,
        J1_rcut_open =  5.0,
        J2_rcut_open = 10.0,
        ),
    )


opt_sections_options = '''
    method cost cycles var_cycles defaults opt_calcs
    '''.split()

opt_sections_required = ['method','cost','cycles','var_cycles','defaults']
opt_sections_defaults = obj(
    v1 = obj(
        method     = 'linear',
        cost       = 'variance',
        cycles     = 12,
        var_cycles = 4,
        defaults   = defaults_version,
        ),
    )


opt_method_defaults = obj(
    linear = obj(
        mm = obj(
            samples           = 128000,  
            warmupsteps       = 25, 
            blocks            = 250,  
            steps             = 1, 
            substeps          = 20, 
            timestep          = 0.5, 
            usedrift          = True, 
            nonlocalpp        = True, 
            usebuffer         = True, 
            minmethod         = 'quartic',
            exp0              = -6,
            bigchange         = 10.0,
            alloweddifference = 1.0e-5, 
            stepsize          = 0.15, 
            nstabilizers      = 1, 
            ),
        yl = obj(
            walkers             = 256, 
            samples             = 655360,  
            warmupsteps         = 1, 
            blocks              = 40, 
            substeps            = 5,     
            stepsbetweensamples = 1, 
            timestep            = 1.0,  
            usedrift            = False, 
            nonlocalpp          = False,
            usebuffer           = False,
            minmethod           = 'quartic',
            gevmethod           = 'mixed',
            minwalkers          = 0.3, 
            maxweight           = 1e9, 
            stepsize            = 0.9, 
            ),
        v1 = obj(
            samples           = 204800,             
            warmupsteps       = 300,                
            blocks            = 100,                
            steps             = 1,                  
            substeps          = 10,                 
            timestep          = 0.3,
            usedrift          = False,                 
            nonlocalpp        = True,                
            usebuffer         = True,                
            minmethod         = 'quartic',            
            exp0              = -6,                 
            bigchange         = 10.0,               
            alloweddifference = 1e-05,              
            stepsize          = 0.15,               
            nstabilizers      = 1,                  
            ),
        ),
    cslinear = obj(
        ls = obj(
            warmupsteps       = 20,
            steps             = 5,
            usedrift          = True,
            timestep          = .8,
            nonlocalpp        = False,
            minmethod         = 'rescale',
            stepsize          = .4,
            beta              = .05,
            gevmethod         = 'mixed',
            alloweddifference = 1e-4,
            bigchange         = 9.,
            exp0              = -16,
            max_its           = 1,
            maxweight         = 1e9,
            minwalkers        = .5,
            nstabilizers      = 3,
            stabilizerscale   = 1,
            usebuffer         = False,
            ),
        jm = obj(
            warmupsteps       = 20,
            usedrift          = True,
            timestep          = .5,
            nonlocalpp        = True,
            minmethod         = 'quartic',
            stepsize          = .4,
            beta              = 0.0,
            gevmethod         = 'mixed',
            alloweddifference = 1.0e-4,
            bigchange         = 9.0,
            exp0              = -16,
            max_its           = 1,
            maxweight         = 1e9,
            minwalkers        = 0.5,
            nstabilizers      = 3,
            stabilizerscale   = 1.0,
            usebuffer         = True,
            )
        ),
    )


vmc_sections_options = [
    'walkers','warmupsteps','blocks','steps',
    'substeps','timestep','checkpoint',
    'J0_warmupsteps','J0_blocks','J0_steps',
    'test_warmupsteps','test_blocks','test_steps',
    ]

dmc_sections_options = [
    'walkers','warmupsteps','blocks','steps',
    'timestep','checkpoint',
    'vmc_samples','vmc_samplesperthread',
    'vmc_walkers','vmc_warmupsteps','vmc_blocks','vmc_steps',
    'vmc_substeps','vmc_timestep','vmc_checkpoint',
    'eq_dmc','eq_warmupsteps','eq_blocks','eq_steps','eq_timestep','eq_checkpoint',
    'J0_warmupsteps','J0_blocks','J0_steps','J0_checkpoint',
    'test_warmupsteps','test_blocks','test_steps',
    'ntimesteps','timestep_factor',
    ]
dmc_sections_required = ['nlmove']



scf_workflow = ['totmag_sys']
scf_required = ['job']
scf_defaults = obj(
    minimal = obj(
        identifier       = 'scf',
        input_type       = 'generic',
        nosym            = True,
        wf_collect       = True,
        use_folded       = True,
        nogamma          = True,
        totmag_sys       = False,
        ),
    v1 = obj(
        identifier       = 'scf',
        input_type       = 'generic',
        diagonalization  = 'david',
        electron_maxstep = 1000,
        conv_thr         = 1e-8,
        mixing_beta      = 0.2,
        occupations      = 'smearing',
        smearing         = 'fermi-dirac',
        degauss          = 0.0001,
        nosym            = True,
        wf_collect       = True,
        use_folded       = True,
        nogamma          = True,
        totmag_sys       = False,
        ),
    )

p2q_workflow = []
p2q_required = ['job']
p2q_defaults = obj(
    minimal = obj(
        identifier = 'p2q',
        write_psir = False,
        ),
    v1 = obj(
        identifier = 'p2q',
        write_psir = False,
        ),
    )

opt_workflow = ['J2_prod','J3_prod','J_defaults']+opt_sections_options
opt_required = ['job']
fixed_defaults = obj(
    J2_prod    = False,
    J3_prod    = False,
    J_defaults = defaults_version,
    )
opt_defaults = obj(
    minimal = obj(
        identifier   = 'opt',
        input_type   = 'basic',
        **fixed_defaults
        ),
    v1 = obj(
        identifier   = 'opt',
        input_type   = 'basic',
        **fixed_defaults
        ),
    )

vmc_workflow = ['J0_prod','J2_prod','J3_prod',
                'J0_test','J2_test','J3_test',
                ] + vmc_sections_options
vmc_required = ['job']
fixed_defaults = obj(
    J0_prod = False,
    J2_prod = False,
    J3_prod = False,
    J0_test = False,
    J2_test = False,
    J3_test = False,
    )
vmc_defaults = obj(
    minimal = obj(
        identifier   = 'vmc',
        input_type   = 'basic',
        **fixed_defaults
        ),
    v1 = obj(
        identifier       = 'vmc',
        input_type       = 'basic',
        walkers          = 1,
        warmupsteps      = 50,
        blocks           = 800,
        steps            = 10,
        substeps         = 3,
        timestep         = 0.3,
        checkpoint       = -1,
        test_warmupsteps = 10,
        test_blocks      = 20,
        test_steps       =  4,
        J0_warmupsteps   = 200,
        J0_blocks        = 800,
        J0_steps         = 100,
        **fixed_defaults
        ),
    )

dmc_workflow = ['J0_prod','J2_prod','J3_prod',
                'J0_test','J2_test','J3_test',
                ]
dmc_required = ['job']
fixed_defaults = obj(
    J0_prod  = False,
    J2_prod  = False,
    J3_prod  = False,
    J0_test  = False,
    J2_test  = False,
    J3_test  = False,
    tmoves   = False,
    locality = False,
    )
dmc_defaults = obj(
    minimal = obj(
        identifier   = 'qmc',
        input_type   = 'basic',
        **fixed_defaults
        ),
    v1 = obj(
        identifier   = 'qmc',
        input_type   = 'basic',
        **fixed_defaults
        ),
    )



def set_null_default(defaults):
    defaults.none = obj()
    defaults[None] = defaults.none
#end def set_null_default

set_null_default(scf_defaults)
set_null_default(p2q_defaults)
set_null_default(opt_defaults)
set_null_default(vmc_defaults)
set_null_default(dmc_defaults)
set_null_default(jastrow_factor_defaults)
set_null_default(opt_sections_defaults)
for method_defaults in opt_method_defaults:
    set_null_default(method_defaults)
#end for








def extract_kwargs(kwargs,required=None,optional=None,defaults=None,default_sources=None,workflow=None,require_empty=False,allow_none=False,contains_defaults=False,encapsulate=True,loc='extract_kwargs'):
    if allow_none and kwargs is None:
        kwargs = {}
    #end if
    if encapsulate:
        kwargs = dict(**kwargs)
    #end if
    kw = obj()
    if contains_defaults:
        if 'defaults' in kwargs:
            defaults = kwargs['defaults']
            del kwargs['defaults']
        #end if
    #end if
    if defaults is not None:
        if isinstance(defaults,str):
            if default_sources is None:
                error('sources for defaults not provided\ndefault set requested: {0}'.format(defaults),loc)
            elif defaults not in default_sources:
                error('invalid default request encountered\ninvalid default set: {0}\nvalid_options_are: {1}'.format(defaults,sorted(default_sources.keys())),loc)
            #end if
            defaults = default_sources[defaults].copy()
        #end if
        defaults = defaults.copy()
        for key,value in defaults.iteritems():
            if key in kwargs:
                value = kwargs[key]
                del kwargs[key]
            #end if
            kw[key] = value
        #end for
    #end if
    if required is not None:
        for key in required:
            if key in kwargs:
                kw[key] = kwargs[key]
                del kwargs[key]
            elif key in kw:
                pass # allow defaults to fill in for required keywords
            else:
                error('a required keyword argument is missing\nmissing keyword argument: {0}'.format(key),loc)
            #end if
        #end for
    #end if
    if optional is not None:
        for key in optional:
            if key in kwargs:
                kw[key] = kwargs[key]
                del kwargs[key]
            #end if
        #end for
    #end if
    if workflow is None:
        return kw
    else:
        wfkw = obj()
        for key in workflow:
            if key in kw:
                wfkw[key] = kw[key]
                del kw[key]
            elif key in kwargs:
                wfkw[key] = kwargs[key]
                del kwargs[key]
            #end if
        #end for
        return kw,wfkw
    #end if
    if require_empty:
        require_empty_kwargs(kwargs,loc)
    #end if
#end def extract_kwargs


def require_empty_kwargs(kwargs,loc='require_empty_kwargs'):
    if len(kwargs)>0:
        error('invalid/unprocessed keywords encountered\ninvalid keywords: {0}'.format(sorted(kwargs.keys())),loc)
    #end if
#end def require_empty_kwargs


def resolve_deps(name,sims,deps,loc='resolve_deps'):
    deplist = []
    missing = []
    for depname,depquant in deps:
        if depname in sims:
            deplist.append((sims[depname],depquant))
        else:
            missing.append(depname)
        #end if
    #end for
    if len(missing)>0:
        keywords = []
        for m in missing:
            if len(m)>=3:
                keywords.append(m[0:3]+'_inputs')
            #end if
        #end for
        error('workflow cannot be run\nsimulation "{0}" depends on other simulations that have not been requested\nmissing simulations: {1}\nthe user needs to provide more detailed input\nthis issue can likely be fixed by providing the following keywords: {2}'.format(name,sorted(missing),sorted(set(keywords))))
    #end if
    return deplist
#end def resolve_deps









def process_jastrow(J,system):
    if isinstance(J,(tuple,list)):
        J = generate_jastrow(*J,system=system)
    #end if
    return J
#end def process_jastrow



def jastrow_factor(**kwargs):
    loc           = kwargs.pop('loc','jastrow_factor')
    defaults      = kwargs.pop('defaults',defaults_version)
    kw = extract_kwargs(
        kwargs          = kwargs,
        defaults        = defaults,
        default_sources = jastrow_factor_defaults,
        require_empty   = True,
        loc             = loc,
        )
    J1           = kw.J1      
    J2           = kw.J2      
    J3           = kw.J3      
    system       = kw.system  
    J1_size      = kw.J1_size 
    J1_rcut      = kw.J1_rcut 
    J2_size      = kw.J2_size 
    J2_rcut      = kw.J2_rcut 
    J2_init      = kw.J2_init 
    J3_isize     = kw.J3_isize
    J3_esize     = kw.J3_esize
    J3_rcut      = kw.J3_rcut 
    J1_rcut_open = kw.J1_rcut_open
    J2_rcut_open = kw.J2_rcut_open
    
    openbc = system.structure.is_open()

    J1 = process_jastrow(J1,system)
    J2 = process_jastrow(J2,system)
    J3 = process_jastrow(J3,system)

    if J1==True:
        if openbc and J1_rcut is None:
            J1_rcut = J1_rcut_open
        #end if
        J1 = generate_jastrow('J1','bspline',J1_size,J1_rcut,system=system)
    #end if
    if J2==True:
        if openbc and J2_rcut is None:
            J2_rcut = J2_rcut_open
        #end if
        J2 = generate_jastrow('J2','bspline',J2_size,J2_rcut,init=J2_init,system=system)
    #end if
    if J3==True:
        J3 = generate_jastrow('J3','polynomial',J3_esize,J3_isize,J3_rcut,system=system)
    #end if

    jastrows = []
    if J1!=False:
        jastrows.append(J1)
    #end if
    if J2!=False:
        jastrows.append(J2)
    #end if
    if J3!=False:
        jastrows.append(J3)
    #end if

    return jastrows
#end def jastrow_factor




def opt_sections(**kwargs):
    if 'opt_calcs' in kwargs:
        return kwargs['opt_calcs']
    #end if
    loc      = kwargs.pop('loc','opt_sections')
    defaults = kwargs.pop('loop_defaults',defaults_version)
    kw = extract_kwargs(
        kwargs          = kwargs,
        required        = opt_sections_required,
        defaults        = defaults,
        default_sources = opt_sections_defaults,
        encapsulate     = False,
        loc             = loc
        )
    method     = kw.method
    cost       = kw.cost
    cycles     = kw.cycles
    var_cycles = kw.var_cycles
    defaults   = kw.defaults
    del kw

    methods = obj(linear=linear,cslinear=cslinear)
    if method not in methods:
        error('invalid optimization method requested\ninvalid method: {0}\nvalid options are: {1}'.format(method,sorted(methods.keys())),loc)
    #end if
    opt = methods[method]

    if len(kwargs)>0:
        valid = set(opt.attributes + opt.parameters)
        invalid = set(kwargs.keys()) - valid
        if len(invalid)>0:
            error('invalid keywords encountered for {0} optimization\ninvalid keywords: {1}\nvalid options are: {2}'.format(method,sorted(invalid),sorted(valid)),loc)
        #end if
    #end if

    opt_inputs = extract_kwargs(
        kwargs          = kwargs,
        defaults        = defaults,
        default_sources = opt_method_defaults[method],
        encapsulate     = False,
        loc             = loc,
        )
    if cost=='variance':
        cost = (0.0,1.0,0.0)
    elif cost=='energy':
        cost = (1.0,0.0,0.0)
    elif isinstance(cost,(tuple,list)) and (len(cost)==2 or len(cost)==3):
        if len(cost)==2:
            cost = (cost[0],0.0,cost[1])
        #end if
    else:
        error('invalid optimization cost function encountered\ninvalid cost fuction: {0}\nvalid options are: variance, energy, (0.95,0.05), etc'.format(cost),loc)
    #end if
    opt_calcs = []
    if abs(cost[0])>1e-6:
        vmin_opt = opt(
            energy               = 0.0,
            unreweightedvariance = 1.0,
            reweightedvariance   = 0.0,
            **opt_inputs
            )
        opt_calcs.append(loop(max=var_cycles,qmc=vmin_opt))
    #end if
    cost_opt = opt(
        energy               = cost[0],
        unreweightedvariance = cost[1],
        reweightedvariance   = cost[2],
        **opt_inputs
        )
    opt_calcs.append(loop(max=cycles,qmc=cost_opt))
    return opt_calcs
#end def opt_sections



def vmc_sections(**kwargs):
    if 'vmc_calcs' in kwargs:
        return kwargs['vmc_calcs']
    #end if
    loc      = kwargs.pop('loc','vmc_sections')
    defaults = kwargs.pop('defaults',defaults_version)
    J0       = kwargs.pop('J0',False)
    test     = kwargs.pop('test',False)
    kw = extract_kwargs(
        kwargs          = kwargs,
        optional        = vmc_sections_options,
        defaults        = defaults,
        default_sources = vmc_defaults,
        require_empty   = True,
        encapsulate     = False,
        )
    if test:
        warmup = kw.test_warmupsteps,
        blocks = kw.test_blocks,
        steps  = kw.test_steps,
    elif J0:
        warmup = kw.J0_warmupsteps
        blocks = kw.J0_blocks
        steps  = kw.J0_steps
    else:
        warmup = kw.warmupsteps
        blocks = kw.blocks
        steps  = kw.steps
    #end if
    vmc_calcs = [
        vmc(
            walkers     = kw.walkers,
            warmupsteps = warmup,
            blocks      = blocks,
            steps       = steps,
            substeps    = kw.substeps,
            timestep    = kw.timestep,
            checkpoint  = kw.checkpoint,
            )
        ]
    return vmc_calcs
#end def vmc_sections



def dmc_sections(**kwargs):
    if 'dmc_calcs' in kwargs:
        return kwargs['vmc_calcs']
    #end if
    loc      = kwargs.pop('loc','dmc_sections')
    defaults = kwargs.pop('defaults',defaults_version)
    J0       = kwargs.pop('J0',False)
    test     = kwargs.pop('test',False)
    kw = extract_kwargs(
        kwargs          = kwargs,
        required        = dmc_sections_required,
        optional        = dmc_sections_options,
        defaults        = defaults,
        default_sources = vmc_defaults,
        require_empty   = True,
        encapsulate     = False,
        )
    if 'vmc_samples' not in kw and 'vmc_samplesperthread' not in kw:
        error('vmc samples (dmc walkers) not specified\nplease provide one of the following keywords: vmc_samples, vmc_samplesperthread',loc)
    #end if
    vsec = vmc(
        walkers     = kw.vmc_walkers,
        warmupsteps = kw.vmc_warmupsteps,
        blocks      = kw.vmc_blocks,
        steps       = kw.vmc_steps,
        substeps    = kw.vmc_substeps,
        timestep    = kw.vmc_timestep,
        checkpoint  = kw.vmc_checkpoint,
        )
    if 'vmc_samples' in kw:
        vsec.samples = kw.vmc_samples
    elif 'vmc_samplesperthread' in kw:
        vsec.samplesperthread = kw.vmc_samplesperthread
    #end if
    dmc_calcs = [vsec]
    if kw.eq_dmc:
        deqsec = dmc(
            walkers 
            )
    #end if

#end def dmc_sections




#  opt_inputs = obj(
#      J2_prod = True,
#      J3_prod = True,
#      )
#  vmc_inputs = obj(
#      J0_test = True,
#      J0_prod = True,
#      J2_prod = True,
#      J3_prod = True,
#      )
#  dmc_inputs = obj(
#      J3_test = True,
#      J3_prod = True,
#      test_calcs = [],  # has a default
#      prod_calcs = [],
#      )


qmcpack_chain_required = ['system','sim_list','dft_pseudos','qmc_pseudos']
qmcpack_chain_defaults = obj(
    scf            = False,
    p2q            = False,
    opt            = False,
    vmc            = False,
    dmc            = False,
    scf_inputs     = None,
    p2q_inputs     = None,
    opt_inputs     = None,
    vmc_inputs     = None,
    dmc_inputs     = None,
    scf_defaults   = defaults_version,
    p2q_defaults   = defaults_version,
    opt_defaults   = defaults_version,
    vmc_defaults   = defaults_version,
    dmc_defaults   = defaults_version,
    orb_source     = None,
    J2_source      = None,
    J3_source      = None,
    processed      = False,
    )


def process_qmcpack_chain_kwargs(kwargs,
                                 defaults      = None,
                                 require_empty = True,
                                 loc           = None,
                                 ):
    if loc is None:
        inloc = None
        loc   = 'process_qmcpack_chain_kwargs'
    else:
        inloc = loc
    #end if
    if defaults is None:
        defaults = qmcpack_chain_defaults,
    else:
        defaults.set_optional(**qmcpack_chain_defaults)
    #end if
    kw = extract_kwargs(
        kwargs        = kwargs,
        required      = qmcpack_chain_required,
        defaults      = defaults,
        require_empty = require_empty,
        loc           = loc,
        )
    kw.scf |= kw.scf_inputs!=None
    kw.p2q |= kw.p2q_inputs!=None
    kw.opt |= kw.opt_inputs!=None
    kw.vmc |= kw.vmc_inputs!=None
    kw.dmc |= kw.dmc_inputs!=None
    
    if kw.scf:
        kw.scf_inputs,kw.scf_workflow = extract_kwargs(
            kwargs            = kw.scf_inputs,
            required          = scf_required,
            defaults          = kw.scf_defaults,
            default_sources   = scf_defaults,
            workflow          = scf_workflow,
            contains_defaults = True,
            allow_none        = True,
            loc               = loc+' scf_inputs',
            )
        kw.scf_inputs.system  = kw.system
        kw.scf_inputs.pseudos = kw.dft_pseudos
    #end if
    if kw.p2q:
        kw.p2q_inputs,kw.p2q_workflow = extract_kwargs(
            kwargs            = kw.p2q_inputs,
            required          = p2q_required,
            defaults          = kw.p2q_defaults,
            default_sources   = p2q_defaults,
            workflow          = p2q_workflow,
            contains_defaults = True,
            allow_none        = True,
            loc               = loc+' p2q_inputs',
            )
    #end if
    if kw.opt:
        kw.opt_inputs,kw.opt_workflow = extract_kwargs(
            kwargs            = kw.opt_inputs,
            required          = opt_required,
            defaults          = kw.opt_defaults,
            default_sources   = opt_defaults,
            workflow          = opt_workflow,
            contains_defaults = True,
            allow_none        = True,
            require_empty     = True,
            loc               = loc+' opt_inputs',
            )
        kw.opt_inputs.system  = kw.system
        kw.opt_inputs.pseudos = kw.qmc_pseudos

        J_defaults = kw.opt_workflow.delete_required('J_defaults')
        jkw = extract_kwargs(
            kwargs          = kw.opt_workflow,
            defaults        = J_defaults,
            default_sources = jastrow_factor_defaults,
            require_empty   = False,
            loc             = loc+' opt_inputs jastrows'
            )
        jkw.system = kw.system
        j2kw = jkw.copy()
        j2kw.set(J1=1,J2=1,J3=0)
        j3kw = jkw.copy()
        j3kw.set(J1=1,J2=1,J3=1)
        kw.J2_inputs = j2kw
        kw.J3_inputs = j3kw
        kw.opt_sec_inputs = extract_kwargs(
            kwargs          = kw.opt_workflow,
            optional        = opt_sections_options,
            require_empty   = False,
            loc             = loc+' opt_inputs opt_methods'
            )
    #end if
    if kw.vmc:
        kw.vmc_inputs,kw.vmc_workflow = extract_kwargs(
            kwargs            = kw.vmc_inputs,
            required          = vmc_required,
            defaults          = kw.vmc_defaults,
            default_sources   = vmc_defaults,
            workflow          = vmc_workflow,
            contains_defaults = True,
            allow_none        = True,
            require_empty     = True,
            loc               = loc+' vmc_inputs',
            )
        kw.vmc_inputs.system  = kw.system
        kw.vmc_inputs.pseudos = kw.qmc_pseudos
        kw.vmc_sec_inputs = extract_kwargs(
            kwargs        = kw.vmc_workflow,
            optional      = vmc_sections_options,
            require_empty = False,
            loc           = loc+' vmc_inputs vmc_methods'
            )
    #end if
    if kw.dmc:
        kw.dmc_inputs,kw.dmc_workflow = extract_kwargs(
            kwargs            = kw.dmc_inputs,
            required          = dmc_required,
            defaults          = kw.dmc_defaults,
            default_sources   = dmc_defaults,
            workflow          = dmc_workflow,
            contains_defaults = True,
            allow_none        = True,
            loc               = loc+' dmc_inputs',
            )
        kw.dmc_inputs.system  = kw.system
        kw.dmc_inputs.pseudos = kw.qmc_pseudos
        kw.dmc_sec_inputs = extract_kwargs(
            kwargs        = kw.dmc_workflow,
            optional      = dmc_sections_options,
            require_empty = False,
            loc           = loc+' dmc_inputs dmc_methods'
            )
    #end if

    del kw.scf_defaults
    del kw.p2q_defaults
    del kw.opt_defaults
    del kw.vmc_defaults
    del kw.dmc_defaults

    if inloc!=None:
        kw.loc = loc
    #end if
    kw.processed = True

    return kw
#end def process_qmcpack_chain_kwargs


def qmcpack_chain(**kwargs):
    loc       = kwargs.pop('loc','qmcpack_chain')
    processed = kwargs.pop('processed',False)
    if processed:
        kw = obj(**kwargs)
    else:
        kw = process_qmcpack_chain_kwargs(kwargs,require_empty=True,loc=loc)
    #end if
    basepath = kw.basepath
    sim_list = kw.sim_list
    sims = obj()

    if kw.orb_source!=None:
        sims.p2q = kw.orb_source
    else:
        if kw.scf:
            scf = generate_pwscf(
                path = os.path.join(basepath,'scf'),
                **kw.scf_inputs
                )
            sims.scf = scf
        #end if

        if kw.p2q:
            deps = resolve_deps('p2q',sims,[('scf','orbitals')],loc)
            p2q = generate_pw2qmcpack(
                path         = os.path.join(basepath,'scf'),
                dependencies = deps,
                **kw.p2q_inputs
                )
            sims.p2q = p2q
        #end if
    #end if
        
    orbdep = [('p2q','orbitals')]
    J2dep  = orbdep + [('optJ2','jastrow')]
    J3dep  = orbdep + [('optJ3','jastrow')]

    if kw.opt:
        if kw.opt_workflow.J2_prod and kw.J2_source is None:
            deps = resolve_deps('optJ2',sims,orbdep,loc)
            optJ2 = generate_qmcpack(
                path         = os.path.join(basepath,'optJ2'),
                jastrows     = jastrow_factor(**kw.J2_inputs),
                calculations = opt_sections(**kw.opt_sec_inputs),
                dependencies = deps,
                **kw.opt_inputs
                )
            sims.optJ2 = optJ2
        #end if
        if kw.opt_workflow.J3_prod and kw.J3_source is None:
            deps = resolve_deps('optJ3',sims,J2dep,loc)
            optJ3 = generate_qmcpack(
                path         = os.path.join(basepath,'optJ3'),
                jastrows     = jastrow_factor(**kw.J3_inputs),
                calculations = opt_sections(**kw.opt_sec_inputs),
                dependencies = deps,
                **kw.opt_inputs
                )
            sims.optJ3 = optJ3
        #end if
    #end if
    if kw.J2_source!=None:
        sims.optJ2 = kw.J2_source
    #end if
    if kw.J3_source!=None:
        sims.optJ3 = kw.J3_source
    #end if

    if kw.vmc:
        if kw.vmc_workflow.J0_test:
            deps = resolve_deps('vmcJ0_test',sims,orbdep,loc)
            vmcJ0_test = generate_qmcpack(
                path         = os.path.join(basepath,'vmcJ0_test'),
                jastrows     = [],
                calculations = vmc_sections(test=1,J0=1,**kw.vmc_sec_inputs),
                dependencies = deps,
                **kw.vmc_inputs
                )
            sims.vmcJ0_test = vmcJ0_test
        #end if
        if kw.vmc_workflow.J0_prod:
            deps = resolve_deps('vmcJ0',sims,orbdep,loc)
            vmcJ0 = generate_qmcpack(
                path         = os.path.join(basepath,'vmcJ0'),
                jastrows     = [],
                calculations = vmc_sections(J0=1,**kw.vmc_sec_inputs),
                dependencies = deps,
                **kw.vmc_inputs
                )
            sims.vmcJ0 = vmcJ0
        #end if
        if kw.vmc_workflow.J2_test:
            deps = resolve_deps('vmcJ2_test',sims,J2dep,loc)
            vmcJ2_test = generate_qmcpack(
                path         = os.path.join(basepath,'vmcJ2_test'),
                jastrows     = [],
                calculations = vmc_sections(test=1,**kw.vmc_sec_inputs),
                dependencies = deps,
                **kw.vmc_inputs
                )
            sims.vmcJ2_test = vmcJ2_test
        #end if
        if kw.vmc_workflow.J2_prod:
            deps = resolve_deps('vmcJ2',sims,J2dep,loc)
            vmcJ2 = generate_qmcpack(
                path         = os.path.join(basepath,'vmcJ2'),
                jastrows     = [],
                calculations = vmc_sections(**kw.vmc_sec_inputs),
                dependencies = deps,
                **kw.vmc_inputs
                )
            sims.vmcJ2 = vmcJ2
        #end if
        if kw.vmc_workflow.J3_test:
            deps = resolve_deps('vmcJ3_test',sims,J3dep,loc)
            vmcJ3_test = generate_qmcpack(
                path         = os.path.join(basepath,'vmcJ3_test'),
                jastrows     = [],
                calculations = vmc_sections(test=1,**kw.vmc_sec_inputs),
                dependencies = deps,
                **kw.vmc_inputs
                )
            sims.vmcJ3_test = vmcJ3_test
        #end if
        if kw.vmc_workflow.J3_prod:
            deps = resolve_deps('vmcJ3',sims,J2dep,loc)
            vmcJ3 = generate_qmcpack(
                path         = os.path.join(basepath,'vmcJ3'),
                jastrows     = [],
                calculations = vmc_sections(**kw.vmc_sec_inputs),
                dependencies = deps,
                **kw.vmc_inputs
                )
            sims.vmcJ3 = vmcJ3
        #end if
    #end if

    if kw.dmc:
        nonloc_labels = {None:'',True:'_tm',False:'_la'}
        nonlocalmoves_default = None
        if system.is_pseudized():
            nonlocalmoves_default = True
        #end if
        nloc_moves = []
        if kw.dmc_workflow.tmoves:
            nloc_moves.append(True)
        #end if
        if kw.dmc_workflow.locality:
            nloc_moves.append(False)
        #end if
        if not kw.dmc_workflow.tmoves and not kw.dmc_workflow.locality:
            nloc_moves.append(nonlocal_moves_default)
        #end if
        for nlmove in nloc_moves:
            nll = nonloc_labels[nlmove]
            if kw.dmc_workflow.J0_test:
                label = 'dmcJ0'+nll+'_test'
                deps = resolve_deps(label,sims,orbdep,loc)
                dmcJ0_test = generate_qmcpack(
                    path         = os.path.join(basepath,label),
                    jastrows     = [],
                    calculations = dmc_sections(nlmove=nlmove,test=1,J0=1,**kw.dmc_sec_inputs),
                    dependencies = deps,
                    **kw.dmc_inputs
                    )
                sims[label] = dmcJ0_test
            #end if
            if kw.dmc_workflow.J0_prod:
                label = 'dmcJ0'+nll
                deps = resolve_deps(label,sims,orbdep,loc)
                dmcJ0 = generate_qmcpack(
                    path         = os.path.join(basepath,label),
                    jastrows     = [],
                    calculations = dmc_sections(nlmove=nlmove,J0=1,**kw.dmc_sec_inputs),
                    dependencies = deps,
                    **kw.dmc_inputs
                    )
                sims[label] = dmcJ0
            #end if
            if kw.dmc_workflow.J2_test:
                label = 'dmcJ2'+nll+'_test'
                deps = resolve_deps(label,sims,J2dep,loc)
                dmcJ2_test = generate_qmcpack(
                    path         = os.path.join(basepath,label),
                    jastrows     = [],
                    calculations = dmc_sections(nlmove=nlmove,test=1,**kw.dmc_sec_inputs),
                    dependencies = deps,
                    **kw.dmc_inputs
                    )
                sims[label] = dmcJ2_test
            #end if
            if kw.dmc_workflow.J2_prod:
                label = 'dmcJ2'+nll
                deps = resolve_deps(label,sims,J2dep,loc)
                dmcJ2 = generate_qmcpack(
                    path         = os.path.join(basepath,label),
                    jastrows     = [],
                    calculations = dmc_sections(nlmove=nlmove,**kw.dmc_sec_inputs),
                    dependencies = deps,
                    **kw.dmc_inputs
                    )
                sims[label] = dmcJ2
            #end if
            if kw.dmc_workflow.J3_test:
                label = 'dmcJ3'+nll+'_test'
                deps = resolve_deps(label,sims,J3dep,loc)
                dmcJ3_test = generate_qmcpack(
                    path         = os.path.join(basepath,label),
                    jastrows     = [],
                    calculations = dmc_sections(nlmove=nlmove,test=1,**kw.dmc_sec_inputs),
                    dependencies = deps,
                    **kw.dmc_inputs
                    )
                sims[label] = dmcJ3_test
            #end if
            if kw.dmc_workflow.J3_prod:
                label = 'dmcJ3'+nll
                deps = resolve_deps(label,sims,J2dep,loc)
                dmcJ3 = generate_qmcpack(
                    path         = os.path.join(basepath,label),
                    jastrows     = [],
                    calculations = dmc_sections(nlmove=nlmove,**kw.dmc_sec_inputs),
                    dependencies = deps,
                    **kw.dmc_inputs
                    )
                sims[label] = dmcJ3
            #end if
        #end if
    #end if

    sim_list.extend(sims.list())

    return sims
#end def qmcpack_chain


ecut_scan_required = ['ecuts','basepath']
ecut_scan_defaults = obj(
    dirname      = 'ecut_scan',
    same_jastrow = True,
    ecut_jastrow = None,
    )
ecut_scan_chain_defaults = obj(
    scf = True,
    p2q = True,
    opt = True,
    vmc = True,
    )
def ecut_scan(**kwargs):
    loc = kwargs.pop('loc','ecut_scan')
    kw = extract_kwargs(
        kwargs   = kwargs,
        required = ecut_scan_required,
        defaults = ecut_scan_defaults,
        loc      = loc
        )
    qckw = process_qmcpack_chain_kwargs(
        kwargs   = kwargs,
        defaults = ecut_scan_chain_defaults,
        loc      = loc
        )
    basepath     = kw.basepath
    dirname      = kw.dirname
    ecuts        = list(kw.ecuts)
    same_jastrow = kw.same_jastrow
    ecut_jastrow = kw.ecut_jastrow
    del kw

    if not qckw.scf:
        error('cannot perform ecut scan, no inputs given for scf calculations',loc)
    #end if

    if same_jastrow:
        if ecut_jastrow is None:
            ecut_jastrow = max(ecuts)
        #end if
        found = False
        n=0
        for ecut in ecuts:
            if abs(ecut_jastrow-ecut)<1e-2:
                found = True
                break
            #end if
            n+=1
        #end for
        if not found:
            error('could not find ecut for fixed jastrow in list\necut searched for: {0}\necut list: {1}'.format(ecut_jastrow,ecuts),loc)
        #end if
        ecuts.pop(n)
        ecuts = [ecut_jastrow]+ecuts
    #end if
    J2_source = None
    J3_source = None
    sims = obj()
    for ecut in ecuts:
        qckw.basepath = os.path.join(basepath,dirname,'ecut_{0}'.format(ecut))
        qckw.scf_inputs.ecutwfc = ecut
        if J2_source is not None:
            qckw.J2_source = J2_source
        #end if
        if J3_source is not None:
            qckw.J3_source = J3_source
        #end if
        qcsims = qmcpack_chain(**qckw)
        if same_jastrow:
            J2_source = qcsims.get_optional('optJ2',None)
            J3_source = qcsims.get_optional('optJ3',None)
        #end if
        sims[ecut] = qcsims
    #end for
    return sims
#end def ecut_scan




if __name__=='__main__':
    print 'simple driver for qmcpack_workflows functions'

    tests = obj(
        opt_sections = 1,
        )

    if tests.opt_sections:
        print
        print 'opt_sections()'
        opt_calcs = opt_sections()
        for calc in opt_calcs:
            print calc
        #end for

        print
        print "opt_sections(defaults=None)"
        opt_calcs = opt_sections(defaults=None)
        for calc in opt_calcs:
            print calc
        #end for

        print
        print "opt_sections(defaults='yl')"
        opt_calcs = opt_sections(defaults='yl')
        for calc in opt_calcs:
            print calc
        #end for

        print
        print "opt_sections(cost='energy')"
        opt_calcs = opt_sections(cost='energy')
        for calc in opt_calcs:
            print calc
        #end for

        print
        print "opt_sections(cost=(0.95,0.05))"
        opt_calcs = opt_sections(cost=(0.95,0.05))
        for calc in opt_calcs:
            print calc
        #end for

        print
        print "opt_sections(cost=(0.90,0.03,0.07),defaults=None)"
        opt_calcs = opt_sections(cost=(0.90,0.03,0.07),defaults=None)
        for calc in opt_calcs:
            print calc
        #end for

        print
        print "opt_sections(cost='energy',cycles=20,var_cycles=10,samples=1000000,defaults=None)"
        opt_calcs = opt_sections(cost='energy',cycles=20,var_cycles=10,samples=1000000,defaults=None)
        for calc in opt_calcs:
            print calc
        #end for

        #print
        #print 'fail test'
        #opt_sections(method='my_opt')
        #opt_sections(defaults='my_defaults')
        #opt_sections(bad_key='this')
        #opt_sections(cost=[0,0,0,0])
    #end if


#end if
