
import os
import itertools
from time import clock
from numpy import ndarray,ceil
from developer import obj,ci,error as dev_error,devlog,DevBase
from physical_system import generate_physical_system
from simulation import Simulation,GenericSimulation,graph_sims
from bundle import bundle as bundle_function
from machines import job,Job
from vasp  import generate_vasp
from pwscf import generate_pwscf
from pwscf_postprocessors import generate_projwfc,generate_pp
from qmcpack_converters import generate_pw2qmcpack
from qmcpack_input import generate_jastrow,loop,linear,cslinear,vmc,dmc
from qmcpack import generate_qmcpack



# remaining dev tasks
#   enable scans over shared parameters (e.g. pseudopotential sets)
#   introduce relax sims (pwscf and vasp) into qmcpack_chain
#   make ppset and add inputs
#   eliminate detected fake sims in project manager
#     probably do traversals to collect fake ones, then iterate to eliminate
#
#   track association between input section and generated sims in qmcpack_chain
#     this is needed in SimRep.__init__
#   set first input of each scan into fake chain prior to generating
#   make "i" the default pattern if scanned variables are not str/int/float
#   enable user provided labels for scan variables
#     enables non-integer directory labeling for non-str/int/float variables
#   consider how to replace deep copies in qmcpack_workflow (marked)
#     test robustness of shallow copy, may be ok
#   enable text print out of workflow tree
#     useful on machines that do not allow plotting
#
#   longer term
#     make way to limit downstream branching action, not just "merge"
#       perhaps at the level of segments based on key matching
#     determine how best to extract output data from scanned simulations
#

#
# enabling scan over system params
#   add system_inputs to qmcpack chain
#     process_qmcpack_chain_kwargs should essentially leave it alone
#   if system_inputs is present
#     generate a system object from inputs (via gen_phys_sys or user function)
#     set kw.system to this
#   if kw.fake then also create a SystemHolder(Simulation) object
#     is fake sim by definition
#     will need to have simlabel='system' for current secrep association
#     make all sims in chain depend on it via 'other'
#     this way, will appear in segmented representation and all variations
#     should it go in the chain sims object or not? it has to
#     SystemHolder's should be eliminated when qmcpack_workflow completes
#       will also need to remove them from workflow sims object
#         perhaps via implementing a 'recurse' function in obj class
#       should probably make temp simlist at qmcpack_workflow start



def error(msg,loc=None,exit=True,trace=True,indent='    ',logfile=devlog):
    header = 'qmcpack_workflows'
    if loc!=None:
        msg+='\nfunction location: {0}'.format(loc)
    #end if
    dev_error(msg,header,exit,trace,indent,logfile)
#end def error



defaults_version = 'v1'


def hashable(v):
    try:
        hash(v)
    except:
        return False
    #end try
    return True
#end def hashable


class Missing:
    def __call__(self,value):
        return isinstance(value,Missing)
    #end def __call__
#end class Missing
missing = Missing()


class Novalue:
    def __call__(self,value):
        return isinstance(value,Novalue)
    #end def __call__

    def __str__(self):
        return ''
    #end def __str__

    def __repr__(self):
        return ''
    #end def __repr__
#end class Novalue
novalue = Novalue()


class CallBase(object):
    strict = True
    defs   = obj()
    location = 'unknown location'
#end class CallBase


class Defaults(CallBase):
    def __call__(self,name,value):
        if not missing(value):
            return value
        elif name in self.defs:
            return self.defs[name]
        elif self.strict:
            error('default value is missing for variable named {0}'.format(name),self.location)
        #end if
    #end def __call__
#end class Defaults
default = Defaults()


class Requirement(CallBase):
    def __call__(self,name,value):
        if missing(value):
            error('an input value has not been provided for required variable named {0}'.format(name),self.location)
        #end if
        return value
    #end def __call__
#end class Requirement
require = Requirement()


class RequireAny(CallBase):
    def __call__(self,*name_values):
        any_present = False
        for name,value in name_values:
            any_present |= not missing(value)
        #end for
        if not any_present:
            error('an input value must be provided for at least one of these variables: {0}'.format([name for (name,value) in name_values]),self.location)
        #end if
    #end def __call__
#end class RequireAny
require_any = RequireAny()


class Assignment(CallBase):
    def __call__(self,o,name,value):
        if not missing(value):
            o[name] = value
        #end if
    #end def __call__
#end class Assignment
assign = Assignment()


class RequireAssignment(CallBase):
    def __call__(self,o,name,value):
        if missing(value):
            error('a value has not been provided for required input variable named "{0}"'.format(name),self.location)
        #end if
        o[name] = value
    #end def __call__
#end class RequireAssignment
assign_require = RequireAssignment()


class DefaultAssignment(CallBase):
    def __call__(self,o,name,value):
        if not missing(value):
            o[name] = value
        elif name in self.defs:
            o[name] = self.defs[name]
        elif self.strict:
            error('default value is missing for variable named "{0}"'.format(name),self.location)
        #end if
    #end def __call__
#end class DefaultAssignment
assign_default = DefaultAssignment()


def set_def_loc(defs,loc,strict=True):
    CallBase.defs     = defs
    CallBase.location = loc
    CallBase.strict   = strict
#end def set_def_loc


def set_loc(loc):
    CallBase.location = loc
#end def set_loc


def assign_defaults(o,defs):
    for name,value in defs.iteritems():
        if name not in o:
            o[name] = value
        #end if
    #end for
#end def assign_defaults


def extract_keywords(o,names,optional=False):
    k = obj()
    if not optional:
        missing = set(names)-set(o.keys())
        if len(missing)>0:
            error('keywords are missing, please provide them\nmissing keywords: {0}\nfunction location: {1}'.format(sorted(missing),CallBase.location))
        #end if
        for name in names:
            k[name] = o[name]
            del o[name]
        #end for
    else:
        for name in names:
            if name in o:
                k[name] = o[name]
                del o[name]
            #end if
        #end for
    #end if
    return k
#end def extract_keywords


def prevent_invalid_input(invalid,valid,loc):
    if len(invalid)>0:
        if isinstance(invalid,(dict,obj)):
            invalid = invalid.keys()
        #end if
        error('invalid input keywords encountered\ninvalid keywords: {0}\nvalid options are: {1}\nfunction location: {2}'.format(sorted(invalid),valid,loc))
    #end if
#end def prevent_invalid_input


def render_parameter_key(v,loc='render_parameter_key'):
    if isinstance(v,list):
        vkey = tuple(v)
    elif isinstance(v,ndarray):
        vkey = tuple(v.ravel())
    else:
        vkey = v
    #end if
    if not hashable(vkey):
        error('parameter value is not hashable\nvalue provided: {0}\nvalue type: {1}\nplease restrict parameter values to basic types such as str,int,float,tuple and combinations of these'.format(v,v.__class__.__name__),loc)
    #end if
    return vkey
#end def render_parameter_key


def render_parameter_label(vkey,loc='render_parameter_label'):
    if isinstance(vkey,(int,float,str)):
        vlabel = str(vkey)
    elif isinstance(vkey,tuple):
        vlabel = str(vkey).replace('(','').replace(')','').replace(' ','').replace(',','_')
    else:
        error('cannot transform parameter value key into a directory name\nvalue key: {0}\ntype: {1}\nplease restrict parameter values to basic types such as str,int,float,tuple and combinations of these'.format(vkey,vkey.__class__.__name__),loc)
    #end if
    return vlabel
#end def render_parameter_label


def capture_inputs(inputs):
    if inputs is None or missing(inputs):
        capture = obj()
    else:
        capture = obj(**inputs)
    #end if
    return capture
#end def capture_inputs


class SimSet(DevBase):
    def __setitem__(self,simname,sim):
        if simname in self:
            self.error('simulation {0} is already set\nthis is a developer error'.format(simname))
        elif not isinstance(sim,(Simulation,SimSet)):
            self.error('only simulation objects are allowed in simset\nreceived type: {0}\nwith name: {1}\nthis is a developer error'.format(sim.__class__.__name__,simname))
        else:
            self.__dict__[simname] = sim
        #end if
    #end def

    def remove_fake(self):
        fake_sims = []
        for k,v in self.iteritems():
            if isinstance(v,Simulation):
                if v.fake():
                    fake_sims.append(k)
                #end if
            elif isinstance(v,SimSet):
                v.remove_fake()
            #end if
        #end for
        for k in fake_sims:
            del self[k]
        #end for
    #end def remove_fake
#end class SimSet


class SectionPlaceHolder(DevBase):
    None
#end class SectionPlaceHolder
placeholder = SectionPlaceHolder


class SimHolder(GenericSimulation):
    cls_simlabel = None
    def __init__(self,*args,**kwargs):
        cls = self.__class__
        Simulation.__init__(self,path='',job=job(app_command='',fake=True),fake_sim=True)
        if cls.cls_simlabel!=None:
            self.simlabel = cls.cls_simlabel
        #end if
    #end def __init__
#end class SimHolder


class SystemHolder(SimHolder):
    cls_simlabel = 'system'
#end class SystemHolder



from memory import resident
tprev = clock()
mprev = resident()
def tpoint(loc,n=0,indent='  '):
    global tprev
    global mprev
    tnow = clock()
    mnow = resident()
    #print '{0} {1} elapsed={2}'.format(n*indent,loc,tnow-tprev)
    print '{0} {1} elapsed={2}  {3:3.2f}'.format(n*indent,loc,tnow-tprev,(mnow-mprev)/1e6)
    tprev = tnow
    mprev = mnow
#end def tpoint
    








jastrow_factor_keys = ['J1','J2','J3',
                       'J1_size','J1_rcut',
                       'J2_size','J2_rcut','J2_init',
                       'J3_isize','J3_esize','J3_rcut',
                       'J1_rcut_open','J2_rcut_open']
jastrow_factor_defaults = obj(
    v1 = obj(
        J1           = True,
        J2           = True,
        J3           = False,
        J1_size      = None,
        J1_rcut      = None,
        J2_size      = None,
        J2_rcut      = None,
        J2_init      = 'zero',
        J3_isize     = 3,
        J3_esize     = 3,
        J3_rcut      = 5.0,
        J1_rcut_open = 5.0,
        J2_rcut_open = 10.0,
        ),
    )



opt_sections_keys = [
    'method','cost','cycles','var_cycles','opt_calcs','blocks',
    'warmupsteps','stepsbetweensamples','timestep','samples',
    'minwalkers','maxweight','usedrift','minmethod','beta',
    'exp0','bigchange','alloweddifference','stepsize',
    'stabilizerscale','nstabilizers','max_its','cgsteps',
    'eigcg','walkers','nonlocalpp','usebuffer','gevmethod',
    'steps','substeps','stabilizermethod','cswarmupsteps',
    'alpha_error','gevsplit','beta_error'
    ]



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


vmc_sections_keys = [
    'walkers','warmupsteps','blocks','steps',
    'substeps','timestep','checkpoint',
    'J0_warmupsteps','J0_blocks','J0_steps',
    'test_warmupsteps','test_blocks','test_steps',
    ]
vmc_sections_defaults = obj(
    v1 = obj(
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
        ),
    )



dmc_sections_keys = [
    'warmupsteps','blocks','steps',
    'timestep','checkpoint',
    'vmc_samples','vmc_samplesperthread',
    'vmc_walkers','vmc_warmupsteps','vmc_blocks','vmc_steps',
    'vmc_substeps','vmc_timestep','vmc_checkpoint',
    'eq_dmc','eq_warmupsteps','eq_blocks','eq_steps','eq_timestep','eq_checkpoint',
    'J0_warmupsteps','J0_blocks','J0_steps','J0_checkpoint',
    'test_warmupsteps','test_blocks','test_steps',
    'ntimesteps','timestep_factor','nlmove',
    ]
dmc_sections_defaults = obj(
    v1 = obj(
        warmupsteps          = 20,
        blocks               = 200,
        steps                = 10,
        timestep             = 0.01,
        checkpoint           = -1,
        vmc_samples          = 2048,
        #vmc_samplesperthread = missing, 
        vmc_walkers          = 1,
        vmc_warmupsteps      = 30,
        vmc_blocks           = 20,
        vmc_steps            = 10,
        vmc_substeps         = 3,
        vmc_timestep         = 0.3,
        vmc_checkpoint       = -1,
        eq_dmc               = False,
        eq_warmupsteps       = 20,
        eq_blocks            = 20,
        eq_steps             = 5,
        eq_timestep          = 0.02,
        eq_checkpoint        = -1,
        J0_warmupsteps       = 40,
        J0_blocks            = 400,
        J0_steps             = 20,
        test_warmupsteps     = 2,
        test_blocks          = 10,
        test_steps           = 2,
        ntimesteps           = 1,
        timestep_factor      = 0.5,    
        nlmove               = None,
        ),
    )




system_workflow_keys = []
system_input_defaults = obj(
    v1 = obj(),
    )



vasp_workflow_keys = []
vasp_input_defaults = obj(
    none = obj(
        ),
    minimal = obj(
        identifier = 'vasp',
        ),
    v1 = obj(
        identifier = 'vasp',
        ),
    )



scf_workflow_keys = [
    'struct_src',
    ]
fixed_defaults = obj(
    struct_src = None,
    )
scf_input_defaults = obj(
    none    = obj(
        **fixed_defaults
        ),
    minimal = obj(
        identifier       = 'scf',
        input_type       = 'generic',
        calculation      = 'scf',
        nosym            = True,
        wf_collect       = True,
        use_folded       = True,
        nogamma          = True,
        **fixed_defaults
        ),
    v1 = obj(
        identifier       = 'scf',
        input_type       = 'generic',
        calculation      = 'scf',
        verbosity        = 'high',
        disk_io          = 'low',
        diagonalization  = 'david',
        electron_maxstep = 1000,
        conv_thr         = 1e-8,
        mixing_beta      = 0.2,
        occupations      = 'smearing',
        smearing         = 'fermi-dirac',
        degauss          = 0.0001,
        nosym            = True, # should be False if nscf is next
        wf_collect       = True, # should be False if nscf is next
        use_folded       = True,
        nogamma          = True,
        **fixed_defaults
        ),
    )

nscf_workflow_keys = [
    'struct_src','dens_src',
    'use_scf_dir',
    ]
fixed_defaults = obj(
    struct_src  = None,
    dens_src    = None,
    use_scf_dir = False,
    )
nscf_input_defaults = obj(
    none    = obj(
        **fixed_defaults
        ),
    minimal = obj(
        identifier       = 'nscf',
        input_type       = 'generic',
        calculation      = 'nscf',
        use_folded       = True,
        nogamma          = True,
        **fixed_defaults
        ),
    v1 = obj(
        identifier       = 'nscf',
        input_type       = 'generic',
        calculation      = 'nscf',
        verbosity        = 'high',
        disk_io          = 'low',
        diagonalization  = 'david',
        electron_maxstep = 1000,
        conv_thr         = 1e-8,
        mixing_beta      = 0.2,
        occupations      = 'smearing',
        smearing         = 'fermi-dirac',
        degauss          = 0.0001,
        use_folded       = True,
        nogamma          = True,
        **fixed_defaults
        ),
    )

p2q_workflow_keys = [
    'orb_src','orb_src_type',
    ]
fixed_defaults = obj(
    orb_src      = None,
    orb_src_type = None,
    )
p2q_input_defaults = obj(
    minimal = obj(
        identifier = 'p2q',
        write_psir = False,
        **fixed_defaults
        ),
    v1 = obj(
        identifier = 'p2q',
        write_psir = False,
        **fixed_defaults
        ),
    )

pwf_workflow_keys = [
    'src',
    ]
fixed_defaults = obj(
    src = None,
    )
pwf_input_defaults = obj(
    minimal = obj(
        identifier = 'pwf',
        **fixed_defaults
        ),
    v1 = obj(
        identifier = 'pwf',
        **fixed_defaults
        ),
    )

pp_workflow_keys = [
    'src',
    ]
fixed_defaults = obj(
    src = None,
    )
pp_input_defaults = obj(
    minimal = obj(
        identifier = 'pp',
        **fixed_defaults
        ),
    v1 = obj(
        identifier = 'pp',
        **fixed_defaults
        ),
    )



opt_workflow_keys = [
    'J2_run','J3_run','J_defaults',
    'struct_src','orb_src','J2_src','J3_src',
    'use_J2',
    ]
fixed_defaults = obj(
    J2_run     = False,
    J3_run     = False,
    struct_src = None,
    orb_src    = None,
    J2_src     = None,
    J3_src     = None,
    use_J2     = False,
    J_defaults = defaults_version,
    )
opt_input_defaults = obj(
    minimal = obj(
        identifier   = 'opt',
        input_type   = 'basic',
        **fixed_defaults
        ),
    v1 = obj(
        identifier     = 'opt',
        input_type     = 'basic',
        spin_polarized = True,
        **fixed_defaults
        ),
    )

vmc_workflow_keys = [
    'J0_run','J2_run','J3_run',
    'J0_test','J2_test','J3_test',
    'struct_src','orb_src','J2_src','J3_src',
    'job','test_job',
    'prior_opt',
    ] 
fixed_defaults = obj(
    J0_run     = False,
    J2_run     = False,
    J3_run     = False,
    J0_test    = False,
    J2_test    = False,
    J3_test    = False,
    struct_src = None,
    orb_src    = None,
    J2_src     = None,
    J3_src     = None,
    job        = None,
    test_job   = None,
    prior_opt  = True,
    )
vmc_input_defaults = obj(
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
        spin_polarized   = True,
        **fixed_defaults
        ),
    )

dmc_workflow_keys = [
    'J0_run','J2_run','J3_run',
    'J0_test','J2_test','J3_test',
    'struct_src','orb_src','J2_src','J3_src',
    'job','test_job',
    'tmoves','locality',
    'prior_opt',
    ]
fixed_defaults = obj(
    J0_run     = False,
    J2_run     = False,
    J3_run     = False,
    J0_test    = False,
    J2_test    = False,
    J3_test    = False,
    struct_src = None,
    orb_src    = None,
    J2_src     = None,
    J3_src     = None,
    job        = None,
    test_job   = None,
    tmoves     = False,
    locality   = False,
    prior_opt  = True,
    )
dmc_input_defaults = obj(
    minimal = obj(
        identifier   = 'qmc',
        input_type   = 'basic',
        **fixed_defaults
        ),
    v1 = obj(
        identifier     = 'qmc',
        input_type     = 'basic',
        spin_polarized = True,
        **fixed_defaults
        ),
    )



def resolve_label(name,ch,loc='resolve_label',inp=False,quant=None,index=None,ind=False):
    ind_req = ind
    del ind
    kw   = ch.kw
    task = ch.task
    wf   = task.workflow
    if index is None:
        index = 1
        # def_srcs tracks the last occurring sim source in the workflow chain
        if name in kw.def_srcs:
            index = kw.def_srcs[name]
        #end if
        if quant is not None:
            if quant=='structure' and 'struct_src' in wf:
                ind = wf.struct_src
            elif quant=='charge_density' and 'dens_src' in wf:
                ind = wf.dens_src
            elif quant=='orbitals' and 'orb_src' in wf:
                ind = wf.orb_src
            elif quant=='jastrow':
                if 'J2' in name and 'J2_src' in wf:
                    ind = wf.J2_src
                elif 'J3' in name and 'J3_src' in wf:
                    ind = wf.J3_src
                #end if
            elif quant=='other' and 'src' in wf:
                ind = wf.src
            #end if
            if ind is not None:
                index = ind
            #end if
        #end if
    #end if
    
    if isinstance(index,Simulation):
        sim = index
        if not inp:
            ret = sim
        else:
            ret = sim,name+'_inputs'
        #end if
    else:
        if index==1:
            sind = ''
        else:
            sind = str(index)
        #end if
        if '_' not in name:
            label  = name+sind
        else:
            name,postfix = name.split('_',1)
            label = name+sind+'_'+postfix
        #end if
        inputs = name+'_inputs'+sind 
        if not inp:
            ret = label
        else:
            ret = label,inputs
        #end if
    #end if
    if ind_req:
        if isinstance(ret,tuple):
            ret = tuple(list(ret)+[index])
        else:
            ret = ret,index
        #end if
    #end if
    return ret
#end def resolve_label


def resolve_path(name,ch,loc='resolve_path'):
    kw   = ch.kw
    task = ch.task
    wf   = task.workflow
    index = None
    if 'scf' in name:
        if task.name=='nscf':
            index = wf.dens_src
        elif task.name=='p2q':
            index = wf.orb_src
        elif task.name=='pwf':
            index = wf.src
        elif task.name=='pp':
            index = wf.src
        #end if
    #end if
    if isinstance(index,Simulation):
        sim = index
        path = sim.locdir
    else:
        label = resolve_label(name,ch,loc,index=index)
        path = os.path.join(wf.basepath,label)
    #end if
    return path
#end def resolve_path


def resolve_deps(name,deps,ch,loc='resolve_deps'):
    kw   = ch.kw
    task = ch.task
    sims = ch.sims
    deplist = []
    missing = []
    missing_inp = []
    for depname,depquant in deps:
        label,inp = resolve_label(depname,ch,loc,inp=True,quant=depquant)
        if isinstance(label,Simulation):
            sim = label
            deplist.append((sim,depquant))
        elif label in sims:
            deplist.append((sims[label],depquant))
        else:
            missing.append(label)
            missing_inp.append(inp)
        #end if
    #end for
    if len(missing)>0:
        error('workflow cannot be run\nsimulation "{0}" depends on other simulations that have not been requested\nmissing simulations: {1}\nthe user needs to provide more detailed input\nthis issue can likely be fixed by providing the following keywords: {2}'.format(name,sorted(missing),sorted(set(missing_inp))))
    #end if
    return deplist
#end def resolve_deps


def get_pseudos(type,shared):
    if not missing(shared.pseudos):
        return shared.pseudos
    elif type=='pwscf':
        return shared.dft_pseudos
    elif type=='qmcpack':
        return shared.qmc_pseudos
    else:
        error('pseudopotential type "{0}" is unrecognized\nthis is a developer error'.format(type),'get_pseudos')
    #end if
#end def get_pseudos






def process_jastrow(J,system):
    if isinstance(J,(tuple,list)):
        J = generate_jastrow(*J,system=system)
    #end if
    return J
#end def process_jastrow



def jastrow_factor(
    J1           = missing,
    J2           = missing,
    J3           = missing,
    system       = missing,
    J1_size      = missing,
    J1_rcut      = missing,
    J2_size      = missing,
    J2_rcut      = missing,
    J2_init      = missing,
    J3_isize     = missing,
    J3_esize     = missing,
    J3_rcut      = missing,
    J1_rcut_open = missing,
    J2_rcut_open = missing,
    defaults     = defaults_version,
    loc          = 'jastrow_factor',
    ):

    set_def_loc(jastrow_factor_defaults[defaults],loc)

    J1           = default('J1'          ,J1          )
    J2           = default('J2'          ,J2          )
    J3           = default('J3'          ,J3          )
    J1_size      = default('J1_size'     ,J1_size     )
    J1_rcut      = default('J1_rcut'     ,J1_rcut     )
    J2_size      = default('J2_size'     ,J2_size     )
    J2_rcut      = default('J2_rcut'     ,J2_rcut     )
    J2_init      = default('J2_init'     ,J2_init     )
    J3_isize     = default('J3_isize'    ,J3_isize    )
    J3_esize     = default('J3_esize'    ,J3_esize    )
    J3_rcut      = default('J3_rcut'     ,J3_rcut     )
    J1_rcut_open = default('J1_rcut_open',J1_rcut_open) 
    J2_rcut_open = default('J2_rcut_open',J2_rcut_open) 
    
    require('system',system)

    if system.structure.units!='B':
        system = system.copy()
        system.structure.change_units('B')
    #end if

    openbc = system.structure.is_open()

    J1 = process_jastrow(J1,system)
    J2 = process_jastrow(J2,system)
    J3 = process_jastrow(J3,system)

    if J1==True:
        if J1_rcut is None:
            if openbc:
                J1_rcut = J1_rcut_open
            else:
                J1_rcut = system.structure.rwigner(1)
            #end if
        #end if
        if J1_size is None:
            J1_size = int(ceil(J1_rcut/0.5))
        #end if
        J1 = generate_jastrow('J1','bspline',J1_size,J1_rcut,system=system)
    #end if
    if J2==True:
        if J2_rcut is None:
            if openbc:
                J2_rcut = J2_rcut_open
            else:
                J2_rcut = system.structure.rwigner(1)
            #end if
        #end if
        if J2_size is None:
            J2_size = int(ceil(J2_rcut/0.5))
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




def opt_sections(
    method              = missing,
    cost                = missing,
    cycles              = missing,
    var_cycles          = missing,
    opt_calcs           = missing,
    # linear/cslinear inputs
    blocks              = missing,
    warmupsteps         = missing,
    stepsbetweensamples = missing,
    timestep            = missing,
    samples             = missing,
    minwalkers          = missing,
    maxweight           = missing,
    usedrift            = missing,
    minmethod           = missing,
    beta                = missing,
    exp0                = missing,
    bigchange           = missing,
    alloweddifference   = missing,
    stepsize            = missing,
    stabilizerscale     = missing,
    nstabilizers        = missing,
    max_its             = missing,
    cgsteps             = missing,
    eigcg               = missing,
    walkers             = missing,
    nonlocalpp          = missing,
    usebuffer           = missing,
    gevmethod           = missing,
    steps               = missing,
    substeps            = missing, 
    stabilizermethod    = missing,
    cswarmupsteps       = missing,
    alpha_error         = missing,
    gevsplit            = missing,
    beta_error          = missing,
    defaults            = defaults_version,
    loc                 = 'opt_sections',
    ):

    if not missing(opt_calcs):
        return opt_calcs
    #end if

    set_def_loc(opt_sections_defaults[defaults],loc)

    method     = default('method'    ,method    )
    cost       = default('cost'      ,cost      )
    cycles     = default('cycles'    ,cycles    )
    var_cycles = default('var_cycles',var_cycles)

    methods = obj(linear=linear,cslinear=cslinear)
    if method not in methods:
        error('invalid optimization method requested\ninvalid method: {0}\nvalid options are: {1}'.format(method,sorted(methods.keys())),loc)
    #end if
    opt = methods[method]

    set_def_loc(opt_method_defaults[method][defaults],loc,strict=False)

    opt_inputs = obj()
    assign_default(opt_inputs,'blocks'             ,blocks             )
    assign_default(opt_inputs,'warmupsteps'        ,warmupsteps        )
    assign_default(opt_inputs,'stepsbetweensamples',stepsbetweensamples)
    assign_default(opt_inputs,'timestep'           ,timestep           )
    assign_default(opt_inputs,'samples'            ,samples            )
    assign_default(opt_inputs,'minwalkers'         ,minwalkers         )
    assign_default(opt_inputs,'maxweight'          ,maxweight          )
    assign_default(opt_inputs,'usedrift'           ,usedrift           )
    assign_default(opt_inputs,'minmethod'          ,minmethod          )
    assign_default(opt_inputs,'beta'               ,beta               )
    assign_default(opt_inputs,'exp0'               ,exp0               )
    assign_default(opt_inputs,'bigchange'          ,bigchange          )
    assign_default(opt_inputs,'alloweddifference'  ,alloweddifference  )
    assign_default(opt_inputs,'stepsize'           ,stepsize           )
    assign_default(opt_inputs,'stabilizerscale'    ,stabilizerscale    )
    assign_default(opt_inputs,'nstabilizers'       ,nstabilizers       )
    assign_default(opt_inputs,'max_its'            ,max_its            )
    assign_default(opt_inputs,'cgsteps'            ,cgsteps            )
    assign_default(opt_inputs,'eigcg'              ,eigcg              )
    assign_default(opt_inputs,'walkers'            ,walkers            )
    assign_default(opt_inputs,'nonlocalpp'         ,nonlocalpp         )
    assign_default(opt_inputs,'usebuffer'          ,usebuffer          )
    assign_default(opt_inputs,'gevmethod'          ,gevmethod          )
    assign_default(opt_inputs,'steps'              ,steps              )
    assign_default(opt_inputs,'substeps'           ,substeps           ) 
    assign_default(opt_inputs,'stabilizermethod'   ,stabilizermethod   )
    assign_default(opt_inputs,'cswarmupsteps'      ,cswarmupsteps      )
    assign_default(opt_inputs,'alpha_error'        ,alpha_error        )
    assign_default(opt_inputs,'gevsplit'           ,gevsplit           )
    assign_default(opt_inputs,'beta_error'         ,beta_error         )

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
    if abs(cost[0])>1e-6 and var_cycles>0:
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



def vmc_sections(
    walkers          = missing,
    warmupsteps      = missing,
    blocks           = missing,
    steps            = missing,
    substeps         = missing,
    timestep         = missing,
    checkpoint       = missing,
    J0_warmupsteps   = missing,
    J0_blocks        = missing,
    J0_steps         = missing,
    test_warmupsteps = missing,
    test_blocks      = missing,
    test_steps       = missing,
    vmc_calcs        = missing,
    J0               = False,
    test             = False,
    defaults         = defaults_version,
    loc              = 'vmc_sections',
    ):

    if not missing(vmc_calcs):
        return vmc_calcs
    #end if

    set_def_loc(vmc_sections_defaults[defaults],loc)

    walkers          = default('walkers'         ,walkers         )
    warmupsteps      = default('warmupsteps'     ,warmupsteps     )
    blocks           = default('blocks'          ,blocks          )
    steps            = default('steps'           ,steps           )
    substeps         = default('substeps'        ,substeps        )
    timestep         = default('timestep'        ,timestep        )
    checkpoint       = default('checkpoint'      ,checkpoint      )
    J0_warmupsteps   = default('J0_warmupsteps'  ,J0_warmupsteps  )
    J0_blocks        = default('J0_blocks'       ,J0_blocks       )
    J0_steps         = default('J0_steps'        ,J0_steps        )
    test_warmupsteps = default('test_warmupsteps',test_warmupsteps)
    test_blocks      = default('test_blocks'     ,test_blocks     )
    test_steps       = default('test_steps'      ,test_steps      )
    
    if test:
        warmup = test_warmupsteps,
        blocks = test_blocks,
        steps  = test_steps,
    elif J0:
        warmup = J0_warmupsteps
        blocks = J0_blocks
        steps  = J0_steps
    else:
        warmup = warmupsteps
        blocks = blocks
        steps  = steps
    #end if
    vmc_calcs = [
        vmc(
            walkers     = walkers,
            warmupsteps = warmup,
            blocks      = blocks,
            steps       = steps,
            substeps    = substeps,
            timestep    = timestep,
            checkpoint  = checkpoint,
            )
        ]
    return vmc_calcs
#end def vmc_sections



def dmc_sections(
    warmupsteps          = missing,
    blocks               = missing,
    steps                = missing,
    timestep             = missing,
    checkpoint           = missing,
    vmc_samples          = missing,
    vmc_samplesperthread = missing, 
    vmc_walkers          = missing,
    vmc_warmupsteps      = missing,
    vmc_blocks           = missing,
    vmc_steps            = missing,
    vmc_substeps         = missing,
    vmc_timestep         = missing,
    vmc_checkpoint       = missing,
    eq_dmc               = missing,
    eq_warmupsteps       = missing,
    eq_blocks            = missing,
    eq_steps             = missing,
    eq_timestep          = missing,
    eq_checkpoint        = missing,
    J0_warmupsteps       = missing,
    J0_blocks            = missing,
    J0_steps             = missing,
    test_warmupsteps     = missing,
    test_blocks          = missing,
    test_steps           = missing,
    ntimesteps           = missing,
    timestep_factor      = missing,    
    nlmove               = missing,
    dmc_calcs            = missing,
    J0                   = False,
    test                 = False,
    defaults             = defaults_version,
    loc                  = 'dmc_sections',
    ):

    if not missing(dmc_calcs):
        return dmc_calcs
    #end if

    set_def_loc(dmc_sections_defaults[defaults],loc)

    warmupsteps      = default('warmupsteps'     ,warmupsteps     )
    blocks           = default('blocks'          ,blocks          )
    steps            = default('steps'           ,steps           )
    timestep         = default('timestep'        ,timestep        )
    checkpoint       = default('checkpoint'      ,checkpoint      )
    vmc_samples      = default('vmc_samples'     ,vmc_samples     )
    vmc_walkers      = default('vmc_walkers'     ,vmc_walkers     )
    vmc_warmupsteps  = default('vmc_warmupsteps' ,vmc_warmupsteps )
    vmc_blocks       = default('vmc_blocks'      ,vmc_blocks      )
    vmc_steps        = default('vmc_steps'       ,vmc_steps       )
    vmc_substeps     = default('vmc_substeps'    ,vmc_substeps    )
    vmc_timestep     = default('vmc_timestep'    ,vmc_timestep    )
    vmc_checkpoint   = default('vmc_checkpoint'  ,vmc_checkpoint  )
    eq_dmc           = default('eq_dmc'          ,eq_dmc          )
    eq_warmupsteps   = default('eq_warmupsteps'  ,eq_warmupsteps  )
    eq_blocks        = default('eq_blocks'       ,eq_blocks       )
    eq_steps         = default('eq_steps'        ,eq_steps        )
    eq_timestep      = default('eq_timestep'     ,eq_timestep     )
    eq_checkpoint    = default('eq_checkpoint'   ,eq_checkpoint   )
    J0_warmupsteps   = default('J0_warmupsteps'  ,J0_warmupsteps  )
    J0_blocks        = default('J0_blocks'       ,J0_blocks       )
    J0_steps         = default('J0_steps'        ,J0_steps        )
    test_warmupsteps = default('test_warmupsteps',test_warmupsteps)
    test_blocks      = default('test_blocks'     ,test_blocks     )
    test_steps       = default('test_steps'      ,test_steps      )
    ntimesteps       = default('ntimesteps'      ,ntimesteps      )
    timestep_factor  = default('timestep_factor' ,timestep_factor )    
    nlmove           = default('nlmove'          ,nlmove          )

    if missing(vmc_samples) and missing(vmc_samplesperthread):
        error('vmc samples (dmc walkers) not specified\nplease provide one of the following keywords: vmc_samples, vmc_samplesperthread',loc)
    #end if

    vsec = vmc(
        walkers     = vmc_walkers,
        warmupsteps = vmc_warmupsteps,
        blocks      = vmc_blocks,
        steps       = vmc_steps,
        substeps    = vmc_substeps,
        timestep    = vmc_timestep,
        checkpoint  = vmc_checkpoint,
        )
    if not missing(vmc_samplesperthread):
        vsec.samplesperthread = vmc_samplesperthread
    elif not missing(vmc_samples):
        vsec.samples = vmc_samples
    #end if

    if J0:
        warmupsteps = J0_warmupsteps
        blocks      = J0_blocks
        steps       = J0_steps
    elif test:
        warmupsteps = test_warmupsteps
        blocks      = test_blocks
        steps       = test_steps
    #end if

    dmc_calcs = [vsec]
    if eq_dmc:
        dmc_calcs.append(
            dmc(
                warmupsteps   = eq_warmupsteps,
                blocks        = eq_blocks,
                steps         = eq_steps,
                timestep      = eq_timestep,
                checkpoint    = eq_checkpoint,
                nonlocalmoves = nlmove,
                )
            )
    #end if
    tfac = 1.0
    for n in range(ntimesteps):
        sfac = 1.0/tfac
        dmc_calcs.append(
            dmc(
                warmupsteps   = int(sfac*warmupsteps),
                blocks        = blocks,
                steps         = int(sfac*steps),
                timestep      = tfac*timestep,
                checkpoint    = checkpoint,
                nonlocalmoves = nlmove,
                )
            )
        tfac *= timestep_factor
    #end for
    
    return dmc_calcs
#end def dmc_sections



def process_system_inputs(inputs,shared,task,loc):
    return inputs
#end def process_system_inputs


def process_vasp_inputs(inputs,shared,task,loc):
    inputs.set(
        system = shared.system,
        )
    return inputs
#end def process_vasp_inputs


def process_scf_inputs(inputs,shared,task,loc):
    if shared.run.nscf:
        if 'nosym' not in task.inputs_in:
            inputs.nosym = False
        #end if
        if 'wf_collect' not in task.inputs_in:
            inputs.wf_collect = False
        #end if
    #end if
    inputs.set(
        system  = shared.system,
        pseudos = get_pseudos('pwscf',shared),
        )
    return inputs
#end def process_scf_inputs


def process_nscf_inputs(inputs,shared,task,loc):
    inputs.set(
        nosym      = True,
        wf_collect = True,
        system     = shared.system,
        pseudos    = get_pseudos('pwscf',shared),
        )
    return inputs
#end def process_nscf_inputs


def process_p2q_inputs(inputs,shared,task,loc):
    return inputs
#end def process_p2q_inputs


def process_pwf_inputs(inputs,shared,task,loc):
    return inputs
#end def process_pwf_inputs


def process_pp_inputs(inputs,shared,task,loc):
    return inputs
#end def process_pp_inputs


def process_opt_inputs(inputs,shared,task,loc):
    inputs.set(
        system  = shared.system,
        pseudos = get_pseudos('qmcpack',shared),
        )

    set_loc(loc+'opt_inputs jastrows')
    jkw = extract_keywords(inputs,jastrow_factor_keys,optional=True)
    jkw.system   = shared.system
    j2kw = jkw.copy()
    j2kw.set(J1=1,J2=1,J3=0)
    j3kw = jkw.copy()
    j3kw.set(J1=1,J2=1,J3=1)
    task.J2_inputs = j2kw
    task.J3_inputs = j3kw

    set_loc(loc+'opt_inputs opt_methods')
    task.sec_inputs = extract_keywords(inputs,opt_sections_keys,optional=True)

    return inputs
#end def process_opt_inputs


def process_vmc_inputs(inputs,shared,task,loc):
    inputs.set(
        system  = shared.system,
        pseudos = get_pseudos('qmcpack',shared),
        )

    set_loc(loc+'vmc_inputs jastrows')
    jkw = extract_keywords(inputs,jastrow_factor_keys,optional=True)
    jkw.system   = shared.system
    j2kw = jkw.copy()
    j2kw.set(J1=1,J2=1,J3=0)
    j3kw = jkw.copy()
    j3kw.set(J1=1,J2=1,J3=1)
    task.J2_inputs = j2kw
    task.J3_inputs = j3kw

    set_loc(loc+'vmc_inputs vmc_methods')
    task.sec_inputs = extract_keywords(inputs,vmc_sections_keys,optional=True)
    return inputs
#end def process_vmc_inputs


def process_dmc_inputs(inputs,shared,task,loc):
    inputs.set(
        system  = shared.system,
        pseudos = get_pseudos('qmcpack',shared),
        )

    set_loc(loc+'dmc_inputs jastrows')
    jkw = extract_keywords(inputs,jastrow_factor_keys,optional=True)
    jkw.system   = shared.system
    j2kw = jkw.copy()
    j2kw.set(J1=1,J2=1,J3=0)
    j3kw = jkw.copy()
    j3kw.set(J1=1,J2=1,J3=1)
    task.J2_inputs = j2kw
    task.J3_inputs = j3kw

    set_loc(loc+'dmc_inputs dmc_methods')
    task.sec_inputs = extract_keywords(inputs,dmc_sections_keys,optional=True)
    return inputs
#end def process_dmc_inputs


sim_input_defaults = obj(
    system = system_input_defaults,
    vasp = vasp_input_defaults,
    scf  = scf_input_defaults,
    nscf = nscf_input_defaults,
    p2q  = p2q_input_defaults,
    pwf  = pwf_input_defaults,
    pp   = pp_input_defaults,
    opt  = opt_input_defaults,
    vmc  = vmc_input_defaults,
    dmc  = dmc_input_defaults,
    )
sim_workflow_keys = obj(
    system = system_workflow_keys,
    vasp = vasp_workflow_keys,
    scf  = scf_workflow_keys,
    nscf = nscf_workflow_keys,
    p2q  = p2q_workflow_keys,
    pwf  = pwf_workflow_keys,
    pp   = pp_workflow_keys,
    opt  = opt_workflow_keys,
    vmc  = vmc_workflow_keys,
    dmc  = dmc_workflow_keys,
    )
process_sim_inp = obj(
    system = process_system_inputs,
    vasp = process_vasp_inputs,
    scf  = process_scf_inputs,
    nscf = process_nscf_inputs,
    p2q  = process_p2q_inputs,
    pwf  = process_pwf_inputs,
    pp   = process_pp_inputs,
    opt  = process_opt_inputs,
    vmc  = process_vmc_inputs,
    dmc  = process_dmc_inputs,
    )

def process_sim_inputs(name,inputs_in,defaults,shared,task,loc):
    set_loc(loc)
    inputs = obj(**inputs_in)
    assign_defaults(inputs,sim_input_defaults[name][defaults])
    inputs = process_sim_inp[name](inputs,shared,task,loc)
    workflow      = extract_keywords(inputs,sim_workflow_keys[name])
    task.inputs   = inputs
    task.workflow = workflow
#end def process_sim_inputs


qmcpack_chain_sim_names = ['system','vasp','scf','nscf','p2q','pwf','pp','opt','vmc','dmc']
def capture_qmcpack_chain_inputs(kwargs):
    kw_deep    = obj()
    kw_shallow = obj()
    inp_shallow = obj()
    shallow_types = tuple([Simulation])#(Job,Simulation)
    for k,v in kwargs.iteritems():
        if '_inputs' not in k or isinstance(v,(tuple,list)):
            if isinstance(v,shallow_types):
                kw_shallow[k] = v
            else:
                kw_deep[k] = v
            #end if
        else:
            shallow = obj()
            deep    = obj()
            for kk,vv in v.iteritems():
                if isinstance(vv,shallow_types):
                    shallow[kk]=vv
                else:
                    deep[kk]=vv
                #end if
            #end for
            inp_shallow[k] = shallow
            kw_deep[k] = deep
        #end if
    #end for
    # deep copy
    kwcap = kw_deep.copy()
    # shallow copy
    kwcap.set(**kw_shallow)
    for k,v in inp_shallow.iteritems():
        kwcap[k].set(**v)
    #end for
    return kwcap
#end def capture_qmcpack_chain_inputs



def process_qmcpack_chain_kwargs(
    basepath        = '',
    orb_source      = None,
    J2_source       = None,
    J3_source       = None,
    system          = missing,
    system_function = missing,
    sim_list        = missing,
    pseudos         = missing,
    dft_pseudos     = missing,
    qmc_pseudos     = missing,
    loc             = 'process_qmcpack_chain_kwargs',
    valid_inputs    = None,
    fake            = False,
    sims            = None,
    **other_kwargs
    ):

    if valid_inputs is None:
        valid_inputs = []
    #end if
    valid_inputs.extend([
            'basepath','orb_source','J2_source','J3_source','system',
            'system_function','sim_list','pseudos','dft_pseudos',
            'qmc_pseudos','sims'
            ])
    for sname in qmcpack_chain_sim_names:
        valid_inputs.append(sname+'_inputs')
    #end for

    if 'system_inputs' in other_kwargs:
        system_inputs = other_kwargs['system_inputs']
    else:
        system_inputs = missing
    #end if

    sim_names = set(qmcpack_chain_sim_names)

    # collect other inputs into tasks
    #   organize all sim_inputs and sim_defaults input data
    invalid = []
    tasks = obj()
    task_defaults = obj()
    for name in sim_names:
        task_defaults[name] = defaults_version
    #end for
    for k,v in other_kwargs.iteritems():
        if '_inputs' in k:
            name,tmp = k.split('_',1)
            if name in sim_names:
                index = tmp.replace('inputs','')
                if len(index)==0:
                    index = 1
                else:
                    try:
                        index = int(index)
                    except:
                        error('inputs index must be an integer\ninput section: {0}\ninvalid index: {1}'.format(k,index),loc)
                    #end try
                #end if
                if not isinstance(v,(dict,obj)):
                    error('{0} must be a dict or obj\ntype received: {1}'.format(k,v.__class__.__name__),loc)
                #end if
                if name not in tasks:
                    tasks[name] = obj()
                #end if
                if index in tasks[name]:
                    error('repeated inputs index encountered\ninput section: {0}\ninvalid index: {1}'.format(k,index),loc)
                #end if
                if isinstance(v,SectionPlaceHolder):
                    tasks[name][index] = placeholder(inputs_in=obj())
                else:
                    tasks[name][index] = obj(inputs_in = obj(**v)) # capture the inputs
                #end if
            else:
                invalid.append(k)
            #end if
        elif k.endswith('_defaults') and k.split('_')[0] in sim_names:
            name = k.split('_')[0]
            if name in sim_names:
                task_defaults[name] = v
            else:
                invalid.append(k)
            #end if
        else:
            invalid.append(k)
        #end if
    #end for
    for name in sim_names:
        if name not in tasks:
            tasks[name] = None
        #end if
    #end for
    prevent_invalid_input(invalid,valid_inputs,loc)

    set_loc(loc)

    kw = obj()

    run = obj()
    for name in sim_names:
        run[name] = tasks[name]!=None
    #end for
    kw.run = run

    kw.basepath   = basepath
    kw.orb_source = orb_source
    kw.J2_source  = J2_source
    kw.J3_source  = J3_source
    kw.pseudos    = pseudos
    kw.fake       = fake
    kw.sims       = sims

    kw.system_inputs_provided = not missing(system_inputs)
    if not missing(system):
        kw.system = system
    elif not missing(system_inputs):
        if missing(system_function):
            system_function = generate_physical_system
        #end if
        system = system_function(**system_inputs)
        kw.system = system
    else:
        assign_require(kw,'system'     ,system       )
    #end if

    assign_require(kw,'sim_list'   ,sim_list     )

    if missing(pseudos):
        if run.scf or run.nscf:
            assign_require(kw,'dft_pseudos',dft_pseudos)
        #end if
        if run.opt or run.vmc or run.dmc:
            assign_require(kw,'qmc_pseudos',qmc_pseudos)
        #end if
    #end if

    for name in qmcpack_chain_sim_names: # order matters
        if run[name]:
            sim_tasks = tasks[name]
            prev_index = None
            for index in sorted(sim_tasks.keys()):
                sim_task = sim_tasks[index]
                if 'vary' in sim_task.inputs_in:
                    vary_index = sim_task.inputs_in.delete('vary')
                else:
                    vary_index = prev_index
                #end if
                if vary_index!=None:
                    inputs_in = obj(**sim_tasks[vary_index].inputs_in)
                    inputs_in.set(**sim_task.inputs_in)
                    sim_task.inputs_in = inputs_in
                #end if
                sim_task.vary = vary_index
                process_sim_inputs(
                    name      = name,
                    inputs_in = sim_task.inputs_in,
                    defaults  = task_defaults[name],
                    shared    = kw,
                    task      = sim_task,
                    loc       = '{0} {1}_inputs'.format(loc,name),
                    )
                if index==1:
                    sind = ''
                else:
                    sind = str(index)
                #end if
                label = name+sind
                sim_task.name  = name
                sim_task.label = label
                sim_task.index = index
                sim_task.workflow.set(
                    basepath    = basepath,
                    path        = os.path.join(basepath,label),
                    vary_index  = vary_index,
                    )
                prev_index = index
            #end for
        #end if
    #end for

    kw.tasks = tasks

    return kw
#end def process_qmcpack_chain_kwargs



def gen_vasp_chain(ch,loc):
    task = ch.task
    wf   = task.workflow
    vasp = generate_vasp(
        path = wf.path,
        **task.inputs
        )
    return vasp
#end def gen_vasp_chain


def gen_scf_chain(ch,loc):
    task = ch.task
    wf   = task.workflow
    scf = generate_pwscf(
        path = wf.path,
        **task.inputs
        )
    return scf
#end def gen_scf_chain


def gen_nscf_chain(ch,loc):
    task = ch.task
    wf   = task.workflow
    deps = resolve_deps('nscf',[('scf','charge_density')],ch,loc)
    # obtain user inputs to scf
    scf_label,scf_index = resolve_label('scf',ch,loc,ind=True)
    scf_inputs    = obj(ch.tasks.scf[scf_index].inputs)
    # pass inputs on to nscf after removing kpoint information
    scf_inputs.delete_optional('kgrid')
    scf_inputs.delete_optional('kshift')
    task.inputs.set_optional(**scf_inputs)
    # generate the nscf sim
    if wf.use_scf_dir:
        path = resolve_path('scf',ch,loc) # added to always run in scf dir
    else:
        path = wf.path
    #end if
    nscf = generate_pwscf(
        path         = path,
        dependencies = deps,
        **task.inputs
        )
    return nscf
#end def gen_nscf_chain


def gen_p2q_chain(ch,loc):
    kw   = ch.kw
    task = ch.task
    wf   = task.workflow
    run  = kw.run
    if run.nscf and wf.orb_src_type!='scf':
        if wf.orb_src is None:
            nscf_label,nscf_index = resolve_label('scf',ch,loc,ind=True)
        else:
            nscf_index = wf.orb_src
        #end if
        use_scf_dir = ch.tasks.nscf[nscf_index].workflow.use_scf_dir
        orb_source = 'nscf'
        if use_scf_dir:
            path_source = 'scf'
        else:
            path_source = 'nscf'
        #end if
    else:
        orb_source  = 'scf'
        path_source = 'scf'
    #end if
    path = resolve_path(path_source,ch,loc)
    #path = resolve_path('scf',ch,loc)  # added to always run in scf dir
    deps = resolve_deps('p2q',[(orb_source,'orbitals')],ch,loc)
    p2q = generate_pw2qmcpack(
        path         = path,
        dependencies = deps,
        **task.inputs
        )
    return p2q
#end def gen_p2q_chain


def gen_pwf_chain(ch,loc):
    kw   = ch.kw
    task = ch.task
    wf   = task.workflow
    run  = kw.run
    if run.nscf:
        source = 'nscf'
    else:
        source = 'scf'
    #end if
    path = resolve_path(source,ch,loc)
    deps = resolve_deps('pwf',[(source,'other')],ch,loc)
    pwf = generate_projwfc(
        path         = path,
        dependencies = deps,
        **task.inputs
        )
    return pwf
#end def gen_pwf_chain


def gen_pp_chain(ch,loc):
    kw   = ch.kw
    task = ch.task
    wf   = task.workflow
    run  = kw.run
    source = 'scf'
    path = resolve_path(source,ch,loc)
    deps = resolve_deps('pp',[(source,'other')],ch,loc)
    pp = generate_pp(
        path         = path,
        dependencies = deps,
        **task.inputs
        )
    return pp
#end def gen_pp_chain


def gen_opt_chain(ch,loc):
    kw   = ch.kw
    task = ch.task
    sims = ch.sims
    wf   = task.workflow
    orbdep = [('p2q','orbitals')]
    J2dep  = orbdep + [(task.label+'_J2','jastrow')]
    if wf.J2_src is None:
        optJ2_dep = orbdep
    else:
        optJ2_dep = orbdep + [('opt_J2','jastrow')]
    #end if
    if wf.J3_src is None:
        if wf.J2_run:
            optJ3_dep = J2dep
        else:
            optJ3_dep = orbdep
        #end if
    elif wf.use_J2:
        optJ3_dep = orbdep + [('opt_J2','jastrow')]
    else:
        optJ3_dep = orbdep + [('opt_J3','jastrow')]
    #end if
    if wf.J2_run:
        label = task.label+'_J2'
        deps = resolve_deps(label,optJ2_dep,ch,loc)
        optJ2 = generate_qmcpack(
            path         = os.path.join(wf.basepath,label),
            jastrows     = jastrow_factor(**task.J2_inputs),
            calculations = opt_sections(**task.sec_inputs),
            dependencies = deps,
            **task.inputs
            )
        sims[label] = optJ2
        kw.def_srcs[task.name+'_J2'] = task.index
    #end if
    if wf.J3_run:
        label = task.label+'_J3'
        deps = resolve_deps(label,optJ3_dep,ch,loc)
        optJ3 = generate_qmcpack(
            path         = os.path.join(wf.basepath,label),
            jastrows     = jastrow_factor(**task.J3_inputs),
            calculations = opt_sections(**task.sec_inputs),
            dependencies = deps,
            **task.inputs
            )
        sims[label] = optJ3
        kw.def_srcs[task.name+'_J3'] = task.index
    #end if
    return None
#end def gen_opt_chain


def gen_vmc(simlabel,ch,depset,J,test=0,loc=''):
    task = ch.task
    sims = ch.sims
    wf   = task.workflow
    deps = resolve_deps(simlabel,depset[J],ch,loc)
    if test:
        qmcjob = wf.test_job
    else:
        qmcjob = wf.job
    #end if
    other_inputs = obj(task.inputs)
    if 'jastrows' in other_inputs:
        jastrows = other_inputs.delete('jastrows')
    elif J=='J0':
        jastrows = []
    elif J=='J2':
        jastrows = jastrow_factor(**task.J2_inputs)
    elif J=='J3':
        jastrows = jastrow_factor(**task.J3_inputs)
    #end if
    qmc = generate_qmcpack(
        path         = os.path.join(wf.basepath,simlabel),
        job          = qmcjob,
        jastrows     = jastrows,
        calculations = vmc_sections(test=test,J0=J=='J0',**task.sec_inputs),
        dependencies = deps,
        **other_inputs
        )
    sims[simlabel] = qmc
#end def gen_vmc


def gen_vmc_chain(ch,loc):
    task = ch.task
    wf   = task.workflow
    orbdep = [('p2q','orbitals')]
    if wf.prior_opt:
        J2dep = orbdep + [('opt_J2','jastrow')]
        J3dep = orbdep + [('opt_J3','jastrow')]
    else:
        J2dep = orbdep
        J3dep = orbdep
    #end if
    depset = obj(
        J0 = orbdep,
        J2 = J2dep,
        J3 = J3dep,
        )
    lab = task.label
    if wf.J0_test:
        gen_vmc(lab+'_J0_test',ch,depset,J='J0',test=1,loc=loc)
    #end if
    if wf.J0_run:
        gen_vmc(lab+'_J0',ch,depset,J='J0',loc=loc)
    #end if
    if wf.J2_test:
        gen_vmc(lab+'_J2_test',ch,depset,J='J2',test=1,loc=loc)
    #end if
    if wf.J2_run:
        gen_vmc(lab+'_J2',ch,depset,J='J2',loc=loc)
    #end if
    if wf.J3_test:
        gen_vmc(lab+'_J3_test',ch,depset,J='J3',test=1,loc=loc)
    #end if
    if wf.J3_run:
        gen_vmc(lab+'_J3',ch,depset,J='J3',loc=loc)
    #end if
    return None
#end def gen_vmc_chain


def gen_dmc(simlabel,ch,depset,J,nlmove=None,test=0,loc=''):
    task = ch.task
    sims = ch.sims
    wf   = task.workflow
    deps = resolve_deps(simlabel,depset[J],ch,loc)
    if test:
        qmcjob = wf.test_job
    else:
        qmcjob = wf.job
    #end if
    other_inputs = obj(task.inputs)
    if 'jastrows' in other_inputs:
        jastrows = other_inputs.delete('jastrows')    
    elif J=='J0':
        jastrows = []
    elif J=='J2':
        jastrows = jastrow_factor(**task.J2_inputs)
    elif J=='J3':
        jastrows = jastrow_factor(**task.J3_inputs)
    #end if
    if 'calculations' not in other_inputs:
        other_inputs.calculations = dmc_sections(nlmove=nlmove,test=test,J0=J=='J0',**task.sec_inputs)
    #end if
    qmc = generate_qmcpack(
        path         = os.path.join(wf.basepath,simlabel),
        job          = qmcjob,
        jastrows     = jastrows,
        dependencies = deps,
        **other_inputs
        )
    sims[simlabel] = qmc
#end def gen_dmc


def gen_dmc_chain(ch,loc):
    kw   = ch.kw
    task = ch.task
    wf   = task.workflow
    orbdep = [('p2q','orbitals')]
    if wf.prior_opt:
        J2dep = orbdep + [('opt_J2','jastrow')]
        J3dep = orbdep + [('opt_J3','jastrow')]
    else:
        J2dep = orbdep
        J3dep = orbdep
    #end if
    depset = obj(
        J0 = orbdep,
        J2 = J2dep,
        J3 = J3dep,
        )
    nonloc_labels = {None:'',True:'_tm',False:'_la'}
    nonlocalmoves_default = None
    if kw.system.pseudized:
        nonlocalmoves_default = True
    #end if
    nloc_moves = []
    if wf.tmoves:
        nloc_moves.append(True)
    #end if
    if wf.locality:
        nloc_moves.append(False)
    #end if
    if not wf.tmoves and not wf.locality:
        nloc_moves.append(nonlocalmoves_default)
    #end if
    lab = task.label
    for nlmove in nloc_moves:
        nll = nonloc_labels[nlmove]
        if wf.J0_test:
            gen_dmc(lab+'_J0'+nll+'_test',ch,depset,J='J0',nlmove=nlmove,test=1,loc=loc)
        #end if
        if wf.J0_run:
            gen_dmc(lab+'_J0'+nll,ch,depset,J='J0',nlmove=nlmove,loc=loc)
        #end if
        if wf.J2_test:
            gen_dmc(lab+'_J2'+nll+'_test',ch,depset,J='J2',nlmove=nlmove,test=1,loc=loc)
        #end if
        if wf.J2_run:
            gen_dmc(lab+'_J2'+nll,ch,depset,J='J2',nlmove=nlmove,loc=loc)
        #end if
        if wf.J3_test:
            gen_dmc(lab+'_J3'+nll+'_test',ch,depset,J='J3',nlmove=nlmove,test=1,loc=loc)
        #end if
        if wf.J3_run:
            gen_dmc(lab+'_J3'+nll,ch,depset,J='J3',nlmove=nlmove,loc=loc)
        #end if
    #end for
    return None
#end def gen_dmc_chain


gen_sim_ch = obj(
    vasp = gen_vasp_chain,
    scf  = gen_scf_chain,
    nscf = gen_nscf_chain,
    p2q  = gen_p2q_chain,
    pwf  = gen_pwf_chain,
    pp   = gen_pp_chain,
    opt  = gen_opt_chain,
    vmc  = gen_vmc_chain,
    dmc  = gen_dmc_chain,
    )

def gen_sim_chain(name,ch,loc):
    kw        = ch.kw
    tasks     = ch.tasks
    sims      = ch.sims
    sim_tasks = tasks[name]
    for index in sorted(sim_tasks.keys()):
        task = sim_tasks[index]
        if not isinstance(task,SectionPlaceHolder):
            ch.set(
                name  = name,
                index = index,
                task  = task,
                )
            sim = gen_sim_ch[name](ch,loc)
            if sim is not None:
                sims[task.label] = sim
            #end if
        #end if
        kw.def_srcs[task.name] = task.index # only the last index will matter
    #end for
#end def gen_sim_chain



def qmcpack_chain(**kwargs):
    loc      = kwargs.pop('loc','qmcpack_chain')

    kw       = process_qmcpack_chain_kwargs(loc=loc,**kwargs)

    sim_list = kw.sim_list
    run      = kw.run
    tasks    = kw.tasks
    
    if kw.fake:
        Simulation.creating_fake_sims = True
    #end if

    if kw.sims is None:
        sims = SimSet()
    else:
        sims = kw.sims
    #end if
    kw.def_srcs = obj()

    ch = obj(
        kw    = kw,
        tasks = tasks,
        sims  = sims,
        )

    if run.vasp:
        gen_sim_chain('vasp',ch,loc)
    #end if

    if kw.orb_source!=None:
        sims.p2q = kw.orb_source
    else:
        if run.scf:
            gen_sim_chain('scf',ch,loc)
        #end if
        if run.nscf:
            gen_sim_chain('nscf',ch,loc)
        #end if
        if run.p2q:
            gen_sim_chain('p2q',ch,loc)
        #end if
        if run.pwf:
            gen_sim_chain('pwf',ch,loc)
        #end if
        if run.pp:
            gen_sim_chain('pp',ch,loc)
        #end if
    #end if

    if run.opt:
        gen_sim_chain('opt',ch,loc)
    #end if
    if kw.J2_source!=None:
        sims.opt_J2 = kw.J2_source
    #end if
    if kw.J3_source!=None:
        sims.opt_J3 = kw.J3_source
    #end if

    if run.vmc:
        gen_sim_chain('vmc',ch,loc)
    #end if

    if run.dmc:
        gen_sim_chain('dmc',ch,loc)
    #end if

    # apply labeling to simulations
    for simlabel,sim in sims.iteritems():
        sim.simlabel = simlabel
    #end for

    if not kw.fake:
        sim_list.extend(sims.list())
    else:
        Simulation.creating_fake_sims = False
    #end if

    # for scanned workflows only
    #   add fake system sim
    #   make all other sims depend on it
    if kw.fake and kw.system_inputs_provided:
        sys = SystemHolder(system=kw.system)
        for sim in sims:
            if len(sim.dependencies)==0:
                sim.depends(sys,'other')
            #end if
        #end for
        sims.system = sys
    #end if

    #print kw.def_srcs
    #exit()

    return sims
#end def qmcpack_chain



def unpack_scan_inputs(scan,loc='unpack_scan_inputs'):
    if not isinstance(scan,(tuple,list)):
        error('parameter "scan" must be a list or tuple\ntype received: {0}'.format(scan.__class__.__name__),loc=loc)
    #end if
    scans = obj()
    for scan_list in scan:
        if not isinstance(scan_list,(tuple,list)):
            error('scan list for parameter "scan" must be a tuple or list\ntype received: {0}'.format(scan_list.__class__.__name__),loc=loc)
        #end if
        misformatted = len(scan_list)<3 
        misformatted|= (len(scan_list)-1)%2!=0 or not isinstance(scan_list[0],str)
        if not misformatted:
            section = scan_list[0]
            vars_in = scan_list[1::2]
            vals_in = scan_list[2::2]
            all_vars = set(vars_in)
            vars = []
            values_in = obj()
            labels_in = obj()
            for var,val_list in zip(vars_in,vals_in):
                misformatted |= not isinstance(var,str) or not isinstance(val_list,(tuple,list,ndarray))
                if isinstance(var,str):
                    if var.endswith('_labels'):
                        vname = var.rsplit('_',1)[0]
                        if vname in all_vars:
                            labels_in[vname] = val_list
                        else:
                            values_in[var] = val_list
                            vars.append(var)
                        #end if
                    else:
                        values_in[var] = val_list
                        vars.append(var)
                    #end if
                #end if
            #end for
            vals = []
            labs = []
            for var in vars_in:
                if var in values_in:
                    vin = values_in[var]
                    vals.append(vin)
                    if var in labels_in:
                        lin = labels_in[var]
                        if len(lin)!=len(vin):
                            error('scan list for parameter "scan" is formatted improperly\nincorrect number of labels provided for scan parameter "{0}"\nnumber of parameter values provided: {1}\nnumber of parameter labels provided: {2}\nparameter labels provided: {3}'.format(var,len(vin),len(lin),lin),loc=loc)
                        #end if
                        labs.append(lin)
                    else:
                        labs.append(len(vin)*[None])
                    #end if
                #end if
            #end for
            del vars_in
            del vals_in
            if not misformatted:
                if len(vars)==1:
                    vals_in = vals[0]
                    labs_in = labs[0]
                    vals = []
                    inds = []
                    labs = []
                    i = 1
                    for v,l in zip(vals_in,labs_in):
                        vals.append((v,))
                        inds.append((i,))
                        labs.append((l,))
                        i+=1
                    #end for
                else:
                    vars = tuple(vars)
                    n=1
                    for vals_in in vals[1:]:
                        if len(vals_in)!=len(vals[0]):
                            error('problem with "scan" input\nall value_lists for section "{0}" must have the same length (they are covarying)\nvar_list "{1}" has length {2}\nvar_list "{3}" has length {4}'.format(section,vars[0],len(vals[0]),vars[n],len(vals[n])),loc=loc)
                        #end if
                        n+=1
                    #end for
                    vals = zip(*vals)
                    labs = zip(*labs)
                    r = range(1,len(vals)+1)
                    inds = zip(*[r for n in range(len(vars))])
                #end if
                if section not in scans:
                    sec = obj(parameters=list(vars),values=vals,indices=inds,labels=labs,values_in=values_in,labels_in=labels_in)
                    scans[section] = sec
                else:
                    sec = scans[section]
                    sec.values_in.transfer_from(values_in)
                    sec.labels_in.transfer_from(labels_in)
                    sec.parameters.extend(vars)
                    values = []
                    for v in sec.values:
                        for v2 in vals:
                            val = tuple(list(v)+list(v2))
                            values.append(val)
                        #end for
                    #end for
                    sec.values = values
                    indices = []
                    for i in sec.indices:
                        for i2 in inds:
                            ind = tuple(list(i)+list(i2))
                            indices.append(ind)
                        #end for
                    #end for
                    sec.indices = indices
                    labels = []
                    for l in sec.labels:
                        for l2 in labs:
                            lab = tuple(list(l)+list(l2))
                            labels.append(lab)
                        #end for
                    #end for
                    sec.labels = labels
                #end if
            #end if
        #end if
        if misformatted:
            error('scan list for parameter "scan" is formatted improperly\nlist must have input section name (string) followed by variable-value_list pairs\nnumber of entries present (should be odd): {0}\nvalue_lists must be of type tuple/list/array\nscan list contents: {1}'.format(len(scan_list),scan_list,loc))
        #end if
    #end for

    return scans
#end def unpack_scan_inputs


def unpack_fix_inputs(fix,loc='unpack_fix_inputs'):
    if not isinstance(fix,(tuple,list)):
        error('parameter "fix" must be a list or tuple\ntype received: {0}'.format(fix.__class__.__name__),loc=loc)
    #end if
    constraints = obj()
    for clist in fix:
        if not isinstance(clist,(tuple,list)):
            error('constraint list for parameter "fix" must be a tuple or list\ntype received: {0}'.format(clist.__class__.__name__),loc=loc)
        #end if
        misformatted = len(clist)<3 
        misformatted|= (len(clist)-1)%2!=0 or not isinstance(clist[0],str)
        if not misformatted:
            sim_name  = clist[0]
            variables = list(clist[1::2])
            selectors = list(clist[2::2])
            for var in variables:
                misformatted |= not isinstance(var,str)
            #end for
            if not misformatted:
                if variables[-1]=='pattern':
                    pattern = list(selectors[-1])
                    variables.pop()
                    selectors.pop()
                    if len(pattern)!=len(variables):
                        error('must have a single pattern element for each variable in contraint list for parameter "fix"\nnumber of pattern elements: {0}\nnumber of variables: {1}\npattern list: {2}\nvariable list: {3}'.format(len(pattern),len(variables),pattern,variables),loc=loc)
                    #end if
                    pattern_misformatted = False
                    pstr = ''
                    prob = ''
                    for v,s,p in zip(variables,selectors,pattern):
                        bad_int = p=='i' and not isinstance(s,int)
                        bad_val = p=='v' and not isinstance(s,(str,int,float))
                        bad_pat = p not in ('v','i')
                        if bad_int or bad_val or bad_pat:
                            pattern_misformatted = True
                            pstr+='({0})'.format(p)
                            prob+='{0} (invalid pattern {1})'.format(v,p)
                        else:
                            pstr+=p
                        #end if
                    #end for
                    if pattern_misformatted:
                        error('constraint values do not match pattern provided\nthe problem variables are: {0}\npattern provided: {1}\nvariables provided: {2}\nfor pattern "v" only an integer,float, or string can be provided\nfor pattern "i" only an integer can be provided'.format(prob,pstr,variables))
                    #end if
                else:
                    pattern = list(len(variables)*'v')
                #end if
            #end if
        #end if
        if misformatted:
            error('constraint list for parameter "fix" is formatted improperly\nlist must have sim name (string) followed by variable-value/index (string/int etc.) pairs\nnumber of entries present (should be odd): {0}\nlist contents: {1}'.format(len(clist),clist,loc))
        #end if
        constraints[sim_name] = obj(
            sim_name  = sim_name,
            variables = variables,
            selectors = selectors,
            pattern   = pattern,
            )
    #end for
    return constraints
#end def unpack_fix_inputs



class WFRep(DevBase):
    @classmethod
    def get_id(cls):
        if not cls.class_has('rep_count'):
            cls.class_set(rep_count=0)
        #end if
        id = cls.rep_count
        cls.rep_count+=1
        return id
    #end def get_id

    def __init__(self,label=None):
        self.identifier   = label
        self.id           = self.__class__.get_id()
        self.label        = label
        self.graphid      = self.id
        self.depth        = 0
        self.dependents   = obj()
        self.dependencies = obj()
    #end def __init__

    def add_dependent(self,rep):
        self.dependents[rep.id]   = rep
        rep.dependencies[self.id] = self
    #end def add_dependent


    def assign_depth(self,depth=None):
        if depth is None:
            depth = self.depth
        #end if
        self.depth = max(self.depth,depth)
        for rep in self.dependents:
            rep.assign_depth(depth+1)
        #end for
    #end def assign_depth


    def assign_graphid_down(self,graphid=None):
        if graphid is None:
            graphid = self.graphid
        else:
            self.graphid = graphid
        #end if
        for rep in self.dependents:
            rep.assign_graphid_down(graphid)
        #end for
    #end def assign_graphid_down


    def assign_graphid_up(self,graphid=None):
        if graphid is None:
            graphid = self.graphid
        else:
            self.graphid = graphid
        #end if
        for rep in self.dependencies:
            rep.assign_graphid_up(graphid)
        #end for
    #end def assign_graphid_up


    def downstream_labels(self,collection=None):
        if collection is None:
            collection = set()
        #end if
        for rep in self.dependents:
            collection.add(rep.label)
        #end for
        for rep in self.dependents:
            rep.downstream_labels(collection)
        #end for
        return collection
    #end def downstream_labels


    def upstream_labels(self,collection=None):
        if collection is None:
            collection = set()
        #end if
        for rep in self.dependencies:
            collection.add(rep.label)
        #end for
        for rep in self.dependencies:
            rep.upstream_labels(collection)
        #end for
        return collection
    #end def upstream_labels

    @property
    def simid(self):
        return self.id
    #end def simid

    @property 
    def simlabel(self):
        return self.label
    #end def simlabel
#end class WFRep


class SimRep(WFRep):
    def __init__(self,sim):
        WFRep.__init__(self)
        self.transfer_from(sim,['identifier'])
        self.id           = sim.simid
        self.label        = sim.simlabel
        self.graphid      = self.simid
        self.secrep       = None
        label = self.simlabel
        if '_' in label:
            label = label.split('_',1)[0]
        #end if
        n=len(label)
        while label[n-1].isdigit():
            n=-1
        #end while
        # get the input section
        #   the association between label and section will not always hold
        #   jtk mark: need to replace this with something better
        #             probably should track this during qmcpack_chain
        self.section = label[:n]+'_inputs'+label[n:]
        #print self.simlabel,self.section
    #end def __init__            
#end class SimRep


class SecRep(WFRep):
    def __init__(self,section,simreps_in):
        WFRep.__init__(self,section)
        self.segrep = None
        sims = obj()
        for sim in simreps_in:
            if sim.section==section:
                sims[sim.simid] = sim
                sim.secrep = self
            #end if
        #end for
        self.simreps = sims
    #end def __init__

    @property
    def section(self):
        return self.label
    #end def section
#end class SecRep


class SegRep(WFRep):
    def __init__(self,label,secreps,invariant=False):
        WFRep.__init__(self,label)
        self.invariant = invariant
        self.secreps = obj()
        for sec in secreps:
            self.secreps[sec.label] = sec
            sec.segrep = self
        #end for
        self.scan_data        = None # non-invariant only
        self.fixed_inputs     = None 
        self.reference_inputs = None
    #end def __init__


    def set_reference(self,active_sections,scans,kwargs,placeholders):
        if not self.invariant:
            self.scan_data = scans[self.label]
        #end if
        fixed = obj()
        for k,v in kwargs.iteritems():
            if k not in active_sections:
                fixed[k] = v
            #end if
        #end for
        self.fixed_inputs = fixed
        ref = obj()
        for section in self.secreps.keys():
            ref[section] = obj(**kwargs[section]) # shallow copy inputs data
            #ref[section].delete_optional('job')
        #end for
        self.placeholders = obj()
        for section,ph in placeholders.iteritems():
            if section not in ref:
                self.placeholders[section] = ph
            #end if
        #end for
        self.reference_inputs = ref
    #end def set_reference


    def check_constraint(self,constraint,tol,first=True,loc='check_constraint'):
        if first:
            constraint.indices = obj()
            constraint.check = obj(
                variables_found = set(),
                matches_found   = set(),
                all_upstream_variables = set(),
                )
        #end if
        sd = self.scan_data
        if sd is not None and len(sd.parameters)>0:
            c = constraint
            for v,s,p in zip(c.variables,c.selectors,c.pattern):
                for k in sd.values_in.keys():
                    c.check.all_upstream_variables.add(k)
                #end for
                if v in sd.values_in:
                    values = sd.values_in[v]
                    c.check.variables_found.add(v)
                    if p=='i':
                        if s>0 and s<len(values)+1:
                            c.check.matches_found.add(v)
                            c.indices[v] = s
                        else:
                            error('problem with "fix" parameter input\nindex provided for constraint variable {0} is out of range\n{1} values provided for {0} in section {2}\nindex must be in the range [1,{1}]\nindex provided: {2}\n'.format(v,len(values),self.label,s),loc=loc)
                        #end if
                    elif p=='v':
                        if isinstance(s,(int,str)):
                            if s in values:
                                c.check.matches_found.add(v)
                                c.indices[v] = values.index(s)+1
                            else:
                                error('problem with "fix" parameter input\nmatch not found for constraint variable "{0}"\nvalue provided: {1}\nvalues present in section {2}: {3}'.format(v,s,self.label,values),loc=loc)
                            #end if
                        elif isinstance(s,float):
                            ftol = min(tol,s*1e-2)
                            found = False
                            index = 1
                            for val in values:
                                if abs(s-val)<ftol:
                                    c.check.matches_found.add(v)
                                    c.indices[v] = index
                                    found=True
                                    break
                                #end if
                                index+=1
                                #end if
                            #end for
                            if not found:
                                error('problem with "fix" parameter input\nmatch not found for constraint variable "{0}"\nvalue provided: {1}\nvalues present in section {2}: {3}\ntolerance provided ("fix_tol"): {4}\ntolerance used: {5}'.format(v,s,self.label,values,tol,ftol),loc=loc)
                            #end if
                        else:
                            error('problem with "fix" parameter input\ninvalid type provided for constraint pattern "v"\nconstraint variable: {0}\ntype provided: {1}\npattern for this type must be "i" (provide an integer instead of type {1})'.format(v,s.__class__.__name__),loc=loc)
                        #end if
                #end if
            #end for
        #end if
        for segrep in self.dependencies:
            segrep.check_constraint(constraint,tol,False,loc)
        #end for
    #end def check_constraint


    def generate_workflow(self,basepath=None,cur_sims=None,cur_inds=None,cur_vals=None,sim_coll=None,system_inputs=None):
        if self.invariant:
            basepath = self.fixed_inputs.basepath
        elif basepath is None:
            self.error('first call must be with invariant segment')
        #end if
        # compose inputs for the current sub-chain/segment
        #  start with reference inputs
        chain_inputs = obj(**self.fixed_inputs)
        for section,inputs in self.reference_inputs.iteritems():
            chain_inputs[section] = obj(**inputs)
        #end for
        chain_inputs.transfer_from(self.placeholders)
        # enable passing of system_inputs
        if system_inputs!=None:
            chain_inputs.system_inputs = system_inputs
        #end if
        new_sims = SimSet()
        if self.invariant:
            sim_coll = SimColl()
            cinds = obj()
            cvals = obj()
            # make invariant part of workflow
            any_present = False
            for sname in qmcpack_chain_sim_names:
                iname = sname+'_inputs'
                if iname in chain_inputs and not isinstance(chain_inputs[iname],SectionPlaceHolder):
                    any_present = True
                #end if
            #end if
            if any_present:
                csims = qmcpack_chain(**chain_inputs)
                new_sims.transfer_from(csims)
                sim_coll.add_sims(csims)
            else:
                csims = SimSet()
            #end if
            # make scanned parts of workflow
            for seg in self.dependents:
                nsims = seg.generate_workflow(basepath,csims,cinds,cvals,sim_coll)
                new_sims[seg.label] = nsims
            #end for
        else:
            # scan over inputs and make current sub-chains
            scan_params  = self.scan_data.parameters
            scan_values  = self.scan_data.values
            scan_indices = self.scan_data.indices
            scan_labels  = self.scan_data.labels
            sec_vary = chain_inputs[self.label]
            cur_sim_keys = set(cur_sims.keys())
            n = 0
            for vset in scan_values:
                iset = scan_indices[n]
                lset = scan_labels[n]
                bpath = basepath
                cinds = obj(cur_inds)
                cvals = obj(cur_vals)
                akey = []
                for k,v,i,l in zip(scan_params,vset,iset,lset):
                    if l is not None:
                        dlabel = l
                        dkey   = l
                        dname  = l
                    elif isinstance(v,(tuple,list,ndarray)):
                        all_basic = True
                        dlabel = ''
                        for vv in v:
                            all_basic &= isinstance(vv,(str,int,float))
                            dlabel += str(vv)+'_'
                        #end if
                        if all_basic and len(v)>0:
                            dkey   = tuple(v)
                            dlabel = dlabel[:-1]
                        else:
                            dkey   = i
                            dlabel = i
                        #end if
                        dname = '{0}_{1}'.format(k,dlabel)
                    else:
                        if isinstance(v,(str,int,float)):
                            dkey = v
                        else:
                            dkey = i
                        #end if
                        dlabel = dkey
                        dname = '{0}_{1}'.format(k,dlabel)
                    #end if
                    bpath = os.path.join(bpath,dname)
                    akey.append(dkey)
                    sec_vary[k] = v
                    cinds[k] = i
                    cvals[k] = v
                #end for
                if len(akey)==1:
                    akey = akey[0]
                else:
                    akey = tuple(akey)
                #end if
                chain_inputs.basepath = bpath
                # make current part of workflow
                qmcpack_chain(sims=cur_sims,**chain_inputs)
                # add new sims to all_sims
                new = set(cur_sims.keys())-cur_sim_keys
                nsims = SimSet()
                for k in new:
                    nsims[k] = cur_sims[k]
                #end for
                sim_coll.add_sims(nsims,cinds,cvals)
                # pass along system_inputs if scanning over it
                if self.label=='system_inputs':
                    system_inputs = sec_vary
                #end if
                # make scanned sub-workflows
                for seg in self.dependents:
                    ns = seg.generate_workflow(bpath,cur_sims,cinds,cvals,sim_coll,system_inputs=system_inputs)
                    nsims[seg.label] = ns
                #end for
                new_sims[akey] = nsims
                # reset current sims to state at beginning of loop
                #  ie remove newly generated sims
                rem = set(cur_sims.keys())-cur_sim_keys
                for k in rem:
                    del cur_sims[k]
                #end for
                n+=1
            #end for
        #end if
        if self.invariant:
            #sim_coll.report()
            return new_sims,sim_coll
        else:
            return new_sims
        #end if
    #end def generate_workflow
#end class SegRep


# full sim collection organized by identifier
#   used to merge simulations together based on "fix" paramater input
class SimColl(DevBase):
    def add_sims(self,sims,indices=None,values=None):
        for sim in sims:
            label = sim.simlabel
            if label not in self:
                self[label] = []
            #end if
            self[label].append((sim,indices,values))
        #end for
    #end def add_sims

    def report(self):
        for label in sorted(self.keys()):
            print label
            for (sim,inds,vals) in self[label]:
                print ' ',inds
                print ' ',vals
            #end for
        #end for
    #end def report
#end class SimColl



def qmcpack_workflow(
    scan           = missing,
    fix            = missing,
    bundle         = missing,
    parameter      = missing,
    values         = missing,
    fix_value      = missing,
    fix_index      = missing,
    fix_tol        = 1e-6,
    graph_workflow = False,
    write_workflow = False,
    sim_list       = missing,
    # mostly for dev work
    graph_simreps  = False,
    write_simreps  = False,
    graph_secreps  = False,
    write_secreps  = False,
    graph_segreps  = False,
    write_segreps  = False,
    **kwargs
      ):
    loc  = kwargs.pop('loc' ,'qmcpack_workflow')
    kwargs['loc'] = loc

    set_loc(loc)

    require('sim_list',sim_list)
    loc_sims = []
    kwargs['sim_list'] = loc_sims

    #tpoint('qmcpack_workflow')

    valid_inputs = ['scan','parameter','values']
    kwargs['valid_inputs'] = valid_inputs

    # do a "chain" workflow if not scanning
    if missing(scan):
        sims = qmcpack_chain(**kwargs)
        if graph_workflow:
            graph_sims(loc_sims,exit=False)
        #end if
        if write_workflow:
            write_sims(loc_sims,exit=False)
        #end if
        sim_list.extend(loc_sims)
        return sims
    #end if


    # handle input for a single parameter scan
    if isinstance(scan,str):
        require('parameter',parameter)
        require('values'   ,values   )
        scan = [(scan,parameter,values)]
    #end if

    # handle input for a single fix/merge request
    if isinstance(fix,str):
        require('parameter',parameter)        
        require_any('fix_value',fix_value,'fix_index',fix_index)
        if not missing(fix_value):
            fix = [(fix,parameter,fix_value)]
        elif not missing(fix_index):
            fix = [(fix,parameter,fix_index,'pattern','i')]
        #end if
    #end if
            
    # unpack scan information
    scans = unpack_scan_inputs(scan,loc=loc)

    # unpack fix/merge information
    if not missing(fix):
        constraints = unpack_fix_inputs(fix,loc=loc)
    #end if

    # handle bundle input
    if not missing(bundle):
        if isinstance(bundle,str):
            bundle = bundle.split()
        elif not isinstance(bundle,(tuple,list)):
            error('problem with "bundle" parameter input\nmust be a list of simulation names to bundle together as a single job\nreceived type: {0}\nwith value: {1}'.format(bundle.__class__.__name__,bundle),loc=loc)
        else:
            for b in bundle:
                if not isinstance(b,str):
                    error('problem with "bundle" parameter input\nmust be a list of simulation names to bundle together as a single job\nsome names provided are not simple strings\nreceived type {0}\nwith value: {1}\nall values provided: {2}'.format(b.__class__.__name__,b,bundle),loc=loc)
                #end if
            #end for
        #end if
        bundle = list(set(bundle))
    #end if

    # make placeholders for input sections to keep track of upstream tasks
    #   when varying ones downstream
    section_placeholders = obj()
    for k,v in kwargs.iteritems():
        if 'inputs' in k and k!='valid_inputs':
            section_placeholders[k] = placeholder()
        #end if
    #end for


    # obtain an example simulation workflow (chain)
    #   convert sims into simpler representation class 
    kwcap = capture_qmcpack_chain_inputs(kwargs)
    sims = qmcpack_chain(fake=True,**kwcap)
    for sim in sims:
        sim.simrep = SimRep(sim)
    #end for
    for sim in sims:
        for dsim in sim.dependents:
            sim.simrep.add_dependent(dsim.simrep)
        #end for
    #end for
    simreps = obj()
    for sim in sims:
        simrep = sim.simrep
        del sim.simrep
        simreps[simrep.simlabel] = simrep
    #end for
    del sims
    for sim in simreps:
        if len(sim.dependencies)==0:
            sim.assign_depth()
        #end if
    #end for

    # determine which input sections are present/active
    active_sections = set()
    for sim in simreps:
        active_sections.add(sim.section)
    #end for

    # coarse grain the sim workflow into an input section representation
    secreps = obj()
    for section in active_sections:
        secreps[section] = SecRep(section,simreps)
    #end for
    for sim in simreps:
        for dsim in sim.dependents:
            if sim.secrep.id!=dsim.secrep.id:
                sim.secrep.add_dependent(dsim.secrep)
            #end if
        #end for
    #end for
    for sec in secreps:
        if len(sec.dependencies)==0:
            sec.assign_depth()
        #end if
    #end for


    # coarse grain the section workflow into a segmented representation
    #   each segment corresponds to scan variation points on the input graph
    #   i.e. input sections within a segment covary w.r.t scan parameters
    remaining_sections = set(active_sections)
    # determine which sections are downstream from others
    downstream_sections = obj()
    for sec in secreps:
        downstream_sections[sec.section] = sec.downstream_labels()
    #end for
    # order section workflow by depth
    depth_list = obj()
    for sec in secreps:
        depth = sec.depth
        if depth not in depth_list:
            depth_list[depth] = []
        #end if
        depth_list[depth].append(sec.section)
    #end for
    # collect sections into segments by moving in reverse depth order
    segreps = obj()
    scan_sections = set(scans.keys())
    for d in reversed(sorted(depth_list.keys())):
        level_sections = depth_list[d]
        for scan_section in scan_sections:
            if scan_section in level_sections:
                ds_sections = downstream_sections[scan_section]
                seg_secreps = [secreps[scan_section]]
                for section in ds_sections:
                    if section in remaining_sections:
                        seg_secreps.append(secreps[section])
                    #end if
                #end for
                seg = SegRep(scan_section,seg_secreps)
                segreps[seg.label] = seg
                remaining_sections -= ds_sections
                if scan_section in remaining_sections:
                    remaining_sections.remove(scan_section)
                #end if
            #end if
        #end for
    #end for
    # make a final segment for the unvarying input sections
    seg_secreps = []
    for section in remaining_sections:
        seg_secreps.append(secreps[section])
    #end for
    segreps.invariant = SegRep('invariant',seg_secreps,invariant=True)
    # assign segment dependents and depth
    for sec in secreps:
        for dsec in sec.dependents:
            if sec.segrep.id!=dsec.segrep.id:
                sec.segrep.add_dependent(dsec.segrep)
            #end if
        #end for
    #end for
    # invariant segment will serve as root node for variations
    #   so make other segments "depend" on it even if they don't
    for seg in segreps:
        if len(seg.dependencies)==0 and seg.label!='invariant':
            segreps.invariant.add_dependent(seg)
        #end if
    #end for
    for seg in segreps:
        if len(seg.dependencies)==0:
            seg.assign_depth()
        #end if
    #end for

    # print input sections comprising each segment
    #for name in sorted(segreps.keys()):
    #    seg = segreps[name]
    #    print name
    #    for section in sorted(seg.secreps.keys()):
    #        print ' ',section
    #    #end for
    ##end for

    #tpoint('coarse rep')

    # expand chain inputs
    #   this allows information to flow e.g. from *_inputs to *_inputs2 
    kwcap = capture_qmcpack_chain_inputs(kwargs)
    kwp = process_qmcpack_chain_kwargs(**kwcap)
    for name,tasklist in kwp.tasks.iteritems():
        #print name
        if tasklist is not None:
            for index,task in tasklist.iteritems():
                if index==1:
                    sec_name = '{0}_inputs'.format(name)
                else:
                    sec_name = '{0}_inputs{1}'.format(name,index)
                #end if
                #print ' ',sec_name
                #print '   ',sorted(kwargs[sec_name].keys())
                #print '   ',sorted(task.inputs_in.keys())
                kwargs[sec_name] = obj(**task.inputs_in)
            #end for
        #end if
    #end for

    #tpoint('expand chain inputs')

    # gather reference and scan inputs into each segment
    for seg in segreps:
        seg.set_reference(active_sections,scans,kwargs,section_placeholders)
    #end for

    #tpoint('set scan reference')


    # make a final check on scan inputs
    if not missing(scan):
        for sec_name,sec in scans.iteritems():
            if sec_name not in secreps:
                error('problem with "scan" parameter input\nno user input was provided for section "{0}"\ninput sections provided: {1}'.format(sec_name,sorted(secreps.keys())),loc=loc)
            #end if
        #end for
    #end if

    # make a final check on fix inputs
    #  and find matching indices for constraints
    if not missing(fix):
        for sim_name,constraint in constraints.iteritems():
            if sim_name not in simreps:
                error('problem with "fix" parameter input\nsimulations with label "{0}" are not being performed\nsimulations being performed: {1}'.format(sim_name,sorted(simreps.keys())))
            #end if
            # walk up the coarse (segmented) DAG and search for constraint variables
            #  upon return, an index for each matching constraint variable is present
            simreps[sim_name].secrep.segrep.check_constraint(constraint,fix_tol,loc=loc)
            c = constraint
            for v,s,p in zip(c.variables,c.selectors,c.pattern):
                if v not in c.check.variables_found:
                    error('problem with "fix" parameter input\nconstraint variable not found in scans\nerroneous constraint variable: {0}\nthere is likely a typo in this user-provided variable name'.format(v),loc=loc)
                #end if
                if v not in c.check.matches_found:
                    error('match not found in scans for constraint variable\nconstraint variable to match: {0}\nvalue to match: {1}'.format(v,s),loc=loc)
                #end if
            #end if
        #end for
    #end if

    #graph_sims(simreps,quants=False,exit=False)
    #graph_sims(secreps,quants=False,exit=False)
    #graph_sims(segreps,quants=False,exit=False)

    # generate the (multi) scanned workflow recursively (heavy)
    sims,sim_coll = segreps.invariant.generate_workflow()


    # impose constraints by merging specified simulations
    if not missing(fix):
        for sim_name,constraint in constraints.iteritems():
            if sim_name not in sim_coll or len(sim_coll[sim_name])==0:
                error('cannot impose constraint on simulation labeled {0}\nthis is likely a developer error'.format(sim_name))
            #end if
            scoll = sim_coll[sim_name]
            if len(scoll)>1:
                fix_vars = sorted(set(constraint.variables))
                fix_inds = constraint.indices.tuple(fix_vars)
                fixed_sims   = obj()
                fixed_simids = set()
                grouped_sims = obj()
                inds0 = scoll[0][1]
                other_vars = sorted(set(inds0.keys())-set(fix_vars))
                for sim,inds,vals in scoll:
                    ifix   = inds.tuple(fix_vars)
                    iother = inds.tuple(other_vars)
                    if ifix==fix_inds:
                        fixed_sims[iother] = sim
                    else:
                        if iother not in grouped_sims:
                            grouped_sims[iother] = []
                        #end if
                        grouped_sims[iother].append(sim)
                    #end if
                #end for
                for iother,fixed_sim in fixed_sims.iteritems():
                    if iother in grouped_sims:
                        for sim in grouped_sims[iother]:
                            fixed_sim.acquire_dependents(sim)
                        #end for
                    #end if
                #end for
            #end if
        #end for
    #end if

    # bundle requested simulation jobs together
    if not missing(bundle) and len(bundle)>0:
        bsims = SimSet()
        sims.bundle = bsims
        for sim_name in bundle:
            if sim_name not in sim_coll:
                error('problem with "bundle" parameter input\nsimulations with requested name are not present\nsimulation name requested for bundling: {0}\nsimulation names present: {1}'.format(sim_name,sorted(sim_coll.keys())),loc=loc)
            #end if
            bsim_list = []
            for scoll in sim_coll[sim_name]:
                sim = scoll[0]
                if not sim.fake():
                    bsim_list.append(sim)
                #end if
            #end for
            bsim = bundle_function(bsim_list)
            loc_sims.append(bsim)
            bsims[sim_name] = bsim
        #end for
    #end if

    # get rid of any fake/temporary simulations
    sims.remove_fake()


    if graph_workflow:
        graph_sims(loc_sims,exit=False)
    #end if
    if write_workflow:
        write_sims(loc_sims,exit=False)
    #end if
            
    if graph_simreps:
        graph_sims(simreps,quants=False,exit=False)
    #end if
    if write_simreps:
        write_sims(simreps,quants=False,exit=False)
    #end if
    if graph_secreps:
        graph_sims(secreps,quants=False,exit=False)
    #end if
    if write_secreps:
        write_sims(secreps,quants=False,exit=False)
    #end if
    if graph_segreps:
        graph_sims(segreps,quants=False,exit=False)
    #end if
    if write_segreps:
        write_sims(segreps,quants=False,exit=False)
    #end if


    sim_list.extend(loc_sims)

    #tpoint('finish workflows')


    #ci()
    #exit()

    return sims
#end def qmcpack_workflow



# deprecated
ecut_scan_chain_defaults = obj(
    scf = True,
    p2q = True,
    opt = True,
    vmc = True,
    )
def ecut_scan(
    ecuts        = missing,
    basepath     = missing,
    dirname      = 'ecut_scan',
    same_jastrow = True,
    ecut_jastrow = None,
    **kwargs
    ):

    loc = 'ecut_scan'

    set_def_loc(obj(),loc)

    require('ecuts'   ,ecuts   )
    require('basepath',basepath)

    ecuts = list(ecuts)

    qckw = process_qmcpack_chain_kwargs(
        defaults = ecut_scan_chain_defaults,
        loc      = loc,
        **kwargs
        )

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



# deprecated
def system_scan(
    basepath     = missing,
    dirname      = 'system_scan',
    systems      = missing,
    sysdirs      = missing,
    syskeys      = None,
    same_jastrow = False, # deprecated
    jastrow_key  = missing,
    J2_source    = None,
    J3_source    = None,
    loc          = 'system_scan',
    **kwargs
    ):
    set_loc(loc)

    require('basepath',basepath)
    require('systems' ,systems )
    require('sysdirs' ,sysdirs )

    if syskeys is None:
        syskeys = sysdirs
    #end if

    if len(systems)==0:
        error('no systems provided',loc)
    #end if
    if len(sysdirs)!=len(systems):
        error('must provide one directory per system via sysdirs keyword\nnumber of directories provided: {0}\nnumber of systems: {1}'.format(len(sysdirs),len(systems)),loc)
    #end if
    if len(syskeys)!=len(systems):
        error('must provide one key per system via syskeys keyword\nnumber of keys provided: {0}\nnumber of systems: {1}'.format(len(syskeys),len(systems)),loc)
    #end if

    same_jastrow = not missing(jastrow_key)
    if same_jastrow:
        if jastrow_key not in set(syskeys):
            error('key used to identify jastrow for use across scan was not found\njastrow key provided: {0}\nsystem keys present: {1}'.format(jastrow_key,sorted(system_keys)),loc)
        #end if
        sys = obj()
        for n in xrange(len(systems)):
            sys[syskeys[n]] = systems[n],sysdirs[n]
        #end for
        key_system,key_sysdir = sys[jastrow_key]
        del sys[jastrow_key]
        systems = [key_system]
        sysdirs = [key_sysdir]
        syskeys = [jastrow_key]
        for syskey in sys.keys():
            system,sysdir = sys[syskey]
            systems.append(system)
            sysdirs.append(sysdir)
            syskeys.append(syskey)
        #end for
    #end if

    sims = obj()
    for n in xrange(len(systems)):
        qckw = process_qmcpack_chain_kwargs(
            defaults = qmcpack_chain_defaults,
            system   = systems[n],
            loc      = loc,
            **kwargs
            )
        qckw.basepath  = os.path.join(basepath,dirname,sysdirs[n])
        qckw.J2_source = J2_source
        qckw.J3_source = J3_source
        qcsims = qmcpack_chain(**qckw)
        if same_jastrow and n==0:
            J2_source = qcsims.get_optional('optJ2',None)
            J3_source = qcsims.get_optional('optJ3',None)
        #end if
        sims[syskeys[n]] = qcsims
    #end for
    return sims
#end def system_scan


# deprecated
def system_parameter_scan(
    basepath  = missing,
    dirname   = 'system_param_scan',
    sysfunc   = missing,
    parameter = missing,
    variable  = missing, # same as parameter, to be deprecated
    values    = missing,
    fixed     = None,
    loc       = 'system_parameter_scan',
    **kwargs
    ):

    set_loc(loc)

    if missing(parameter):
        parameter = variable
    #end if

    require('basepath' ,basepath )
    require('sysfunc'  ,sysfunc  )
    require('parameter',parameter)
    require('values'   ,values   )

    systems = []
    sysdirs = []
    syskeys = []
    for v in values:
        params = obj()
        params[parameter] = v
        if fixed!=None:
            params.set(**fixed)
        #end if
        system = sysfunc(**params)
        vkey   = render_parameter_key(v,loc)
        vlabel = render_parameter_label(vkey,loc)
        sysdir = '{0}_{1}'.format(parameter,vlabel)
        systems.append(system)
        sysdirs.append(sysdir)
        syskeys.append(vkey)
    #end for

    sp_sims = system_scan(
        basepath = basepath,
        dirname  = dirname,
        systems  = systems,
        sysdirs  = sysdirs,
        syskeys  = syskeys,
        loc      = loc,
        **kwargs
        )

    return sp_sims
#end def system_parameter_scan




# deprecated
def input_parameter_scan(
    basepath     = missing,
    dirname      = 'input_param_scan',
    section      = missing,
    parameter    = missing,
    parameterset = missing,
    variable     = missing, # same as parameter, to be deprecated
    values       = missing,
    tags         = missing,
    jastrow_key  = missing,
    J2_source    = None,
    J3_source    = None,
    loc          = 'input_parameter_scan',
    **kwargs
    ):

    set_loc(loc)

    if missing(parameter):
        parameter = variable
    #end if

    require('basepath' ,basepath )
    require('section'  ,section  )
    if missing(parameterset):
        require('parameter',parameter)
        require('values'   ,values   ) 
        if not missing(tags) and len(tags)!=len(values):
            error('must provide one tag (directory label) per parameter value\nnumber of values: {0}\nnumber of tags: {1}\nvalues provided: {2}\ntags provided: {3}'.format(len(values),len(tags),values,tags),loc)
        #end if
    else:
        require('tags',tags)
        for pname,pvalues in parameterset.iteritems():
            if len(tags)!=len(pvalues):
                error('must provide one tag (directory label) per parameter value\nparameter name: {0}\nnumber of values: {1}\nnumber of tags: {2}\nvalues provided: {3}\ntags provided: {4}'.format(pname,len(pvalues),len(tags),pvalues,tags),loc)
            #end if
        #end for
    #end if


    paramdirs = []
    paramkeys = []
    if missing(parameterset):
        n = 0
        for v in values:
            if not missing(tags):
                vkey = tags[n]
            else:
                vkey = render_parameter_key(v,loc)
            #end if
            vlabel = render_parameter_label(vkey,loc)
            paramdir = '{0}_{1}'.format(parameter,vlabel)
            paramdirs.append(paramdir)
            paramkeys.append(vkey)
            n+=1
        #end for
    else:
        paramdirs.extend(tags)
        paramkeys.extend(tags)
    #end if

    same_jastrow = not missing(jastrow_key)
    if same_jastrow:
        jastrow_key = render_parameter_key(jastrow_key,loc)
        if not jastrow_key in set(paramkeys):
            error('input parameter value for fixed Jastrow not found\njastrow key provided: {0}\nparameter keys present: {1}'.format(jastrow_key,paramkeys),loc)
        #end if
        index = paramkeys.index(jastrow_key)
        paramdirs.insert(0,paramdirs.pop(index))
        paramkeys.insert(0,paramkeys.pop(index))
        if missing(parameterset):
            values = list(values)
            values.insert(0,values.pop(index))
        else:
            for pname,pvalues in parameterset.iteritems():
                pvalues = list(pvalues)
                pvalues.insert(0,pvalues.pop(index))
                parameterset[pname] = pvalues
            #end for
        #end if
    #end if

    ip_sims = obj()
    for n in xrange(len(paramkeys)):
        qckw = process_qmcpack_chain_kwargs(
            defaults = qmcpack_chain_defaults,
            loc      = loc,
            **kwargs
            )
        # just do the simplest thing for now: duplicate full workflows
        # later should branch workflow only downstream from varying sim
        if 'inputs' not in section:
            error('section must be a name of the form *_inputs\nsection provided: {0}\ninput sections present: {1}'.format(section,[s for s in sorted(qckw.keys()) if 'inputs' in s]),loc)
        elif section not in qckw:
            error('input section not found\ninput section provided: {0}\ninput sections present: {1}'.format(section,[s for s in sorted(qckw.keys()) if 'inputs' in s]),loc)
        #end if
        if missing(parameterset):
            qckw[section][parameter] = values[n]
        else:
            sec = qckw[section]
            for pname,pvalues in parameterset.iteritems():
                pvalue = pvalues[n]
                if not novalue(pvalue):
                    sec[pname] = pvalue
                #end if
            #end for
        #end if

        qckw.basepath  = os.path.join(basepath,dirname,paramdirs[n])
        qckw.J2_source = J2_source
        qckw.J3_source = J3_source
        qcsims = qmcpack_chain(**qckw)
        if same_jastrow and n==0:
            J2_source = qcsims.get_optional('optJ2',None)
            J3_source = qcsims.get_optional('optJ3',None)
        #end if
        ip_sims[paramkeys[n]] = qcsims
    #end for

    return ip_sims
#end def input_parameter_scan




# deprecated
def input_parameter_scan2(
    basepath     = missing,
    dirname      = 'input_param_scan',
    section      = missing,
    parameter    = missing,
    parameterset = missing,
    variable     = missing, # same as parameter, to be deprecated
    values       = missing,
    tags         = missing,
    jastrow_key  = missing,
    J2_source    = None,
    J3_source    = None,
    loc          = 'input_parameter_scan',
    **kwargs
    ):

    set_loc(loc)

    if missing(parameter):
        parameter = variable
    #end if

    require('basepath' ,basepath )
    require('section'  ,section  )
    if missing(parameterset):
        require('parameter',parameter)
        require('values'   ,values   ) 
        if not missing(tags) and len(tags)!=len(values):
            error('must provide one tag (directory label) per parameter value\nnumber of values: {0}\nnumber of tags: {1}\nvalues provided: {2}\ntags provided: {3}'.format(len(values),len(tags),values,tags),loc)
        #end if
    else:
        require('tags',tags)
        for pname,pvalues in parameterset.iteritems():
            if len(tags)!=len(pvalues):
                error('must provide one tag (directory label) per parameter value\nparameter name: {0}\nnumber of values: {1}\nnumber of tags: {2}\nvalues provided: {3}\ntags provided: {4}'.format(pname,len(pvalues),len(tags),pvalues,tags),loc)
            #end if
        #end for
    #end if


    paramdirs = []
    paramkeys = []
    if missing(parameterset):
        n = 0
        for v in values:
            if not missing(tags):
                vkey = tags[n]
            else:
                vkey = render_parameter_key(v,loc)
            #end if
            vlabel = render_parameter_label(vkey,loc)
            paramdir = '{0}_{1}'.format(parameter,vlabel)
            paramdirs.append(paramdir)
            paramkeys.append(vkey)
            n+=1
        #end for
        parameterset = obj()
        parameterset[parameter] = values
    else:
        paramdirs.extend(tags)
        paramkeys.extend(tags)
    #end if

    same_jastrow = not missing(jastrow_key)
    if same_jastrow:
        jastrow_key = render_key(jastrow_key,loc)
        if not jastrow_key in set(paramkeys):
            error('input parameter value for fixed Jastrow not found\njastrow key provided: {0}\nparameter keys present: {1}'.format(jastrow_key,paramkeys),loc)
        #end if
        index = paramkeys.index(jastrow_key)
        paramdirs.insert(0,paramdirs.pop(index))
        paramkeys.insert(0,paramkeys.pop(index))
        for pname,pvalues in parameterset.iteritems():
            pvalues = list(pvalues)
            pvalues.insert(0,pvalues.pop(index))
            parameterset[pname] = pvalues
        #end for
    #end if

    varying = obj()
    for name in qmcpack_chain_sim_names:
        varying[name] = section==name+'_inputs'
    #end for
    varying.nscf |= varying.scf
    varying.p2q  |= varying.nscf
    varying.opt  |= varying.p2q
    varying.vmc  |= varying.opt
    varying.dmc  |= varying.opt

    suppress_J0 = False
    suppress_J2 = False
    for pname in parameterset.keys():
        if pname.startswith('J3'):
            suppress_J0 = True
            suppress_J2 = True
        elif pname.startswith('J2'):
            suppress_J0 = True
        #end if
    #end for

    ip_sims = obj()
    for n in xrange(len(paramkeys)):

        kwcap = capture_qmcpack_chain_inputs(**kwargs)

        if 'inputs' not in section:
            error('section must be a name of the form *_inputs\nsection provided: {0}\ninput sections present: {1}'.format(section,[s for s in sorted(kwcap.keys()) if 'inputs' in s]),loc)
        elif section not in kwcap:
            error('input section not found\ninput section provided: {0}\ninput sections present: {1}'.format(section,[s for s in sorted(kwcap.keys()) if 'inputs' in s]),loc)
        #end if
        sec = kwcap[section]
        for pname,pvalues in parameterset.iteritems():
            pvalue = pvalues[n]
            if not novalue(pvalue):
                sec[pname] = pvalue
            #end if
        #end for

        qckw = process_qmcpack_chain_kwargs(
            defaults = qmcpack_chain_defaults,
            loc      = loc,
            **kwcap
            )

        if n>0:
            for name,vary in varying.iteritems():
                if not vary:
                    qckw[name]           = False
                    qckw[name+'_inputs'] = None
                #end if
            #end for
            if 'vmc_inputs' in qckw:
                if suppress_J0:
                    qckw.vmc_inputs.J0_run = 0
                    qckw.vmc_inputs.J0_test = 0
                #end if
                if suppress_J2:
                    qckw.vmc_inputs.J2_run = 0
                    qckw.vmc_inputs.J2_test = 0
                #end if
            #end if
            if 'opt_inputs' in qckw:
                if suppress_J2:
                    qckw.opt_inputs.J2_run = 0
                    qckw.opt_inputs.J2_test = 0
                #end if
            #end if
            if 'dmc_inputs' in qckw:
                if suppress_J0:
                    qckw.dmc_inputs.J0_run = 0
                    qckw.dmc_inputs.J0_test = 0
                #end if
                if suppress_J2:
                    qckw.dmc_inputs.J2_run = 0
                    qckw.dmc_inputs.J2_test = 0
                #end if
            #end if
        #end if

        qckw.basepath  = os.path.join(basepath,dirname,paramdirs[n])
        qckw.J2_source = J2_source
        qckw.J3_source = J3_source
        qcsims = qmcpack_chain(**qckw)
        if n==0:
            if not varying.p2q:
                orb_source = qcsims.get_optional('p2q',None)
            elif not varying.opt:
                J2_source = qcsims.get_optional('optJ2',None)
                J3_source = qcsims.get_optional('optJ3',None)
            elif varying.opt and suppress_J2:
                J3_source = qcsims.get_optional('optJ3',None)
            elif same_jastrow:
                J2_source = qcsims.get_optional('optJ2',None)
                J3_source = qcsims.get_optional('optJ3',None)
            #end if
        #end if
        ip_sims[paramkeys[n]] = qcsims
    #end for

    return ip_sims
#end def input_parameter_scan2




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
