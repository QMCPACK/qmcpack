##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  qmcpack_converters.py                                             #
#    Nexus interfaces for orbital converter tools of QMCPACK:        #
#    pw2qmcpack and convert4qmc                                      #
#                                                                    #
#  Content summary:                                                  #
#    generate_pw2qmcpack                                             #
#      User-facing function to create pw2qmcpack simulation objects. #
#                                                                    #
#    generate_pw2qmcpack_input                                       #
#      User-facing funcion to create input for pw2qmcpack.           #
#                                                                    #
#    Pw2qmcpack                                                      #
#      Simulation class for pw2qmcpack.                              #
#                                                                    #
#    Pw2qmcpackInput                                                 #
#      SimulationInput class for pw2qmcpack.                         #
#                                                                    #
#    Pw2qmcpackAnalyzer                                              #
#      SimulationAnalyzer class for pw2qmcpack.                      #
#                                                                    #
#                                                                    #
#    Convert4qmcInput                                                #
#      Class representing command line interface of convert4qmc.     #
#                                                                    #
#    Convert4qmcAnalyzer                                             #
#      Placeholder class for output analysis.                        #
#                                                                    #
#    Convert4qmc                                                     #
#      Class representing convert4qmc instance.                      #
#                                                                    #
#    generate_convert4qmc_input                                      #
#      Function to generate arbitrary convert4qmc input.             #
#                                                                    #
#    generate_convert4qmc                                            #
#      Function to generate Convert4qmc simulation object.           #
#                                                                    #
#====================================================================#



import os
from generic import obj
from simulation import Simulation,SimulationInput,SimulationAnalyzer
from pwscf import Pwscf
from gamess import Gamess
from pyscf_sim import Pyscf
from quantum_package import QuantumPackage


# read/write functions associated with pw2qmcpack only
def read_str(sv):
    return sv.strip('"').strip("'")
#end def read_str

def read_int(sv):
    return int(sv)
#end def read_int

def read_float(sv):
    return float(sv.replace('d','e').replace('D','e'))
#end def read_float

bconv = {'.true.':True,'.false.':False}
def read_bool(sv):
    return bconv[sv]
#end def read_bool

def write_str(val):
    return "'"+val+"'"
#end def write_str

def write_int(val):
    return str(val)
#end def write_int

def write_float(val):
    return str(val)
#end def write_float

def write_bool(val):
    return '.'+str(val).lower()+'.'
#end def write_bool

readval={str:read_str,int:read_int,float:read_float,bool:read_bool}
writeval={str:write_str,int:write_int,float:write_float,bool:write_bool}


class Pw2qmcpackInput(SimulationInput):
    ints   = []
    floats = []
    strs   = ['outdir','prefix']
    bools  = ['write_psir']

    var_types = dict()
    for v in ints:
        var_types[v]=int
    #end for
    for v in floats:
        var_types[v]=float
    #end for
    for v in strs:
        var_types[v]=str
    #end for
    for v in bools:
        var_types[v]=bool
    #end for

    allowed = set(ints+floats+strs+bools)

    def read_text(self,contents,filepath=None):
        lines = contents.split('\n')
        inside = False
        for l in lines:
            if inside:
                tokens = l.split(',')
                for t in tokens:
                    ts = t.strip()
                    if ts!='' and ts!='/':
                        name,value = ts.split('=')
                        name = name.strip()
                        value= value.strip()
                        if name in self.allowed:
                            vtype = self.var_types[name]
                            value = readval[vtype](value)
                            sobj[name]=value
                        else:
                            self.error('encountered unknown variable: '+name)
                        #end if
                    #end if
                #end for
            #end if
            if '&' in l:
                inside=True
                section = l.lstrip('&').lower()
                sobj = obj()
                self[section]=sobj
            elif l.strip()=='/':
                inside=False
            #end if
        #end for
    #end def read_text


    def write_text(self,filepath=None):
        contents = ''
        for sname,section in self.items():
            contents+='&'+sname+'\n'
            for name,value in section.items():
                vtype = type(value)
                contents += '  '+name+' = '+writeval[vtype](value)+'\n'
            #end for
            contents+='/\n'
        #end for
        return contents
    #end def write_text


    def __init__(self,filepath=None,**vars):
        if filepath!=None:
            self.read(filepath)
        else:
            inputpp = obj()
            for name,value in vars.items():
                inputpp[name] = value
            #end for
            self.inputpp = inputpp
        #end if
    #end def __init__
#end class Pw2qmcpackInput


def generate_pw2qmcpack_input(prefix='pwscf',outdir='pwscf_output',write_psir=False):
    pw = Pw2qmcpackInput(
        prefix     = prefix,
        outdir     = outdir,
        write_psir = write_psir
        )
    return pw
#end def generate_pw2qmcpack_input


def read_eshdf_eig_data(filename, Ef_list):
    import numpy as np
    from numpy import array,pi
    from numpy.linalg import inv
    from unit_converter import convert
    from hdfreader import read_hdf
    from developer import error

    def h5int(i):
        return array(i,dtype=int)[0]
    #end def h5int

    h        = read_hdf(filename,view=True)
    axes     = array(h.supercell.primitive_vectors)
    kaxes    = 2*pi*inv(axes).T
    nk       = h5int(h.electrons.number_of_kpoints)
    ns       = h5int(h.electrons.number_of_spins)
    if (len(Ef_list) == 1 and ns == 2):
        # Using the same E_fermi for up and down electrons
        E_fermi = Ef_list[0]
        Ef_list = np.array([E_fermi, E_fermi])
    elif len(Ef_list) != ns:
        msg = 'Ef "%s" must have same length as nspin=%d' % (str(Ef_list), ns)
        error(msg)
    data     = obj()
    for k in range(nk):
        kp = h.electrons['kpoint_'+str(k)]
        for s, Ef in zip(range(ns), Ef_list):
            E_fermi = Ef+1e-8
            eig_s = []
            path = 'electrons/kpoint_{0}/spin_{1}'.format(k,s)
            spin = h.get_path(path)
            eig = convert(array(spin.eigenvalues),'Ha','eV')
            nst = h5int(spin.number_of_states)
            for st in range(nst):
                e = eig[st]
                if e<E_fermi:
                    eig_s.append(e)
                #end if
            #end for
            data[k,s] = obj(
                kpoint = array(kp.reduced_k),
                eig    = array(eig_s),
                )
        #end for
    #end for
    res = obj(
        orbfile  = filename,
        axes     = axes,
        kaxes    = kaxes,
        nkpoints = nk,
        nspins   = ns,
        data     = data,
        )
    return res
#end def read_eshdf_eig_data


def gcta_occupation(wfh5, ntwist):
  nspin = wfh5.nspins
  nk = wfh5.nkpoints
  nprim = nk//ntwist
  assert nprim*ntwist == nk
  nelecs_at_twist = []
  for itwist in range(ntwist):
    iks = range(itwist, nk, ntwist)
    # calculate nelec for each spin
    nelecs = []
    for ispin in range(nspin):
      nl = [len(wfh5.data[ik, ispin].eig) for ik in iks]
      nup = sum(nl)
      nelecs.append(nup)
      if nspin == 1:
        nelecs.append(nup)
    nelecs_at_twist.append(nelecs)
  return nelecs_at_twist
#end gcta_occupation


class Pw2qmcpackAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0):
        if isinstance(arg0,Simulation):
            sim = arg0
            self.infile = sim.infile
            prefix,outdir = sim.input.inputpp.tuple('prefix','outdir')
            self.dir = sim.locdir
            self.h5file = os.path.join(sim.locdir,outdir,prefix+'.pwscf.h5')
        else:
            self.infile = arg0
        #end if
    #end def __init__

    def analyze(self, Ef_list=None):
      if Ef_list is not None:
        self.wfh5 = read_eshdf_eig_data(self.h5file, Ef_list)
      #end if
    #end def analyze

    def get_result(self,result_name):
        self.not_implemented()
    #end def get_result
#end class Pw2qmcpackAnalyzer


class Pw2qmcpack(Simulation):
    input_type = Pw2qmcpackInput
    analyzer_type = Pw2qmcpackAnalyzer
    generic_identifier = 'pw2qmcpack'
    application = 'pw2qmcpack.x'
    application_properties = set(['serial'])
    application_results    = set(['orbitals','gc_occupation'])

    def check_result(self,result_name,sim):
        calculating_result = False
        if result_name=='orbitals':
            calculating_result = True
        elif result_name=='gc_occupation':
            calculating_result = True
        #end if        
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj()
        inputpp = self.input.inputpp
        prefix = 'pwscf'
        outdir = './'
        if 'prefix' in inputpp:
            prefix = inputpp.prefix
        #end if
        if 'outdir' in inputpp:
            outdir = inputpp.outdir
        #end if
        if outdir.startswith('./'):
            outdir = outdir[2:]
        #end if
        if result_name=='orbitals':
            result.h5file   = os.path.join(self.locdir,outdir,prefix+'.pwscf.h5')
            result.ptcl_xml = os.path.join(self.locdir,outdir,prefix+'.ptcl.xml')
            result.wfs_xml  = os.path.join(self.locdir,outdir,prefix+'.wfs.xml')
        elif result_name=='gc_occupation':
            pass  # defer to Qmcpack.incorporate_result
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if        
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        implemented = True
        if result_name=='orbitals':
            if isinstance(sim,Pwscf):
                pwin = sim.input.control
                p2in = self.input.inputpp
                pwprefix = 'pwscf'
                p2prefix = 'pwscf'
                pwoutdir = './'
                p2outdir = './'
                if 'prefix' in pwin:
                    pwprefix = pwin.prefix
                #end if
                if 'prefix' in p2in:
                    p2prefix = p2in.prefix
                #end if
                if 'outdir' in pwin:
                    pwoutdir = pwin.outdir
                #end if
                if 'outdir' in p2in:
                    p2outdir = p2in.outdir
                #end if
                if pwoutdir.startswith('./'):
                    pwoutdir = pwoutdir[2:]
                #end if
                if p2outdir.startswith('./'):
                    p2outdir = p2outdir[2:]
                #end if
                pwdir = os.path.abspath(os.path.join(sim.locdir ,pwoutdir))
                p2dir = os.path.abspath(os.path.join(self.locdir,p2outdir))
                errors = False
                if pwdir!=p2dir:
                    self.error('to use orbitals, '+self.generic_identifier+' must have the same outdir as pwscf\n  pwscf outdir: '+pwdir+'\n  '+self.generic_identifier+' outdir: '+p2dir,exit=False)
                    errors = True
                #end if
                if pwprefix!=p2prefix:
                    self.error('to use orbitals, '+self.generic_identifier+' must have the same prefix as pwscf\n  pwscf prefix: '+pwprefix+'\n  '+self.generic_identifier+' prefix: '+p2prefix,exit=False)
                    errors = True
                #end if
                if errors:
                    self.error(self.generic_identifier+' cannot use orbitals from pwscf')
                #end if
            else:
                implemented = False
            #end if
        else:
            implemented = False
        #end if
        if not implemented:
            self.error('ability to incorporate result "{0}" from {1} has not been implemented'.format(result_name,sim.__class__.__name__))
        #end if                
    #end def incorporate_result


    def check_sim_status(self):
        outfile = os.path.join(self.locdir,self.outfile)
        fobj = open(outfile,'r')
        output = fobj.read()
        fobj.close()
        inputpp = self.input.inputpp
        prefix = 'pwscf'
        outdir = './'
        if 'prefix' in inputpp:
            prefix = inputpp.prefix
        #end if
        if 'outdir' in inputpp:
            outdir = inputpp.outdir
        #end if
        if outdir.startswith('./'):
            outdir = outdir[2:]
        #end if
        h5file   = os.path.join(self.locdir,outdir,prefix+'.pwscf.h5')
        ptcl_xml = os.path.join(self.locdir,outdir,prefix+'.ptcl.xml')
        wfs_xml  = os.path.join(self.locdir,outdir,prefix+'.wfs.xml')
        must_exist = [h5file,ptcl_xml,wfs_xml]

        files_exist = True
        for file in must_exist:
            files_exist = files_exist and os.path.exists(file)
        #end for
        outfin = True
        #outfin = outfin and 'esh5 create' in output
        #outfin = outfin and 'Creating electrons' in output
        outfin = outfin and 'npw=' in output
        outfin = outfin and 'ik=' in output

        outfin = outfin or 'JOB DONE' in output

        success = files_exist and outfin

        #self.finished = success and self.job.finished

        # pw2qmcpack has too many variants to assess completion based on log output
        #   assume (optimistically) that job completion indicates success
        self.finished = files_exist and self.job.finished
    #end def check_sim_status


    def get_output_files(self):
        output_files = []
        return output_files
    #end def get_output_files


    def app_command(self):
        return self.app_name+'<'+self.infile
    #end def app_command
#end class Pw2qmcpack




def generate_pw2qmcpack(**kwargs):
    sim_args,inp_args = Simulation.separate_inputs(kwargs)

    if not 'input' in sim_args:
        sim_args.input = generate_pw2qmcpack_input(**inp_args)
    #end if
    pw2qmcpack = Pw2qmcpack(**sim_args)

    return pw2qmcpack
#end def generate_pw2qmcpack










class Convert4qmcInput(SimulationInput):

    input_codes = '''
        pyscf              
        qp                 
        gaussian           
        casino             
        vsvb               
        gamess             
        gamess_ascii       
        gamess_fmo         
        gamess_xml         
        '''.split()

    input_order = input_codes + '''
        prefix             
        hdf5               
        add_cusp           
        psi_tag            
        ion_tag            
        no_jastrow         
        production         
        orbitals
        multidet
        gridtype
        first
        last
        size
        ci                 
        read_initial_guess 
        target_state       
        natural_orbitals   
        threshold          
        opt_det_coeffs
        zero_ci            
        add_3body_J
        '''.split()

    input_aliases = obj(
        pyscf              = 'pyscf',
        qp                 = 'QP',
        gaussian           = 'gaussian',
        casino             = 'casino',
        vsvb               = 'VSVB',
        gamess             = 'gamess',
        gamess_ascii       = 'gamess',
        gamess_fmo         = 'gamessFMO',
        gamess_xml         = 'gamesxml', # not a typo
        prefix             = 'prefix',
        hdf5               = 'hdf5',
        add_cusp           = 'addCusp',
        psi_tag            = 'psi_tag',
        ion_tag            = 'ion_tag',
        no_jastrow         = 'nojastrow',
        production         = 'production',
        orbitals           = 'orbitals',
        multidet           = 'multidet',
        gridtype           = 'gridtype',
        first              = 'first',
        last               = 'last',
        size               = 'size',
        ci                 = 'ci',
        read_initial_guess = 'readInitialGuess',
        target_state       = 'TargetState',
        natural_orbitals   = 'NaturalOrbitals',
        threshold          = 'threshold',
        opt_det_coeffs     = 'optDetCoeffs',
        zero_ci            = 'zeroCi',
        add_3body_J        = 'add3BodyJ',
        )

    input_types = obj(
        app_name           = str, # executable name
        pyscf              = str, # file path
        qp                 = str, # file path
        gaussian           = str, # file path
        casino             = str, # file path
        vsvb               = str, # file path
        gamess             = str, # file path
        gamess_ascii       = str, # file path
        gamess_fmo         = str, # file path
        gamess_xml         = str, # file path
        prefix             = str, # any name
        hdf5               = bool,
        add_cusp           = bool,
        psi_tag            = str, # wavefunction tag
        ion_tag            = str, # particeset tag
        no_jastrow         = bool,
        production         = bool,
        orbitals           = str,
        multidet           = str,
        gridtype           = str,
        first              = float,
        last               = float,
        size               = int,
        ci                 = str, # file path
        read_initial_guess = int,
        target_state       = int,
        natural_orbitals   = int,
        threshold          = float,
        opt_det_coeffs     = bool,
        zero_ci            = bool,
        add_3body_J        = bool,
        )

    input_defaults = obj(
        app_name           = 'convert4qmc',
        pyscf              = None, # input codes
        qp                 = None,
        gaussian           = None,
        casino             = None,
        vsvb               = None,
        gamess             = None, 
        gamess_ascii       = None,
        gamess_fmo         = None,
        gamess_xml         = None,
        prefix             = None, # general options
        hdf5               = False,
        add_cusp           = False,
        psi_tag            = None,
        ion_tag            = None,
        no_jastrow         = False,
        production         = False,
        orbitals           = None,
        multidet           = None,
        gridtype           = None,
        first              = None,
        last               = None,
        size               = None,
        ci                 = None, # gamess specific below
        read_initial_guess = None,
        target_state       = None,
        natural_orbitals   = None,
        threshold          = None,
        opt_det_coeffs     = False,
        zero_ci            = False,
        add_3body_J        = False,# deprecated
        )


    def __init__(self,**kwargs):
        # check that only allowed keyword inputs are provided
        invalid = set(kwargs.keys())-set(self.input_types.keys())
        if len(invalid)>0:
            self.error('invalid inputs encountered\ninvalid keywords: {0}\nvalid keyword inputs are: {1}'.format(sorted(invalid),sorted(self.input_types.keys())))
        #end if

        # assign inputs
        self.set(**kwargs)

        # assign default values
        self.set_optional(**self.input_defaults)

        # check that all keyword inputs are valid
        self.check_valid()
    #end def __init__


    def check_valid(self,exit=True):
        valid = True
        # check that all inputs have valid types and assign them
        for k,v in self.items():
            if v is not None and not isinstance(v,self.input_types[k]):
                valid = False
                if exit:
                    self.error('keyword input {0} must be of type {1}\nyou provided a value of type {2}\nplease revise your input and try again'.format(k,self.input_types[k].__name__),v.__class__.__name__)
                #end if
                break
            #end if
        #end for
        return valid
    #end def check_valid


    def set_app_name(self,app_name):
        self.app_name = app_name
    #end def set_app_name


    def input_code(self):
        input_code = None
        for k in self.input_codes:
            if k in self and self[k] is not None:
                if input_code is not None:
                    input_code = None
                    break
                else:
                    input_code = self[k]
                #end if
            #end if
        #end for
        return input_code
    #end def input_code


    def has_input_code(self):
        return self.input_code() is not None
    #end def has_input_code


    def app_command(self):
        self.check_valid()
        c = self.app_name
        for k in self.input_order:
            if k in self:
                v = self[k]
                n = self.input_aliases[k]
                if isinstance(v,bool):
                    if v:
                        c += ' -{0}'.format(n)
                    #end if
                elif v is not None:
                    c += ' -{0} {1}'.format(n,str(v))
                #end if
            #end if
        #end for
        return c
    #end def app_command


    def read(self,filepath):
        None
    #end def read


    def write_text(self,filepath=None):
        return self.app_command()
    #end def write_text


    def output_files(self):
        prefix = 'sample'
        if self.prefix!=None:
            prefix = self.prefix
        #end if
        wfn_file  = prefix+'.Gaussian-G2.xml'
        ptcl_file = prefix+'.Gaussian-G2.ptcl.xml'
        return wfn_file,ptcl_file
    #end def output_files
#end class Convert4qmcInput



def generate_convert4qmc_input(**kwargs):
    return Convert4qmcInput(**kwargs)
#end def generate_convert4qmc_input



class Convert4qmcAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0):
        if isinstance(arg0,Simulation):
            self.infile = arg0.infile
        else:
            self.infile = arg0
        #end if
    #end def __init__

    def analyze(self):
        None
    #end def analyze
#end class Convert4qmcAnalyzer



class Convert4qmc(Simulation):
    input_type             = Convert4qmcInput
    analyzer_type          = Convert4qmcAnalyzer
    generic_identifier     = 'convert4qmc'
    application            = 'convert4qmc'
    application_properties = set(['serial'])
    application_results    = set(['orbitals','particles'])
    renew_app_command      = True

    def __init__(self,*args,**kwargs):
        Simulation.__init__(self,*args,**kwargs)
        self.input_code = None
    #end def __init__


    def set_app_name(self,app_name):
        self.app_name = app_name
        self.input.set_app_name(app_name)
    #end def set_app_name


    def propagate_identifier(self):
        None
        #self.input.prefix = self.identifier
    #end def propagate_identifier


    def get_prefix(self):
        input = self.input
        prefix = 'sample'
        if input.prefix is not None:
            prefix = input.prefix
        #end if
        return prefix
    #end def get_prefix


    def list_output_files(self):
        # try to support both pre and post v3.3.0 convert4qmc
        prefix = self.get_prefix()
        wfn_file  = prefix+'.Gaussian-G2.xml'
        ptcl_file = prefix+'.Gaussian-G2.ptcl.xml'
        if not os.path.exists(os.path.join(self.locdir,ptcl_file)):
            if self.input.no_jastrow:
                wfn_file  = prefix+'.wfnoj.xml'
            else:
                wfn_file  = prefix+'.wfj.xml'
            #end if
            ptcl_file = prefix+'.structure.xml'
        #end if
        return wfn_file,ptcl_file
    #end def list_output_files


    def check_result(self,result_name,sim):
        calculating_result = False
        if result_name=='orbitals':
            calculating_result = True
        elif result_name=='particles':
            calculating_result = True
        #end if        
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj()
        input = self.input
        wfn_file,ptcl_file = self.list_output_files()
        if result_name=='orbitals':
            result.location = os.path.join(self.locdir,wfn_file)
            orbfile = self.get_prefix()+'.orbs.h5'
            result.orbfile = os.path.join(self.locdir,orbfile)
        elif result_name=='particles':
            result.location = os.path.join(self.locdir,ptcl_file)
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if        
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        implemented = True
        input = self.input
        if isinstance(sim,Gamess):
            self.input_code = 'gamess'
            if result_name=='orbitals':
                orbpath = os.path.relpath(result.location,self.locdir)
                if result.scftyp=='mcscf':
                    input.gamess_ascii = orbpath
                    input.ci           = orbpath
                elif result.scftyp=='none': # cisd, etc
                    input.gamess_ascii = orbpath
                    input.ci           = orbpath
                    if result.mos>0:
                        input.read_initial_guess = result.mos
                    elif result.norbitals>0:
                        input.read_initial_guess = result.norbitals
                    #end if
                else:
                    input.gamess_ascii = orbpath
                #end if
                self.job.app_command = input.app_command()
            else:
                implemented = False
            #end if
        elif isinstance(sim,Pyscf):
            self.input_code = 'pyscf'
            if result_name=='orbitals':
                orbpath = os.path.relpath(result.h5_file,self.locdir)
                input.orbitals = orbpath
            else:
                implemented = False
            #end if
        elif isinstance(sim,QuantumPackage):
            self.input_code = 'qp'
            if result_name=='orbitals':
                orbpath = os.path.relpath(result.outfile,self.locdir)
                input.orbitals = orbpath
            else:
                implemented = False
            #end if
        else:
            implemented = False
        #end if
        if not implemented:
            self.error('ability to incorporate result "{0}" from {1} has not been implemented'.format(result_name,sim.__class__.__name__))
        #end if
    #end def incorporate_result


    def check_sim_status(self):
        output = open(os.path.join(self.locdir,self.outfile),'r').read()
        #errors = open(os.path.join(self.locdir,self.errfile),'r').read()

        # Recent versions of convert4qmc no longer produce the orbs.h5 file.
        # Instead, the file produced directly by e.g. Pyscf is used instead.
        # Therefore, make a symlink to the previously produced file in 
        # place of the orbs.h5 file.
        orbs     = self.input.orbitals
        finished = self.job.finished
        h5_orbs  = orbs is not None and orbs.endswith('.h5')
        if finished and h5_orbs:
            orbfile     = self.get_prefix()+'.orbs.h5'
            orbfilepath = os.path.join(self.locdir,orbfile)
            h5_orbs_missing = not os.path.exists(orbfilepath)
            if h5_orbs_missing:
                cwd = os.getcwd()
                os.chdir(self.locdir)
                os.system('ln -s {} {}'.format(orbs,orbfile))
                os.chdir(cwd)
            #end if
        #end if

        success = 'QMCGaussianParserBase::dump' in output
        for filename in self.list_output_files():
            success &= os.path.exists(os.path.join(self.locdir,filename))
        #end for

        self.failed = not success
        self.finished = self.job.finished
    #end def check_sim_status


    def get_output_files(self):
        output_files = []
        return output_files
    #end def get_output_files


    def app_command(self):
        return self.input.app_command()
    #end def app_command
#end class Convert4qmc



def generate_convert4qmc(**kwargs):
    sim_args,inp_args = Simulation.separate_inputs(kwargs)
    if 'identifier' in sim_args and not 'prefix' in inp_args:
        inp_args.prefix = sim_args.identifier
    #end if

    if not 'input' in sim_args:
        sim_args.input = generate_convert4qmc_input(**inp_args)
    #end if
    convert4qmc = Convert4qmc(**sim_args)

    return convert4qmc
#end def generate_convert4qmc






class PyscfToAfqmcInput(SimulationInput):

    input_order = '''
        help
        input
        output
        wavefunction
        qmcpack_input
        cholesky_threshold
        kpoint
        gdf
        ao
        cas
        disable_ham
        num_dets
        real_ham
        verbose
        '''.split()

    input_flags = obj(
        help               = 'h',
        input              = 'i',
        output             = 'o',
        wavefunction       = 'w',
        qmcpack_input      = 'q',
        cholesky_threshold = 't',
        kpoint             = 'k',
        gdf                = 'g',
        ao                 = 'a',
        cas                = 'c',
        disable_ham        = 'd',
        num_dets           = 'n',
        real_ham           = 'r',
        verbose            = 'v',
        )

    input_types = obj(
        app_name           = str,
        help               = bool,
        input              = str,
        output             = str,
        wavefunction       = str,
        qmcpack_input      = str,
        cholesky_threshold = float,
        kpoint             = bool,
        gdf                = bool,
        ao                 = bool,
        cas                = tuple,
        disable_ham        = bool,
        num_dets           = int,
        real_ham           = int,
        verbose            = bool,
        )

    input_defaults = obj(
        app_name           = 'pyscf_to_afqmc.py',
        help               = False,
        input              = None,
        output             = None,
        wavefunction       = None,
        qmcpack_input      = None,
        cholesky_threshold = None,
        kpoint             = False,
        gdf                = False,
        ao                 = False,
        cas                = None,
        disable_ham        = False,
        num_dets           = None,
        real_ham           = None,
        verbose            = False,
        )


    def __init__(self,**kwargs):
        # reassign inputs provided via short flag names
        for k,v in PyscfToAfqmcInput.input_flags.items():
            if v in kwargs:
                kwargs[k] = kwargs.pop(v)
            #end if
        #end for

        # check that only allowed keyword inputs are provided
        invalid = set(kwargs.keys())-set(self.input_types.keys())
        if len(invalid)>0:
            self.error('invalid inputs encountered\ninvalid keywords: {0}\nvalid keyword inputs are: {1}'.format(sorted(invalid),sorted(self.input_types.keys())))
        #end if

        # assign inputs
        self.set(**kwargs)

        # assign default values
        self.set_optional(**self.input_defaults)

        # check that all keyword inputs are valid
        self.check_valid()
    #end def __init__


    def check_valid(self,exit=True):
        valid = True
        # check that all inputs have valid types and assign them
        for k,v in self.items():
            if v is not None and not isinstance(v,self.input_types[k]):
                valid = False
                if exit:
                    self.error('keyword input "{0}" must be of type "{1}"\nyou provided a value of type "{2}"\nplease revise your input and try again'.format(k,self.input_types[k].__name__),v.__class__.__name__)
                #end if
                break
            #end if
        #end for
        if 'cas' in self and self.cas is not None:
            if len(self.cas)!=2:
                valid = False
                if exit:
                    self.error('keyword input "cas" must contain only two elements\nnumber of elements provided: {}\nvalue provided: {}'.format(len(self.cas),self.cas))
                #end if
            #end if
            noninteger = False
            for v in self.cas:
                noninteger |= not isinstance(v,int)
            #end for
            if noninteger:
                valid = False
                if exit:
                    self.error('keyword input "cas" must contain two integers\nvalue provided: {}'.format(self.cas))
                #end if
            #end if
        #end if
        return valid
    #end def check_valid


    def is_valid(self):
        return self.check_valid(exit=False)
    #end def is_valid


    def set_app_name(self,app_name):
        self.app_name = app_name
    #end def set_app_name


    def app_command(self):
        self.check_valid()
        c = self.app_name
        for k in self.input_order:
            if k in self:
                v = self[k]
                n = self.input_flags[k]
                if isinstance(v,bool):
                    if v:
                        c += ' -{0}'.format(n)
                    #end if
                elif isinstance(v,tuple):
                    vs = ''
                    for tv in v:
                        vs += '{},'.format(tv)
                    #end for
                    c += ' -{0} {1}'.format(n,vs[:-1])
                elif v is not None:
                    c += ' -{0} {1}'.format(n,str(v))
                #end if
            #end if
        #end for
        return c
    #end def app_command


    def read(self,filepath):
        None
    #end def read


    def write_text(self,filepath=None):
        return self.app_command()
    #end def write_text
#end class PyscfToAfqmcInput



def generate_pyscf_to_afqmc_input(**kwargs):
    return PyscfToAfqmcInput(**kwargs)
#end def generate_pyscf_to_afqmc_input



class PyscfToAfqmcAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0):
        if isinstance(arg0,Simulation):
            self.infile = arg0.infile
        else:
            self.infile = arg0
        #end if
    #end def __init__

    def analyze(self):
        None
    #end def analyze
#end class PyscfToAfqmcAnalyzer



class PyscfToAfqmc(Simulation):
    input_type             = PyscfToAfqmcInput
    analyzer_type          = PyscfToAfqmcAnalyzer
    generic_identifier     = 'pyscf2afqmc'
    application            = 'pyscf_to_afqmc.py'
    application_properties = set(['serial'])
    application_results    = set(['wavefunction','hamiltonian'])
    renew_app_command      = True


    def set_app_name(self,app_name):
        self.app_name = app_name
        self.input.set_app_name(app_name)
    #end def set_app_name


    def check_result(self,result_name,sim):
        calculating_result = False
        if result_name=='wavefunction':
            calculating_result = self.input.output is not None
        elif result_name=='hamiltonian':
            calculating_result = self.input.output is not None
        #end if        
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj()
        input = self.input
        if result_name in ('wavefunction','hamiltonian'):
            result.h5_file = os.path.join(self.locdir,input.output)
            if input.qmcpack_input is not None:
                result.xml = os.path.join(self.locdir,input.qmcpack_input)
            #end if
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if        
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        implemented = True
        input = self.input
        if isinstance(sim,Pyscf):
            if result_name=='wavefunction':
                chkfile = os.path.relpath(result.chkfile,self.locdir)
                input.input = chkfile
            else:
                implemented = False
            #end if
        else:
            implemented = False
        #end if
        if not implemented:
            self.error('ability to incorporate result "{0}" from {1} has not been implemented'.format(result_name,sim.__class__.__name__))
        #end if
    #end def incorporate_result       


    def check_sim_status(self):
        output = open(os.path.join(self.locdir,self.outfile),'r').read()

        success = '# Finished.' in output
        success &= os.path.exists(os.path.join(self.locdir,self.input.output))

        self.failed = not success
        self.finished = self.job.finished
    #end def check_sim_status


    def get_output_files(self):
        output_files = []
        return output_files
    #end def get_output_files


    def app_command(self):
        return self.input.app_command()
    #end def app_command
#end class PyscfToAfqmc



def generate_pyscf_to_afqmc(**kwargs):
    sim_args,inp_args = Simulation.separate_inputs(kwargs)
    if 'identifier' in sim_args:
        if 'output' not in inp_args:
            inp_args.output = '{}.afqmc.h5'.format(sim_args.identifier)
        #end if
        if 'qmcpack_input' not in inp_args:
            inp_args.qmcpack_input = '{}.afqmc.xml'.format(sim_args.identifier)
        #end if
    #end if

    if not 'input' in sim_args:
        sim_args.input = generate_pyscf_to_afqmc_input(**inp_args)
    #end if
    sim = PyscfToAfqmc(**sim_args)

    return sim
#end def generate_pyscf_to_afqmc
