##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  qmcpack_converters.py                                             #
#    Nexus interfaces for orbital converter tools of QMCPACK:        #
#    pw2qmcpack, wfconvert, and convert4qmc                          #
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
#    WfconvertInput                                                  #
#      SimulationInput class for wfconvert.                          #     
#                                                                    #
#    WfconvertAnalyzer                                               #
#      SimulationAnalyzer class for wfconvert.                       #        
#                                                                    #
#    Wfconvert                                                       #
#      Simulation class for wfconvert.                               #
#                                                                    #
#    generate_wfconvert                                              #
#      User-facing function to generate wfconvert simulation objects.#
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
    for v in floats:
        var_types[v]=float
    for v in strs:
        var_types[v]=str
    for v in bools:
        var_types[v]=bool

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
        for sname,section in self.iteritems():
            contents+='&'+sname+'\n'
            for name,value in section.iteritems():
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
            for name,value in vars.iteritems():
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

    def analyze(self):
        if False:
            import h5py
            self.log('Fixing h5 file',n=5)

            path = os.path.split(self.h5file)[0]
            print os.getcwd()
            print os.listdir('./')
            if os.path.exists(path):
                print os.listdir(path)
            #end if
            print self.h5file

            h = h5py.File(self.h5file)
            if 'electrons' in h:
                elec = h['electrons']
                nkpoints = 0
                for name,val in elec.iteritems():
                    if name.startswith('kpoint'):
                        nkpoints+=1
                    #end for
                #end if
                nkold = elec['number_of_kpoints'][0] 
                self.log('Were',nkold,'kpoints, now',nkpoints,'kpoints',n=6)
                elec['number_of_kpoints'][0] = nkpoints
            #end for        
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
    application_results    = set(['orbitals'])

    def check_result(self,result_name,sim):
        calculating_result = False
        inputpp = self.input.inputpp
        if result_name=='orbitals':
            calculating_result = True
        else:
            calculating_result = False
            self.error('ability to check for result '+result_name+' has not been implemented')
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
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if        
        return result
    #end def get_result

    def incorporate_result(self,result_name,result,sim):
        if result_name=='orbitals':
            pwin = sim.input.control
            p2in = self.input.inputpp
            pwprefix = 'pwscf'
            p2prefix = 'pwscf'
            pwoutdir = './'
            p2outdir = './'
            if 'prefix' in pwin:
                pwprefix = pwin.prefix
            if 'prefix' in p2in:
                p2prefix = p2in.prefix
            if 'outdir' in pwin:
                pwoutdir = pwin.outdir
            if 'outdir' in p2in:
                p2outdir = p2in.outdir
            if pwoutdir.startswith('./'):
                pwoutdir = pwoutdir[2:]
            if p2outdir.startswith('./'):
                p2outdir = p2outdir[2:]
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
            self.error('ability to incorporate result '+result_name+' has not been implemented')
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















class WfconvertInput(SimulationInput):
    def __init__(self,app_name='wfconvert',h5in='MISSING.h5',h5out='wfconvert.h5',spline=False,format='eshdf',factor=None):
        self.app_name = app_name
        self.h5in = h5in
        self.h5out= h5out
        self.spline = spline
        self.format = format
        self.factor = factor
    #end def __init__

#wfconvert --nospline --eshdf diamond.h5 out/diamond.pwscf.h5 >& diamond-wfconvert.out 
    def set_app_name(self,app_name):
        self.app_name = app_name
    #end def set_app_name

    def app_command(self):
        c = self.app_name+' '
        if not self.spline:
            c+= '--nospline '
        #end if
        c+='--'+self.format+' '+self.h5out+' '+self.h5in
        return c
    #end def app_command
        

    def read(self,filepath):
        None
    #end def read

    def write_text(self,filepath=None):
        return self.app_command()
    #end def write_text
#end class WfconvertInput


def generate_wfconvert_input(app_name='wfconvert',h5in='MISSING.h5',h5out='wfconvert.h5',spline=False,format='eshdf',factor=None):
    wi = WfconvertInput(
        app_name = app_name,
        h5in   = h5in,
        h5out  = h5out,
        spline = spline,
        format = format,
        factor = factor
        )
    return wi
#end def generate_wfconvert_input


class WfconvertAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0):
        if isinstance(arg0,Simulation):
            sim = arg0
            self.infile = sim.infile
            self.dir    = sim.locdir
            self.h5file = os.path.join(sim.locdir,sim.input.h5out)
        else:
            self.infile = arg0
        #end if
    #end def __init__

    def analyze(self):
        if False:
            import h5py
            self.log('Fixing h5 file',n=5)
            h = h5py.File(self.h5file)
            if 'electrons' in h:
                elec = h['electrons']
                nkpoints = 0
                for name,val in elec.iteritems():
                    if name.startswith('kpoint'):
                        nkpoints+=1
                    #end for
                #end if
                nkold = elec['number_of_kpoints'][0] 
                self.log('Were',nkold,'kpoints, now',nkpoints,'kpoints',n=6)
                elec['number_of_kpoints'][0] = nkpoints
            #end for        
        #end if
    #end def analyze
#end class WfconvertAnalyzer



class Wfconvert(Simulation):
    input_type             = WfconvertInput
    analyzer_type          = WfconvertAnalyzer
    generic_identifier     = 'wfconvert'
    application            = 'wfconvert'
    application_properties = set(['serial'])
    application_results    = set(['orbitals'])

    def set_app_name(self,app_name):
        self.app_name = app_name
        self.input.set_app_name(app_name)
    #end def set_app_name

    def check_result(self,result_name,sim):
        calculating_result = False
        if result_name=='orbitals':
            calculating_result = True
        else:
            calculating_result = False
            self.error('ability to check for result '+result_name+' has not been implemented')
        #end if        
        return calculating_result
    #end def check_result

    def get_result(self,result_name,sim):
        result = obj()
        if result_name=='orbitals':
            result.h5file   = os.path.join(self.locdir,self.input.h5out)
            result.outfile  = os.path.join(self.locdir,self.outfile)
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if        
        return result
    #end def get_result

    def incorporate_result(self,result_name,result,sim):
        if result_name=='orbitals':
            self.input.h5in = os.path.relpath(result.h5file,self.locdir)
            self.job.app_command = self.input.app_command()
        else:
            self.error('ability to incorporate result '+result_name+' has not been implemented')
        #end if                
    #end def incorporate_result

    def check_sim_status(self):
        outfile = os.path.join(self.locdir,self.outfile)
        errfile = os.path.join(self.locdir,self.errfile)
        fobj = open(outfile,'r')
        output = fobj.read()
        fobj.close()
        fobj = open(errfile,'r')
        errors = fobj.read()
        fobj.close()
        h5file = os.path.join(self.locdir,self.input.h5out)
        file_exists = os.path.exists(h5file)
        outfin = 'Successfully read' in errors and 'numSpins' in errors
        outfin = outfin and 'Writing laplacians' in output

        success = file_exists and outfin

        self.finished = success
    #end def check_sim_status

    def get_output_files(self):
        output_files = []
        return output_files
    #end def get_output_files

    def app_command(self):
        # app_name is passed along in post_init
        return self.input.app_command()
    #end def app_command
#end class Wfconvert




def generate_wfconvert(**kwargs):
    sim_args,inp_args = Simulation.separate_inputs(kwargs)

    if not 'input' in sim_args:
        sim_args.input = generate_wfconvert_input(**inp_args)
    #end if
    wfconvert = Wfconvert(**sim_args)

    return wfconvert
#end def generate_wfconvert









class Convert4qmcInput(SimulationInput):
    def __init__(self,
                 app_name           = 'convert4qmc',
                 prefix             = None,
                 gamess_ascii       = None,
                 ci                 = None,
                 read_initial_guess = None,
                 natural_orbitals   = None,
                 threshold          = None,
                 zero_ci            = False,
                 add_3body_J        = False
                 ):
        self.prefix             = prefix
        self.app_name           = app_name
        self.gamess_ascii       = gamess_ascii      
        self.ci                 = ci                
        self.read_initial_guess = read_initial_guess
        self.natural_orbitals   = natural_orbitals  
        self.threshold          = threshold         
        self.zero_ci            = zero_ci           
        self.add_3body_J        = add_3body_J       
    #end def __init__

    def set_app_name(self,app_name):
        self.app_name = app_name
    #end def set_app_name

    def app_command(self):
        c = self.app_name
        if self.prefix!=None:
            c += ' -prefix '+self.prefix
        #end if
        if self.gamess_ascii!=None:
            c += ' -gamessAscii '+self.gamess_ascii
        #end if
        if self.ci!=None:
            c += ' -ci '+self.ci
        #end if
        if self.threshold!=None:
            c += ' -threshold '+str(self.threshold)
        #end if
        if self.zero_ci:
            c += ' -zeroCi'
        #end if
        if self.read_initial_guess!=None:
            c += ' -readInitialGuess '+str(self.read_initial_guess)
        #end if
        if self.natural_orbitals!=None:
            c += ' -NaturalOrbitals '+str(self.natural_orbitals)
        #end if
        if self.add_3body_J:
            c += ' -add3BodyJ'
        #end if
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

    def set_app_name(self,app_name):
        self.app_name = app_name
        self.input.set_app_name(app_name)
    #end def set_app_name

    def propagate_identifier(self):
        None
        #self.input.prefix = self.identifier
    #end def propagate_identifier

    def check_result(self,result_name,sim):
        calculating_result = False
        if result_name=='orbitals':
            calculating_result = True
        elif result_name=='particles':
            calculating_result = True
        else:
            calculating_result = False
            self.error('ability to check for result '+result_name+' has not been implemented')
        #end if        
        return calculating_result
    #end def check_result

    def get_result(self,result_name,sim):
        result = obj()
        input = self.input
        if result_name=='orbitals':
            wfn_file,ptcl_file = input.output_files()
            result.location = os.path.join(self.locdir,wfn_file)
        elif result_name=='particles':
            wfn_file,ptcl_file = input.output_files()
            result.location = os.path.join(self.locdir,ptcl_file)
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if        
        return result
    #end def get_result

    def incorporate_result(self,result_name,result,sim):
        if result_name=='orbitals':
            orbpath = os.path.relpath(result.location,self.locdir)
            if result.scftyp=='mcscf':
                self.input.gamess_ascii = orbpath
                self.input.ci           = orbpath
            elif result.scftyp=='none': # cisd, etc
                self.input.gamess_ascii = orbpath
                self.input.ci           = orbpath
                if result.mos>0:
                    self.input.read_initial_guess = result.mos
                elif result.norbitals>0:
                    self.input.read_initial_guess = result.norbitals
                #end if
            else:
                self.input.gamess_ascii = orbpath
            #end if
            self.job.app_command = self.input.app_command()
        else:
            self.error('ability to incorporate result '+result_name+' has not been implemented')
        #end if                
    #end def incorporate_result

    def check_sim_status(self):
        output = open(os.path.join(self.locdir,self.outfile),'r').read()
        #errors = open(os.path.join(self.locdir,self.errfile),'r').read()

        success = 'QMCGaussianParserBase::dump' in output
        for filename in self.input.output_files():
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
