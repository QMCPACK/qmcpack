
import os
from generic import obj
from simulation import Simulation,SimulationInput,SimulationAnalyzer


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

    def read_contents(self,contents):
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
    #end def read_contents

    def write_contents(self):
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
    #end def write_contents


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


def generate_pw2qmcpack_input(prefix='pwscf',outdir='pwscf_output',write_psir=True):
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


