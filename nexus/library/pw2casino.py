##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  pw2casino.py                                                      #
#    Nexus interface for the pw2casino tool distributed with PWSCF.  #
#    pseudopotentials between file formats.                          #
#                                                                    #
#  Content summary:                                                  #
#    generate_pw2casino                                              #
#      User-facing function to create pw2casino simulation objects.  #
#                                                                    #
#    generate_pw2casino_input                                        #
#      User-facing funcion to create input for pw2casino.            #
#                                                                    #
#    Pw2casino                                                       #
#      Simulation class for pw2casino.                               #
#                                                                    #
#    Pw2casinoInput                                                  #
#      SimulationInput class for pw2casino.                          #
#                                                                    #
#    Pw2casinoAnalyzer                                               #
#      SimulationAnalyzer class for pw2casino.                       #
#      Reads data from pw2casino log output into numeric form.       #
#                                                                    #
#====================================================================#


import os
from generic import obj
from unit_converter import convert
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


class Pw2casinoInput(SimulationInput):
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
                            print name in self.allowed
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
#end class Pw2casinoInput


def generate_pw2casino_input(prefix='pwscf',outdir='pwscf_output'):
    pw = Pw2casinoInput(
        prefix     = prefix,
        outdir     = outdir
        )
    return pw
#end def generate_pw2casino_input


 
class Pw2casinoAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0=None,arg1=None):
        if isinstance(arg0,Simulation):
            sim = arg0
            self.dir = sim.locdir
            self.infile = sim.infile
            self.outfile = sim.outfile
            prefix,outdir = sim.input.inputpp.tuple('prefix','outdir')
        elif arg0 is not None:
            if arg1 is None:
                self.error('must specify output file')
            #end if
            self.dir,self.infile = os.path.split(arg0)
            self.outfile = os.path.split(arg1)[1]
        #end if
    #end def __init__


    def analyze(self):
        outfile = os.path.join(self.dir,self.outfile)
        if not os.path.exists(outfile):
            self.error('output file does not exist\n  file path: '+outfile)
        #end if
        fobj = open(outfile,'r')
        lines = fobj.read().split('\n')
        fobj.close()

        self.units = 'Ha'

        energies = obj()
        for line in lines:
            ls = line.strip()
            name = None
            if ls.startswith('Kinetic energy'):
                name = 'kinetic'
            elif ls.startswith('Local energy'):
                name = 'local'
            elif ls.startswith('Non-Local energy'):
                name = 'nonlocal'
            elif ls.startswith('Ewald energy'):
                name = 'ewald'
            elif ls.startswith('xc contribution'):
                name = 'xc'
            elif ls.startswith('hartree energy'):
                name = 'hartree'
            elif ls.startswith('Total energy'):
                name = 'total'
            #end if
            if name!=None:
                energies[name] = float(ls.split()[2])
            #end if
        #end for
        self.energies = energies
    #end def analyze


    def change_units(self,new_unit):
        old_unit = self.units
        for name in self.energy.keys():
            self.energy[name] = convert(self.energy[name],old_unit,new_unit)
        #end for
        self.units = new_unit
    #end def change_units
        

    def get_result(self,result_name):
        self.not_implemented()
    #end def get_result
#end class Pw2casinoAnalyzer


class Pw2casino(Simulation):
    input_type = Pw2casinoInput
    analyzer_type = Pw2casinoAnalyzer
    generic_identifier = 'pw2casino'
    application = 'pw2casino.x'
    application_properties = set(['serial'])
    application_results    = set(['orbitals'])

    def check_result(self,result_name,sim):
        calculating_result = False
        inputpp = self.input.inputpp
        if result_name=='kinetic_energy':
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
        if result_name=='':
            None
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

        success = 'Total energy' in output

        self.finished = success and self.job.finished
    #end def check_sim_status


    def get_output_files(self):
        output_files = []
        return output_files
    #end def get_output_files


    def app_command(self):
        return self.app_name+'<'+self.infile
    #end def app_command
#end class Pw2casino




def generate_pw2casino(**kwargs):
    kw = set(kwargs.keys())
    sim_kw = kw & Simulation.allowed_inputs
    inp_kw = (kw - sim_kw)
    sim_args = dict()
    inp_args  = dict()
    for kw in sim_kw:
        sim_args[kw] = kwargs[kw]
    #end for
    for kw in inp_kw:
        inp_args[kw] = kwargs[kw]
    #end for    

    sim_args['input'] = generate_pw2casino_input(**inp_args)
    pw2casino = Pw2casino(**sim_args)

    return pw2casino
#end def generate_pw2casino


