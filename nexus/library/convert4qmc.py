

import os
from generic import obj
from simulation import Simulation,SimulationInput,SimulationAnalyzer
from debug import *


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

    def write_contents(self):
        return self.app_command()
    #end def write_contents

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
    application_results    = set(['orbitals'])

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
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if        
        return result
    #end def get_result

    def incorporate_result(self,result_name,result,sim):
        if result_name=='orbitals':
            self.input.gamess_ascii = os.path.relpath(result.location,self.locdir)
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

    if not 'input' in sim_args:
        sim_args.input = generate_convert4qmc_input(**inp_args)
    #end if
    convert4qmc = Convert4qmc(**sim_args)

    return convert4qmc
#end def generate_convert4qmc
