
import os
from numpy import array,ndarray,abs
from generic import obj
from developer import DevBase
from debug import *
from simulation import Simulation
from gamess_input import GamessInput,generate_gamess_input,FormattedGroup,KeywordGroup,GuessGroup
from gamess_analyzer import GamessAnalyzer



class Gamess(Simulation):
    input_type         = GamessInput
    analyzer_type      = GamessAnalyzer
    generic_identifier = 'gamess'
    application        = 'gamess.x' 
    infile_extension   = '.inp'
    application_properties = set(['serial','mpi'])
    application_results    = set(['orbitals'])

    ericfmt = None
    mcppath = None

    @staticmethod
    def settings(ericfmt=None,mcppath=None):
        Gamess.ericfmt = ericfmt
        Gamess.mcppath = mcppath
    #end def settings


    def post_init(self):
        # gamess seems to need lots of environment variables to run properly
        # nearly all of these are names of output/work files
        # setup the environment to run gamess
        if not isinstance(self.ericfmt,str):
            self.error('you must set ericfmt with settings() or Gamess.settings()')
        #end if
        env = obj()
        for file,unit in GamessInput.file_units.iteritems():
            env[file] = '{0}.F{1}'.format(self.identifier,str(unit).zfill(2))
        #end for
        env.INPUT   = self.infile
        env.ERICFMT = self.ericfmt
        env.MCPPATH = self.mcppath
        self.job.set_environment(**env)
    #end def post_init


    def check_result(self,result_name,sim):
        input = self.input 
        if result_name=='orbitals':
            calculating_result = 'contrl' in input and 'scftyp' in input.contrl and input.contrl.scftyp in ('rhf','rohf','uhf')
        else:
            calculating_result = False
        #end if
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj()
        analyzer = self.load_analyzer_image()
        if result_name=='orbitals':
            result.location = os.path.join(self.locdir,self.outfile)
            result.vec = None
            result.norbitals = 0
            if 'punch' in analyzer and 'vec' in analyzer.punch:
                result.norbitals = analyzer.punch.norbitals
                result.vec = analyzer.punch.vec
            #end if
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        input = self.input
        if result_name=='orbitals':
            if result.vec is None or result.norbitals<1:
                self.error('could not obtain orbitals from previous GAMESS run')
            #end if
            if not 'guess' in input:
                input.guess = GuessGroup()
            #end if
            input.guess.set(
                guess = 'moread',
                norb  = result.norbitals
                )
            input.vec = FormattedGroup(result.vec)
        else:
            self.error('ability to incorporate result '+result_name+' has not been implemented')
        #end if
    #end def incorporate_result


    def app_command(self):
        return self.app_name+' '+self.infile.replace('.inp','')      
    #end def app_command


    def check_sim_status(self):
        output = open(os.path.join(self.locdir,self.outfile),'r').read()
        #errors = open(os.path.join(self.locdir,self.errfile),'r').read()
        
        self.failed = 'EXECUTION OF GAMESS TERMINATED -ABNORMALLY-' in output
        self.finished = self.failed or 'EXECUTION OF GAMESS TERMINATED NORMALLY' in output
    #end def check_sim_status


    def get_output_files(self):
        output_files = []
        return output_files
    #end def get_output_files
#end class Gamess



def generate_gamess(**kwargs):
    sim_args,inp_args = Simulation.separate_inputs(kwargs,copy_pseudos=False)

    if not 'input' in sim_args:
        sim_args.input = generate_gamess_input(**inp_args)
    #end if
    gamess = Gamess(**sim_args)

    return gamess
#end def generate_gamess














