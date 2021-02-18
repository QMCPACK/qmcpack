##################################################################
##  (c) Copyright 2020-  by Jaron T. Krogel                     ##
##################################################################


from simulation import Simulation
from rmg_input import RmgInput,generate_rmg_input
from rmg_analyzer import RmgAnalyzer



class Rmg(Simulation):
    input_type             = RmgInput
    analyzer_type          = RmgAnalyzer
    generic_identifier     = 'rmg'
    application            = 'rmg-cpu' 
    application_properties = set(['serial','mpi'])
    application_results    = set([''])


    def check_result(self,result_name,sim):
        calculating_result = False
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = None
        self.error('Ability to get result '+result_name+' has not been implemented.')
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        self.error('ability to incorporate result '+result_name+' has not been implemented')
    #end def incorporate_result


    def app_command(self):
        return self.app_name+' '+self.infile
    #end def app_command


    def check_sim_status(self):
        # assume all is well
        self.succeeded = True
        self.failed    = False
        self.finished  = self.job.finished
    #end def check_sim_status


    def get_output_files(self): # returns list of output files to save
        return []
    #end def get_output_files
#end class Rmg




def generate_rmg(**kwargs):
    sim_args,inp_args = Rmg.separate_inputs(kwargs)

    if not 'input' in sim_args:
        sim_args.input = generate_rmg_input(**inp_args)
    #end if
    rmg = Rmg(**sim_args)

    return rmg
#end def generate_rmg
