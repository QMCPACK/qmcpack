##################################################################
##  (c) Copyright 2018-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  pyscf_sim.py                                                      #
#    Nexus interface for the PySCF simulation framework.             #
#                                                                    #
#  Content summary:                                                  #
#    Pyscf                                                           #
#      Simulation class for PySCF                                    #
#                                                                    #
#    generate_pyscf                                                  #
#      User-facing function to generate Pyscf simulation objects.    #
#====================================================================#


import os
from generic import obj
from execute import execute
from simulation import Simulation
from pyscf_input import PyscfInput,generate_pyscf_input
from pyscf_analyzer import PyscfAnalyzer
from developer import ci



class Pyscf(Simulation):
    input_type         = PyscfInput
    analyzer_type      = PyscfAnalyzer
    generic_identifier = 'pyscf'
    infile_extension   = '.py'
    application        = 'python'
    application_properties = set(['serial','mpi'])
    application_results    = set([]) 


    def check_result(self,result_name,sim):
        return False
    #end def check_result


    def get_result(self,result_name,sim):
        self.not_implemented()
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        not_implemented = False
        if not_implemented:
            self.error('ability to incorporate result '+result_name+' has not been implemented')
        #end if
    #end def incorporate_result


    def check_sim_status(self):
        # success of a generic pyscf script is too hard to assess
        # burden of when to initiate dependent simulations left to user
        self.failed   = False
        self.finished = self.job.finished
    #end def check_sim_status


    def get_output_files(self):
        output_files = []
        return output_files
    #end def get_output_files


    def app_command(self):
        app_command = self.app_name+' '+self.infile
        return app_command
    #end def app_command

#end class Pyscf



def generate_pyscf(**kwargs):
    sim_args,inp_args = Pyscf.separate_inputs(kwargs)

    if not 'input' in sim_args:
        if 'input_type' in inp_args:
            input_type = inp_args.input_type
            del inp_args.input_type
        #end if
        sim_args.input = generate_pyscf_input(**inp_args)
    #end if
    py = Pyscf(**sim_args)

    return py
#end def generate_pyscf

