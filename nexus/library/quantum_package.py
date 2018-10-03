##################################################################
##  (c) Copyright 2018-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  quantum_package.py                                                #
#    Nexus interface for the Quantum Package simulation code.        #
#                                                                    #
#  Content summary:                                                  #
#    QuantumPackage                                                  #
#      Simulation class for Quantum Package.                         #
#                                                                    #
#    generate_quantum_package                                        #
#      User-facing function to generate Quantum Package simulation   #
#      objects.                                                      #
#====================================================================#


import os
from execute import execute
from simulation import Simulation
from quantum_package_input import QuantumPackageInput,generate_quantum_package_input
from quantum_package_analyzer import QuantumPackageAnalyzer



class QuantumPackage(Simulation):
    input_type         = QuantumPackageInput
    analyzer_type      = QuantumPackageAnalyzer
    generic_identifier = 'qp'
    application        = 'quantum_package' 
    application_properties = set(['serial','mpi'])
    application_results    = set([]) 

    qprc = None

    @staticmethod
    def settings(qprc=None):
        # path to quantum_package.rc file
        QuantumPackage.qprc = qprc
        if qprc is not None:
            if not isinstance(qprc,str):
                QuantumPackage.class_error('settings input "qprc" must be a path\nreceived type: {0}\nwith value: {1}'.format(qprc.__class__.__name__,qprc))
            elif not os.path.exists(qprc):
                QuantumPackage.class_error('quantum_package.rc file does not exist\nfile path provided via "qprc" in settings\nfile path: {0}'.format(qprc))
            #end if
            execute('source '+qprc)
        #end if
    #end def settings

    @staticmethod
    def restore_default_settings():
        QuantumPackage.qprc = None
    #end def restore_default_settings


    def check_result(self,result_name,sim):
        return False
    #end def check_result


    def get_result(self,result_name,sim):
        self.not_implemented()
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        self.not_implemented()
    #end def incorporate_result


    def check_sim_status(self):
        success = True
        self.finished = success and self.job.finished
    #end def check_sim_status


    def get_output_files(self):
        output_files = []
        return output_files
    #end def get_output_files


    def app_command(self):
        return ''
    #end def app_command

#end class QuantumPackage



def generate_quantum_package(**kwargs):
    sim_args,inp_args = QuantumPackage.separate_inputs(kwargs)

    if not 'input' in sim_args:
        input_type = inp_args.input_type
        del inp_args.input_type
        sim_args.input = generate_quantum_package_input(input_type,**inp_args)
    #end if
    qp = QuantumPackage(**sim_args)

    return qp
#end def generate_quantum_package
    
