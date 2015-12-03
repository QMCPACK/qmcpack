##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  sqd.py                                                            #
#    Nexus interface with the SQD atomic Hartree-Fock code shipped   #
#    with QMCPACK.                                                   #
#    Essentially an early fork of qmcpack_input.py.                  #
#                                                                    #
#  Content summary:                                                  #
#    Sqd                                                             #
#      Simulation class for SQD.                                     #
#      Handles passing orbital data to QMCPACK.                      #
#                                                                    #
#    generate_spd                                                    #
#      User-facing function to generate SQD simulation objects.      #
#                                                                    #
#                                                                    #
#====================================================================#


import os
from generic import obj
from physical_system import PhysicalSystem
from simulation import Simulation
from sqd_input import SqdInput,generate_sqd_input,hunds_rule_filling
from sqd_analyzer import SqdAnalyzer


class Sqd(Simulation):
    input_type             = SqdInput
    analyzer_type          = SqdAnalyzer
    generic_identifier     = 'sqd'
    infile_extension       = '.in.xml'
    application            = 'sqd'
    application_properties = set(['serial'])
    application_results    = set(['orbitals','jastrow'])


    def propagate_identifier(self):
        #self.input.simulation.project.id = self.identifier
        self.input.simulation.project.id = self.system.structure.elem[0]
    #end def propagate_identifier


    def check_result(self,result_name,sim):
        calculating_result = False
        if result_name=='orbitals':
            calculating_result = True
        elif result_name=='jastrow':
            calculating_result = True
        else:
            self.error('ability to check for result '+result_name+' has not been implemented')
        #end if        
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj()
        analyzer = self.load_analyzer_image()
        if result_name=='orbitals':
            outfiles = self.input.get_output_info(list=False)
            result.set(
                dir     = self.locdir,
                h5file  = outfiles.h5,
                qmcfile = outfiles.qmc
                )
        elif result_name=='jastrow':
            result.set(
                rcut = analyzer.find_rcut(qcut=1e-3),
                B    = 1./analyzer.moment(n=1)
                )
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if        
        del analyzer
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        self.error('ability to incorporate result '+result_name+' has not been implemented')
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

        outfiles = self.input.get_output_info(list=False)
        files_exist = True
        for file in outfiles:
            file_loc = os.path.join(self.locdir,file)
            files_exist = files_exist and os.path.exists(file_loc)
        #end for
        ran_to_end = False
        if files_exist:
            fobj = open(os.path.join(self.locdir,outfiles.log),'r')
            log = fobj.read()
            fobj.close()
            ran_to_end =  'E_tot' in log
        #end if

        aborted = 'Fatal Error' in errors

        self.failed   = aborted
        self.finished = files_exist and ran_to_end and self.job.finished and not aborted 
    #end def check_sim_status


    def get_output_files(self):
        output_files = self.input.get_output_info()
        return output_files
    #end def get_output_files


    def app_command(self):
        if self.job.app_name is None:
            app_name = self.app_name
        else:
            app_name = self.job.app_name
        #end if
        return app_name+' '+self.infile      
    #end def app_command
#end class Sqd







def generate_sqd(**kwargs):
    overlapping_kw = set(['system'])
    kw = set(kwargs.keys())
    sim_kw = kw & Simulation.allowed_inputs
    inp_kw = (kw - sim_kw) | (kw & overlapping_kw)    
    sim_args = dict()
    inp_args  = dict()
    for kw in sim_kw:
        sim_args[kw] = kwargs[kw]
    #end for
    for kw in inp_kw:
        inp_args[kw] = kwargs[kw]
    #end for    
    if 'system' in inp_args:
        sys = inp_args['system']
        if isinstance(sys,str):
            inp_args['system'] = str(sys)
        elif isinstance(sys,PhysicalSystem):
            inp_args['system'] = sys.copy()
        #end if
    #end if

    sim_args['input'] = generate_sqd_input(**inp_args)
    sqd = Sqd(**sim_args)

    return sqd
#end def generate_sqd
