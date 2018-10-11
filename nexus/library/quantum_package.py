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
from generic import obj
from execute import execute
from simulation import Simulation
from quantum_package_input import QuantumPackageInput,generate_quantum_package_input
from quantum_package_analyzer import QuantumPackageAnalyzer



class QuantumPackage(Simulation):
    input_type         = QuantumPackageInput
    analyzer_type      = QuantumPackageAnalyzer
    generic_identifier = 'qp'
    infile_extension   = '.ezfio'
    application        = 'qp_run'
    application_properties = set(['serial','mpi'])
    application_results    = set([]) 

    qprc = None

    slave_partners = obj(
        SCF     = 'qp_ao_ints',
        fci_zmq = 'selection_davidson_slave',
        )

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
        #end if
    #end def settings

    @staticmethod
    def restore_default_settings():
        QuantumPackage.qprc = None
    #end def restore_default_settings


    def post_init(self):
        qprc = QuantumPackage.qprc
        if qprc is None:
            self.error('cannot run quantum package\nplease provide path to quantum_package.rc in settings via argument "qprc"')
        #end if
        self.job.presub += '\nsource {0}\n'.format(os.path.abspath(qprc))
    #end def post_init


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


    def get_slave(self):
        rc = self.input.run_control
        sp = QuantumPackage.slave_partners
        slave = None
        if 'slave' in rc:
            slave = rc.slave
        elif rc.run_type in sp:
            slave = sp[rc.run_type]
        #end if
        return slave
    #end def get_slave


    def app_command(self):
        run_type = self.input.run_control.run_type
        app_command = self.app_name+' '+run_type+' '+self.infile
        job = self.job
        slave = self.get_slave()
        split_nodes  = job.nodes>1 and job.full_command is None
        split_nodes &= slave is not None
        if split_nodes:
            slave_command = self.app_name+' -slave {0} {1}'.format(slave,self.infile)
            outfile = self.outfile
            errfile = self.errfile
            prefix,ext = outfile.split('.',1)
            slave_outfile = prefix+'_slave.'+ext
            prefix,ext = errfile.split('.',1)
            slave_errfile = prefix+'_slave.'+ext

            job.divert_out_err()
            job1,job2 = job.split_nodes(1)
            job1.app_command = app_command
            job2.app_command = slave_command
            s  = ''
            s += job1.run_command()+' >{0} 2>{1}&\n'.format(outfile,errfile)
            s += 'sleep {0}\n'.format(self.input.run_control.sleep)
            s += job2.run_command()+' >{0} 2>{1}\n'.format(slave_outfile,slave_errfile)
            job.full_command = s
        #end if

        return app_command
    #end def app_command

#end class QuantumPackage



def generate_quantum_package(**kwargs):
    sim_args,inp_args = QuantumPackage.separate_inputs(kwargs)

    if not 'input' in sim_args:
        if 'input_type' in inp_args:
            input_type = inp_args.input_type
            del inp_args.input_type
        #end if
        sim_args.input = generate_quantum_package_input(**inp_args)
    #end if
    qp = QuantumPackage(**sim_args)

    return qp
#end def generate_quantum_package
    
