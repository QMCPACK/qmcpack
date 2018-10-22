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
from developer import ci



class QuantumPackage(Simulation):
    input_type         = QuantumPackageInput
    analyzer_type      = QuantumPackageAnalyzer
    generic_identifier = 'qp'
    infile_extension   = '.ezfio'
    application        = 'qp_run'
    application_properties = set(['serial','mpi'])
    application_results    = set([]) 

    allow_overlapping_files = True

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

    def pre_init(self):
        prefix = self.input.run_control.prefix
        self.infile = prefix + self.infile_extension
    #end def pre_init


    def post_init(self):
        qprc = QuantumPackage.qprc
        if qprc is None:
            self.error('cannot run quantum package\nplease provide path to quantum_package.rc in settings via argument "qprc"')
        #end if
        self.job.presub += '\nsource {0}\n'.format(os.path.abspath(qprc))
    #end def post_init


    def write_prep(self):
        # write an ascii representation of the input changes
        infile = self.identifier+'.in'
        infile = os.path.join(self.locdir,infile)
        f = open(infile,'w')
        s = self.input.delete_optional('structure',None)
        f.write(str(self.input))
        if s is not None:
            self.input.structure = s
        #end if
        f.close()

        # copy ezfio directory from dependencies
        qp_dirs = []
        for dep in self.dependencies:
            dsim = dep.sim
            if isinstance(dsim,QuantumPackage):
                d_ezfio = os.path.join(dsim.locdir,dsim.infile)
                s_ezfio = os.path.join(self.locdir,self.infile)
                d_ezfio = os.path.abspath(d_ezfio)
                s_ezfio = os.path.abspath(s_ezfio)
                sync_record = os.path.join(self.locdir,self.identifier+'.sync_record')
                if s_ezfio!=d_ezfio:
                    qp_dirs.append(d_ezfio)
                    if not os.path.exists(sync_record):
                        if not os.path.exists(s_ezfio):
                            os.makedirs(s_ezfio)
                        #end if
                        command = 'rsync -av {0}/ {1}/'.format(d_ezfio,s_ezfio)
                        out,err,rc = execute(command)
                        if rc>0:
                            self.warn('rsync of ezfio directory failed\nall runs depending on this one will be blocked\nsimulation identifier: {0}\nlocal directory: {1}\nattempted rsync command: {2}'.format(self.identifier,self.locdir,command))
                            self.failed = True
                            self.block_dependents()
                        else:
                            f = open(sync_record,'w')
                            f.write(command+'\n')
                            f.close()
                        #end if
                    #end if
                #end if
            #end if
        #end for
        if len(qp_dirs)>1:
            qpd = ''
            for d in qp_dirs:
                qpd += d+'\n'
            #end for
            self.error('quantum package run depends on multiple others with distinct ezfio directories\ncannot determine which run to copy ezfio directory from\nezfio directories from prior runs:\n{0}'.format(qpd))
        #end if
    #end def write_prep


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

        # get run controls
        input = self.input
        rc = input.run_control

        # make the basic app command, no splitting etc
        run_type = rc.run_type
        app_command = self.app_name+' '+run_type+' '+self.infile

        # prepare local vars in case splitting or other tricks are needed
        fc = ''
        job = self.job

        # add cis-loop runs if requested
        if 'cis_loop' in rc:
            nloop = 0
            if isinstance(rc.cis_loop,bool) and rc.cis_loop:
                nloop = 2
            else:
                nloop = rc.cis_loop
            #end for
            if nloop>0:
                jloop = job.clone()
                fc+='\n'
                for n in range(nloop):
                    jloop.app_command = self.app_name+' cis '+self.infile
                    fc += jloop.run_command()+' >{0}_{1}.out 2>{0}_{1}.err\n'.format(self.identifier,n)
                    jloop.app_command = self.app_name+' save_natorb '+self.infile
                    fc += jloop.run_command()+'\n'
                #end for
                fc+='\n'
                integrals = [
                    'integrals_monoelec/disk_access_ao_one_integrals',
                    'integrals_monoelec/disk_access_mo_one_integrals',
                    'integrals_bielec/disk_access_ao_integrals',
                    'integrals_bielec/disk_access_mo_integrals',
                    ]
                cl = ''
                for integral in integrals:
                    isec,ivar = integral.split('/')
                    if input.present(ivar):
                        val = input.delete(ivar)
                        cl += 'echo "{0}" > {1}/{2}\n'.format(val,self.infile,integral)
                    #end if
                #end for
                if len(cl)>0:
                    fc+=cl+'\n'
                #end if
            #end if
        #end if

        # perform master-slave job splitting if necessary
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
            fc += job1.run_command()+' >{0} 2>{1}&\n'.format(outfile,errfile)
            fc += 'sleep {0}\n'.format(self.input.run_control.sleep)
            fc += job2.run_command()+' >{0} 2>{1}\n'.format(slave_outfile,slave_errfile)

            if 'davidson' in slave and not input.present('distributed_davidson'):
                input.set(distributed_davidson=True)
            #end if
        elif len(fc)>0:
            job.divert_out_err()
            job.app_command = app_command
            fc += job.run_command()+' >{0} 2>{1}\n'.format(self.outfile,self.errfile)
        #end if

        if len(fc)>0:
            job.full_command = fc
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
        if 'prefix' not in inp_args and 'identifier' in sim_args:
            inp_args['prefix'] = sim_args['identifier']
        #end if
        sim_args.input = generate_quantum_package_input(**inp_args)
    #end if
    qp = QuantumPackage(**sim_args)

    return qp
#end def generate_quantum_package
    
