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
from nexus_base import nexus_core
from simulation import Simulation
from quantum_package_input import QuantumPackageInput,generate_quantum_package_input,read_qp_value
from quantum_package_analyzer import QuantumPackageAnalyzer
from gamess import Gamess
from developer import ci



class QuantumPackage(Simulation):
    input_type         = QuantumPackageInput
    analyzer_type      = QuantumPackageAnalyzer
    generic_identifier = 'qp'
    infile_extension   = '.ezfio'
    application        = 'qp_run'
    application_properties = set(['serial','mpi'])
    application_results    = set(['orbitals']) 

    allow_overlapping_files = True

    qprc = None

    slave_partners = obj(
        scf = 'scf',
        fci = 'fci',
        )

    @staticmethod
    def settings(qprc=None):
        # path to quantum_package.rc file
        QuantumPackage.qprc = qprc
        if qprc is not None and not nexus_core.status_only:
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
                        if rc!=0:
                            self.warn('rsync of ezfio directory failed\nall runs depending on this one will be blocked\nsimulation identifier: {0}\nlocal directory: {1}\nattempted rsync command: {2}'.format(self.identifier,self.locdir,command))
                            self.failed = True
                            self.block_dependents()
                        else:
                            f = open(sync_record,'w')
                            f.write(command+'\n')
                            f.close()
                            execute('qp_edit -c {0}'.format(d_ezfio))
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
        calculating_result = False
        rc = self.input.run_control
        if result_name=='orbitals':
            calculating_result  = rc.run_type=='save_for_qmcpack'
            calculating_result |= rc.save_for_qmcpack
        #end if
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj()
        rc = self.input.run_control
        if result_name=='orbitals':
            if rc.run_type=='save_for_qmcpack':
                result.outfile = os.path.join(self.locdir,self.outfile)
            elif rc.save_for_qmcpack:
                result.outfile = os.path.join(self.locdir,'{0}_savewf.out'.format(self.identifier))
            else:
                self.error("cannot get orbitals\ntracking of save_for_qmcpack is somehow corrupted\nthis is a developer error")
            #end if
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        not_implemented = False
        if isinstance(sim,Gamess):
            if result_name=='orbitals':
                loc_file = self.input.run_control.prefix
                loc_out = os.path.join(self.locdir,loc_file)
                gms_out = result.outfile
                command = 'cp {0} {1}'.format(gms_out,loc_out)
                out,err,rc = execute(command)
                if rc!=0:
                    self.warn('copying GAMESS output failed\nall runs depending on this one will be blocked\nsimulation identifier: {0}\nlocal directory: {1}\nattempted command: {2}'.format(self.identifier,self.locdir,command))
                    self.failed = True
                    self.block_dependents()
                #end if
                command = 'qp_convert_output_to_ezfio '+loc_file
                cwd = os.getcwd()
                os.chdir(self.locdir)
                out,err,rc = execute(command)
                os.chdir(cwd)
                if rc!=0:
                    self.warn('creation of ezfio file from GAMESS output failed\nall runs depending on this one will be blocked\nsimulation identifier: {0}\nlocal directory: {1}\nattempted command: {2}'.format(self.identifier,self.locdir,command))
                    self.failed = True
                    self.block_dependents()
                #end if
            else:
                not_implemented = True
            #end if
        else:
            not_implemented = True
        #end if
        if not_implemented:
            self.error('ability to incorporate result "{}" from {} has not been implemented',result_name,sim.__class__.__name__)
        #end if
    #end def incorporate_result


    def attempt_files(self):
        return (self.outfile,self.errfile)
    #end def attempt_files


    def check_sim_status(self):
        # get the run type
        input = self.input
        rc = self.input.run_control
        scf    = rc.run_type=='scf'
        sel_ci = rc.run_type=='fci'

        # assess successful completion of the run
        #   currently a check only exists for HF/SCF runs
        #   more sophisticated checks can be added in the future
        failed = False
        if scf:
            outfile = os.path.join(self.locdir,self.outfile)
            f = open(outfile,'r')
            output = f.read()
            f.close()
            hf_not_converged = '* SCF energy' not in output
            failed |= hf_not_converged
        #end if
        self.failed = failed
        self.finished = self.job.finished

        # check to see if the job needs to be restarted
        conv_dets = 'converge_dets' in rc and rc.converge_dets
        n_det_max = input.get('n_det_max')
        if sel_ci and conv_dets and n_det_max is not None:
            n_det = None
            n_det_path = os.path.join(self.locdir,self.infile,'determinants/n_det')
            if os.path.exists(n_det_path):
                n_det = read_qp_value(n_det_path)
                if isinstance(n_det,int) and n_det<n_det_max:
                    self.save_attempt()
                    input.set(read_wf=True)
                    self.reset_indicators()
                #end if
            #end if
        #end if
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
                    'ao_one_e_ints/io_ao_one_e_integrals',
                    'mo_one_e_ints/io_mo_one_e_integrals',
                    'ao_two_e_ints/io_ao_two_e_integrals',
                    'mo_two_e_ints/io_mo_two_e_integrals',
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

        # check for post-processing operations and save the job in current state
        postprocessors = ['save_natorb',
                          'four_idx_transform',
                          'save_for_qmcpack']
        postprocess = obj()
        jpost = None
        for pp in postprocessors:
            if pp in rc and rc[pp]:
                postprocess[pp] = True
                if jpost is None:
                    jpost = job.clone()
                #end if
            else:
                postprocess[pp] = False
            #end if
        #end for

        # perform master-slave job splitting if necessary
        slave = self.get_slave()
        split_nodes  = job.nodes is not None and job.nodes>1 and job.full_command is None
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

            if 'fci' in slave and not input.present('distributed_davidson'):
                input.set(distributed_davidson=True)
            #end if
        elif len(fc)>0 or jpost is not None:
            job.divert_out_err()
            job.app_command = app_command
            fc += job.run_command()+' >{0} 2>{1}\n'.format(self.outfile,self.errfile)
        #end if

        if postprocess.save_natorb:
            jno = jpost.serial_clone()
            fc += '\n'
            jno.app_command = self.app_name+' save_natorb '+self.infile
            fc += jno.run_command()+' >{0}_natorb.out 2>{0}_natorb.err\n'.format(self.identifier)
        #end if

        if postprocess.four_idx_transform:
            jfit = jpost.serial_clone()
            fc += '\n'
            fc += 'echo "Write" > {}/mo_two_e_ints/io_mo_two_e_integrals\n'.format(self.infile)
            jfit.app_command = self.app_name+' four_idx_transform '+self.infile
            fc += jfit.run_command()+' >{0}_fit.out 2>{0}_fit.err\n'.format(self.identifier)
        #end if

        if postprocess.save_for_qmcpack:
            jsq = jpost.serial_clone()
            fc += '\n'
            jsq.app_command = self.app_name+' save_for_qmcpack '+self.infile
            fc += jsq.run_command()+' >{0}_savewf.out 2>{0}_savewf.err\n'.format(self.identifier)
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
    
