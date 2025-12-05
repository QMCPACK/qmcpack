##################################################################
##  (c) Copyright 2020-  by Jaron T. Krogel                     ##
##################################################################


import os
from simulation import Simulation
from generic import obj
from rmg_input import RmgInput,generate_rmg_input
from rmg_analyzer import RmgAnalyzer
from execute import execute



class Rmg(Simulation):
    input_type             = RmgInput
    analyzer_type          = RmgAnalyzer
    generic_identifier     = 'rmg'
    application            = 'rmg-cpu' 
    application_properties = set(['serial','mpi'])
    application_results    = set(['orbitals','wavefunctions'])


    def __init__(self,**sim_args):
        sync_from_scf = sim_args.pop('sync_from_scf',True)
        Simulation.__init__(self,**sim_args)
        self.sync_from_scf = False
        #calc = self.input.ontrol.get('calculation_mode',None)
        if self.get('calculation_mode')=='NSCF':
            self.sync_from_scf = sync_from_scf
        #end if
    #end def __init__




    def check_result(self,result_name,sim):
        calculating_result = False
        input =self.input
        if result_name=='orbitals': 
            conv_requested  = self.input.write_qmcpack_restart
            calculating_result = conv_requested 
        if result_name=='wavefunctions' :
            calculating_result = True 
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj() 
        input = self.input
        outdir = 'Waves'
        result.locdir   = self.locdir
        result.outdir   = os.path.join(self.locdir,outdir)
        if result_name=='orbitals':
            h5file = 'wave.out.h5'
            result.h5file = os.path.join(self.locdir,outdir,h5file)
        if result_name=='wavefunctions':
            result.location = os.path.join(self.locdir,outdir,'wave.out')
        #self.error('Ability to get result '+result_name+' has not been implemented.')
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        #self.error('ability to incorporate result '+result_name+' has not been implemented')
        if result_name=='wavefunctions':
            res_path = os.path.abspath(result.locdir)
            loc_path = os.path.abspath(self.locdir)
            if res_path==loc_path:
                None # don't need to do anything if in same directory
            elif self.sync_from_scf: # rsync output into nscf dir
                outdir = os.path.join(self.locdir,'Waves')
                #outdir = os.path.join(self.locdir,c.outdir)
                command = 'rsync -av {0}/* {1}/'.format(result.outdir,outdir)
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                #end if
                sync_record = os.path.join(outdir,'nexus_sync_record')
                if not os.path.exists(sync_record):
                    execute(command)
                    f = open(sync_record,'w')
                    f.write('\n')
                    f.close()
                #end if
            #end if
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

    if 'input' not in sim_args:
        sim_args.input = generate_rmg_input(**inp_args)
    #end if
    rmg = Rmg(**sim_args)

    return rmg
#end def generate_rmg
