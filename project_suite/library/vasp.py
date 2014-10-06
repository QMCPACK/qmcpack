
import os
from generic import obj
from simulation import Simulation
from vasp_input import VaspInput,generate_vasp_input
from vasp_analyzer import VaspAnalyzer




class Vasp(Simulation):
    input_type         = VaspInput
    analyzer_type      = VaspAnalyzer
    generic_identifier = 'vasp'
    application        = 'vasp' 
    application_properties = set(['serial','mpi'])
    application_results    = set([]) 

    allow_overlapping_files = True

    vasp_save_files = 'INCAR KPOINTS POSCAR CONTCAR DOSCAR EIGENVAL IBZKPT OSZICAR OUTCAR PCDAT XDATCAR vasprun.xml'.split()

    def set_files(self):
        self.infile  = 'INCAR'
        self.outfile = self.identifier + self.outfile_extension
        self.errfile = self.identifier + self.errfile_extension
    #end def set_files


    def check_result(self,result_name,sim):
        return False
    #end def check_result


    def get_result(self,result_name,sim):
        self.not_implemented()
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        self.not_implemented()
    #end def incorporate_result


    def app_command(self):
        return self.app_name
    #end def app_command


    def check_sim_status(self):
        success = False
        outpath = os.path.join(self.locdir,self.identifier+'.OUTCAR')
        exists  = os.path.exists(outpath)
        if not exists:
            outpath = os.path.join(self.locdir,'OUTCAR')
            exists  = os.path.exists(outpath)
        #end if
        if exists:
            outcar = open(outpath,'r').read()
            success = 'General timing and accounting' in outcar
        #end if
        self.finished = success
    #end def check_sim_status


    def get_output_files(self):
        output_files = []
        for file in self.vasp_save_files:
            native_file = os.path.join(self.locdir,file)
            save_file   = os.path.join(self.locdir,self.identifier+'.'+file)
            if os.path.exists(native_file):
                os.system('cp {0} {1}'.format(native_file,save_file))
                output_files.append(file)
            #end if
        #end for
        return output_files
    #end def get_output_files
#end class Vasp



def generate_vasp(**kwargs):
    sim_args,inp_args = Simulation.separate_inputs(kwargs)

    sim_args.input = generate_vasp_input(**inp_args)
    vasp = Vasp(**sim_args)

    return vasp
#end def generate_vasp
