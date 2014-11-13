
import os
from generic import obj
from simulation import Simulation,SimulationInput,SimulationAnalyzer
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
    sim_args,inp_args = Simulation.separate_inputs(kwargs,copy_pseudos=False)

    sim_args.input = generate_vasp_input(**inp_args)
    vasp = Vasp(**sim_args)

    return vasp
#end def generate_vasp





# VASP HT (hand template) classes and functions






class VaspHTAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0=None,xml=False,analyze=False):
        self.info = obj(xml=xml)
        prefix = None
        if isinstance(arg0,Simulation):
            sim = arg0
            infile = sim.infile
            path   = sim.locdir
        elif arg0!=None:
            path,infile = os.path.split(arg0)
            if infile=='':
                infile = None
            #end if
            if infile!=None:
                if not infile.endswith('INCAR'):
                    self.error('please provide the path to an INCAR file')
                #end if
                prefix = infile.replace('INCAR','').strip()
                if prefix=='':
                    prefix=None
                #end if
            #end if
        else:
            self.info.xml = False
            return
        #end if
        self.info.set(
            path   = path,
            infile = infile,
            prefix = prefix
            )
        if analyze:
            self.analyze()
        #end if
    #end def __init__

    def analyze(self):
        if self.info.xml:
            xmlfile = 'vasprun.xml'
            if self.info.prefix!=None:
                xmlfile = self.info.prefix+xmlfile
            #end if
            self.xmldata = read_vxml(os.path.join(self.info.path,xmlfile))
        #end if
    #end def analyze
#end class VaspHTAnalyzer




class VaspHT(Simulation):
    input_type         = SimulationInput
    analyzer_type      = VaspHTAnalyzer
    generic_identifier = 'vasp_ht'
    infile_extension   = None
    application        = 'vasp'
    application_properties = set(['serial','mpi'])
    application_results    = set([])

    allow_overlapping_files = True

    vasp_save_files = 'INCAR KPOINTS POSCAR CONTCAR DOSCAR EIGENVAL IBZKPT OSZICAR OUTCAR PCDAT XDATCAR vasprun.xml'.split()

    all_inputs  = 'INCAR KPOINTS POSCAR POTCAR'.split()
    all_outputs = 'CHG CHGCAR CONTCAR DOSCAR EIGENVAL IBZKPT OSZICAR OUTCAR PCDAT PROCAR WAVECAR XDATCAR'.split()

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
        command_line_args = ''
        return self.app_name + command_line_args
    #end def app_command


    def check_sim_status(self):
        success = False
        outpath = os.path.join(self.locdir,'OUTCAR')
        if os.path.exists(outpath):
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
#end class VaspHT



def generate_vasp_ht(**kwargs):
    sim_args,inp_args = Simulation.separate_inputs(kwargs,copy_pseudos=False)

    if not 'input' in sim_args:
        VaspHT.class_error('input keyword is required','generate_vasp_ht')
    #end if
    vasp_ht = VaspHT(**sim_args)

    return vasp_ht
#end def generate_vasp_ht
