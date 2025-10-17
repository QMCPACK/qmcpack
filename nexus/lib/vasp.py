##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  vasp.py                                                           #
#    Nexus interface for the VASP simulation code.                   #
#                                                                    #
#  Content summary:                                                  #
#    Vasp                                                            #
#      Simulation class for VASP.                                    #
#      Handles passing of structure dependencies for relax and NEB.  #
#                                                                    #
#    generate_vasp                                                   #
#      User-facing function to generate VASP simulation objects.     #
#                                                                    #
#    VaspHT                                                          #
#      Original VASP interface supporting only template inputs.      #
#                                                                    #
#    generate_vasp_ht                                                #
#      User interface to VaspHT.                                     #
#                                                                    #
#====================================================================#


import os
from generic import obj
from simulation import Simulation,SimulationInput,SimulationAnalyzer
from vasp_input import VaspInput,generate_vasp_input,generate_poscar,Poscar
from vasp_analyzer import VaspAnalyzer,read_vxml
from structure import Structure
from debug import *



class Vasp(Simulation):
    input_type         = VaspInput
    analyzer_type      = VaspAnalyzer
    generic_identifier = 'vasp'
    application        = 'vasp' 
    application_properties = set(['serial','mpi'])
    application_results    = set(['structure']) 

    allow_overlapping_files = True

    vasp_save_files = 'INCAR KPOINTS POSCAR CONTCAR DOSCAR EIGENVAL IBZKPT OSZICAR OUTCAR PCDAT XDATCAR vasprun.xml'.split()

    def set_files(self):
        self.infile  = 'INCAR'
        self.outfile = self.identifier + self.outfile_extension
        self.errfile = self.identifier + self.errfile_extension
    #end def set_files


    def check_result(self,result_name,sim):
        input = self.input
        if result_name=='structure':
            calculating_result = input.producing_structure()
        else:
            calculating_result = False
        #end if
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj()
        input = self.input
        if result_name=='structure':
            # get structure from CONTCAR
            ccfile = os.path.join(self.locdir,self.identifier+'.CONTCAR')
            if not os.path.exists(ccfile):
                self.error('CONTCAR file does not exist for relax simulation at '+self.locdir)
            #end if
            contcar = Poscar(ccfile)
            structure = Structure()
            if contcar.elem is not None:
                structure.read_poscar(ccfile)
            else:
                elem,elem_count = self.system.structure.order_by_species()
                structure.read_poscar(ccfile,elem=elem)
            #end if
            if input.poscar.dynamic is not None:
                structure.freeze(
                    input.poscar.dynamic,
                    negate = True
                    )
            #end if
            result.structure = structure
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        input = self.input
        if result_name=='structure':
            if input.performing_neb():
                if 'neb_structures' not in self:
                    self.neb_structures = []
                #end if
                neb_structures = self.neb_structures
                if len(neb_structures)>1:
                    self.error('NEB simulation at {0} depends on more than two structures\n  please check your inputs'.format(self.locdir))
                #end if
                neb_structures.append(result.structure.copy())
                if len(neb_structures)==2:
                    input.setup_neb(*neb_structures,images=input.incar.images)
                #end if
            else:
                input.poscar = generate_poscar(result.structure)
            #end if
        else:
            self.error('ability to incorporate result '+result_name+' has not been implemented')
        #end if  
    #end def incorporate_result


    def app_command(self):
        return self.app_name
    #end def app_command


    def check_sim_status(self):
        outpaths = []
        if not self.input.performing_neb():
            outpath = os.path.join(self.locdir,self.identifier+'.OUTCAR')
            exists  = os.path.exists(outpath)
            if not exists:
                outpath = os.path.join(self.locdir,'OUTCAR')
            #end if
            outpaths.append(outpath)
        else:
            for i in range(self.input.incar.images):
                outpaths.append(os.path.join(self.locdir,str(i+1).zfill(2),'OUTCAR'))
            #end for
        #end if
        success = False
        all_exist = True
        for outpath in outpaths:
            all_exist &= os.path.exists(outpath)
        #end for
        if all_exist:
            success = True
            for outpath in outpaths:
                outcar = open(outpath,'r').read()
                success &= 'General timing and accounting' in outcar
            #end for
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
    sim_args,inp_args = Vasp.separate_inputs(kwargs,copy_pseudos=False)

    sim_args.input = generate_vasp_input(**inp_args)
    vasp = Vasp(**sim_args)

    return vasp
#end def generate_vasp


