import os
from numpy import array
from generic import obj
from physical_system import PhysicalSystem
from converters import Pw2qmcpack
from simulation import Simulation
from pwscf_input import PwscfInput,generate_pwscf_input
from pwscf_analyzer import PwscfAnalyzer


class Pwscf(Simulation):
    input_type = PwscfInput
    analyzer_type = PwscfAnalyzer
    generic_identifier = 'pwscf'
    application = 'pw.x'
    application_properties = set(['serial','mpi'])
    application_results    = set(['charge_density','orbitals','structure'])

    #def propagate_identifier(self):
    #    self.input.control.prefix = self.identifier        
    ##end def propagate_identifier


    def __init__(self,**sim_args):
        has_group_atoms = 'group_atoms' in sim_args
        group_atoms = True
        if has_group_atoms:
            group_atoms = sim_args['group_atoms']
            del sim_args['group_atoms']
        #end if
        Simulation.__init__(self,**sim_args)
        if group_atoms and isinstance(self.system,PhysicalSystem):
            self.system.structure.group_atoms()
        #end if
    #end def post_init


    def write_prep(self):
        #make sure the output directory exists
        outdir = os.path.join(self.locdir,self.input.control.outdir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        #end if
    #end def write_prep


    def check_result(self,result_name,sim):
        input = self.input
        control = input.control
        if result_name=='charge_density':
            calculating_result = True
        elif result_name=='orbitals':
            calculating_result = 'wf_collect' in control and control.wf_collect
            #if isinstance(sim,Pw2qmcpack):
            #    sim.input.inputpp.prefix = self.identifier
            ##end if
        elif result_name=='structure':
            calculating_result = control.calculation.lower() == 'relax'
        else:
            calculating_result = False
            self.error('ability to check for result '+result_name+' has not been implemented')
        #end if
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj()        
        input = self.input
        control = input.control
        prefix = 'pwscf'
        outdir = './'
        if 'prefix' in control:
            prefix = control.prefix
        #end if
        if 'outdir' in control:
            outdir = control.outdir
        #end if
        if outdir.startswith('./'):
            outdir = outdir[2:]
        #end if
        if result_name=='charge_density':
            result.location = os.path.join(self.locdir,outdir,prefix+'.save','charge-density.dat')
        elif result_name=='orbitals':
            result.location = os.path.join(self.locdir,outdir,prefix+'.wfc1')
        elif result_name=='structure':
            pa = self.load_analyzer_image()
            structs = pa.structures
            pos,atoms = structs[len(structs)-1].tuple('positions','atoms')
            scale = self.input.system['celldm(1)']
            pos   = scale*array(pos)
            atoms = array(atoms)
            result.structure = obj(
                positions = pos,
                atoms     = atoms
                )
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if
        return result
    #end def get_result

    def incorporate_result(self,result_name,result,sim):
        if result_name=='charge_density':
            c = self.input.control
            link_loc = os.path.join(self.locdir,c.outdir,c.prefix+'.save')
            cd_loc = result.location
            cd_rel = os.path.relpath(cd_loc,link_loc)
            cwd = os.getcwd()
            if not os.path.exists(link_loc):
                os.makedirs(link_loc)
            #end if
            os.chdir(link_loc)
            os.system('ln -s '+cd_rel+' charge-density.dat')
            os.chdir(cwd)
        elif result_name=='structure':
            structure = self.system.structure
            structure.change_units('B')
            relstruct = result.structure
            structure.set(
                pos   = relstruct.positions,
                atoms = relstruct.atoms
                )
            input = self.input
            preserve_kp = 'k_points' in input and 'specifier' in input.k_points and input.k_points.specifier=='automatic'
            if preserve_kp:
                kp = input.k_points.copy()
            #end if
            input.incorporate_system(self.system)
            if preserve_kp:
                input.k_points = kp
            #end if
        else:
            self.error('ability to incorporate result '+result_name+' has not been implemented')
        #end if        
    #end def incorporate_result

    def check_sim_status(self):
        outfile = os.path.join(self.locdir,self.outfile)
        fobj = open(outfile,'r')
        output = fobj.read()
        fobj.close()
        self.finished = 'JOB DONE' in output
    #end def check_sim_status

    def get_output_files(self):
        output_files = []
        return output_files
    #end def get_output_files

    def app_command(self):
        #ac = self.app_name+' -in '+self.infile
        ac = self.app_name+' -input '+self.infile
        return ac
    #end def app_command
#end class Pwscf



def generate_pwscf(**kwargs):
    has_input = 'input_type' in kwargs
    if has_input:
        input_type = kwargs['input_type']
        del kwargs['input_type']
    #end if
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
    if 'pseudos' in inp_args:
        if 'files' in sim_args:
            sim_args['files'] = list(sim_args['files'])
        else:
            sim_args['files'] = list()
        #end if
        sim_args['files'].extend(list(inp_args['pseudos']))
    #end if
    if 'system' in inp_args and isinstance(inp_args['system'],PhysicalSystem):
        inp_args['system'] = inp_args['system'].copy()
    #end if

    sim_args['input'] = generate_pwscf_input(input_type,**inp_args)
    pwscf = Pwscf(**sim_args)

    return pwscf
#end def generate_pwscf
