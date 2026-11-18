##################################################################
##  (c) Copyright 2020-  by Jaron T. Krogel                     ##
##################################################################


import os

from .developer import obj
from .execute import execute
from .pseudoset import PseudoSet
from .rmg_analyzer import RmgAnalyzer
from .rmg_input import RmgInput, generate_rmg_input
from .simulation import Simulation


class Rmg(Simulation):
    input_type             = RmgInput
    analyzer_type          = RmgAnalyzer
    generic_identifier     = 'rmg'
    application            = 'rmg-cpu'
    application_properties = frozenset({'serial','mpi'})
    application_results    = frozenset({'charge_density','orbitals'})

    default_wave_function_file = 'Waves/wave.out'


    def __init__(self,**sim_args):
        sync_from_scf = sim_args.pop('sync_from_scf',True)
        Simulation.__init__(self,**sim_args)
        self.sync_from_scf = sync_from_scf and self.input.run_mode=='nscf'
    #end def __init__


    def wave_function_file(self,kind):
        """Return the input or output restart file name of this run (relative to locdir)."""
        key = kind+'_wave_function_file'
        if key in self.input:
            return self.input[key]
        #end if
        return self.default_wave_function_file
    #end def wave_function_file


    def check_result(self,result_name,sim):
        input = self.input
        if result_name=='charge_density':
            calculating_result = input.run_mode in (None,'scf')
        elif result_name=='orbitals':
            calculating_result = 'write_qmcpack_restart' in input and bool(input.write_qmcpack_restart)
        else:
            calculating_result = False
        #end if
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj()
        outfile = self.wave_function_file('output')
        if result_name=='charge_density':
            result.locdir   = self.locdir
            result.location = os.path.join(self.locdir,outfile)
            result.outdir   = os.path.dirname(result.location)
            for name in ('wavefunction_grid','processor_grid'):
                result[name] = self.input[name] if name in self.input else None
            #end for
        elif result_name=='orbitals':
            # RMG writes the ESHDF file for QMCPACK next to its restart files
            result.h5file = os.path.join(self.locdir,outfile+'.h5')
        else:
            msg = 'ability to get result '+result_name+' has not been implemented'
            raise NotImplementedError(msg)
        #end if
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        if result_name=='charge_density':
            input = self.input
            if input.run_mode!='nscf':
                msg = (
                    'RMG can only use the charge density of another run in an NSCF calculation.\n'
                    f'Please set calculation_mode="NSCF" for the run in {self.locdir}.'
                    )
                raise ValueError(msg)
            #end if
            # NSCF restarts read the per-processor files of the SCF run, so the
            # real-space grid and its processor decomposition must be the same
            for name in ('wavefunction_grid','processor_grid'):
                scf_val = result[name]
                val = input[name] if name in input else None
                if scf_val is not None and val is not None and tuple(scf_val)!=tuple(val):
                    msg = (
                        f'{name} of the NSCF run does not match the SCF run it restarts from.\n'
                        f'  SCF  ({result.locdir}): {tuple(scf_val)}\n'
                        f'  NSCF ({self.locdir}): {tuple(val)}'
                        )
                    raise ValueError(msg)
                #end if
            #end for
            res_path = os.path.abspath(result.locdir)
            loc_path = os.path.abspath(self.locdir)
            infile = self.wave_function_file('input')
            if res_path==loc_path:
                pass # don't need to do anything if in same directory
            elif self.sync_from_scf: # rsync the scf restart files into the nscf dir
                if os.path.basename(infile)!=os.path.basename(result.location):
                    msg = (
                        'input_wave_function_file of the NSCF run must have the same file name as\n'
                        'output_wave_function_file of the SCF run when sync_from_scf is used.\n'
                        f'  SCF output: {result.location}\n'
                        f'  NSCF input: {infile}'
                        )
                    raise ValueError(msg)
                #end if
                outdir = os.path.dirname(os.path.join(self.locdir,infile))
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                #end if
                # skip the SCF's QMCPACK file so a failed NSCF cannot leave it behind as its own
                command = f"rsync -av --exclude='*.h5' {result.outdir}/ {outdir}/"
                sync_record = os.path.join(outdir,'nexus_sync_record')
                if not os.path.exists(sync_record):
                    print(f'    Running rsync for the {outdir} directory. This might take a while.')
                    execute(command)
                    print(f'    Completed rsync for the {outdir} directory.')
                    with open(sync_record,'w') as f:
                        f.write('\n')
                    #end with
                #end if
            else: # read the scf restart files in place
                input.input_wave_function_file = os.path.relpath(result.location,self.locdir)
            #end if
        else:
            msg = 'ability to incorporate result '+result_name+' has not been implemented'
            raise NotImplementedError(msg)
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
    pseudos = kwargs.get('pseudos',None)
    if pseudos is not None:
        system = kwargs.get('system',None)
        pseudos = PseudoSet.get_pseudos(
            pseudos = pseudos,
            system = system,
            code = 'rmg',
            )
        kwargs['pseudos'] = pseudos
        kwargs['files'] = list(kwargs.get('files',[])) + list(pseudos.values())
    #end if

    sim_args,inp_args = Rmg.separate_inputs(kwargs,sim_kw=['sync_from_scf'])

    if 'input' not in sim_args:
        sim_args.input = generate_rmg_input(**inp_args)
    #end if
    rmg = Rmg(**sim_args)

    return rmg
#end def generate_rmg
