##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  gamess.py                                                         #
#    Nexus interface to the GAMESS simulation code.                  #
#                                                                    #
#  Content summary:                                                  #
#    GamessInput                                                     #
#      Input class for the GAMESS code.                              #
#      Capable of reading/writing arbitrary GAMESS input files.      #
#                                                                    #
#    generate_gamess_input                                           #
#      User function to create arbitrary GAMESS input.               #
#                                                                    #
#    KeywordGroup                                                    #
#      Represents an arbitrary keyword group in the input file.       #
#                                                                    #
#    KeywordSpecGroup                                                #
#      Base class for specialized keyword groups.                    #
#      Derived classes enforce the keyword specification.            #
#      See ContrlGroup, SystemGroup, GuessGroup, ScfGroup,           #
#        McscfGroup, DftGroup, GugdiaGroup, DrtGroup, CidrtGroup,    #
#        and DetGroup                                                #
#                                                                    #
#    FormattedGroup                                                  #
#      Represents strict machine-formatted input groups.             #
#                                                                    #
#====================================================================#


import os
import numpy as np
from numpy import array,ndarray,abs
from generic import obj
from developer import DevBase
from debug import *
from simulation import Simulation
from gamess_input import GamessInput,generate_gamess_input,FormattedGroup,KeywordGroup,GuessGroup,GIarray
from gamess_analyzer import GamessAnalyzer



class Gamess(Simulation):
    input_type         = GamessInput
    analyzer_type      = GamessAnalyzer
    generic_identifier = 'gamess'
    application        = 'gamess.x' 
    infile_extension   = '.inp'
    application_properties = set(['serial','mpi'])
    application_results    = set(['orbitals'])

    ericfmt = None
    mcppath = None

    @staticmethod
    def settings(ericfmt=None,mcppath=None):
        Gamess.ericfmt = ericfmt
        Gamess.mcppath = mcppath
    #end def settings

    @staticmethod
    def restore_default_settings():
        Gamess.ericfmt = None
        Gamess.mcppath = None
    #end def restore_default_settings


    def __init__(self,**kwargs):
        self.mo_reorder = None
        mo_reorder = kwargs.pop('mo_reorder',None)
        if mo_reorder is not None:
            self.mo_reorder = [s.lower() for s in mo_reorder]
        #end if
        Simulation.__init__(self,**kwargs)
    #end def __init__


    def init_job_extra(self):
        # gamess seems to need lots of environment variables to run properly
        # nearly all of these are names of output/work files
        # setup the environment to run gamess
        if not isinstance(self.ericfmt,str):
            self.error('you must set ericfmt with settings() or Gamess.settings()')
        #end if
        env = obj()
        for file,unit in GamessInput.file_units.items():
            env[file] = '{0}.F{1}'.format(self.identifier,str(unit).zfill(2))
        #end for
        env.INPUT   = self.infile
        env.ERICFMT = self.ericfmt
        env.MCPPATH = self.mcppath
        self.job.set_environment(**env)
    #end def init_job_extra


    def check_result(self,result_name,sim):
        input = self.input 
        if result_name=='orbitals':
            calculating_result = 'contrl' in input and 'scftyp' in input.contrl and input.contrl.scftyp.lower() in ('rhf','rohf','uhf','mcscf','none')
        else:
            calculating_result = False
        #end if
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj()
        input    = self.input
        analyzer = self.load_analyzer_image()
        if result_name=='orbitals':
            result.location  = os.path.join(self.locdir,self.outfile)
            result.outfile   = result.location
            result.vec       = None # vec from punch
            result.norbitals = 0    # orbital count in punch
            result.mos       = 0    # orbital count (MO's) from log file
            result.scftyp    = input.contrl.scftyp.lower()
            if 'counts' in analyzer and 'mos' in analyzer.counts:
                result.mos = analyzer.counts.mos
            #end if
            if 'punch' in analyzer and 'vec' in analyzer.punch:
                result.norbitals = analyzer.punch.norbitals
                result.vec       = analyzer.punch.vec
            #end if
            if 'orbitals' in analyzer:
                result.orbitals = analyzer.orbitals
            #end if
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        input = self.input
        if result_name=='orbitals':
            if result.vec is None or result.norbitals<1:
                self.error('could not obtain orbitals from previous GAMESS run')
            #end if
            if not 'guess' in input:
                input.guess = GuessGroup()
            #end if
            if 'norb' in input.guess: # user provided norb
                norb = input.guess.norb
            else:
                norb = result.norbitals
            #end if
            if 'norder' not in input.guess:
                input.guess.clear()
            #end if
            if self.mo_reorder is not None:
                if 'orbitals' not in result:
                    self.error('Orbital information from prior calculation "{}" located at {} cannot be found. You requested orbital reordering via the  "mo_reorder" input keyword.  Due to missing information, this operation cannot be performed.  The current simulation "{}" is located at {}.'.format(sim.identifier,sim.locdir,self.identifier,self.locdir))
                    self.block()
                #end if
                guess_inputs = obj()
                ecounts = self.system.particles.electron_counts()
                orbs = result.orbitals
                order_map = obj(up='iorder',down='jorder')
                nelec_map = obj(up=ecounts[0],down=ecounts[1])
                for spin,vname in order_map.items():
                    nelec = nelec_map[spin]
                    if len(self.mo_reorder)<nelec:
                        self.error('Too few symmetries provided in "mo_reorder" for spin "{0}".\nNumber of electrons with spin "{0}": {1}\nNumber of entries in "mo_reorder": {2}\nContents of "mo_reorder": {3}\nSimulation identifier: {4}\nSimulation location: {5}'.format(spin,nelec,len(self.mo_reorder),self.mo_reorder,self.identifier,self.locdir))
                    #end if
                    symmetries = [s.lower() for s in orbs[spin].symmetry]
                    missing = set(self.mo_reorder)-set(symmetries)
                    if len(missing)>0:
                        self.error('Symmetries provided by "mo_reorder" keyword are not found in the outputted MOs.\nSet of symmetries provided in "mo_reorder": {}\nSet of symmetries present in MOs: {}\nContents of "mo_reorder": {}\nSimulation identifier: {}\nSimulation location: {}'.format(sorted(set(self.mo_reorder)),sorted(set(symmetries)),self.mo_reorder,self.identifier,self.locdir))
                    #end if
                    occ  = np.zeros(len(symmetries),dtype=bool)
                    for symm in self.mo_reorder[:nelec]:
                        for i,(s,o) in enumerate(zip(symmetries,occ)):
                            if not o and symm==s:
                                occ[i] = True
                                break
                            #end if
                        #end for
                    #end for
                    if occ.sum()<nelec:
                        self.error('Too few orbitals occupied based on "mo_reorder" request.\nNumber of orbitals occupied: {}\nNumber of spin "{}" electrons: {}\nContents of "mo_reorder": {}\nSimulation identifier: {}\nSimulation location: {}'.format(occ.sum(),spin,nelec,self.mo_reorder,self.identifier,self.locdir))
                    #end if
                    indices = np.arange(len(symmetries),dtype=int)[occ]+1
                    start = 0
                    found = False
                    for i in range(len(indices)):
                        start = i
                        if indices[i]!=i+1:
                            found = True
                            break
                        #end if
                    #end if
                    if found:
                        reduced_indices = indices[start:]
                        start+=1
                    else:
                        reduced_indices = []
                    #end if
                    if len(reduced_indices)>0:
                        guess_inputs.norder = 1
                        guess_inputs[vname] = GIarray({start:reduced_indices})
                    #end if
                #end for
                input.guess.set(**guess_inputs)
            #end if
            input.guess.set(
                guess = 'moread',
                norb  = norb,
                prtmo = True,
                )
            input.vec = FormattedGroup(result.vec)
        else:
            self.error('ability to incorporate result '+result_name+' has not been implemented')
        #end if
    #end def incorporate_result


    def app_command(self):
        if self.app_name == 'rungms':
            return 'rungms '+self.infile
        else:
          return self.app_name+' '+self.infile.replace('.inp','')      
        #end if
    #end def app_command


    def check_sim_status(self):
        output = open(os.path.join(self.locdir,self.outfile),'r').read()
        #errors = open(os.path.join(self.locdir,self.errfile),'r').read()
        
        self.failed = 'EXECUTION OF GAMESS TERMINATED -ABNORMALLY-' in output
        self.finished = self.failed or 'EXECUTION OF GAMESS TERMINATED NORMALLY' in output
    #end def check_sim_status


    def get_output_files(self):
        output_files = []
        return output_files
    #end def get_output_files


    def output_filename(self,name):
        name = name.upper()
        if name not in GamessInput.file_units:
            self.error('gamess does not produce a file matching the requested description: {0}'.format(name))
        #end if
        unit = GamessInput.file_units[name]
        filename = '{0}.F{1}'.format(self.identifier,str(unit).zfill(2))
        return filename
    #end def output_filename


    def output_filepath(self,name):
        filename = self.output_filename(name)
        filepath = os.path.join(self.locdir,filename)
        filepath = os.path.abspath(filepath)
        return filepath
    #end def
#end class Gamess



def generate_gamess(**kwargs):
    sim_args,inp_args = Gamess.separate_inputs(kwargs,copy_pseudos=False,sim_kw=['mo_reorder'])

    if not 'input' in sim_args:
        sim_args.input = generate_gamess_input(**inp_args)
    #end if
    gamess = Gamess(**sim_args)

    return gamess
#end def generate_gamess














