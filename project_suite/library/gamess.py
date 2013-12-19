
import os
from generic import obj
from simulation import Simulation,SimulationInput,SimulationAnalyzer


# PLEASE READ THIS
#   
#   depending on what you want to do with a simulation
#     you will have to implement different functions below
#     here are a few use cases and the functions required
#     details about each function are given within the classes
#
#   use cases
#     1) standalone simulation
#         project suite drives this simulation in isolation of others
#           i.e., one performs parameter scans to drive several independent gamess runs
#         in this setting, a gamess simulation does not provide information to 
#           other simulations (e.g. orbitals to qmcpack) and does not accept 
#           information from prior simulations (e.g. structure from pwscf or gamess)
#      
#         the input file will be read from a template file
#           and modified to obtain the desired inputs
#         one could also provide the input longhand in python 
#           in a form GamessInput understands (this depends on your implementation)
#
#         required functions to be implemented:
#           GamessInput: read_contents, write_contents
#           Gamess:      app_command, check_sim_status
#
#     2) generated standalone simulation
#          as above, but with fully generated input files
#            generate functions provide a short-hand of minimal vars for input
#            structure information for the input is extracted from 
#              a standard PhysicalSystem object
#
#          required functions to be implemented:
#           GamessInput: read_contents, write_contents, incorporate_system
#           Gamess:      app_command, check_sim_status
#           generate_gamess_input
#           
#     3) simulation that provides information to subsequent chained simulations
#          as above (with or without #2)
#            other simulations can request and get information about 
#              results produced by this simulation
#            (e.g. relaxed structure data, location of orbital files, etc.)
#            this information is used by the others to populate input files
#
#         required functions to be implemented:
#           GamessInput: read_contents, write_contents
#           Gamess:      app_command,check_sim_status, 
#                        check_result, get_result
#
#           if required to get needed output information:
#           GamessAnalyzer: analyze
#
#     3) simulation that provides/receives info to/from other simulations
#          as above (with or without #2)
#            this simulation can request info from other sims 
#            info is used to populate own input file
#
#         required functions to be implemented:
#           GamessInput: read_contents, write_contents
#           Gamess:      app_command,check_sim_status, 
#                        check_result, get_result,
#                        incorporate_result
#
#           if required to get needed output information:
#           GamessAnalyzer: analyze
# 



class GamessInput(SimulationInput):
    def __init__(self,filepath=None):
        # optional
        #  below is a convenient default
        #  but it can be changed to anything desired
        if filepath!=None:
            self.read(filepath)
        #end if
    #end def __init__


    def read_contents(self,contents):
        # required
        #  the string 'contents' contains the text of an input file
        #  translate contents into an internal representation of the input
        # for example, an input file looking like:
        #   section_a
        #     var_a  1
        #     var_b  word
        #   section_b
        #     var_c  1e4
        #     var_d  T
        # might be translated into:
        #   sec = obj()
        #   sec['var_a'] = 1
        #   sec['var_b'] = 'word'
        #   self['section_a'] = sec
        #   sec = obj()
        #   sec['var_c'] = 1e4
        #   sec['var_d'] = True
        #   self['section_b'] = sec
        # with some line parsing:
        #   for line in contents.splitlines():
        #      # parse lines here
        #   #end for
        # the resulting input object can then be worked with in a simple way:
        #    >>> input = DemoInput()
        #    >>> input
        #      section_a             obj
        #      section_b             obj
        #    
        #    >>> print input
        #      section_a
        #        var_a           = 1
        #        var_b           = word
        #      end section_a
        #      section_b
        #        var_c           = 10000.0
        #        var_d           = True
        #      end section_b
        #    
        #    >>> input.section_b.var_c = 25.0
        None
    #end def read_contents


    def write_contents(self):
        # required
        #  translate the internal representation of input into a string
        # for the above example, this might look like:
        #   contents = ''
        #   for secname,sec in self.iteritems():
        #       contents += secname+'\n'
        #       for val,val in sec.iteritems():
        #          contents += '  '+var+' = '+val2str(val)
        #       #end for
        #   #end for
        contents = ''
        return contents
    #end def write_contents


    def incorporate_system(self,system):
        # optional
        #  only necessary if you want to populate atomic positions, etc
        #  from a PhysicalSystem object
        # if you don't want to implement it, no action is required
        self.not_implemented()
    #end def incorporate_system
#end class GamessInput



def generate_gamess_input(
    # probably lots of keyword arguments
    # kw1    = default_val1,
    # kw2    = default_val2,
    # system = None,
    # ...
    ):
    # optional
    #  only necessary if you want to make gamess input files
    #  with the fewest relevant variables
    # if you don't want to implement it, uncomment the following line
    #exit()

    gi = GamessInput()

    # use keyword inputs to complete gamess input

    # a common feature is to read information from a PhysicalSystem object
    #  (#electrons, atom positions, etc)
    # if system!=None
    #     gi.incorporate_system(system)
    # #end if

    return gi
#end def generate_gamess_input



class GamessAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0=None):
        # optional
        #  only necessary if you want to use results from output files
        #   to inform the inputs of subsequent simulations
        #   or if you want to have a general purpose class to scrape
        #   and process simulation data
        # below is a reasonable default implementation
        # if you don't want to implement it, just uncomment the following line
        #return

        self.path  = None
        self.input = None
        infile = None
        if isinstance(arg0,Simulation):
            sim = arg0
            infile = os.path.join(sim.locdir,sim.infile)
        else:
            infile = arg0
        #end if
        if infile!=None:
            self.path = os.path.dirname(infile)
            self.input = GamessInput(infile)
        #end if
    #end def __init__


    def analyze(self):
        # optional
        #  only necessary if you want to use results from output files
        #   to inform the inputs of subsequent simulations
        #   or if you want to have a general purpose class to scrape
        #   and process simulation data
        # if you don't want to implement it, no action is required
        None
    #end def analyze
#end class GamessAnalyzer



class Gamess(Simulation):
    input_type         = GamessInput
    analyzer_type      = GamessAnalyzer
    generic_identifier = 'gamess'
    application        = 'gamess_exe' #replace with default name of gamess executable
    application_properties = set(['serial','mpi'])
    application_results    = set(['orbitals']) #what gamess produces that other simulations can use

    def check_result(self,result_name,sim):
        # optional
        #  only necessary if another simulation depends on this one
        #  e.g.
        #  other_sim.depends(gamess_sim,'orbitals')  or similar
        # if you don't want to implement it, uncomment the line below
        #return False
        calculating_result = False
        input = self.input # a GamessInput object
        # check the input to see if result is being calculated
        #  (e.g. result_name='orbitals')
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        # optional
        #  only necessary if another simulation depends on this one
        #  e.g.
        #  other_sim.depends(gamess_sim,'orbitals')  or similar
        # if you don't want to implement it, uncomment the line below
        #self.not_implemented()

        result = obj()
        input    = self.input
        #analyzer = self.load_analyzer_image()

        # package information about a result/product in the result object
        # for example, if orbitals are requested, 
        # the path to the orbital file might be provided:
        # result.orbital_file = '/path/to/orbital/file'
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        # optional
        #  only necessary if this gamess sim depends on another sim
        #  e.g.
        #  gamess_sim.depends(other_sim,'structure')  or similar
        # if you don't want to implement it no action is required
        self.not_implemented()
    #end def incorporate_result


    def app_command(self):
        # required
        #  specify command line arguments to the executable, such as the input file
        #    e.g. command_line_args = ' '+self.infile
        command_line_args = ''
        return self.app_name + command_line_args
    #end def app_command


    def check_sim_status(self):
        # required
        #  read output/error files to check whether simulation has
        #    completed successfully
        #  one could also check whether all output files exist
        output = open(os.path.join(self.locdir,self.outfile),'r').read()
        errors = open(os.path.join(self.locdir,self.errfile),'r').read()
        
        success = False
        # check output and errors
        #  set success=True if run completed successfully

        self.finished = success and self.job.finished
    #end def check_sim_status


    def get_output_files(self):
        # optional
        #  if provided, the listed output files will be copied to the results directory
        # if you don't want to implement it, no action is required
        output_files = []
        return output_files
    #end def get_output_files
#end class Gamess



def generate_gamess(**kwargs):
    # optional
    #  the following code should work provided
    #  generate_gamess_input is suitably defined
    overlapping_kw = set(['system'])
    kw       = set(kwargs.keys())
    sim_kw   = kw & Simulation.allowed_inputs
    inp_kw   = (kw - sim_kw) | (kw & overlapping_kw)    
    sim_args = dict()
    inp_args = dict()
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

    sim_args['input'] = generate_gamess_input(**inp_args)
    gamess = Gamess(**sim_args)

    return gamess
#end def generate_gamess
