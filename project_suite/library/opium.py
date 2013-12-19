
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
#           i.e., one performs parameter scans to drive several independent opium runs
#         in this setting, a opium simulation does not provide information to 
#           other simulations (e.g. pseudopotentials to qmcpack) 
#      
#         the input file will be read from a template file
#           and modified to obtain the desired inputs
#         one could also provide the input longhand in python 
#           in a form OpiumInput understands (this depends on your implementation)
#
#         required functions to be implemented:
#           OpiumInput: read_contents, write_contents
#           Opium:      app_command, check_sim_status
#
#     2) generated standalone simulation
#          as above, but with fully generated input files
#            generate functions provide a short-hand of minimal vars for input
#
#          required functions to be implemented:
#           OpiumInput: read_contents, write_contents
#           Opium:      app_command, check_sim_status
#           generate_opium_input
#           
#     3) simulation that provides information to subsequent chained simulations
#          as above (with or without #2)
#            other simulations can request and get information about 
#              results produced by this simulation
#            (e.g. relaxed structure data, location of orbital files, etc.)
#            this information is used by the others to populate input files
#
#         required functions to be implemented:
#           OpiumInput: read_contents, write_contents
#           Opium:      app_command,check_sim_status, 
#                        check_result, get_result
#
#           if required to get needed output information:
#           OpiumAnalyzer: analyze



class OpiumInput(SimulationInput):
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
#end class OpiumInput



def generate_opium_input(
    # probably lots of keyword arguments
    # kw1    = default_val1,
    # kw2    = default_val2,
    # system = None,
    # ...
    ):
    # optional
    #  only necessary if you want to make opium input files
    #  with the fewest relevant variables
    # if you don't want to implement it, uncomment the following line
    #exit()

    gi = OpiumInput()

    # use keyword inputs to complete opium input

    # a common feature is to read information from a PhysicalSystem object
    #  (#electrons, atom positions, etc)
    # if system!=None
    #     gi.incorporate_system(system)
    # #end if

    return gi
#end def generate_opium_input



class OpiumAnalyzer(SimulationAnalyzer):
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
            self.input = OpiumInput(infile)
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
#end class OpiumAnalyzer



class Opium(Simulation):
    input_type         = OpiumInput
    analyzer_type      = OpiumAnalyzer
    generic_identifier = 'opium'
    application        = 'opium_exe' #replace with default name of opium executable
    application_properties = set(['serial','mpi'])
    application_results    = set(['orbitals']) #what opium produces that other simulations can use

    def check_result(self,result_name,sim):
        # optional
        #  only necessary if another simulation depends on this one
        #  e.g.
        #  other_sim.depends(opium_sim,'pseudopotential')  or similar
        # if you don't want to implement it, uncomment the line below
        #return False
        calculating_result = False
        input = self.input # a OpiumInput object
        # check the input to see if result is being calculated
        #  (e.g. result_name='pseudopotential')
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        # optional
        #  only necessary if another simulation depends on this one
        #  e.g.
        #  other_sim.depends(opium_sim,'pseudopotential')  or similar
        # if you don't want to implement it, uncomment the line below
        #self.not_implemented()

        result = obj()
        input    = self.input
        #analyzer = self.load_analyzer_image()

        # package information about a result/product in the result object
        # for example, if pseudopotentials are requested, 
        # the path to the pseudopotential file might be provided:
        # result.pseudopotential_file = '/path/to/pseudopotential/file'
        return result
    #end def get_result


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
#end class Opium



def generate_opium(**kwargs):
    # optional
    #  the following code should work provided
    #  generate_opium_input is suitably defined
    kw = set(kwargs.keys())
    sim_kw = kw & Simulation.allowed_inputs
    inp_kw = (kw - sim_kw)
    sim_args = dict()
    inp_args  = dict()
    for kw in sim_kw:
        sim_args[kw] = kwargs[kw]
    #end for
    for kw in inp_kw:
        inp_args[kw] = kwargs[kw]
    #end for    

    sim_args['input'] = generate_opium_input(**inp_args)
    opium = Opium(**sim_args)

    return opium
#end def generate_opium
