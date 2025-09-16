
#====================================================================#
#  custom_simulation.py                                              #
#    Provides a generic custom simulation class for users to        #
#    develop custom simulation objects that can interact with       #
#    existing Nexus simulations. Supports array manipulation,       #
#    bash scripts, and other custom tasks.                          #
#                                                                    #
#  Content summary:                                                  #
#    CustomSimulationInput                                           #
#      Input class for custom simulations. Handles script content   #
#      and dependency data incorporation.                            #
#                                                                    #
#    CustomSimulationAnalyzer                                        #
#      Analyzer class for custom simulations. Checks completion.    #
#                                                                    #
#    CustomSimulation                                                #
#      Main simulation class for custom tasks.                      #
#                                                                    #
#    generate_custom_simulation                                      #
#      User-facing function to create custom simulation objects.    #
#                                                                    #
#====================================================================#


import os
from simulation import Simulation, SimulationInput, SimulationAnalyzer
from generic import obj


class CustomSimulationInput(SimulationInput):
    def __init__(self, script_content=None, completion_file=None, error_file=None, **kwargs):
        super().__init__()
        self.script_content = script_content or ""
        self.parameters = kwargs
        self.dependency_data = {}  # Will store lists of results for each result name
        self.completion_file = completion_file or "output.txt"
        self.error_file = error_file or "error.txt"

    def write_text(self, filepath=None):
        # Replace placeholders in script with actual dependency data
        script = self.script_content
        # Replace any {result_name} placeholders with directory paths
        for result_name, results in self.dependency_data.items():
            if isinstance(results, list):
                # Multiple dependencies with same result name
                for i, result in enumerate(results):
                    placeholder = f'{{{result_name}_{i}}}' if i > 0 else f'{{{result_name}}}'
                    
                    if hasattr(result, 'outdir'):
                        script = script.replace(placeholder, os.path.abspath(result.outdir))
                    elif hasattr(result, 'locdir'):
                        script = script.replace(placeholder, os.path.abspath(result.locdir))
                    else:
                        script = script.replace(placeholder, os.path.abspath(os.path.dirname(result.location)))
            else:
                # Single dependency
                if hasattr(results, 'outdir'):
                    script = script.replace(f'{{{result_name}}}', os.path.abspath(results.outdir))
                elif hasattr(results, 'locdir'):
                    script = script.replace(f'{{{result_name}}}', os.path.abspath(results.locdir))
                else:
                    script = script.replace(f'{{{result_name}}}', os.path.abspath(os.path.dirname(results.location)))
        
        return script
    
    def incorporate_result(self, result_name, result, sim):  # pylint: disable=unused-argument
        # Store dependency data for use in script
        if result_name not in self.dependency_data:
            self.dependency_data[result_name] = result
        else:
            # Convert to list if we have multiple dependencies with same result name
            if not isinstance(self.dependency_data[result_name], list):
                self.dependency_data[result_name] = [self.dependency_data[result_name]]
            self.dependency_data[result_name].append(result)


class CustomSimulationAnalyzer(SimulationAnalyzer):
    def __init__(self, sim):
        self.sim = sim
        self.results = {}
    
    def analyze(self):
        completion_file = self.sim.input.completion_file
        error_file = self.sim.input.error_file
        self.results['completed'] = os.path.exists(os.path.join(self.sim.locdir, completion_file))
        self.results['error'] = os.path.exists(os.path.join(self.sim.locdir, error_file))
    
    def save(self, filepath):
        # Override save to avoid pickling issues
        import pickle
        with open(filepath, 'wb') as f:
            pickle.dump(self.results, f)
    
    def load(self, filepath):
        # Override load to match save
        import pickle
        with open(filepath, 'rb') as f:
            self.results = pickle.load(f)


class CustomSimulation(Simulation):
    input_type = CustomSimulationInput
    analyzer_type = CustomSimulationAnalyzer
    generic_identifier = 'custom'
    application = 'custom'
    application_properties = set(['serial'])
    application_results = set(['output'])
    
    
    def app_command(self):
        """
        Generate the command to execute the custom simulation script.
        
        This method provides flexibility for different script types by reading the 
        executable from the job object. The job object can specify different 
        executables (bash, python3, perl, etc.) via the 'app' property in Job object.
        After the job object is initialized, the app_name property is set to the executable.
        
        How it works:
        1. Check if the job has an 'app_name' property set (e.g., 'python3', 'bash', 'perl')
        2. If yes, use that executable + the input file name (e.g., "python3 script.py")
        3. If no, fall back to bash as the default (e.g., "bash script.sh")
        
        Examples:
        - job(serial=True, app='python3') -> "python3 script.py"
        - job(serial=True, app='perl') -> "perl script.pl" 
        - job(serial=True) -> "bash script.sh" (default)
        
        The input file (self.infile) contains the actual script content that was
        written by the CustomSimulationInput.write_text() method.
        """
        # Check if job specifies a custom executable via app_name
        if hasattr(self.job, 'app_name') and self.job.app_name is not None:
            # Use the specified executable + input file
            return f"{self.job.app_name} {self.infile}"
        else:
            # Fall back to bash as default for shell scripts
            return f"bash {self.infile}"
    
    def check_sim_status(self):
        completion_file = self.input.completion_file
        completion_file_exists = os.path.exists(os.path.join(self.locdir, completion_file))

        error_file = self.input.error_file
        error_file_exists = os.path.exists(os.path.join(self.locdir, error_file))

        # Check if simulation has failed
        self.failed = error_file_exists
        
        # If error file exists, mark job as finished (but failed)
        if error_file_exists and not self.job.finished:
            self.job.finished = True
            self.job.status = self.job.states.finished
            print(f"Simulation {self.identifier} failed - error file found: {error_file}")
        
        # If completion file exists, mark job as finished (success)
        elif completion_file_exists and not self.job.finished:
            self.job.finished = True
            self.job.status = self.job.states.finished
            print(f"Simulation {self.identifier} completed successfully")
        
        # Simulation is finished if either completion or error file exists
        self.finished = self.job.finished and (completion_file_exists or error_file_exists)
    
    def check_result(self, result_name, sim):
        # Check if the simulation is calculating the requested result
        if result_name in self.application_results:
            return True  # Custom simulations always calculate their declared results
        return False
    
    def get_result(self, result_name, sim):
        # Return the simulation directory as the result
        if result_name in self.application_results:
            result = obj()
            result.location = self.locdir
            result.outdir = self.locdir  # Use locdir as outdir for custom simulations
            result.locdir = self.locdir
            return result
        return None
    
    def incorporate_result(self, result_name, result, sim):
        # Store dependency results for use in placeholder replacement
        if not hasattr(self.input, 'dependency_data'):
            self.input.dependency_data = {}
        
        if result_name not in self.input.dependency_data:
            self.input.dependency_data[result_name] = result
        else:
            # Convert to list if we have multiple dependencies with same result name
            if not isinstance(self.input.dependency_data[result_name], list):
                self.input.dependency_data[result_name] = [self.input.dependency_data[result_name]]
            self.input.dependency_data[result_name].append(result)
    
    def get_output_files(self):
        return [self.input.completion_file]


def generate_custom_simulation(
    identifier='custom',
    path='./',
    job=None,
    script_content="",
    completion_file="output.txt",
    error_file="error.txt",
    dependencies=None,
    **kwargs
):
    """
    Generate a custom simulation that can run any type of script.
    
    Parameters:
    -----------
    identifier : str
        Unique name for the simulation
    path : str  
        Directory path where simulation will run
    job : Job object
        Job configuration. Use job(serial=True, app_name='executable') to specify
        the script interpreter (e.g., 'python3', 'bash', 'perl', 'matlab', etc.)
    script_content : str
        The actual script content to execute (can be bash, python, perl, etc.)
    completion_file : str
        File that indicates successfulsimulation completion (default: "output.txt")
    error_file : str
        File that indicates unsuccessful simulation completion (default: "error.txt")
    dependencies : list
        List of (simulation, result_name) tuples for dependency chaining
        
    Examples:
    ---------
    # Python script
    job(serial=True, app='python3')
    
    # Bash script  
    job(serial=True, app='bash')
    
    # Perl script
    job(serial=True, app='perl')
    
    # Default (bash)
    job(serial=True)
    """
    sim_args, inp_args = Simulation.separate_inputs(kwargs)
    
    inp_args['script_content'] = script_content
    inp_args['completion_file'] = completion_file
    inp_args['error_file'] = error_file
    sim_args['input'] = CustomSimulationInput(**inp_args)
    
    sim_args['identifier'] = identifier
    sim_args['path'] = path
    
    sim_args['job'] = job
    
    sim = CustomSimulation(**sim_args)
    
    if dependencies:
        for dep_sim, result_name in dependencies:
            sim.depends(dep_sim, result_name)
    
    return sim
