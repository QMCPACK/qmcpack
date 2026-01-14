# Nexus Readme

## Contents
  1) About Nexus 
  2) Installation instructions
  3) Summary of library files 


## 1) About Nexus
  
  Nexus was written by Jaron T. Krogel starting in 2012 at UIUC.
  Most of the code was developed at ORNL after Dec. 2012.



## 2) Installation instructions (Linux only)

### a. Ensure that Python and NumPy are available


The vast majority Linux distributions come with Python pre-installed.

To check for Python, type "python3 --version" at the command line. You should see something like "Python 3.6.9".  Nexus is only compatible with Python 3.x (2.x is no longer supported). If Python 3 is not present on your system, you can install it by following the instructions at [python-guide.org](https://docs.python-guide.org/starting/installation/)

With Python present, check for pip.  Type "pip --version" at the command line.  You should see something like:

    pip 9.0.1 from /usr/lib/python3/dist-packages (python 3.6)

For basic (minimum) functioning of Nexus, `numpy` is required. Optional dependencies include `scipy`, `matplotlib`, `h5py`, `pydot`, `spglib`, `pycifrw`, `cif2cell` and `seekpath`. If not already present, these can all be installed via pip:

    pip install --user numpy
    pip install --user scipy
    pip install --user matplotlib
    pip install --user h5py
    pip install --user pydot
    pip install --user spglib
    pip install --user PyCifRW
    pip install --user cif2cell
    pip install --user seekpath

Please see the Nexus user guide (PDF) for more information about installing dependencies and how they are used by Nexus.


### b. "Install" Nexus

All that should be needed to get Nexus working on your system is to add the path to its library files to the PYTHONPATH  environment variable. 
   
If the path to the directory this README file sits in is:
   
     /your/path/to/nexus/
   
Then you can add one line to the .bashrc file in your home directory:
   
    export PYTHONPATH=/your/path/to/nexus/:$PYTHONPATH
   
After adding this line, type `source $HOME/.bashrc`.  Now check that Nexus is available by typing `python` then `import nexus`. You should be able to use Nexus if the import proceeds without incident (i.e. no messages displayed).

The executables packaged with Nexus can be used once they are added to the PATH environment variable.  Do this by adding the following line to your .bashrc file:
   
    export PATH=/your/path/to/nexus/bin/:$PATH

Both of these steps can alternately be performed by the installer packaged with Nexus. To use it, simply type the following at the command line:

    >/your/path/to/nexus/install

Check the end of your `.bashrc` (or equivalent) to make sure the 
`$PATH` and `$PYTHONPATH` variables have been set properly.

If you want the binaries to be copied to a location outside 
the downloaded Nexus distribution, then just provide that as 
a path to the installer:

    >/your/path/to/nexus/install /some/other/location


## 3) Summary of important library files

**The information below can assist in navigating the Nexus program. See headers of individual Python files for more detail.**

### a. Core facilities


- `generic.py`, `developer.py`

    Abstract base classes for all Nexus classes.  Defines abilities such as querying, printing, etc.  Also provides developer interface functions such as error reporting and redirection for unimplemented base class functions.  Void class protects imports and halts execution if items from unavailable modules are encountered.
 
- `project_manager.py`

    Designed for communicating with computing environments, (HPC, workstation) and scheduling jobs.

- `structure.py`

    Contains the Structure class and structure generator functions.
 
- `physical_system.py`

    Contains PhysicalSystem class and particle information.
 
- `simulation.py`

    Contains Simulation, SimulationInput, and SimulationAnalyzer base classes.  All core functionality of Simulation objects is defined here.  Also contains derived simulation classes (sim, input, analyzer) for generic simulation codes.  Enables driving of virtually any single or multi- input file simulation code by using input file templates.
 
- `machine.py`

    Contains Job, Machine, Workstation, and Supercomputer classes. Also contains all derived Supercomputer classes (Titan, Edison, Mira, etc.)
 
- `project_manager.py`

    Contains ProjectManager class.  Code controlling workflow management is generally shared between ProjectManager and Simulation.
 
- `bundle.py`
    
    Contains a facility (the "bundle" function) for bundling many simulation jobs into a single, larger batch job.


### b. Derived classes for specific simulation codes

- QMCPACK
    - Simulation:
        - `qmcpack.py`
    - SimulationInput:
        - `qmcpack_input.py`
    - SimulationAnalyzer:
        - `qmcpack_analyzer.py`

    Supporting analysis classes:
    - `qmcpack_analyzer_base.py`
    - `qmcpack_method_analyzers.py`
    - `qmcpack_property_analyzers.py`
    - `qmcpack_quantity_analyzers.py`
    - `qmcpack_result_analyzers.py`
 
    Orbital conversion tools for QMCPACK:
    - `convert4qmc.py`
    - `pw2qmcpack.py`
    - `wfconvert.py`
 
- PWSCF
    - Simulation:
        - `pwscf.py`
    - SimulationInput:
        - `pwscf_input.py` (read/write any input file)
    - SimulationAnalyzer:
        - `pwscf_analyzer.py` (parse log output and   data~file.xml)

- VASP
    - Simulation:
        - `vasp.py`
    - SimulationInput:
        - `vasp_input.py` (read/write any input file)
    - SimulationAnalyzer:
        - `vasp_analyzer.py` (parse OUTCAR and vasprun.xml)
 
- GAMESS
    - Simulation:
        - `gamess.py`
    - SimulationInput:
        - `gamess_input.py` (read/write any input file)
    - SimulationAnalyzer:
        - `gamess_analyzer.py` (parse log output)


 c. Miscellaneous facilities
 -------------------------------------------------

   debug.py
   ~~~~~~~~
     Adds support for interactive code break points.
 
   periodic_table.py
   ~~~~~~~~~~~~~~~~~
     Access to periodic table data in a structured object format.
 
   numerics.py, extended_numpy.py
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Statistics and curve fitting support.
 
   unit_converter.py
   ~~~~~~~~~~~~~~~~~
     Unit conversion of scalars and arrays.
 
   xmlreader.py
   ~~~~~~~~~~~~
     Class to convert XML file into structured object format.
 
   hdfreader.py
   ~~~~~~~~~~~~
     Class to convert HDF5 file into structured object containing 
     numpy arrays.
 
   fileio.py
   ~~~~~~~~~
     Interface to XSF files.  Other will go here later.
 
   pseudopotential.py
   ~~~~~~~~~~~~~~~~~~
     Classes representing pseudopotentials.  Functions supporting 
     pseudopotential conversion.
 
   template_simulation.py
   ~~~~~~~~~~~~~~~~~~~~~~
     Example and explanation for developers interested in 
     implementing new derived Simulation classes for codes 
     not yet supported by Nexus.


