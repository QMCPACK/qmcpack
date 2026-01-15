# Nexus Readme

## Contents
  1) About Nexus 
  2) Installation instructions
  3) Summary of library files 


## 1) About Nexus
  
  Nexus was written by Jaron T. Krogel starting in 2012 at UIUC.
  Most of the code was developed at ORNL after Dec. 2012.



## 2) Installation instructions (Linux only)

### a. Ensure that Python is available


The vast majority Linux distributions come with Python pre-installed.

To check for Python, type `python3 --version` at the command line. You should see something like "Python 3.6.9".  Nexus is only compatible with Python >3.6 (2.x is no longer supported). If Python 3 is not present on your system, you can install it by following the instructions at [python-guide.org](https://docs.python-guide.org/starting/installation/)

With Python present, check for pip.  Type `pip --version` at the command line (older versions of Python may require using `python3 -m pip --version`).  You should see something like:
```
pip 9.0.1 from /usr/lib/python3/dist-packages (python 3.6)
```
For legacy installation, this is fine, however for direct installation via a package manager such as `uv` or `pip`, you will either need to [install uv](https://docs.astral.sh/uv/getting-started/installation/) or, provided you have `pip` you should run the command
```bash
> python3 -m pip install -U pip
```
which will upgrade your `pip` installation to the latest version.

### **Nexus requires `numpy>=1.13.0` to function. If you are installing via the legacy install instructions, you must download `numpy` manually through your package manager.**

### b. Nexus Install Instructions

The latest versions of Nexus include packaging as a Python project that enable the installation of Nexus directly with a package manager such as `pip` or `uv`. For installation with `pip` you will need at least version 10.0.0. 

---
#### `pip` installation (one line)

```bash
> python3 -m pip install "nexus@git+https://github.com/QMCPACK/qmcpack.git@develop#subdirectory=nexus"
```
---
#### `uv` installation (with venv)

```bash
> uv venv .venv
> uv pip install "nexus@git+https://github.com/QMCPACK/qmcpack.git@develop#subdirectory=nexus"
```
---

Installing via a package manager will automatically install the required dependency (NumPy), however it is strongly encouraged to install the optional dependencies as well. The current optional dependencies are in two groups, one is for quality of life improvements and near-complete access to Nexus's capabilities, which can be installed by adding `[qol]` to the install command, and another group is for full functionality of Nexus, which is denoted `[full]`. These can either be added individually or combined by adding `[qol,full]` in `nexus@git`, which would look like
```bash
> <python3 -m | uv> pip install "nexus[qol,full]@git+https://github.com/QMCPACK/qmcpack.git@develop#subdirectory=nexus"
```

The `[qol]` dependencies are
- `scipy`
- `matplotlib`
- `h5py`
- `spglib`

The `[full]` dependencies are
- `pydot`
    - This requires GraphViz to be installed, which can be done through your Linux distribution's package manager
        - For Debian/Ubuntu, `> apt-get install graphviz`
        - For Fedora, `> dnf install graphviz`
- `cif2cell`
    - `cif2cell` requires `pycifrw` to function, however `pip` and `uv` will install `pycifrw` by default when you install `cif2cell`
- `seekpath`

If you would like to install these packages manually, this can be done via `pip` or `uv`:

    <python3 -m | uv> pip install --user scipy
    <python3 -m | uv> pip install --user matplotlib
    <python3 -m | uv> pip install --user h5py
    <python3 -m | uv> pip install --user pydot
    <python3 -m | uv> pip install --user spglib
    <python3 -m | uv> pip install --user cif2cell
    <python3 -m | uv> pip install --user seekpath

Please see the Nexus user guide (PDF) for more information about installing dependencies and how they are used by Nexus.

### c. Legacy Installation Instructions

All that should be needed to get Nexus working on your system is to add the path to its library files to the PYTHONPATH  environment variable. 
   
If the path to the directory this README file sits in is:
```bash
/your/path/to/nexus/
```
Then you can add one line to the .bashrc file in your home directory:
```bash
export PYTHONPATH=/your/path/to/nexus/:$PYTHONPATH
```
After adding this line, type `source $HOME/.bashrc`.  Now check that Nexus is available by typing `python` then `import nexus`. You should be able to use Nexus if the import proceeds without incident (i.e. no messages displayed).

The executables packaged with Nexus can be used once they are added to the PATH environment variable.  Do this by adding the following line to your .bashrc file:
```bash
export PATH=/your/path/to/nexus/bin/:$PATH
```
Both of these steps can alternately be performed by the installer packaged with Nexus. To use it, simply type the following at the command line:
```bash
> /your/path/to/nexus/install
```
Check the end of your `.bashrc` (or equivalent) to make sure the `$PATH` and `$PYTHONPATH` variables have been set properly.

If you want the binaries to be copied to a location outside the downloaded Nexus distribution, then just provide that as a path to the installer:
```bash
> /your/path/to/nexus/install /some/other/location
```

## 3) Summary of important library files

**The information below can assist in navigating the Nexus program. See headers of individual Python files for more detail.**

### a. Core facilities


- `generic.py`, `developer.py`

    Abstract base classes for all Nexus classes.  Defines abilities such as querying, printing, etc.  Also provides developer interface functions such as error reporting and redirection for unimplemented base class functions.  Void class protects imports and halts execution if items from unavailable modules are encountered.
 
- `project_manager.py`

    Designed for communicating with computing environments, (HPC, workstation) and scheduling jobs.

- `structure.py`

    Contains the `Structure` class and structure generator functions.
 
- `physical_system.py`

    Contains `PhysicalSystem` class and particle information.
 
- `simulation.py`

    Contains `Simulation`, `SimulationInput`, and `SimulationAnalyzer` base classes.  All core functionality of `Simulation` objects is defined here.  Also contains derived simulation classes (sim, input, analyzer) for generic simulation codes.  Enables driving of virtually any single or multi- input file simulation code by using input file templates.
 
- `machine.py`

    Contains `Job`, `Machine`, `Workstation`, and `Supercomputer` classes. Also contains all derived `Supercomputer` classes (`Titan`, `Edison`, `Mira`, etc.)
 
- `project_manager.py`

    Contains `ProjectManager` class.  Code controlling workflow management is generally shared between `ProjectManager` and `Simulation`.
 
- `bundle.py`
    
    Contains a facility (the `bundle` function) for bundling many simulation jobs into a single, larger batch job.


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
        - `pwscf_input.py`
            - read/write any input file
    - SimulationAnalyzer:
        - `pwscf_analyzer.py`
            - parse log output and   data~file.xml

- VASP
    - Simulation:
        - `vasp.py`
    - SimulationInput:
        - `vasp_input.py`
            - read/write any input file
    - SimulationAnalyzer:
        - `vasp_analyzer.py`
            - parse OUTCAR and vasprun.xml
 
- GAMESS
    - Simulation:
        - `gamess.py`
    - SimulationInput:
        - `gamess_input.py`
            - read/write any input file
    - SimulationAnalyzer:
        - `gamess_analyzer.py`
            - parse log output


### c. Miscellaneous facilities

- `debug.py`
    - Adds support for interactive code break points.
 
- `periodic_table.py`
    - Access to periodic table data in a structured object format.
 
- `numerics.py`
    - Statistics and curve fitting support.
 
- `unit_converter.py`
    - Unit conversion of scalars and arrays.
 
- `xmlreader.py`
    - Class to convert XML file into structured object format.
 
- `hdfreader.py`
    - Class to convert HDF5 file into structured object containing numpy arrays.
 
- `fileio.py`
    - Interface to XSF files.  Other will go here later.
 
- `pseudopotential.py`
    - Classes representing pseudopotentials.
    - Functions supporting pseudopotential conversion.
