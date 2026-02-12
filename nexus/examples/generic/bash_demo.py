#!/usr/bin/env python3

"""
Bash Demo: Demonstration of GenericSimulation class with bash commands.
This shows how to use GenericSimulation with Python scripts that execute bash commands.

Features demonstrated:
- Python data generation
- Python script using subprocess to execute bash commands
- Template-based dependency handling
- Output file tracking for completion status
"""

import os
from nexus import settings, job, run_project
from nexus import generate_simulation, input_template

# Configure Nexus settings
settings(
    runs       = 'runs/generic_bash',
    results    = '',
    sleep      = 1,
    machine    = 'ws16',
    )

print("Bash Demo: GenericSimulation with bash commands...")
print("This shows: Python data generation -> Python script using subprocess for bash commands")
print()

# First simulation: Data generator (Python) using GenericSimulation
# Read the script from file
script_path = os.path.join('scripts', 'data_generator.py')

# Create data generator using GenericSimulation
dg_outfiles = ["data/"+f+".txt" for f in 'matrix statistics x_values y_values'.split()] + ["data_generation_complete.txt"]
# data_generation_complete.txt file is produced by data_generator.py script
# Tracking a completion file is optional, if not provided, 
#   the simulation will be considered finished after the script execution is complete.
data_generator = generate_simulation(
    identifier = 'data_generator',
    path       = 'data_generator',
    job        = job(serial=True, app='python3'),
    input      = script_path,  # Pass script file path
    outfiles   = dg_outfiles,  # Specify completion files
    )

# Second simulation: Bash command executor with python dependency
lister_script_path = os.path.join('scripts', 'list_directory.sh')

input_dl = input_template(filepath=lister_script_path)
input_dl.assign(output=os.path.abspath(data_generator.locdir))

bash_executor = generate_simulation(
    identifier = 'bash_executor',
    path       = 'bash_executor',
    job        = job(serial=True, app='bash'),
    input      = input_dl,  # Pass script file path
    outfiles   = ['list_directory_complete.txt'],  # Specify completion file
    dependencies = [(data_generator, 'other')],
    )

# Run the project
run_project()
