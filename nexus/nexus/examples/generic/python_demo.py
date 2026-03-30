#!/usr/bin/env python3

"""
Python Demo: Demonstration of GenericSimulation class with numpy data processing.
This shows how to use GenericSimulation with Python scripts for data manipulation.

"""

import os
from nexus import settings, job, run_project
from nexus import generate_simulation, input_template

# Configure Nexus settings
settings(
    results    = '',
    runs       = 'runs/generic_numpy',
    sleep      = 1,
    machine    = 'ws16',
    )

print("Python Demo: GenericSimulation with numpy data processing...")
print("This will show how dependencies work with real data manipulation")
print()

# First simulation: Data generator - creates simple numpy arrays using GenericSimulation
# Read the script from file
script_path = os.path.join('scripts', 'data_generator.py')

# Create data generator using GenericSimulation
dg_outfiles = ["data/"+f+".txt" for f in 'matrix statistics x_values y_values'.split()] + ["data_generation_complete.txt"]
data_generator = generate_simulation(
    identifier = 'data_generator',
    path       = 'data_generator',
    job        = job(serial=True, app='python3'),
    input      = script_path,  # Pass script file path
    outfiles   = dg_outfiles,  # Specify completion files
    )

processor_script_path = os.path.join('scripts', 'data_processor.py')

input_dp = input_template(filepath=processor_script_path)
input_dp.assign(output=os.path.abspath(data_generator.locdir))
# Create data processor using GenericSimulation
data_processor = generate_simulation(
    identifier = 'data_processor',
    path       = 'data_processor',
    job        = job(serial=True, app='python3'),
    input      = input_dp,  # Pass script file path
    outfiles   = ['data_processing_complete.txt'],  # Specify completion file
    dependencies = [(data_generator, 'other')],
    )

# Run the project
run_project()
