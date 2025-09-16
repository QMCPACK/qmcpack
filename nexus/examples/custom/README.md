# Custom Simulations in Nexus

This directory contains examples and documentation for creating custom simulation objects in the Nexus framework. Custom simulations allow you to run any type of script (Python, Bash, Perl, etc.) and integrate them into the Nexus dependency system.

## Overview

Custom simulations provide a flexible way to:
- Run arbitrary scripts and executables
- Chain simulations with dependencies
- Process data from other simulations
- Integrate with existing Nexus workflows

## Quick Start

```python
from nexus import generate_custom_simulation, job, run_project

# Create a simple Python simulation
python_sim = generate_custom_simulation(
    identifier='my_python_sim',
    path='./my_sim',
    job=job(serial=True, app='python3'),
    script_content='''
#!/usr/bin/env python3
import sys
import os

try:
    # Your Python code starts here
    print("Hello from custom simulation!")
    # Your Python code ends here
    
    # REQUIRED: Create completion file
    with open('completion_file.txt', 'w') as f:
        f.write("Simulation completed successfully\\n")
    print("Simulation completed successfully")
    
except Exception as e:
    # REQUIRED: Create error file
    with open('error_file.txt', 'w') as f:
        f.write(f"Simulation failed: {str(e)}\\n")
    print(f"Simulation failed: {str(e)}")
    sys.exit(1)
''',
    completion_file="completion_file.txt",
    error_file="error_file.txt"
)

# Create a bash simulation that depends on the Python one
bash_sim = generate_custom_simulation(
    identifier='my_bash_sim',
    path='./my_bash_sim',
    job=job(serial=True, app='bash'),
    script_content='''
#!/bin/bash
set -e

create_completion_file() {
    echo "Simulation completed successfully" > completion_file.txt
}

create_error_file() {
    echo "Simulation failed: $1" > error_file.txt
}

trap 'create_error_file "Script failed at line $LINENO"' ERR

# Your bash code starts here
echo "Hello from bash simulation!"
echo "Dependency data from: {output}"
# Your bash code ends here

create_completion_file
echo "Simulation completed successfully"
''',
    completion_file="completion_file.txt",
    error_file="error_file.txt",
    dependencies=[(python_sim, 'output')]
)

# Run the project
run_project()
```

## How Dependencies Work

The dependency system uses placeholder replacement to inject actual file paths into your scripts. This allows simulations to access data from other simulations automatically.

### Single Dependency

When you have one dependency:

```python
bash_sim = generate_custom_simulation(
    dependencies=[(python_sim, 'output')]
)
```

**What happens:**
1. `python_sim` produces an `'output'` result (its directory path)
2. In `bash_sim`'s script, `{output}` gets replaced with the absolute path to `python_sim`'s directory
3. The replacement happens automatically when the script is written

**Example:**
```bash
# In your bash script
echo "Data from dependency: {output}"
# Gets replaced with:
echo "Data from dependency: /absolute/path/to/python_sim/directory"
```

### Multiple Dependencies with Same Result Name

When multiple simulations provide the same result name, the system creates numbered placeholders:

```python
# Two simulations both provide 'output'
sim1 = generate_custom_simulation(identifier='sim1', ...)
sim2 = generate_custom_simulation(identifier='sim2', ...)

# A third simulation depends on both
bash_sim = generate_custom_simulation(
    dependencies=[
        (sim1, 'output'),  # First 'output'
        (sim2, 'output')   # Second 'output'
    ]
)
```

**Placeholder naming:**
- `{output}` → First dependency (sim1's directory)
- `{output_1}` → Second dependency (sim2's directory)  
- `{output_2}` → Third dependency (if there was one)
- And so on...

Starting with python3.7 ordering in lists are always preserved, therefore this will work in the same order in the dependencies object. 
However, using same result names is not particularly useful in practice since they are only used as placeholders, therefore a unique name could be used for each dependency simulation instead without loss of functionality. 

**Example script:**
```bash
#!/bin/bash
echo "First dependency: {output}"      # /path/to/sim1
echo "Second dependency: {output_1}"   # /path/to/sim2

# Process data from both
ls -la "{output}/data"
ls -la "{output_1}/data"

# Compare results
diff "{output}/results.txt" "{output_1}/results.txt"
```

### Different Result Names

You can use different result names for different types of data:

```python
# Simulations providing different result types
data_sim = generate_custom_simulation(identifier='data', ...)
config_sim = generate_custom_simulation(identifier='config', ...)

# Analysis simulation depends on both
analysis_sim = generate_custom_simulation(
    dependencies=[
        (data_sim, 'data'),      # Provides 'data' result
        (config_sim, 'config')   # Provides 'config' result
    ]
)
```

**Script with different result types:**
```bash
#!/bin/bash
echo "Data directory: {data}"        # /path/to/data_sim
echo "Config directory: {config}"    # /path/to/config_sim

# Load configuration
source "{config}/settings.conf"

# Process data
python process_data.py --input "{data}/raw_data.txt" --config "{config}/params.json"
```
