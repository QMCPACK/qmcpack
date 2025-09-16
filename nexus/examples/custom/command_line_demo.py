#!/usr/bin/env python3

"""
Simple demonstration of custom simulation dependencies using command line executables.
This shows a basic example: Python data generation followed by bash ls command.

TEMPLATE-BASED ERROR HANDLING APPROACH:
=======================================

This example follows the template-based error handling approach for custom simulations.
Each script follows this pattern:

1. Use appropriate shebang for the script type
2. Wrap all user code with error handling wrappers
3. Create completion file at the end of successful execution

Python template structure:
```python
#!/usr/bin/env python3
import sys
import os

# Completion and Error handling wrappers
def create_completion_file():
    with open('completion_file.txt', 'w') as f:
        f.write('Simulation completed successfully')

def create_error_file(error_msg):
    with open('error_file.txt', 'w') as f:
        f.write(f'Simulation failed: {error_msg}')
    sys.exit(1)

# completion_file.txt and error_file.txt should match the completion_file and error_file 
# parameters in the generate_custom_simulation function

# Try/except block to wrap all user code
try:
    # Your Python code goes here
    
    # Completion handling wrappers
    create_completion_file()
    print("Simulation completed successfully")
    
except Exception as e:
    create_error_file(str(e))
```

Bash template structure:
```bash
#!/bin/bash
set -e  # Exit on any error

# completion_file.txt and error_file.txt should match the completion_file and error_file 
# parameters in the generate_custom_simulation function

# Completion and Error handling wrappers
create_completion_file() {
    echo "Simulation completed successfully" > completion_file.txt
}

create_error_file() {
    echo "Simulation failed: $1" > error_file.txt
}

# Start trapping errors before the command is executed
trap 'create_error_file "Script failed at line $LINENO"' ERR

# Your command line operations go here

# Completion handling wrappers
create_completion_file
echo "Simulation completed successfully"
```

Key requirements:
- File names in the simulation object must match completion_file.txt and error_file.txt 
- Use meaningful error messages and exit codes for trapping errors
- Create completion file only on success
"""

from nexus import settings, job, run_project
from nexus import generate_custom_simulation

# Configure Nexus settings
settings(
    runs       = 'runs/cmd_line',
    results    = '',
    sleep      = 1,
    machine    = 'ws16',
    )

print("Demonstrating simple custom simulation dependencies using command line executables...")
print("This shows: Python data generation -> Bash ls command")
print()

# First simulation: Data generator (Python) - same as simple_numpy_demo
data_generator_script = """#!/usr/bin/env python3
import sys
import os
import numpy as np

# Auto-generated error handling wrapper
def create_completion_file():
    with open('data_generation_complete.txt', 'w') as f:
        f.write('Data generation completed successfully')

def create_error_file(error_msg):
    with open('data_generation_error.txt', 'w') as f:
        f.write(f'Data generation failed: {error_msg}')
    sys.exit(1)

try:
    # ==========================================
    # USER CODE STARTS HERE
    # ==========================================
    
    print("=== Data Generator Simulation ===")
    print("Generating sample data using numpy")
    
    # Create output directory
    os.makedirs('data', exist_ok=True)
    
    # Generate sample data
    np.random.seed(42)  # For reproducible results
    data = np.random.randn(100, 3)  # 100 points, 3 features
    
    # Save data to file
    np.savetxt('data/sample_data.txt', data, header='x y z', comments='')
    
    # Generate some statistics
    stats = {
        'mean': np.mean(data, axis=0),
        'std': np.std(data, axis=0),
        'min': np.min(data, axis=0),
        'max': np.max(data, axis=0)
    }
    
    # Save statistics
    with open('data/statistics.txt', 'w') as f:
        f.write("Data Statistics:\\n")
        f.write(f"Mean: {stats['mean']}\\n")
        f.write(f"Std:  {stats['std']}\\n")
        f.write(f"Min:  {stats['min']}\\n")
        f.write(f"Max:  {stats['max']}\\n")
    
    print(f"Generated data with shape: {data.shape}")
    print(f"Mean values: {stats['mean']}")
    print("Data generation completed successfully!")
    
    # ==========================================
    # USER CODE ENDS HERE
    # ==========================================
    
    # Auto-generated completion
    create_completion_file()
    print("Data generation completed successfully")
    
except Exception as e:
    create_error_file(str(e))
"""

data_generator = generate_custom_simulation(
    identifier = 'data_generator',
    path       = 'demo/data_generator',
    job        = job(serial=True, app='python3'),
    script_content = data_generator_script,
    completion_file = "data_generation_complete.txt",
    error_file = "data_generation_error.txt",
    )

# Second simulation: Directory lister (Bash) - simple ls command
directory_lister_script = """#!/bin/bash
set -e  # Exit on any error

# Auto-generated error handling wrapper
create_completion_file() {
    echo "Directory listing completed successfully" > cmd_line_complete.txt
}

create_error_file() {
    echo "Directory listing failed: $1" > cmd_line_error.txt
}

# Trap errors
trap 'create_error_file "Script failed at line $LINENO"' ERR

# ==========================================
# USER CODE STARTS HERE
# ==========================================

echo "=== Directory Lister Simulation ==="
echo "This simulation lists files in the data generator directory"
echo "It depends on data from: {output}"
echo "The placeholder {output} has been replaced with the actual dependency path"

# Check if dependency data exists
if [ -d "{output}" ]; then
    echo "Found dependency data at: {output}"
    echo "Contents of dependency directory:"
    ls -la "{output}"
    
    # List files in the data subdirectory
    data_dir="{output}/data"
    if [ -d "$data_dir" ]; then
        echo "Found data directory, listing contents:"
        ls -la "$data_dir"
        
        # Show file sizes and types
        echo "File details:"
        file "$data_dir"/*
        
        # Count lines in text files
        echo "Line counts:"
        for file in "$data_dir"/*.txt; do
            if [ -f "$file" ]; then
                echo "$(basename "$file"): $(wc -l < "$file") lines"
            fi
        done
        
        echo "Directory listing completed successfully!"
        
    else
        echo "Error: Data directory not found in dependency"
        exit 1
    fi
else
    echo "Error: Dependency directory not found"
    exit 1
fi

# ==========================================
# USER CODE ENDS HERE
# ==========================================

# Auto-generated completion
create_completion_file
echo "Directory listing completed successfully"
"""

directory_lister = generate_custom_simulation(
    identifier = 'directory_lister',
    path       = 'demo/directory_lister',
    job        = job(serial=True, app='bash'),
    script_content = directory_lister_script,
    completion_file = "cmd_line_complete.txt",
    error_file = "cmd_line_error.txt",
    dependencies = [(data_generator, 'output')],
    )

# Run the project
run_project()
