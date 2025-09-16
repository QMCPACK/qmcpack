#!/usr/bin/env python3

"""
Simple demonstration of custom simulation dependencies using numpy arrays.
This shows how dependencies work with real data manipulation.

TEMPLATE-BASED ERROR HANDLING APPROACH:
=======================================

This example follows the template-based error handling approach for custom simulations.
Each python script follows this pattern:

1. Import sys and os at the top
2. Wrap all user code in a try/except block
3. Create completion file at the end of successful execution
4. Create error file in the except block
5. Use sys.exit(1) for failures

Template structure:
```python
#!/usr/bin/env python3
import sys
import os

try:    
    # Your simulation code goes here
    
    # REQUIRED: Create completion file
    with open('completion_file.txt', 'w') as f:
        f.write("Simulation completed successfully\\n")
    print("Simulation completed successfully")
    
except Exception as e:
    # REQUIRED: Create error file
    error_msg = str(e)
    with open('error_file.txt', 'w') as f:
        f.write(f"Simulation failed: {error_msg}\\n")
    print(f"Simulation failed: {error_msg}")
    sys.exit(1)
```

Key requirements:
- File names must match completion_file and error_file parameters
- Use meaningful error messages
- Always call sys.exit(1) on failure
- Create completion file only on success
"""

from nexus import settings, job, run_project
from nexus import generate_custom_simulation

# Configure Nexus settings
settings(
    results    = '',
    runs       = 'runs/numpy',
    sleep      = 1,
    machine    = 'ws16',
    )

print("Demonstrating simple custom simulation dependencies with numpy arrays...")
print("This will show how dependencies work with real data manipulation")
print()

# First simulation: Data generator - creates simple numpy arrays
# Following the template-based error handling approach
data_generator_script = """#!/usr/bin/env python3
import sys
import os

# Template-based error handling for custom simulations
try:
    # ==========================================
    # USER CODE STARTS HERE
    # ==========================================
    
    import numpy as np

    print("=== Data Generator Simulation ===")
    print("Creating simple numpy arrays for other simulations to use")

    # Create output directory
    os.makedirs('data', exist_ok=True)

    # Generate simple test data
    print("Generating test data...")

    # 1. Simple 1D array
    x = np.linspace(0, 10, 100)
    y = np.sin(x) + 0.1 * np.random.random(100)
    np.savetxt('data/x_values.txt', x)
    np.savetxt('data/y_values.txt', y)
    print(f"Generated x_values.txt with {len(x)} points")
    print(f"Generated y_values.txt with {len(y)} points")

    # 2. Simple 2D array
    matrix = np.random.random((10, 10))
    np.savetxt('data/matrix.txt', matrix)
    print(f"Generated matrix.txt with shape {matrix.shape}")

    # 3. Simple statistics
    stats = {
        'x_min': x.min(),
        'x_max': x.max(),
        'y_min': y.min(),
        'y_max': y.max(),
        'matrix_mean': matrix.mean(),
        'matrix_std': matrix.std()
    }

    # Save statistics
    with open('data/statistics.txt', 'w') as f:
        for key, value in stats.items():
            f.write(f"{key}: {value:.6f}\\n")

    print("Statistics saved to statistics.txt")
    print("Data generation completed successfully")
    
    # ==========================================
    # USER CODE ENDS HERE
    # ==========================================
    
    # REQUIRED: Create completion file
    with open('data_generation_complete.txt', 'w') as f:
        f.write("Data generation completed successfully\\n")
    print("Simulation completed successfully")
    
except Exception as e:
    # REQUIRED: Create error file
    error_msg = str(e)
    with open('data_generation_error.txt', 'w') as f:
        f.write(f"Data generation failed: {error_msg}\\n")
    print(f"Simulation failed: {error_msg}")
    sys.exit(1)
"""

data_generator = generate_custom_simulation(
    identifier = 'data_generator',
    path       = 'demo/data_generator',
    job        = job(serial=True, app='python3'),
    script_content = data_generator_script,
    completion_file = "data_generation_complete.txt",
    error_file = "data_generation_error.txt",
    )

# Second simulation: Data processor (depends on data_generator)
# Following the template-based error handling approach
data_processor_script = """#!/usr/bin/env python3
import sys
import os

# Template-based error handling for custom simulations
try:
    # ==========================================
    # USER CODE STARTS HERE
    # ==========================================
    
    import numpy as np

    print("=== Data Processor Simulation ===")
    print("Processing data from the data generator")
    print("This simulation depends on data from: {output}")
    print("The placeholder {output} has been replaced with the actual dependency path")

    # Check if dependency data exists
    if os.path.exists("{output}"):
        print("Found dependency data at: {output}")
        print("Contents of dependency directory:")
        for item in os.listdir("{output}"):
            print(f"  {item}")
        
        # Process the data from data_generator
        data_dir = os.path.join("{output}", "data")
        if os.path.exists(data_dir):
            print("Found data directory, processing...")
            
            # Load the data
            x = np.loadtxt(os.path.join(data_dir, "x_values.txt"))
            y = np.loadtxt(os.path.join(data_dir, "y_values.txt"))
            matrix = np.loadtxt(os.path.join(data_dir, "matrix.txt"))
            
            print(f"Loaded x: {x.shape}, y: {y.shape}, matrix: {matrix.shape}")
            
            # Perform some analysis
            print("\\nPerforming analysis...")
            
            # 1. Basic statistics
            x_stats = {
                'mean': np.mean(x),
                'std': np.std(x),
                'min': np.min(x),
                'max': np.max(x)
            }
            
            y_stats = {
                'mean': np.mean(y),
                'std': np.std(y),
                'min': np.min(y),
                'max': np.max(y)
            }
            
            matrix_stats = {
                'mean': np.mean(matrix),
                'std': np.std(matrix),
                'min': np.min(matrix),
                'max': np.max(matrix)
            }
            
            # 2. Simple calculations
            correlation = np.corrcoef(x, y)[0, 1]
            matrix_det = np.linalg.det(matrix)
            
            # 3. Save results
            print("Saving processed results...")
            
            # Save statistics
            with open('processed_statistics.txt', 'w') as f:
                f.write("X Statistics:\\n")
                for key, value in x_stats.items():
                    f.write(f"  {key}: {value:.6f}\\n")
                
                f.write("\\nY Statistics:\\n")
                for key, value in y_stats.items():
                    f.write(f"  {key}: {value:.6f}\\n")
                
                f.write("\\nMatrix Statistics:\\n")
                for key, value in matrix_stats.items():
                    f.write(f"  {key}: {value:.6f}\\n")
                
                f.write(f"\\nCorrelation (x,y): {correlation:.6f}\\n")
                f.write(f"Matrix determinant: {matrix_det:.6f}\\n")
            
            # Save processed arrays
            np.savetxt('processed_x.txt', x)
            np.savetxt('processed_y.txt', y)
            np.savetxt('processed_matrix.txt', matrix)
            
            print("Data processing completed successfully!")
            
        else:
            raise FileNotFoundError("Data directory not found in dependency")
    else:
        raise FileNotFoundError("Dependency directory not found")
    
    # ==========================================
    # USER CODE ENDS HERE
    # ==========================================
    
    # REQUIRED: Create completion file
    with open('data_processing_complete.txt', 'w') as f:
        f.write("Data processing completed successfully\\n")
    print("Simulation completed successfully")
    
except Exception as e:
    # REQUIRED: Create error file
    error_msg = str(e)
    with open('data_processing_error.txt', 'w') as f:
        f.write(f"Data processing failed: {error_msg}\\n")
    print(f"Simulation failed: {error_msg}")
    sys.exit(1)
"""

data_processor = generate_custom_simulation(
    identifier = 'data_processor',
    path       = 'demo/data_processor',
    job        = job(serial=True, app='python3'),
    script_content = data_processor_script,
    completion_file = "data_processing_complete.txt",
    
    dependencies = [(data_generator, 'output')],
    )

# Run the project
run_project()
