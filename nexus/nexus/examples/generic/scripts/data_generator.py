#!/usr/bin/env python3
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
            f.write(f"{key}: {value:.6f}\n")

    print("Statistics saved to statistics.txt")
    print("Data generation completed successfully")
    
    # ==========================================
    # USER CODE ENDS HERE
    # ==========================================
    
    # REQUIRED: Create completion file
    with open('data_generation_complete.txt', 'w') as f:
        f.write("Data generation completed successfully\n")
    print("Simulation completed successfully")
    
except Exception as e:
    # REQUIRED: Create error file
    error_msg = str(e)
    with open('data_generation_error.txt', 'w') as f:
        f.write(f"Data generation failed: {error_msg}\n")
    print(f"Simulation failed: {error_msg}")
    sys.exit(1)
