#!/usr/bin/env python3
import sys
import os

# Template-based error handling for custom simulations
try:
    # ==========================================
    # USER CODE STARTS HERE
    # ==========================================
    
    import numpy as np
    output = '${output}'

    print("=== Data Processor Simulation ===")
    print("Processing data from the data generator")
    print("This simulation depends on data from: ", output)
    print("The placeholder, ", output, " has been replaced with the actual dependency path")

    # Check if dependency data exists
    if os.path.exists(output):
        print("Found dependency data at: ", output)
        print("Contents of dependency directory:")
        for item in os.listdir(output):
            print(f"  {item}")
        
        # Process the data from data_generator
        data_dir = os.path.join(output, "data")
        if os.path.exists(data_dir):
            print("Found data directory, processing...")
            
            # Load the data
            x = np.loadtxt(os.path.join(data_dir, "x_values.txt"))
            y = np.loadtxt(os.path.join(data_dir, "y_values.txt"))
            matrix = np.loadtxt(os.path.join(data_dir, "matrix.txt"))
            
            print(f"Loaded x: {x.shape}, y: {y.shape}, matrix: {matrix.shape}")
            
            # Perform some analysis
            print("\nPerforming analysis...")
            
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
                f.write("X Statistics:\n")
                for key, value in x_stats.items():
                    f.write(f"  {key}: {value:.6f}\n")
                
                f.write("\nY Statistics:\n")
                for key, value in y_stats.items():
                    f.write(f"  {key}: {value:.6f}\n")
                
                f.write("\nMatrix Statistics:\n")
                for key, value in matrix_stats.items():
                    f.write(f"  {key}: {value:.6f}\n")
                
                f.write(f"\nCorrelation (x,y): {correlation:.6f}\n")
                f.write(f"Matrix determinant: {matrix_det:.6f}\n")
            
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
        f.write("Data processing completed successfully\n")
    print("Simulation completed successfully")
    
except Exception as e:
    # REQUIRED: Create error file
    error_msg = str(e)
    with open('data_processing_error.txt', 'w') as f:
        f.write(f"Data processing failed: {error_msg}\n")
    print(f"Simulation failed: {error_msg}")
    sys.exit(1)
