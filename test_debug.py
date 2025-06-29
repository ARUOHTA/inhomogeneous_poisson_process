#!/usr/bin/env python3
"""
Debug script to test NW estimation step by step
"""

import sys
import os
import gc
import psutil

# Add path
sys.path.append('/home/ohta/dev/bayesian_statistics')

from bayesian_statistics.model3_config import Model3Config, Model3Pipeline

def get_memory_usage():
    """Get current memory usage"""
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    return memory_info.rss / 1024 / 1024  # MB

def test_step_by_step():
    """Test NW estimation step by step"""
    
    print(f"Initial memory: {get_memory_usage():.1f} MB")
    
    # Create config
    config = Model3Config(
        data_dir="/home/ohta/dev/bayesian_statistics/data/",
        x_min=138,
        x_max=141, 
        y_min=34,
        y_max=37,
        nw_sigma=500,
        nw_sigma_for_sites=0.1,
    )
    
    # Initialize pipeline
    pipeline = Model3Pipeline(config)
    
    # Run preprocessing
    print("Starting preprocessing...")
    preprocessor = pipeline.run_preprocessing()
    print(f"After preprocessing: {get_memory_usage():.1f} MB")
    
    # Test each step manually
    from bayesian_statistics.model3_nadaraya_watson import NadarayaWatsonEstimator
    
    nw_estimator = NadarayaWatsonEstimator(
        sigma=config.nw_sigma,
        sigma_for_sites=config.nw_sigma_for_sites,
    )
    
    print("Testing step 1: create_explanatory_variables")
    W_grids, W_sites = preprocessor.create_explanatory_variables(config.nw_variable_names)
    print(f"W_grids shape: {W_grids.shape}, dtype: {W_grids.dtype}")
    print(f"W_sites shape: {W_sites.shape}, dtype: {W_sites.dtype}")
    print(f"After explanatory variables: {get_memory_usage():.1f} MB")
    
    print("\nTesting step 2: load_tobler_distances")
    distances = preprocessor.load_tobler_distances()
    print(f"distances shape: {distances.shape}, dtype: {distances.dtype}")
    print(f"After distances: {get_memory_usage():.1f} MB")
    
    print("\nTesting step 3: calculate_kernel_weights")
    weights = nw_estimator.calculate_kernel_weights(distances, nw_estimator.sigma)
    print(f"weights shape: {weights.shape}, dtype: {weights.dtype}")
    print(f"After weights: {get_memory_usage():.1f} MB")
    
    print("\nTesting step 4: calculate_distance_W")
    print(f"Memory usage estimate: {W_grids.shape[0] * W_sites.shape[0] * W_grids.shape[1] * 8 / 1e9:.2f} GB")
    
    try:
        distance_W = nw_estimator.calculate_distance_W(W_grids, W_sites)
        print(f"distance_W shape: {distance_W.shape}, dtype: {distance_W.dtype}")
        print(f"After distance_W: {get_memory_usage():.1f} MB")
        
        print("\nTesting step 5: K function and prod")
        weights_W = nw_estimator.K(distance_W, nw_estimator.sigma).prod(axis=2)
        print(f"weights_W shape: {weights_W.shape}, dtype: {weights_W.dtype}")
        print(f"After weights_W: {get_memory_usage():.1f} MB")
        
        print("✅ All steps completed successfully!")
        
    except Exception as e:
        print(f"❌ Error at step 4 or 5: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_step_by_step()