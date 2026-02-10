
import sys
import os
import numpy as np

# Add python directory to path
sys.path.append(r"c:\Users\daybr\Git_Works\CMOST_experiment\python")

def verify_module(module_name):
    try:
        mod = __import__(f"settings.{module_name}", fromlist=['settings'])
        settings = mod.settings
        
        print(f"Successfully imported {module_name}")
        print(f"  Keys: {len(settings)} keys found.")
        
        # Check specific expected keys based on inspection
        expected_keys = ['AdvRisk', 'Benchmarks', 'Calibration', 'DirectCancerRate']
        for k in expected_keys:
            if k in settings:
                val = settings[k]
                if isinstance(val, (list, np.ndarray)):
                    print(f"  {k}: Array/List with length {len(val)}")
                elif isinstance(val, dict):
                     print(f"  {k}: Dict with keys {list(val.keys())}")
                else:
                    print(f"  {k}: {type(val)}")
            else:
                print(f"  WARNING: Expected key '{k}' not found.")
                
        return True
    except Exception as e:
        print(f"Failed to verify {module_name}: {e}")
        return False

if __name__ == "__main__":
    print("Verifying Settings Modules...")
    verify_module("CMOST13")
    verify_module("CMOST19")
    verify_module("CMOST8")
