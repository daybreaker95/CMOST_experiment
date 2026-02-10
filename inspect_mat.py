
import scipy.io
import os
import numpy as np

def print_structure(data, indent=0):
    indent_str = " " * indent
    if isinstance(data, dict):
        for key in data:
            if key.startswith('__'): continue
            print(f"{indent_str}{key}:")
            print_structure(data[key], indent + 2)
    elif isinstance(data, np.ndarray):
        print(f"{indent_str}Array shape: {data.shape}, dtype: {data.dtype}")
        if data.dtype.names:
            for name in data.dtype.names:
                print(f"{indent_str}Field: {name}")
                # logical_indices = [slice(None)] * (len(data.shape) or 1) # Not used
                # Just take the first element if it's an array of structs
                if data.size > 0:
                    val = data[0] if data.ndim > 0 else data.item()
                    # If it's a structured array, data[0] allows access by field name
                    try:
                         print_structure(val[name], indent + 2)
                    except:
                        pass
        else:
             if data.size < 10:
                 print(f"{indent_str}Values: {data}")
    elif isinstance(data, np.void): # For struct elements
        if data.dtype.names:
             for name in data.dtype.names:
                print(f"{indent_str}Field: {name}")
                print_structure(data[name], indent + 2)
    else:
        print(f"{indent_str}Type: {type(data)} Value: {str(data)[:100]}")

base_dir = r"c:\Users\daybr\Git_Works\CMOST_experiment\Settings"
files = [f for f in os.listdir(base_dir) if f.endswith('.mat')]

for f in files:
    print(f"\n--- Structure of {f} ---")
    try:
        data = scipy.io.loadmat(os.path.join(base_dir, f))
        print_structure(data)
    except Exception as e:
        print(f"Error reading {f}: {e}")
