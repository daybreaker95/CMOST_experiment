
import scipy.io
import os
import numpy as np
import pprint

def clean_data(data):
    """
    Recursively cleans data from scipy.io.loadmat.
    - Handles nested structs (void objects).
    - Flattens 1x1 arrays.
    - Converts object arrays to lists or keeps as arrays if multidimensional.
    """
    if isinstance(data, dict):
        return {k: clean_data(v) for k, v in data.items() if not k.startswith('__')}
    elif isinstance(data, np.ndarray):
        # Handle 0-d arrays
        if data.size == 0:
            return []
        
        # Handle structured arrays (fields)
        if data.dtype.names:
            # If it's a structured array, we might want to convert it to a list of dicts
            # or a single dict if it's 1x1
            if data.size == 1:
                # Treat as a single struct -> dict
                # data.item() returns the void object (tuple-like)
                # But to get fields we can iterate over dtype.names
                # data[0] is the first element
                element = data.item() # This is a numpy.void object
                # But element doesn't easy convert to dict. 
                # Better to use data[name] slicing which preserves structure
                clean_dict = {}
                for name in data.dtype.names:
                    clean_dict[name] = clean_data(data[name])
                return clean_dict
            else:
                # List of structs
                # This is trickier because data['Field'] returns an array of that field
                # We want a list of dicts: [ {Field: val1}, {Field: val2} ]
                # But MATLAB usually stores "struct array" as a single object where fields are arrays.
                # Let's check the shape.
                # Actually, usually in these files, the top level is 1x1.
                # If we have a real array of structs:
                result_list = []
                for i in range(data.size):
                    # We need to flatten the index
                    # This is getting complicated for generic cases.
                    # Let's try to convert per-field.
                    pass
                
                # Simplified approach for now: return as dict of arrays
                clean_dict = {}
                for name in data.dtype.names:
                    clean_dict[name] = clean_data(data[name])
                return clean_dict

        # Handle object arrays (cells)
        if data.dtype == 'O':
            # Flatten 1x1 object array
            if data.size == 1:
                return clean_data(data.item())
            # Convert 1D object array to list
            elif data.ndim == 1:
                return [clean_data(x) for x in data]
            # Convert 2D object array to list of lists? 
            # Or keep as numpy array of objects if useful?
            # Let's try to keep as native python list structure for cleaner file output
            else:
                return [clean_data(x) for x in data.flatten()] # Simple flattening for now, can improve if needed
        
        # Handle numeric arrays
        # Flatten 1x1 numeric array
        if data.size == 1:
            return data.item()
        
        # Squeeze dimensions of 1 (e.g. 5x1 -> 5,)
        data = np.squeeze(data)
        if data.ndim == 0: # It became a scalar
             return data.item()
             
        return data.tolist() # Convert to list for nicer writing to file (avoid array(...))

    elif isinstance(data, np.void):
        # Struct element
        if data.dtype.names:
             return {name: clean_data(data[name]) for name in data.dtype.names}
        return str(data)
    else:
        return data

def convert_mat_to_py(mat_file, py_file):
    try:
        mat_data = scipy.io.loadmat(mat_file)
        clean_settings = clean_data(mat_data)
        
        # Extract the likely root key (e.g. 'Calibration' or just all keys)
        # The inspection showed keys like 'temp' containing the data?
        # Let's look at the inspection output again.
        # CMOST13.mat -> 'temp'
        # CMOST8.mat -> 'temp'
        # It seems 'temp' is the main variable.
        
        content = {}
        if 'temp' in clean_settings:
            content = clean_settings['temp']
        else:
            content = clean_settings

        with open(py_file, 'w', encoding='utf-8') as f:
            f.write("# Auto-generated from {}\n".format(os.path.basename(mat_file)))
            f.write("import numpy as np\n")
            f.write("nan = np.nan\n")
            f.write("inf = np.inf\n\n")
            f.write("settings = ")
            f.write(pprint.pformat(content, indent=4, width=120))
            f.write("\n")
            
        print(f"Converted {mat_file} -> {py_file}")
        return True
    except Exception as e:
        print(f"Failed to convert {mat_file}: {e}")
        return False

if __name__ == "__main__":
    base_dir = r"c:\Users\daybr\Git_Works\CMOST_experiment\Settings"
    target_dir = r"c:\Users\daybr\Git_Works\CMOST_experiment\python\settings"
    
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    files = [f for f in os.listdir(base_dir) if f.endswith('.mat')]
    
    for f in files:
        mat_path = os.path.join(base_dir, f)
        py_name = os.path.splitext(f)[0] + ".py"
        py_path = os.path.join(target_dir, py_name)
        convert_mat_to_py(mat_path, py_path)
