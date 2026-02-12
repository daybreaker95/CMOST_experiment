import os
import numpy as np
import pandas as pd

directory = r'C:\Users\daybr\Git_Works\CMOST_experiment\python\Results'

for filename in os.listdir(directory):
    if filename.endswith('.npz'):
        file_path = os.path.join(directory, filename)
        with np.load(file_path, allow_pickle=True) as data:
            # Define pairs to merge
            pairs = {
                'Variable': 'Var_Legend',
                'Benchmark': 'BM_Description'
            }
            
            processed_keys = set()
            
            # Process pairs
            for val_key, desc_key in pairs.items():
                if val_key in data.files and desc_key in data.files:
                    try:
                        # Handle 0-d arrays for values
                        values = data[val_key]
                        if values.ndim == 0:
                            values = [values.item()]
                            
                        # Handle 0-d arrays for descriptions
                        descriptions = data[desc_key]
                        if descriptions.ndim == 0:
                            descriptions = [descriptions.item()]
                        
                        # Create combined DataFrame
                        # Ensure lengths match or handle mismatch
                        min_len = min(len(values), len(descriptions))
                        df = pd.DataFrame({
                            'Description': descriptions[:min_len],
                            'Value': values[:min_len]
                        })
                        
                        output_name = f"{os.path.splitext(filename)[0]}_{val_key}_Combined.csv"
                        df.to_csv(os.path.join(directory, output_name), index=False)
                        
                        processed_keys.add(val_key)
                        processed_keys.add(desc_key)
                        print(f"Created combined file: {output_name}")
                    except Exception as e:
                        print(f"Failed to combine {val_key} and {desc_key}: {e}")

            # Process remaining keys
            for key in data.files:
                if key not in processed_keys:
                    output_name = f"{os.path.splitext(filename)[0]}_{key}.csv" if len(data.files) > 1 else f"{os.path.splitext(filename)[0]}.csv"
                    content = data[key]
                    if content.ndim == 0:
                        content = [content.item()]
                    pd.DataFrame(content).to_csv(os.path.join(directory, output_name), index=False)
