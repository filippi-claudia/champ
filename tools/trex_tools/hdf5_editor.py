import h5py
import numpy as np
import sys
import os
import argparse

def interactive_hdf5_editor(file_path: str) -> None:
    """
    Provides a CLI interface to explore groups, view datasets with full precision,
    and edit values in-place within an HDF5 file.
    """
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' not found.")
        return

    # Open the file in 'r+' (read/write) mode. It must already exist.
    with h5py.File(file_path, 'r+') as h5f:
        
        while True:
            # 1. Show all top-level groups
            groups = [k for k, v in h5f.items() if isinstance(v, h5py.Group)]
            if not groups:
                print("No groups found in the root of the file.")
                return
                
            print("\n--- Available Groups ---")
            for i, g in enumerate(groups):
                print(f"[{i}] {g}")
                
            g_input = input("\nSelect a group number (or 'q' to quit): ").strip()
            if g_input.lower() == 'q':
                return
                
            try:
                g_idx = int(g_input)
                selected_group = h5f[groups[g_idx]]
            except (ValueError, IndexError):
                print("Invalid selection. Try again.")
                continue

            while True:
                # 2. Show datasets within the selected group
                datasets = [k for k, v in selected_group.items() if isinstance(v, h5py.Dataset)]
                if not datasets:
                    print("No datasets found in this group.")
                    break
                    
                print(f"\n--- Datasets in '{groups[g_idx]}' ---")
                for i, d in enumerate(datasets):
                    print(f"[{i}] {d}")
                    
                d_input = input("\nSelect a dataset number ('b' to go back, 'q' to quit): ").strip()
                if d_input.lower() == 'q':
                    return
                elif d_input.lower() == 'b':
                    break
                    
                try:
                    d_idx = int(d_input)
                    dset_name = datasets[d_idx]
                    dset = selected_group[dset_name]
                except (ValueError, IndexError):
                    print("Invalid selection. Try again.")
                    continue

                # 3. Print dataset with full precision
                print(f"\n--- Dataset: {dset.name} ---")
                print(f"Shape: {dset.shape}, Dtype: {dset.dtype}")
                
                # Guard against printing massive datasets directly to the terminal
                total_elements = np.prod(dset.shape)
                print_data = True
                if total_elements > 10000:
                    print(f"\nWarning: Dataset contains {total_elements} elements.")
                    proceed = input("Printing this may flood your terminal. Print anyway? (y/n): ")
                    if proceed.lower() != 'y':
                        print_data = False

                if print_data:
                    print("\nData values (Indexed View):")
                    data = dset[...]
                    if data.ndim == 0:
                        val = data.item() if hasattr(data, 'item') else data
                        if isinstance(val, bytes):
                            val = val.decode('utf-8', errors='replace')
                        print(f"[]: {repr(val)}")
                    else:
                        for idx, val in np.ndenumerate(data):
                            idx_str = ", ".join(map(str, idx)) if len(idx) > 1 else str(idx[0])
                            # Use item() to extract native Python type to ensure repr() shows full precision
                            val_py = val.item() if hasattr(val, 'item') else val
                            if isinstance(val_py, bytes):
                                val_py = val_py.decode('utf-8', errors='replace')
                            print(f"[{idx_str}]: {repr(val_py)}")

                # 4. In-place Editing Loop
                print("\n--- Edit Mode ---")
                print("Enter the index of the element to edit.")
                print("  - For 1D arrays: enter a single number (e.g., '5')")
                print("  - For multi-dimensional arrays: enter comma-separated numbers (e.g., '0, 1')")
                print("  - To replace all occurrences of a specific value, type 'replaceAll'")
                print("Type 'b' to go back, 'q' to quit.")
                
                quit_script = False
                while True:
                    idx_str = input("\nIndex to edit (or 'replaceAll', 'b' to go back, 'q' to quit): ").strip()
                    if idx_str.lower() == 'q':
                        quit_script = True
                        break
                    elif idx_str.lower() == 'b':
                        break
                        
                    if idx_str.lower() == 'replaceall':
                        try:
                            old_val_str = input("Enter the value to be replaced: ")
                            new_val_str = input("Enter the new value: ")
                            
                            old_val = np.array(old_val_str, dtype=dset.dtype)
                            new_val = np.array(new_val_str, dtype=dset.dtype)
                            
                            # Compute mask
                            data = dset[...]
                            # Use isclose for floats to avoid tiny precision issues during parsing, otherwise exact match
                            if np.issubdtype(dset.dtype, np.floating):
                                mask = np.isclose(data, old_val, rtol=1e-15, atol=1e-15)
                            else:
                                mask = (data == old_val)
                                
                            count = np.sum(mask)
                            
                            if count > 0:
                                data[mask] = new_val
                                dset[...] = data
                                print(f"Successfully replaced {count} occurrences of {old_val} with {new_val}")
                            else:
                                print(f"Value {old_val} not found in the dataset.")
                        except Exception as e:
                            print(f"\nError: {e}")
                        continue
                        
                    try:
                        # Parse the input string into a tuple of integers for multi-dimensional indexing
                        if ',' in idx_str:
                            idx = tuple(int(x.strip()) for x in idx_str.split(','))
                        else:
                            idx = int(idx_str)
                            
                        # Fetch and display the exact current value at that index
                        current_val = dset[idx]
                        current_val_py = current_val.item() if hasattr(current_val, 'item') else current_val
                        if isinstance(current_val_py, bytes):
                            current_val_py = current_val_py.decode('utf-8', errors='replace')
                        print(f"Current value at {idx}: {repr(current_val_py)}")
                        
                        val_str = input("Enter new value: ")
                        
                        # Cast the user's string input strictly to the dataset's existing data type
                        new_val = np.array(val_str, dtype=dset.dtype)
                        
                        # Perform the in-place modification directly on disk via h5py's slicing translation
                        dset[idx] = new_val
                        
                        print(f"Successfully updated {idx} to {new_val}")
                        
                    except Exception as e:
                        print(f"\nError: {e}")
                        print(f"Ensure the index falls within shape {dset.shape} and value matches dtype {dset.dtype}.")
                        
                if quit_script:
                    return

if __name__ == "__main__":
    description = """
Interactive HDF5 Editor
=======================
This script provides an interactive command-line interface to explore, view,
and edit the contents of HDF5 files. 

Features:
  - Hierarchical Navigation: Browse through top-level groups and select datasets.
  - Indexed and Decoded Viewing: See full precision float arrays alongside their exact 
    indices. Bytes-encoded text datasets are safely decoded to UTF-8 for readability.
  - Value Modification: Replace specific indices within 1D and multi-dimensional matrices. 
  - Bulk Replacement: Enter the 'replaceAll' command to globally swap a specific value 
    with a new one within a dataset without needing to update elements one-by-one!

How to Use:
  1. Run this script passing the target HDF5 file path as an argument.
  2. Enter the number corresponding to the target group and then dataset.
  3. Inside the internal Edit Mode:
      - Type an index (e.g., '5' for 1D or '0, 1' for 2D) to update a single value.
      - Type 'replaceAll' to identify a value and swap it universally across the dataset.
      - Type 'b' to go back one menu step (e.g., returning to dataset selection).
      - Type 'q' to completely exit the script safely.

Note: Data types are rigidly enforced. Inserted text or numbers will be cast directly into the parent dataset's native dtype. Keep backups of key files!
"""
    parser = argparse.ArgumentParser(
        description=description, 
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "file_path", 
        type=str, 
        help="Path to the HDF5 (.h5, .hdf5) file to edit."
    )
    
    args = parser.parse_args()
    interactive_hdf5_editor(args.file_path)