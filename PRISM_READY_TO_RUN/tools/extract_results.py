#!/usr/bin/env python3
"""
Simple utility script to extract loop models and CSV reports from subfolders.

This script searches the current directory for folders, then copies any
files containing "LOOP" in their name or ending in ".csv" into a new
folder named "<original_folder>_processed".
"""

import os
import shutil
import pandas as pd



numbers_pdbs_to_copy = 100
modeling_results_directory = "modeling_results"
destination_folder = "modeling_results_processed"

# Get the modeling folder:
def find_modeling_folder():
    for f in os.listdir('.'):
        if os.path.isdir(f) and f == modeling_results_directory:
            return f
    print(f"Not found '{modeling_results_directory}'")
    return None


# Get the csv
def find_csv(modeling_folder):
    csv_files = []
    for file in os.listdir(modeling_folder):
        if file.endswith('.csv'):
            csv_files.append(os.path.join(modeling_folder, file))
    if len(csv_files) == 1:
        return csv_files[0]
    elif len(csv_files) > 1:
        print(f"Warning: Found {len(csv_files)} CSV files in {modeling_folder}. Please copy the correct one manually.")
    else:
        print(f"No CSV file found in {modeling_folder} to copy.")
    

# Get the pdbs to copy
def get_pdbs_to_copy(modeling_folder, csv_file_path):
    paths_pdbs_to_copy = []
    df = pd.read_csv(csv_file_path)
    pdbs = df['Model Name'].head(numbers_pdbs_to_copy).tolist()
    for pdb in pdbs:
        pdb_path = os.path.join(modeling_folder, pdb)
        paths_pdbs_to_copy.append(pdb_path)
        print(f"Selected: {pdb}")
    
    if len(paths_pdbs_to_copy) == 0:
        print(f"No PDB files found in {modeling_folder} to copy.")
    else:
        print(f"Selected {len(paths_pdbs_to_copy)} PDB files to copy.") 

    return paths_pdbs_to_copy


# Copy the pdbs to the destination folder
def copy_data(paths_pdbs_to_copy, destination_folder, csv_file_path):
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)
    for pdb_path in paths_pdbs_to_copy:
        shutil.copy(pdb_path, destination_folder)
    print(f"Copied {len(paths_pdbs_to_copy)} PDB files to {destination_folder}")
    shutil.copy(csv_file_path, destination_folder)
    print(f"Copied {csv_file_path} to {destination_folder}")

# Add the rank index to the pdb file name
def add_rank_to_pdbs(destination_folder):
    renames = 0
    for file in os.listdir(destination_folder):
        if file.endswith('.csv'):
            df = pd.read_csv(os.path.join(destination_folder, file))
            df_top = df.head(numbers_pdbs_to_copy)
            name_to_rank = dict(zip(df_top['Model Name'], df_top['Rank']))
            for old_name in os.listdir(destination_folder):
                if not old_name.endswith('.pdb'):
                    continue
                if old_name in name_to_rank:
                    rank = name_to_rank[old_name]
                    new_name = f"RANK_{rank}_{old_name}"
                    old_path = os.path.join(destination_folder, old_name)
                    new_path = os.path.join(destination_folder, new_name)
                    os.rename(old_path, new_path)
                    print(f"Renamed {old_name} to {new_name}")
                    renames += 1
    print(f"Renamed {renames} PDB files")



def main():
    modeling_folder = find_modeling_folder()
    csv_file_path = find_csv(modeling_folder)
    paths_pdbs_to_copy = get_pdbs_to_copy(modeling_folder, csv_file_path)
    copy_data(paths_pdbs_to_copy, destination_folder, csv_file_path)
    add_rank_to_pdbs(destination_folder)
    
if __name__ == "__main__":
    main()


    
