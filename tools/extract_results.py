#!/usr/bin/env python3
"""
Simple utility script to extract loop models and CSV reports from subfolders.

This script searches the current directory for folders, then copies any
files containing "LOOP" in their name or ending in ".csv" into a new
folder named "<original_folder>_processed".
"""

import os
import shutil

# Get a list of all directories in the current path
all_folders = [f for f in os.listdir('.') if os.path.isdir(f)]
# Exclude any folders that already end in '_processed'
model_folders = [f for f in all_folders if not f.endswith('_processed')]

for folder in model_folders:
    print(f"Processing folder: {folder}")

    models_to_copy = []
    csv_to_copy = []

    # Find models and CSVs inside the folder
    for file in os.listdir(folder):
        file_path = os.path.join(folder, file)
        if "LOOP" in file:
            models_to_copy.append(file_path)
        elif file.endswith('.csv'):
            csv_to_copy.append(file_path)

    # Create destination folder
    destination_folder = f"{folder}_processed"
    os.makedirs(destination_folder, exist_ok=True)

    # Copy CSV (only if one is found)
    if len(csv_to_copy) == 1:
        shutil.copy(csv_to_copy[0], os.path.join(destination_folder, os.path.basename(csv_to_copy[0])))
        print(f"Copied CSV file to {destination_folder}")
    elif len(csv_to_copy) > 1:
        print(f"Warning: Found {len(csv_to_copy)} CSV files in {folder}. Please copy the correct one manually.")
    else:
        print(f"No CSV file found in {folder} to copy.")

    # Copy models
    for model_path in models_to_copy:
        dest_path = os.path.join(destination_folder, os.path.basename(model_path))
        try:
            if os.path.isdir(model_path):
                shutil.copytree(model_path, dest_path, dirs_exist_ok=True)
            else:
                shutil.copy(model_path, dest_path)
        except Exception as e:
            print(f"  Error copying {model_path} to {dest_path}: {e}")

    print(f"Copied {len(models_to_copy)} models to {destination_folder}\n")