#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM (ui_utils)

'''
Utility functions for the PRISM Streamlit dashboard.
Handles Nextflow execution, PDB visualization, and file management.
'''

import os
import subprocess
import pandas as pd
import py3Dmol
import streamlit as st
import json
import sys
from pathlib import Path
from typing import List, Dict, Any, Union, Optional, Generator
from PRISM.config import settings

def run_nextflow() -> subprocess.Popen:
    '''
    Launches the Nextflow pipeline.
    '''
    cmd = ["nextflow", "run", "orchestrator.nf"]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return process

def get_nextflow_progress(process: subprocess.Popen) -> Generator[str, None, None]:
    '''
    Generator to read Nextflow output and yield progress updates.
    '''
    for line in process.stdout:
        yield line

def visualize_pdb(pdb_path: Union[str, Path], style: str = "cartoon", color: str = "spectrum") -> str:
    '''
    Visualizes a PDB file using py3Dmol with configurable styles.
    '''
    with open(pdb_path, 'r') as f:
        pdb_data = f.read()
    
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_data, 'pdb')
    
    if style == "cartoon":
        view.setStyle({'cartoon': {'color': color}})
    elif style == "sphere":
        view.setStyle({'sphere': {'color': color}})
    elif style == "stick":
        view.setStyle({'stick': {'color': color}})
    elif style == "line":
        view.setStyle({'line': {'color': color}})
    elif style == "cross":
        view.setStyle({'cross': {'color': color}})

    view.zoomTo()
    return view._make_html()

def list_files_in_dir(directory: Union[str, Path]) -> List[Dict[str, str]]:
    '''
    Lists files in a directory with metadata.
    '''
    path = Path(directory)
    if not path.exists():
        return []
    files = []
    for f in path.iterdir():
        if f.is_file():
            files.append({
                "name": f.name,
                "size": f"{f.stat().st_size / 1024:.1f} KB",
                "modified": pd.to_datetime(f.stat().st_mtime, unit='s').strftime('%Y-%m-%d %H:%M:%S')
            })
    return files

def run_tool(script_name: str, args: Optional[List[str]] = None) -> str:
    '''
    Executes a tool from the tools/ directory.
    '''
    tool_path = Path("tools") / script_name
    if not tool_path.exists():
        return f"Error: Tool {script_name} not found."
    
    cmd = [sys.executable, str(tool_path)]
    if args:
        cmd.extend(args)
    
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return result.stdout

def get_score_distribution_data() -> Optional[pd.DataFrame]:
    '''
    Reads the final ranking CSV and returns data for plotting.
    '''
    csv_path = settings.FINAL_RANKING_CSV
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path)
        return df
    return None

def load_ranking_csv() -> Optional[pd.DataFrame]:
    '''
    Loads the final ranking CSV if it exists.
    '''
    csv_path = settings.FINAL_RANKING_CSV
    if os.path.exists(csv_path):
        return pd.read_csv(csv_path)
    return None

def delete_file(directory: Union[str, Path], filename: str) -> bool:
    '''
    Deletes a file from the specified directory.
    '''
    file_path = Path(directory) / filename
    if file_path.exists():
        os.remove(file_path)
        return True
    return False

def save_uploaded_file(uploaded_file: Any, target_dir: Union[str, Path]) -> str:
    '''
    Saves an uploaded file to the target directory.
    '''
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    file_path = os.path.join(target_dir, uploaded_file.name)
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return file_path
