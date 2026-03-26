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
    Output is redirected to the 'output_tools' folder.
    '''
    output_dir = Path("output_tools")
    output_dir.mkdir(exist_ok=True)

    tool_path = Path("tools").absolute() / script_name
    if not tool_path.exists():
        return f"Error: Tool {script_name} not found."
    
    processed_args = []
    if args:
        for arg in args:
            if (arg.startswith("./") or "/" in arg or os.path.exists(arg)) and not arg.startswith(("-", "--")):
                processed_args.append(str(Path(arg).absolute()))
            else:
                processed_args.append(arg)
    
    cmd = [sys.executable, str(tool_path)]
    cmd.extend(processed_args)
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, cwd=str(output_dir))
        return result.stdout
    except subprocess.CalledProcessError as e:
        error_msg = f"Error executing tool: {script_name}\n"
        error_msg += f"Return code: {e.returncode}\n"
        if e.stdout:
            error_msg += f"\n--- Standard Output ---\n{e.stdout}"
        if e.stderr:
            error_msg += f"\n--- Standard Error ---\n{e.stderr}"
        return error_msg
    except Exception as e:
        return f"Unexpected error running tool: {e}"

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

def delete_path(path: Union[str, Path]) -> bool:
    '''
    Deletes a file or directory.
    '''
    p = Path(path)
    if not p.exists():
        return False
    
    try:
        if p.is_file():
            p.unlink()
        elif p.is_dir():
            import shutil
            shutil.rmtree(p)
        return True
    except Exception:
        return False

def create_directory(path: Union[str, Path]) -> bool:
    '''
    Creates a directory if it doesn't exist.
    '''
    try:
        Path(path).mkdir(parents=True, exist_ok=True)
        return True
    except Exception:
        return False

def move_path(src: Union[str, Path], dst: Union[str, Path]) -> bool:
    '''
    Moves a file or directory to a new location.
    '''
    import shutil
    try:
        shutil.move(str(src), str(dst))
        return True
    except Exception:
        return False

def copy_path(src: Union[str, Path], dst: Union[str, Path]) -> bool:
    '''
    Copies a file or directory to a new location.
    '''
    import shutil
    try:
        src_path = Path(src)
        if src_path.is_dir():
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)
        return True
    except Exception:
        return False

def list_root_dirs() -> List[str]:
    '''
    Lists directories in the project root, excluding hidden ones.
    '''
    dirs = [d for d in os.listdir(".") if os.path.isdir(d) and not d.startswith(".")]
    return sorted(dirs)

def get_all_project_files() -> List[str]:
    '''
    Recursively lists all files in the project for autocompletion, 
    excluding hidden and internal directories.
    '''
    all_files = []
    exclude_dirs = {".pixi", ".git", "__pycache__", ".snakemake", ".nextflow", "work"}
    for root, dirs, files in os.walk("."):
        # Filter out excluded directories in-place to prevent walking into them
        dirs[:] = [d for d in dirs if d not in exclude_dirs and not d.startswith(".")]
        
        for f in files:
            if f.startswith("."): continue
            p = os.path.relpath(os.path.join(root, f), ".")
            all_files.append(p)
    return sorted(all_files)

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
