#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM Orchestrator (Python Version)
# Replaces run_prism.sh

"""
PRISM Workflow Orchestrator.

This script submits the full PRISM pipeline to the Slurm scheduler.

Usage:
    python3 orchestrator.py
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path


current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

try:
    from PRISM import config
except ImportError:
    print("[ORCHESTRATOR][FATAL] Could not import PRISM.config. Run this script from the project root.")
    sys.exit(1)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Define specific paths
PROJECT_ROOT = current_dir
RESULTS_DIR = os.path.join(PROJECT_ROOT, "modeling_results")
LOG_DIR = os.path.join(PROJECT_ROOT, "logs")

# Define the MODELLER execution command
PYTHON_EXE = shutil.which("python3")
MOD_PYTHON_CMD = f"/home/richard/bin/modeller10.7/bin/modpy.sh {PYTHON_EXE}"

# Ensure directories exist
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

# =============================================================================
# HELPER: JOB SUBMISSION
# =============================================================================

def submit_slurm_job(job_name, command_str, dep_job_id=None, array_range=None, 
                     cpus=1, mem="4G", output_prefix=None):
    """
    Submits a job to Slurm using sbatch and returns the Job ID.
    """
    if output_prefix is None:
        output_prefix = job_name.lower()

    # Base sbatch command
    sbatch_args = [
        "sbatch",
        "--parsable",  # Returns only the job ID
        f"--job-name={job_name}",
        f"--output={LOG_DIR}/{output_prefix}_%A_%a.log" if array_range else f"--output={LOG_DIR}/{output_prefix}_%j.log",
        f"--error={LOG_DIR}/{output_prefix}_%A_%a.err" if array_range else f"--error={LOG_DIR}/{output_prefix}_%j.err",
        f"--cpus-per-task={cpus}",
        f"--mem={mem}"
    ]

    if dep_job_id:
        sbatch_args.append(f"--dependency=afterok:{dep_job_id}")

    if array_range:
        sbatch_args.append(f"--array={array_range}")

    wrapped_cmd = f"export PYTHONPATH={PROJECT_ROOT}:${{PYTHONPATH}}; {command_str}"
    sbatch_args.append(f"--wrap={wrapped_cmd}")

    print(f"\n[ORCHESTRATOR] Submitting {job_name}...")
    try:
        result = subprocess.run(sbatch_args, capture_output=True, text=True, check=True)
        job_id = result.stdout.strip()
        print(f" > [ORCHESTRATOR] Job ID: {job_id}")
        if array_range:
            print(f" > [ORCHESTRATOR] Array: {array_range}")
        return job_id
    except subprocess.CalledProcessError as e:
        print(f"[ORCHESTRATOR][FATAL] Submission failed: {e.stderr}")
        sys.exit(1)

# =============================================================================
# MAIN WORKFLOW
# =============================================================================

def main():
    print("="*60)
    print("      PRISM PIPELINE: WORKFLOW ORCHESTRATION")
    print("="*60)
    print(f"Project Root : {PROJECT_ROOT}")
    print(f"Python Config: Loaded successfully")
    print(f"AutoModel    : {config.NUM_SPLITTED_PROCESSES} array tasks")
    print(f"LoopModel    : {config.NUM_MODELS_TO_REFINE} array tasks")
    print("-" * 60)

    # --- PRECALCULATION? ---
    if config.RSR_INI_PRECALCULATION:
        print("\n[ORCHESTRATOR] Running in PRECALCULATION mode.")    

    # --- STAGE 1: Prereq-CDE ---
    # Generates the alignment files
    cmd_s1 = f"{MOD_PYTHON_CMD} -m PRISM.controller --stage prereq-cde"
    
    id_s1 = submit_slurm_job(
        job_name="PRISM_S1_CDE",
        command_str=cmd_s1,
        cpus=1,
        mem="4G",
        output_prefix="S1_cde"
    )

    # --- STAGE 2: AutoModel (Job Array) ---
    # Generates initial models in parallel

    # ----- RSR_INI_PRECALCULATION MODE -----
    if config.RSR_INI_PRECALCULATION:
        print("\n[ORCHESTRATOR] Precalculation active: Submitting S2 as a single task.")
        cmd_s2 = (
            f'TASK_DIR="{RESULTS_DIR}/precalc_task"; '
            f'mkdir -p ${{TASK_DIR}}; '
            f'cd ${{TASK_DIR}}; '
            f'{MOD_PYTHON_CMD} -m PRISM.controller --stage automodel --job-id 1; '
            f'cd {RESULTS_DIR}; '
            f'rmdir ${{TASK_DIR}} 2>/dev/null || true'
        )

        id_s2 = submit_slurm_job(
            job_name="PRISM_S2_Precalc",
            command_str=cmd_s2,
            dep_job_id=id_s1,
            array_range=None,
            cpus=1,
            mem="8G",
            output_prefix="S2_precalc"
        )

        sys.exit(0)
    
    # ----- DEFAULT MODE ----- 
    else:
        cmd_s2 = (
            f'TASK_DIR="{RESULTS_DIR}/task_${{SLURM_ARRAY_TASK_ID}}"; '
            f'mkdir -p ${{TASK_DIR}}; '
            f'cd ${{TASK_DIR}}; '
            f'{MOD_PYTHON_CMD} -m PRISM.controller --stage automodel --job-id ${{SLURM_ARRAY_TASK_ID}}; '
        f'mv *.pdb {RESULTS_DIR}/ 2>/dev/null; '
        f'cd {RESULTS_DIR}; '
        f'rmdir ${{TASK_DIR}} 2>/dev/null || true'
    )

    id_s2 = submit_slurm_job(
        job_name="PRISM_S2_Auto",
        command_str=cmd_s2,
        dep_job_id=id_s1,
        array_range=f"1-{config.NUM_SPLITTED_PROCESSES}",
        cpus=10,
        mem="80G",
        output_prefix="S2_auto"
    )

    # --- STAGE 3: Rank-AutoModel ---
    # Ranks and renames the output of S2
    cmd_s3 = (
        f'cd {RESULTS_DIR}; '
        f'{MOD_PYTHON_CMD} -m PRISM.controller --stage rank-automodel'
    )

    id_s3 = submit_slurm_job(
        job_name="PRISM_S3_Rank",
        command_str=cmd_s3,
        dep_job_id=id_s2,
        cpus=1,
        mem="4G",
        output_prefix="S3_rank"
    )

    # --- STAGE 4: Loop Refinement (Job Array) ---
    # Refines the best models from S3
    cmd_s4 = (
        f'MODEL_NAME="AUTO_${{SLURM_ARRAY_TASK_ID}}.pdb"; '
        f'TASK_DIR="{RESULTS_DIR}/loop_${{SLURM_ARRAY_TASK_ID}}"; '
        f'mkdir -p ${{TASK_DIR}}; '
        f'cp {RESULTS_DIR}/${{MODEL_NAME}} ${{TASK_DIR}}/; '
        f'cd ${{TASK_DIR}}; '
        f'{MOD_PYTHON_CMD} -m PRISM.controller --stage loopmodel --input-model ${{MODEL_NAME}}; '
        f'mv *.pdb {RESULTS_DIR}/ 2>/dev/null; '
        f'cd {RESULTS_DIR}; '
        f'rm -r ${{TASK_DIR}}'
    )

    id_s4 = submit_slurm_job(
        job_name="PRISM_S4_Loop",
        command_str=cmd_s4,
        dep_job_id=id_s3,
        array_range=f"1-{config.NUM_MODELS_TO_REFINE}",
        cpus=10,
        mem="80G",
        output_prefix="S4_loop"
    )

    # --- STAGE 5: Final Ranking ---
    # Ranks the refined loop models
    cmd_s5 = (
        f'cd {RESULTS_DIR}; '
        f'{MOD_PYTHON_CMD} -m PRISM.controller --stage final-rank'
    )

    id_s5 = submit_slurm_job(
        job_name="PRISM_S5_Final",
        command_str=cmd_s5,
        dep_job_id=id_s4,
        cpus=1,
        mem="4G",
        output_prefix="S5_final"
    )

    print("="*60)
    print(" Workflow submitted successfully.")
    print(" Monitor with: squeue -u $USER")
    print("="*60)

if __name__ == "__main__":
    main()