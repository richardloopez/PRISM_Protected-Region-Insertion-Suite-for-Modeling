#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM (Protected-Region Insertion Suite for Modeling)
#
# Citation:
# If you use this software in your research, please cite:
# Lopez-Corbalan, R.

"""
Homology modeling module with parallel execution support.

This module handles the initial homology modeling phase using multiple
templates while keeping experimental regions completely fixed.
"""

import os
from typing import List, Set

from modeller import *
from modeller.parallel import Job
from modeller.automodel import *

from . import config
from .config import (
    ALIGN_CODES_TEMPLATES,
    ALIGN_CODE_SEQUENCE,
    NUM_MODELS_AUTO,
    NUM_MODELS_TO_REFINE,
    CHAIN_ID
)
from .custom_models import FixedRegionAutoModel


def run_automodel(env: Environ, align_file: str, job: Job,
                  experimental_residues: Set[int]) -> List[str]:
    """
    Execute AutoModel using all templates to generate initial homology models.

    Generates models with gap filling, ranks them by DOPE-HR score, renames them
    with traceable nomenclature, and returns the top N models for loop refinement.

    Args:
        env: MODELLER environment
        align_file: Path to PIR alignment file
        job: Job object for parallel processing
        experimental_residues: Set of residue numbers that must NOT be optimized

    Returns:
        List of PDB filenames for the best models selected for loop refinement

    Note:
        Models are renamed to AUTO_1.pdb, AUTO_2.pdb, etc. in order of quality
    """
    initial_models_for_loop_names: List[str] = []

    knowns_tuple = tuple(ALIGN_CODES_TEMPLATES)

    print(f"\n{'='*80}")
    print(f"[STEP 3] Starting AutoModel (Gap Filling) with {NUM_MODELS_AUTO} models...")
    print(f"[STEP 3] Using {len(knowns_tuple)} templates: {knowns_tuple}")
    print(f"[STEP 3] 🔒 FIXED MODE: {len(experimental_residues)} experimental residues will NOT be optimized")
    print(f"{'='*80}\n")

    a = FixedRegionAutoModel(env,
                             experimental_residues=experimental_residues,
                             chain_id=CHAIN_ID,
                             alnfile=align_file,
                             knowns=knowns_tuple,
                             sequence=ALIGN_CODE_SEQUENCE,
                             assess_methods=(assess.DOPEHR, assess.GA341))

    a.use_parallel_job(job)
    a.starting_model = 1
    a.ending_model = NUM_MODELS_AUTO
    a.library_schedule = autosched.slow
    a.max_var_iterations = 1000
    a.make()

    results_auto = a.outputs
    if not results_auto:
        print("[ERROR] AutoModel failed to generate models.")
        return []

    sorted_auto_models = sorted(results_auto, key=lambda x: x.get('DOPE-HR score', 9999999.0))

    print(f"\n[STEP 3.1] Renaming {len(sorted_auto_models)} AutoModel models with traceable nomenclature.")

    for model_rank, model_info in enumerate(sorted_auto_models):
        old_name = model_info['name']
        new_name = f'AUTO_{model_rank+1}.pdb'

        try:
            os.rename(old_name, new_name)
            if (model_rank + 1) <= NUM_MODELS_TO_REFINE:
                initial_models_for_loop_names.append(new_name)
        except Exception as e:
            print(f"    [ERROR] Could not rename {old_name} to {new_name}. Error: {e}")

    print(f"\n[STEP 3.2] ✓ AutoModel completed. {len(initial_models_for_loop_names)} models selected for loop refinement.")

    return initial_models_for_loop_names
