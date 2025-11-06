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
Loop refinement module with parallel execution support.

This module handles sequential loop refinement for flexible regions while
maintaining experimental coordinates completely fixed.
"""

import os
from typing import List, Tuple, Set

from modeller import *
from modeller.automodel import *
from modeller.parallel import Job
from modeller.automodel import *

from . import config
from .config import ALIGN_CODE_SEQUENCE, CHAIN_ID, NUM_MODELS_LOOP
from .custom_models import FixedRegionLoopRefiner


def run_loop_refinement(env: Environ, job: Job, initial_models_names: List[str],
                       loop_ranges: List[Tuple[int, int]],
                       experimental_residues: Set[int]):
    """
    Execute sequential loop refinement with DOPE-HR assessment for initial models.

    Automatically excludes loops that contain any experimental residues from
    the main template to ensure coordinate integrity.

    Args:
        env: MODELLER environment
        job: Job object for parallel processing
        initial_models_names: List of initial model PDB filenames
        loop_ranges: List of (start, end) tuples defining loops to refine
        experimental_residues: Set of residue numbers that must NOT be refined

    Note:
        Models are renamed with traceable nomenclature:
        AUTO_1_LOOP1_R1.pdb = Auto model 1, Loop 1, Refinement 1
    """

    if not loop_ranges:
        print("\n[STEP 4] Skipping loop refinement: No flexible loops defined.")
        return

    print(f"\n{'='*80}")
    print(f"[STEP 4] Filtering loops that contain experimental residues from MAIN_PDB...")
    print(f"{'='*80}\n")

    valid_loop_ranges = [r for r in loop_ranges if 4 <= (r[1] - r[0] + 1) <= 1000]

    if not valid_loop_ranges:
        print("\n[STEP 4] Skipping loop refinement: No loops meet length requirements (4-1000 residues).")
        return

    print(f"\n{'='*80}")
    print(f"[STEP 4.1] Starting targeted refinement for {len(valid_loop_ranges)} valid segments...")
    print(f"[STEP 4.1] 🔒 GUARANTEE: No loop contains experimental residues from MAIN_PDB")
    print(f"{'='*80}\n")

    for model_index, initial_pdb_file in enumerate(initial_models_names):

        base_name = initial_pdb_file.replace('.pdb', '')

        print(f"\n --- Processing Base Model #{model_index+1}: {initial_pdb_file} ({base_name}) ---")

        current_best_pdb_for_thread = initial_pdb_file
        current_base_name_for_refinement = base_name

        for j, (start, end) in enumerate(valid_loop_ranges):

            print(f"  > Refining Loop {j+1}/{len(valid_loop_ranges)}: Residues {start} to {end}")

            try:
                ml = FixedRegionLoopRefiner(env,
                                           inimodel=current_best_pdb_for_thread,
                                           sequence=ALIGN_CODE_SEQUENCE,
                                           loop_start=start,
                                           loop_end=end,
                                           chain_id=CHAIN_ID,
                                           experimental_residues=experimental_residues)

                ml.use_parallel_job(job)
                ml.loop.starting_model = 1
                ml.loop.ending_model = NUM_MODELS_LOOP
                ml.loop.md_level = refine.slow_large
                ml.md_level = None
                ml.loop.assess_methods = (assess.DOPEHR, assess.GA341)
                ml.max_var_iterations = 1000

                ml.make()

                loop_models_of_this_step = ml.loop.outputs

                if loop_models_of_this_step:
                    sorted_loop_outputs_by_loop_dopeHR = sorted(loop_models_of_this_step,
                                                                key=lambda x: x.get('DOPE-HR score', 9999999.0))

                    for m, model_info in enumerate(sorted_loop_outputs_by_loop_dopeHR):
                        old_name = model_info['name']
                        new_loop_name = f'{current_base_name_for_refinement}_LOOP{j+1}_R{m+1}.pdb'

                        try:
                            os.rename(old_name, new_loop_name)
                        except Exception as e:
                            print(f"    [ERROR] Could not rename {old_name} to {new_loop_name}. Error: {e}")

                    best_of_this_loop = sorted_loop_outputs_by_loop_dopeHR[0]
                    current_best_pdb_for_thread = f'{current_base_name_for_refinement}_LOOP{j+1}_R1.pdb'
                    current_base_name_for_refinement = f'{current_base_name_for_refinement}_LOOP{j+1}_R1'

                else:
                    print(f"    -> Warning: FixedRegionLoopRefiner did not generate valid results for {start}-{end}. Using previous best model.")

            except ValueError as ve:
                print(f"     > VALIDATION ERROR in Loop {start}-{end}: {ve}")
                print(f"     > This loop will be SKIPPED for safety.")
                continue
            except Exception as e:
                print(f"     > FATAL ERROR in FixedRegionLoopRefiner for {start}-{end} (Base Model #{model_index+1}): {e}")
                continue

        print(f"\n[STEP 4.2] ✓ Loop refinement completed for base model #{model_index + 1}.")
