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
Main controller for the PRISM modeling pipeline.

This module orchestrates the complete protein modeling workflow:
1. Setup and alignment preparation
2. Experimental residue identification
3. Loop detection from secondary structure
4. Homology modeling with fixed regions
5. Loop refinement
6. Final model evaluation and ranking
"""

import sys
import os
from modeller import *
from modeller.parallel import Job, LocalWorker

from . import config
from .config import (
    ALIGNMENT_FILE, ALIGNMENT_CDE_FILE, NUM_PROCESSORS, USE_MANUAL_ALIGNMENT,
    MODELING_RESULTS_DIR, INPUT_DIR,
    SS2_FILE_PATH, sequence_full,
    MAIN_ALIGN_CODE_TEMPLATE, ALIGN_CODE_SEQUENCE,
    EXPERIMENTAL_FLANK_SIZE,
    MANUAL_ALIGNMENT_FILE, MANUAL_ALIGNMENT_CDE_FILE
)

from . import utils
from . import homology_modeling
from . import loop_refinement
from .fixed_region_utils import identify_experimental_residues
from .utils import (
    get_coil_residues, group_ranges, read_sequences_from_ali_temp,
    generate_pir_files_AUTO, extract_cde_line,
    add_cde_line_to_pir
)


def main_workflow():
    """
    Execute the complete PRISM modeling pipeline with fixed experimental regions.

    Pipeline Steps:
        1. Setup directories and MODELLER environment
        2. Prepare alignment (manual or automatic)
        3. Identify experimental residues from main template
        4. Detect loops from secondary structure (coil + experimental flanks)
        5. Run homology modeling with fixed core
        6. Run loop refinement on flexible regions
        7. Evaluate and rank all generated models

    The pipeline ensures experimental coordinates remain completely fixed
    while allowing targeted refinement of flexible loop regions.
    """

    print(f"\n[SETUP] Ensuring results directory exists: {MODELING_RESULTS_DIR}")
    os.makedirs(MODELING_RESULTS_DIR, exist_ok=True)

    print(f"[SETUP] Changing working directory to: {MODELING_RESULTS_DIR}")
    os.chdir(MODELING_RESULTS_DIR)

    print("\n" + "="*80)
    print("PRISM: PROTECTED-REGION INSERTION SUITE FOR MODELING")
    print("Homology Modeling Pipeline with Fixed Experimental Regions")
    print(f"Working Directory: {os.getcwd()}")
    print("="*80 + "\n")

    env = Environ()

    env.io.atom_files_directory = ['.', INPUT_DIR, '../atom_files']
    env.io.hetatm = True

    try:
        log.verbose()
    except Exception as e:
        print(f"[WARNING] Failed to enable MODELLER logging. Execution will continue. Error: {e}")

    env.jobs = config.NUM_PROCESSORS
    job = Job()

    print(f"[PARALLEL] Configuring {NUM_PROCESSORS} local workers.")
    for _ in range(NUM_PROCESSORS):
        job.append(LocalWorker())
    job.start()

    print(f"\n{'='*80}")
    print(f"[STEP 1] Alignment Preparation")
    print(f"{'='*80}\n")

    ali_file_for_automodel: str = ""
    ali_file_for_analysis: str = ""

    try:
        if USE_MANUAL_ALIGNMENT:
            print(f"[STEP 1] Using Manual Alignment (USE_MANUAL_ALIGNMENT = True)")

            ali_file_for_automodel = MANUAL_ALIGNMENT_FILE
            ali_file_for_analysis = MANUAL_ALIGNMENT_CDE_FILE

            if not os.path.exists(ali_file_for_automodel):
                raise FileNotFoundError(f"Manual alignment file NOT found: {ali_file_for_automodel}")

            print(f"  > File for AutoModel (clean): {os.path.basename(ali_file_for_automodel)}")
            print(f"  > Generating Analysis file (with CDE) at: {os.path.basename(ali_file_for_analysis)}")

            utils.add_cde_line_to_pir(
                clean_ali_path=ali_file_for_automodel,
                cde_ali_path=ali_file_for_analysis,
                ss2_file_path=SS2_FILE_PATH,
                target_seq_full=sequence_full,
                align_code_sequence=ALIGN_CODE_SEQUENCE
            )

            print(f"  > Analysis file (CDE) generated: {os.path.basename(ali_file_for_analysis)}")

        else:
            print(f"[STEP 1] Using Automatic Alignment (USE_MANUAL_ALIGNMENT = False)")
            ali_file_for_automodel = ALIGNMENT_FILE
            ali_file_for_analysis = ALIGNMENT_CDE_FILE

            utils.generate_pir_files_AUTO(
                env, ali_file_for_automodel, ali_file_for_analysis
            )

        print(f"\n[STEP 1.1] Reading sequence map from: {os.path.basename(ali_file_for_analysis)}")
        all_aligned_seqs_map = utils.read_sequences_from_ali_temp(ali_file_for_analysis)

        if MAIN_ALIGN_CODE_TEMPLATE not in all_aligned_seqs_map:
            raise ValueError(f"Could not find main template '{MAIN_ALIGN_CODE_TEMPLATE}' in {ali_file_for_analysis}")
        if ALIGN_CODE_SEQUENCE not in all_aligned_seqs_map:
            raise ValueError(f"Could not find target sequence '{ALIGN_CODE_SEQUENCE}' in {ali_file_for_analysis}")

        main_template_aligned_seq = all_aligned_seqs_map[MAIN_ALIGN_CODE_TEMPLATE]
        target_aligned_seq = all_aligned_seqs_map[ALIGN_CODE_SEQUENCE]

        cde_line = utils.extract_cde_line(ali_file_for_analysis)

    except Exception as e:
        print(f"\n[FATAL ERROR] Failed to prepare PIR files. Terminating. Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    print(f"\n{'='*80}")
    print(f"[STEP 2] Experimental Residue Identification (from MAIN_PDB)")
    print(f"{'='*80}\n")

    experimental_residues = identify_experimental_residues(main_template_aligned_seq, target_aligned_seq)

    print(f"\n{'='*80}")
    print(f"[STEP 3] Loop Detection (Coil + Experimental Flanks 'C')")
    print(f"{'='*80}\n")

    try:
        all_coil_residues = get_coil_residues(SS2_FILE_PATH, sequence_full)
        print(f"[STEP 3] Total 'C' (Coil) residues detected: {len(all_coil_residues)}")

        experimental_ranges = group_ranges(list(experimental_residues))
        print(f"[STEP 3] Continuous experimental ranges detected: {experimental_ranges}")

        flank_residues = set()
        N = EXPERIMENTAL_FLANK_SIZE
        for start, end in experimental_ranges:
            for i in range(N):
                if start + i <= end:
                    flank_residues.add(start + i)
            for i in range(N):
                if end - i >= start:
                    flank_residues.add(end - i)

        print(f"[STEP 3] Flank size (per edge) to consider: {N} residues.")
        print(f"[STEP 3] Flank residues (potentially refinable) identified: {len(flank_residues)}")

        truly_fixed_residues = experimental_residues - flank_residues
        print(f"[STEP 3] Core residues (total): {len(experimental_residues)}")
        print(f"[STEP 3] Deep core residues (TRULY fixed): {len(truly_fixed_residues)}")

        refinable_residues = all_coil_residues - truly_fixed_residues
        refinable_flank = flank_residues.intersection(all_coil_residues)

        print(f"[STEP 3] {len(refinable_residues)} final residues selected for refinement (All 'C' - Deep Core).")
        print(f"         (This includes {len(refinable_flank)} 'C' residues from experimental flanks)")

        loop_ranges_to_refine = group_ranges(list(refinable_residues))

        range_strings = [f"[{s}-{e}]" if s != e else f"[{s}]" for s, e in loop_ranges_to_refine]
        print(f"\n[STEP 3] Final loop ranges to refine: {range_strings}")

    except Exception as e:
        print(f"\n[FATAL ERROR] Failed to detect loops. Terminating. Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    print(f"\n{'='*80}")
    print(f"[STEP 4] Homology Modeling (AutoModel with Fixed Regions)")
    print(f"{'='*80}\n")

    initial_models_names = homology_modeling.run_automodel(
        env,
        ali_file_for_automodel,
        job,
        experimental_residues
    )

    print(f"\n{'='*80}")
    print(f"[STEP 5] Loop Refinement (Loops 'C' + Flanks 'C')")
    print(f"{'='*80}\n")

    if initial_models_names:
        loop_refinement.run_loop_refinement(
            env, job, initial_models_names, loop_ranges_to_refine, truly_fixed_residues
        )

    print("[PARALLEL] All MODELLER processes completed.")

    print(f"\n{'='*80}")
    print(f"[STEP 6] Final Evaluation and Ranking")
    print(f"{'='*80}\n")

    final_ranking, best_final_model = utils.final_evaluation_and_ranking(env)

    if best_final_model:
        print(f"\n{'='*80}")
        print(f"FINAL RESULT")
        print(f"{'='*80}")
        print(f"\nThe highest quality model (most negative DOPE-HR) is:")
        print(f"  Name: {best_final_model['name']}")
        print(f"  Z-score: {best_final_model['DOPEHR Z-score']:.3f}")
        print(f"\n⚠️  IMPORTANT: This model maintains the FIXED experimental deep core.")
        print(f"              'C' loops and 'C' flanks (flank={N} res) were refined.")
        print(f"{'='*80}\n")


if __name__ == '__main__':
    main_workflow()
