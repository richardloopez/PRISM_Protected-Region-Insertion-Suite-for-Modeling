#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM (controller)

"""
Main Controller for the PRISM pipeline.

Orchestrates the protein modeling workflow by delegating tasks to specialized
modules. Designed for external parallelization (e.g., Slurm arrays).

Stages:
- prereq-cde:     Generate CDE alignment file.
- automodel:      Generate homology models (parallel split).
- rank-automodel: Evaluate and rename initial homology models.
- loop-model:    Refine loops on a specific input model.
- final-rank:     Final evaluation of all refined models.
"""

import sys
import argparse
import math 
from pathlib import Path

from . import config
from . import utils 
from . import psipred_client
from . import modeling_engine

def main():
    """Main function for the PRISM pipeline."""
    parser = argparse.ArgumentParser(description="PRISM Pipeline Controller")

    parser.add_argument(
        '--stage',
        type=str,
        required=True,
        choices=['prereq-cde', 'automodel', 'rank-automodel', 'loopmodel', 'final-rank'],
        help="Pipeline stage to execute"
    )

    parser.add_argument(
        '--job-id',
        type=int,
        help="Job ID for 'automodel' splitting (e.g., Slurm Array Task ID)."
    )

    parser.add_argument(
        '--input-model',
        type=str,
        help="Input PDB for 'loop-model' stage (e.g., 'AUTO_1.pdb')."
    )

    args = parser.parse_args()


    env, job = utils.setup_environment()

    if config.PERFORM_PSIPRED_PREDICTION:
        if not Path(config.SS2_FILE_PATH).exists():
            print(f"\n[ENVIRONMENT][main] PSIPRED prediction enabled and SS2 file missing. Running client...")
            try:
                psipred_client.run_psipred_request()
            except Exception as e:
                print(f"\n[ENVIRONMENT][main] PSIPRED prediction failed: {e}")
                sys.exit(1) 
        else:
            print(f"\n[ENVIRONMENT][main] PSIPRED prediction enabled, but SS2 file already exists. Skipping...")
    
    try:
        if args.stage == 'prereq-cde':
            utils.run_prereq_cde(env)
        
        elif args.stage == 'automodel':
            if not args.job_id:
                parser.error("--job-id is required for 'automodel' stage.") 
            models_per_job = math.ceil(config.NUM_MODELS_AUTO / config.NUM_SPLITTED_PROCESSES)
            start = (args.job_id - 1) * models_per_job + 1
            end = min(args.job_id * models_per_job, config.NUM_MODELS_AUTO)

            ali_file, _, _, fixed_res = utils.run_prerequisites(env)

            modeling_engine.run_automodel(env, ali_file, job, fixed_res, start, end)
        
        elif args.stage == 'rank-automodel':
            utils.run_rank_automodel_models(env)
        
        elif args.stage == 'loopmodel':
            if not args.input_model:
                parser.error("--input-model is required for 'loopmodel' stage.")
            
            _, loop_ranges, truly_fixed, _ = utils.run_prerequisites(env)
            
            modeling_engine.run_loop_model(
                env,
                job,
                [args.input_model],
                loop_ranges,
                truly_fixed
            )
        
        elif args.stage == 'final-rank':
            utils.final_evaluation_and_ranking(env)

    except Exception as e:
        print(f"\n[ENVIRONMENT][main] Pipeline stage '{args.stage}' failed.")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    print(f"\n[ENVIRONMENT][main] Pipeline stage '{args.stage}' completed successfully.")

if __name__ == '__main__':
    main()

            