#!/usr/bin/env python3
"""Main controller for the PRISM pipeline.

Orchestrates the protein modeling workflow by delegating tasks to specialised
modules.  Designed for external parallelisation (e.g. Slurm arrays).

Stages:
    prereq-cde:      Generate CDE alignment file.
    automodel:        Generate homology models (parallel split).
    rank-automodel:   Evaluate and rename initial homology models.
    loopmodel:        Refine loops on a specific input model.
    final-rank:       Final evaluation of all refined models.
"""

from __future__ import annotations

import argparse
import logging
import math
from pathlib import Path

from . import config, modeling_engine, psipred_client, utils

logger = logging.getLogger("PRISM.controller")


def main() -> None:  # noqa: PLR0912
    """Parse CLI arguments and dispatch the requested pipeline stage."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    parser = argparse.ArgumentParser(description="PRISM Pipeline Controller")

    parser.add_argument(
        "--stage",
        type=str,
        required=True,
        choices=[
            "prereq-cde",
            "automodel",
            "rank-automodel",
            "loopmodel",
            "final-rank",
        ],
        help="Pipeline stage to execute.",
    )

    parser.add_argument(
        "--job-id",
        type=int,
        help="Job ID for 'automodel' splitting (e.g., Slurm Array Task ID).",
    )

    parser.add_argument(
        "--input-model",
        type=str,
        help="Input PDB for 'loopmodel' stage (e.g., 'AUTO_1.pdb').",
    )

    parser.add_argument(
        "--phase",
        type=str,
        choices=["precalculation", "precomputed"],
        help="Phase for 'replicator' stage (e.g., 'precalculation', 'precomputed').",
    )

    parser.add_argument(
        "--input-mode",
        type=str,
        choices=["precalculation", "precomputed", "normal"],
        help="Input mode for 'automodel' stage.",
    )

    args = parser.parse_args()

    env, job = utils.setup_environment()

    if config.PERFORM_PSIPRED_PREDICTION:
        if not Path(config.SS2_FILE_PATH).exists():
            logger.info(
                "PSIPRED prediction enabled and SS2 file missing. Running client..."
            )
            psipred_client.run_psipred_request()
        else:
            logger.info(
                "PSIPRED prediction enabled, but SS2 file already exists. Skipping..."
            )

    if args.stage == "prereq-cde":
        utils.run_prereq_cde(env)

    elif args.stage == "automodel":
        if args.job_id is None:
            parser.error("--job-id is required for 'automodel' stage.")
        if not args.input_mode:
            parser.error("--input-mode is required for 'automodel' stage.")

        models_per_job = math.ceil(
            config.TOTAL_HOMOLOGY_MODELS / config.TOTAL_PARALLEL_JOBS
        )
        start = (args.job_id - 1) * models_per_job + 1
        end = min(args.job_id * models_per_job, config.TOTAL_HOMOLOGY_MODELS)
        input_mode = args.input_mode

        if config.EXECUTION_PARADIGM in (
            "prism-power",
            "precalculation",
            "precomputed",
        ):
            knowns, ali_file = utils.prepare_prism_power_files(phase=input_mode)
            _, _, _, fixed_res = utils.run_prerequisites(env)
        elif config.EXECUTION_PARADIGM == "normal":
            ali_file, _, _, fixed_res = utils.run_prerequisites(env)
            knowns = config.PDB_TEMPLATE_FILES_NAMES
        else:
            logger.error(
                "Unknown EXECUTION_PARADIGM: '%s'. "
                "Expected one of: prism-power, precalculation, precomputed, normal.",
                config.EXECUTION_PARADIGM,
            )
            return

        modeling_engine.run_automodel(
            env, ali_file, job, fixed_res, start, end, knowns, input_mode
        )

    elif args.stage == "rank-automodel":
        utils.run_rank_automodel_models(env)

    elif args.stage == "loopmodel":
        if not args.input_model:
            parser.error("--input-model is required for 'loopmodel' stage.")

        _, loop_ranges, truly_fixed, _ = utils.run_prerequisites(env)

        modeling_engine.run_loop_model(
            env,
            job,
            [args.input_model],
            loop_ranges,
            truly_fixed,
        )

    elif args.stage == "final-rank":
        utils.final_evaluation_and_ranking(env)

    logger.info("Pipeline stage '%s' completed successfully.", args.stage)


if __name__ == "__main__":
    main()
