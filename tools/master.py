#!/usr/bin/env python3
"""PRISM Master Pipeline Tool.

Orchestrates the full PRISM pipeline for a single system, from PDB preparation
to template unification, automating configuration updates and alignment
post-processing.
"""

from __future__ import annotations

import argparse
import logging
import re
import shutil
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT_BOOTSTRAP = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT_BOOTSTRAP) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT_BOOTSTRAP))

from PRISM import config  # noqa: E402
from tools import pdb_utils  # noqa: E402
from tools.pdb_utils import INPUT_DIR, OUTPUT_TOOLS_DIR, PROJECT_ROOT  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("tools.master")


def run_command(cmd: list[str], description: str) -> str:
    """Run a shell command and return its combined stdout and stderr.

    Args:
        cmd: List of command arguments.
        description: Human-readable description for logging.

    Returns:
        The combined standard output and error of the command.
    """
    logger.info("Running: %s", description)
    logger.debug("Command: %s", " ".join(cmd))
    result = subprocess.run(  # noqa: S603
        cmd,
        capture_output=True,
        text=True,
        check=True,
        cwd=str(PROJECT_ROOT),
    )
    return result.stdout + result.stderr


def get_suggested_radius(output: str) -> float | None:
    """Extract the suggested BLOCK_REPULSION_RADIUS from tool output.

    Args:
        output: The combined output from calc_block_distance.py.

    Returns:
        The suggested radius as a float, or None if not found.
    """
    match = re.search(r"SUGGESTED BLOCK_REPULSION_RADIUS:\s+([\d.]+)", output)
    return float(match.group(1)) if match else None


def fix_alignment_file(
    ali_path: Path, target_codes: list[str], ligand_res_count: int
) -> None:
    """Append ligand points to specific sequences in a PIR alignment file.

    Args:
        ali_path: Path to the .ali file.
        target_codes: List of alignment codes (e.g. ['9TKV_prism_prep.pdb', 'FullSeq'])
                      that should receive the ligand points.
        ligand_res_count: Number of points ('.') to add.
    """
    if ligand_res_count <= 0:
        return
    if not ali_path.exists():
        raise FileNotFoundError(f"Alignment file not found: {ali_path}")

    lines = ali_path.read_text().splitlines()
    new_lines = []
    current_code = None
    modifying = False

    for line in lines:
        if line.startswith(">P1;"):
            current_code = line[4:].strip()
            modifying = current_code in target_codes
            new_lines.append(line)
        elif modifying and line.endswith("*"):
            base_seq = line.rstrip("*")
            fixed_seq = f"{base_seq}/" + ("." * ligand_res_count) + "*"
            new_lines.append(fixed_seq)
            modifying = False
        else:
            new_lines.append(line)

    ali_path.write_text("\n".join(new_lines) + "\n")


def setup_input_files(exp_path: Path, pred_path: Path) -> Path:
    """Ensure input PDBs are in the input directory.

    Returns:
        The path to the input directory.
    """
    INPUT_DIR.mkdir(exist_ok=True)

    if not (INPUT_DIR / exp_path.name).exists():
        shutil.copy(exp_path, INPUT_DIR / exp_path.name)
    if not (INPUT_DIR / pred_path.name).exists():
        shutil.copy(pred_path, INPUT_DIR / pred_path.name)

    return INPUT_DIR


def step1_initial_config(exp_name: str, pred_name: str, overlap: int) -> str:
    """Step 1: Perform initial config.yaml update."""
    exp_stem = Path(exp_name).stem
    suffix = config.settings.ALIGN_CODE_SEQUENCE
    ali_name = f"{exp_stem}_prism_prep.pdb_{suffix}.ali"

    logger.info("Step 1: Initial config.yaml update")
    config.settings.MANUAL_ALIGNMENT_BASENAME = ali_name
    config.settings.PDB_TEMPLATE_FILES_NAMES = [exp_name, pred_name]
    config.settings.MOBILE_FLANK_RESIDUES = overlap
    config.settings.save_settings()
    return ali_name


def step2_prep_pdbs(
    input_dir: Path, exp_name: str, prot: str, lig: str, ptm_args: list[str]
) -> Path:
    """Step 2: Run prep_prism_pdb.py prep."""
    logger.info("Step 2: Preparing PDBs")
    prep_cmd = [
        sys.executable,
        "tools/prep_prism_pdb.py",
        "prep",
        str(input_dir / exp_name),
        prot,
        lig,
        *ptm_args,
    ]
    run_command(prep_cmd, "prep_prism_pdb.py prep")
    return input_dir / f"{Path(exp_name).stem}_prism_prep.pdb"


def step3_4_radius(prep_exp_pdb: Path) -> None:
    """Steps 3 and 4: Calculate and update block repulsion radius."""
    logger.info("Step 3: Calculating block distances")
    dist_cmd = [
        sys.executable,
        "tools/calc_block_distance.py",
        str(prep_exp_pdb),
        "--protein_chain",
        "A",
        "--blk_chain",
        "B",
    ]
    dist_output = run_command(dist_cmd, "calc_block_distance.py")
    suggested_radius = get_suggested_radius(dist_output)

    if suggested_radius is not None:
        logger.info("Step 4: Updating BLOCK_REPULSION_RADIUS to %s", suggested_radius)
        config.settings.BLOCK_REPULSION_RADIUS = suggested_radius
        config.settings.save_settings()
    else:
        logger.warning("Could not find suggested radius in output.")


def step5_alignment(prep_exp_pdb: Path, prep_pred_pdb: Path) -> None:
    """Step 5: Run run_alignment.py with temporary ligand removal."""
    logger.info("Step 5: Running alignment with temporary ligand removal")

    prep_exp_name = prep_exp_pdb.name
    prep_pred_name = prep_pred_pdb.name
    config.settings.PDB_TEMPLATE_FILES_NAMES = [prep_exp_name, prep_pred_name]
    config.settings.PRISM_POWER_SETTINGS.PRECALCULATION = {
        prep_exp_name: 1000,
        prep_pred_name: 1,
    }
    config.settings.PRISM_POWER_SETTINGS.PRECOMPUTED = {
        prep_exp_name: 1,
        prep_pred_name: 1,
    }
    config.settings.save_settings()

    backup_pdb = prep_exp_pdb.with_suffix(".pdb.bak")
    shutil.copy(prep_exp_pdb, backup_pdb)

    try:
        pdb_utils.remove_chain(backup_pdb, prep_exp_pdb, "B")
        run_command([sys.executable, "tools/run_alignment.py"], "run_alignment.py")
    finally:
        shutil.move(backup_pdb, prep_exp_pdb)


def main() -> None:
    """Execute the master PRISM pipeline."""
    parser = argparse.ArgumentParser(description="PRISM Master Pipeline Tool")
    parser.add_argument("--experimental_pdb", required=True, help="Experimental PDB")
    parser.add_argument("--prediction_pdb", required=True, help="Prediction PDB")
    parser.add_argument("--protein_chains", required=True, help="e.g., 'A'")
    parser.add_argument("--ligand_chains", required=True, help="e.g., 'B'")
    parser.add_argument("--ptm_args", nargs="*", default=[], help="Optional PTMs")
    parser.add_argument("--overlap", type=int, default=3, help="Overlap (default: 3)")

    args = parser.parse_args()
    exp_path, pred_path = Path(args.experimental_pdb), Path(args.prediction_pdb)

    # Initialize configuration settings
    config.settings = config.load_settings()

    # 0. Setup
    input_dir = setup_input_files(exp_path, pred_path)
    exp_name, pred_name = exp_path.name, pred_path.name

    # 1. Initial Config
    step1_initial_config(exp_name, pred_name, args.overlap)

    # 2. Prep PDBs
    prep_exp_pdb = step2_prep_pdbs(
        input_dir, exp_name, args.protein_chains, args.ligand_chains, args.ptm_args
    )
    prep_pred_pdb = step2_prep_pdbs(input_dir, pred_name, "A", "None", [])

    # 3-4. Distance & Radius
    step3_4_radius(prep_exp_pdb)

    # 5. Alignment
    step5_alignment(prep_exp_pdb, prep_pred_pdb)

    # 6. Post-process alignment
    logger.info("Step 6: Restoring ligand residues in alignment sequence")
    ligand_res_count = pdb_utils.count_residues(prep_exp_pdb, "B")

    exp_stem = exp_path.stem
    pattern = f"{exp_stem}_prism_prep.pdb_*.ali"
    matches = list(OUTPUT_TOOLS_DIR.glob(pattern))

    if not matches:
        raise FileNotFoundError(
            f"Could not find alignment file matching {pattern} in {OUTPUT_TOOLS_DIR}"
        )

    ali_path = matches[0]
    target_codes = [prep_exp_pdb.name, config.settings.ALIGN_CODE_SEQUENCE]
    fix_alignment_file(ali_path, target_codes, ligand_res_count)

    # 7. Unify templates
    logger.info("Step 7: Unifying templates")
    unify_cmd = [
        sys.executable,
        "tools/unify_templates.py",
        str(ali_path),
        "--overlap",
        str(args.overlap),
    ]
    run_command(unify_cmd, "unify_templates.py")

    # 8. Final Config update
    suffix = config.settings.ALIGN_CODE_SEQUENCE
    unified_ali = f"{exp_stem}_prism_prep.pdb_{suffix}_unified.ali"
    logger.info("Step 8: Setting final alignment to %s", unified_ali)
    prep_unified_exp_pdb = prep_exp_pdb.stem + "_unified.pdb"
    prep_unified_pred_pdb = prep_pred_pdb.stem + "_unified.pdb"
    config.settings.MANUAL_ALIGNMENT_BASENAME = unified_ali
    config.settings.PDB_TEMPLATE_FILES_NAMES = [
        prep_unified_exp_pdb,
        prep_unified_pred_pdb,
    ]
    config.settings.PRISM_POWER_SETTINGS.PRECALCULATION = {
        prep_unified_exp_pdb: 1000,
        prep_unified_pred_pdb: 1,
    }
    config.settings.PRISM_POWER_SETTINGS.PRECOMPUTED = {
        prep_unified_exp_pdb: 1,
        prep_unified_pred_pdb: 1,
    }
    config.settings.save_settings()

    logger.info("Master pipeline completed successfully!")


if __name__ == "__main__":
    main()
