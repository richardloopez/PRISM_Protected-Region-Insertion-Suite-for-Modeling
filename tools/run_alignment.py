#!/usr/bin/env python3
"""Alignment execution tool for the PRISM pipeline.

Automates the generation of sequence alignments (PIR format) using
MODELLER's salign function. This script handles the temporary override
of the ``USE_MANUAL_ALIGNMENT`` configuration flag, ensures PSIPRED
predictions are available if needed, and standardises output PIR
headers for robust modeling.

Execution:
    Run from the project root:
    ``python3 tools/run_alignment.py``

Notes:
    This tool depends on PRISM.config for paths and flags. It will
    temporarily set ``USE_MANUAL_ALIGNMENT`` to False in config.yaml and
    True to ensure automatic alignment execution, restoring it afterwards.
"""

import logging
import shutil
import sys
from pathlib import Path

# Ensure project root is in sys.path for PRISM imports
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Ensure tools directory is in sys.path for pdb_utils
TOOLS_DIR = PROJECT_ROOT / "tools"
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("tools.run_alignment")

try:
    import pdb_utils

    from PRISM import config, psipred_client, utils
except ImportError:
    logger.exception(
        "Could not import PRISM modules. Ensure you are running from the "
        "project root or tools/ directory."
    )
    sys.exit(1)


CONFIG_YAML_PATH = PROJECT_ROOT / "config.yaml"


def bullet_proof_ali(ali_file: Path) -> None:
    """Standardise PIR headers in an alignment file.

    Replaces residue ranges with ``FIRST:@:END:@`` to allow automatic
    residue detection by MODELLER, making the alignment more robust to
    PDB numbering variants.

    Args:
        ali_file: Path to the ``.ali`` file to bullet-proof.
    """
    if not ali_file.exists():
        logger.error("Alignment file not found for bullet-proofing: %s", ali_file)
        return

    lines = ali_file.read_text().splitlines()
    new_lines: list[str] = []
    current_code: str | None = None

    for line in lines:
        stripped = line.strip()
        if stripped.startswith(">P1;"):
            current_code = stripped.split(";")[1].strip()
            new_lines.append(stripped)
        elif current_code and stripped.startswith("structure"):
            new_lines.append(f"structure:{current_code}:FIRST:@:END:@::::")
        elif current_code and stripped.startswith("sequence"):
            new_lines.append(f"sequence:{current_code}:FIRST:@:END:@::::")
        else:
            new_lines.append(line)

    ali_file.write_text("\n".join(new_lines) + "\n")


def print_alignment_warning() -> None:
    """Print a warning regarding BLK residues and alignment artifacts."""
    logger.warning("\n" + "!" * 80)
    logger.warning("ALIGNMENT POST-PROCESSING ADVISORY")
    logger.warning("!" * 80)
    logger.warning(
        "When using BLK residues in Chain B (standard PRISM practice), MODELLER may"
    )
    logger.warning("introduce redundant gaps ('-')")
    logger.warning("")
    logger.warning(
        "If these artifacts are present in your generated .ali files, please:"
    )
    logger.warning("Manually remove the extra '-' from the sequences.")
    logger.warning("!" * 80 + "\n")


def main() -> None:
    """Entry point for the standalone alignment tool."""
    original_manual_val = config.USE_MANUAL_ALIGNMENT

    try:
        if original_manual_val:
            logger.info(
                "USE_MANUAL_ALIGNMENT is True — temporarily setting "
                "to False in config.yaml and in-memory."
            )
            config.settings.USE_MANUAL_ALIGNMENT = False
            config.settings.save_settings()
            config.USE_MANUAL_ALIGNMENT = False

        env, _job = utils.setup_environment()
        if config.PERFORM_PSIPRED_PREDICTION:
            if not Path(config.SS2_FILE_PATH).exists():
                logger.info(
                    "PSIPRED prediction enabled and SS2 file missing. Running client..."
                )
                psipred_client.run_psipred_request()
            else:
                logger.info(
                    "PSIPRED prediction enabled, but SS2 file already exists. "
                    "Skipping..."
                )

        utils.run_prereq_cde(env)

        results_dir = Path(config.MODELING_RESULTS_DIR)
        output_tools_dir = pdb_utils.OUTPUT_TOOLS_DIR
        output_tools_dir.mkdir(parents=True, exist_ok=True)

        for ali_file in results_dir.glob("*.ali"):
            dest = output_tools_dir / ali_file.name
            shutil.move(str(ali_file), str(dest))
            logger.info("Moved %s -> %s", ali_file.name, dest)

            bullet_proof_ali(dest)

        logger.info("Alignment and bullet-proofing completed successfully.")

    finally:
        if original_manual_val:
            logger.info("Restoring USE_MANUAL_ALIGNMENT to True in config.yaml.")
            config.settings.USE_MANUAL_ALIGNMENT = True
            config.settings.save_settings()
            config.USE_MANUAL_ALIGNMENT = True


if __name__ == "__main__":
    try:
        main()
    finally:
        print_alignment_warning()
