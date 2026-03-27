#!/usr/bin/env python3
"""
Standalone Alignment Tool (PREREQ_CDE only).

Runs only the alignment step of the PRISM pipeline by directly importing
and calling the existing PRISM functions. No alignment logic is redefined
here — everything is delegated to PRISM.utils and PRISM.config.

Usage (from project root):
    python3 tools/run_alignment.py
"""

import sys
import shutil
import logging
import yaml
from pathlib import Path
from typing import Union, List

project_root = Path(__file__).resolve().parent.parent
tools_dir = Path(__file__).resolve().parent
config_yaml_path = project_root / "config.yaml"
sys.path.insert(0, str(project_root))

logger = logging.getLogger("tools.run_alignment")

def main() -> None:
    with open(config_yaml_path, "r") as f:
        raw_config = yaml.safe_load(f)

    original_value = raw_config.get("USE_MANUAL_ALIGNMENT", False)

    try:
        if original_value:
            logger.info("USE_MANUAL_ALIGNMENT is True — temporarily setting to False in config.yaml.")
            raw_config["USE_MANUAL_ALIGNMENT"] = False
            with open(config_yaml_path, "w") as f:
                yaml.dump(raw_config, f, default_flow_style=False, sort_keys=False)

        from PRISM import config
        from PRISM import psipred_client
        from PRISM.utils import setup_environment, run_prereq_cde

        env, _job = setup_environment()

        if config.PERFORM_PSIPRED_PREDICTION:
            if not Path(config.SS2_FILE_PATH).exists():
                logger.info("PSIPRED prediction enabled and SS2 file missing. Running client...")
                psipred_client.run_psipred_request()
            else:
                logger.info("PSIPRED prediction enabled, but SS2 file already exists. Skipping...")

        run_prereq_cde(env)

        results_dir = Path(config.MODELING_RESULTS_DIR)
        for ali_file in results_dir.glob("*.ali"):
            dest = tools_dir / ali_file.name
            shutil.move(str(ali_file), str(dest))
            logger.info(f"Moved {ali_file.name} -> {dest}")
            
            bullet_proof_ali(dest)

        logger.info("Alignment and bullet-proofing completed successfully.")

    finally:
        if original_value:
            logger.info("Restoring USE_MANUAL_ALIGNMENT to True in config.yaml.")
            raw_config["USE_MANUAL_ALIGNMENT"] = True
            with open(config_yaml_path, "w") as f:
                yaml.dump(raw_config, f, default_flow_style=False, sort_keys=False)

def bullet_proof_ali(ali_file: Path) -> None:
    """
    Standardizes PIR headers in an alignment file to use Modeller's universal shortcuts.
    
    Replaces residue ranges with 'FIRST:@:END:@' to allow automatic residue 
    detection by Modeller, making the alignment more robust to PDB numbering variants.
    """
    if not ali_file.exists():
        logger.error(f"Alignment file not found for bullet-proofing: {ali_file}")
        return

    lines = ali_file.read_text().splitlines()
    new_lines = []
    current_code = None

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
    """Prints a professional warning regarding BLK residues and alignment artifacts."""
    logger.warning("\n" + "!" * 80)
    logger.warning("ALIGNMENT POST-PROCESSING ADVISORY")
    logger.warning("!" * 80)
    logger.warning("When using BLK residues in Chain B (standard PRISM practice), MODELLER may")
    logger.warning("introduce redundant chain-break artifacts ('/') in the template sequence")
    logger.warning("and corresponding gaps ('-') in the target sequence.")
    logger.warning("")
    logger.warning("If these artifacts are present in your generated .ali files, please:")
    logger.warning("1. Manually remove the '/' from the main template sequence.")
    logger.warning("2. Manually remove the extra '-' from the other sequences.")
    logger.warning("!" * 80 + "\n")

if __name__ == "__main__":
    main()
    print_alignment_warning()
