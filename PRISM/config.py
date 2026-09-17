#!/usr/bin/env python3
"""Configuration module for the PRISM pipeline.

Defines the ``PrismConfig`` Pydantic model that loads and validates all
user-facing parameters from ``config.yaml``, computes derived paths, and
exposes them as module-level attributes via a lazy ``__getattr__`` hook.

Design note — module-level injection:
    This module uses a ``__getattr__`` override so that downstream code can
    write ``from PRISM import config; config.CHAIN_ID`` without explicitly
    calling ``load_settings()``.  On first attribute access the YAML
    file is loaded once, and every validated field is promoted to a module
    global.  While unconventional, this pattern keeps the ergonomic benefit
    of attribute-style access while guaranteeing Pydantic validation.
"""

from __future__ import annotations

import logging
import re
import threading
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import (
    BaseModel,
    EmailStr,
    Field,
    computed_field,
    field_validator,
    model_validator,
)

logger = logging.getLogger("PRISM.config")

# ============================================================================
# PYDANTIC CONFIGURATION MODELS
# ============================================================================


class PrismPowerConfig(BaseModel):
    """Per-phase template replica weights for the *prism-power* paradigm."""

    PRECALCULATION: dict[str, int]
    PRECOMPUTED: dict[str, int]


class PrismConfig(BaseModel):
    """Central configuration model for the PRISM pipeline.

    All parameters are loaded from ``config.yaml`` and validated here.
    Computed fields derive full paths from basenames + directory constants.
    """

    # -- 1. Sequence & Identity -----------------------------------------------
    ALIGN_CODE_SEQUENCE: str
    CHAIN_ID: str
    BLK_CHAIN_ID: str | None = None

    # -- 2. Input Files (basenames only) --------------------------------------
    FASTA_FILE_BASENAME: str
    SS2_FILE_BASENAME: str
    MANUAL_ALIGNMENT_BASENAME: str
    CUSTOM_INIFILE_BASENAME: str
    CUSTOM_RSRFILE_BASENAME: str
    PDB_TEMPLATE_FILES_NAMES: list[str]

    # -- 3. Modeling Parameters -----------------------------------------------
    TOTAL_HOMOLOGY_MODELS: int = Field(gt=0)
    TOP_MODELS_FOR_REFINEMENT: int = Field(ge=0)
    LOOP_MODELS_PER_TARGET: int = Field(gt=0)
    NUM_BEST_FINAL_MODELS: int | float

    # -- 4. Experimental Region Control ---------------------------------------
    BLOCK_REPULSION_RADIUS: float
    MOBILE_FLANK_RESIDUES: int
    REFINE_FLANKS_DURING_AUTOMODEL: bool
    USE_MANUAL_OPTIMIZATION_SELECTION: bool
    MANUAL_OPTIMIZATION_RESIDUES: list[int]
    USE_MANUAL_FIXATION_SELECTION: bool
    MANUAL_FIXATION_RESIDUES: list[int]

    # -- 5. Execution Engine --------------------------------------------------
    EXECUTION_PARADIGM: Literal[
        "normal", "precalculation", "precomputed", "prism-power"
    ]
    MODELLER_CORES: int = Field(gt=0)
    TOTAL_PARALLEL_JOBS: int = Field(gt=0)
    PRISM_POWER_SETTINGS: PrismPowerConfig | None = None

    # -- 6. External Services -------------------------------------------------
    PERFORM_PSIPRED_PREDICTION: bool
    PSIPRED_EMAIL: EmailStr
    PSIPRED_POLL_INTERVAL: int
    USE_MANUAL_ALIGNMENT: bool

    # -- 7. Directories (internal constants) ----------------------------------
    INPUT_DIR_NAME: str = "input"
    MODELING_RESULTS_DIR_NAME: str = "modeling_results"
    PSIPRED_RESULTS_DIR_NAME: str = "psipred_results"

    # ---- Validators ---------------------------------------------------------

    @model_validator(mode="after")
    def check_num_models_to_refine(self) -> PrismConfig:
        """Ensure refinement count does not exceed total models."""
        if self.TOP_MODELS_FOR_REFINEMENT > self.TOTAL_HOMOLOGY_MODELS:
            raise ValueError("Refinement count exceeds total models.")
        return self

    @field_validator("NUM_BEST_FINAL_MODELS", mode="before")
    @classmethod
    def parse_num_best_final_models(cls, val: Any) -> int | float:
        """Parse NUM_BEST_FINAL_MODELS, converting 'inf' to float('inf')."""
        if isinstance(val, str) and val.lower() == "inf":
            return float("inf")
        try:
            int_val = int(val)
        except (ValueError, TypeError) as err:
            raise ValueError(
                "NUM_BEST_FINAL_MODELS must be a positive int or 'inf'"
            ) from err
        else:
            if int_val < 1:
                raise ValueError(
                    "NUM_BEST_FINAL_MODELS must be a positive int or 'inf'"
                )
            return int_val

    @field_validator(
        "MANUAL_OPTIMIZATION_RESIDUES", "MANUAL_FIXATION_RESIDUES", mode="before"
    )
    @classmethod
    def expand_residue_ranges(cls, v: Any) -> list[int]:
        """Expand range notation (e.g. ``'10-20'``) into individual residue indices.

        Args:
            v: Raw value from YAML — may contain integers or ``'start-end'``
                strings.

        Returns:
            Flat list of integer residue indices.
        """
        if not isinstance(v, list):
            return v

        expanded: list[int] = []
        for item in v:
            if isinstance(item, str) and "-" in item:
                start_str, end_str = item.split("-", 1)
                start, end = int(start_str.strip()), int(end_str.strip())
                if start > end:
                    logger.warning(
                        "Residue range '%s' is specified in reverse order "
                        "(%d > %d); swapping to %d-%d.",
                        item,
                        start,
                        end,
                        end,
                        start,
                    )
                    start, end = end, start
                expanded.extend(range(start, end + 1))
            else:
                expanded.append(int(item))
        return expanded

    # ---- Computed fields (derived paths) ------------------------------------

    @computed_field
    @property
    def PROJECT_ROOT(self) -> str:  # noqa: N802
        """Absolute path to the repository root directory."""
        return str(Path(__file__).resolve().parent.parent)

    @computed_field
    @property
    def INPUT_DIR(self) -> str:  # noqa: N802
        """Absolute path to the input directory."""
        return str(Path(self.PROJECT_ROOT) / self.INPUT_DIR_NAME)

    @computed_field
    @property
    def MODELING_RESULTS_DIR(self) -> str:  # noqa: N802
        """Absolute path to the modeling results directory."""
        return str(Path(self.PROJECT_ROOT) / self.MODELING_RESULTS_DIR_NAME)

    @computed_field
    @property
    def PSIPRED_RESULTS_DIR(self) -> str:  # noqa: N802
        """Absolute path to the PSIPRED results directory."""
        return str(Path(self.PROJECT_ROOT) / self.PSIPRED_RESULTS_DIR_NAME)

    @computed_field
    @property
    def CUSTOM_INIFILE_PATH(self) -> str:  # noqa: N802
        """Full path to the precomputed initial structure file."""
        return str(Path(self.INPUT_DIR) / self.CUSTOM_INIFILE_BASENAME)

    @computed_field
    @property
    def CUSTOM_RSRFILE_PATH(self) -> str:  # noqa: N802
        """Full path to the precomputed restraint file."""
        return str(Path(self.INPUT_DIR) / self.CUSTOM_RSRFILE_BASENAME)

    @computed_field
    @property
    def FASTA_FILE_PATH(self) -> str:  # noqa: N802
        """Full path to the target FASTA sequence file."""
        return str(Path(self.INPUT_DIR) / self.FASTA_FILE_BASENAME)

    @computed_field
    @property
    def SS2_FILE_PATH(self) -> str:  # noqa: N802
        """Full path to the PSIPRED secondary structure prediction file."""
        return str(Path(self.INPUT_DIR) / self.SS2_FILE_BASENAME)

    @computed_field
    @property
    def MANUAL_ALIGNMENT_CDE_BASENAME(self) -> str:  # noqa: N802
        """Basename for the CDE-annotated manual alignment file."""
        p = Path(self.MANUAL_ALIGNMENT_BASENAME)
        return f"{p.stem}_cde{p.suffix}"

    @computed_field
    @property
    def MANUAL_ALIGNMENT_FILE(self) -> str:  # noqa: N802
        """Full path to the manual alignment file."""
        return str(Path(self.INPUT_DIR) / self.MANUAL_ALIGNMENT_BASENAME)

    @computed_field
    @property
    def MANUAL_ALIGNMENT_CDE_FILE(self) -> str:  # noqa: N802
        """Full path to the CDE-annotated manual alignment file."""
        return str(Path(self.INPUT_DIR) / self.MANUAL_ALIGNMENT_CDE_BASENAME)

    @computed_field
    @property
    def PDB_TEMPLATE_FILES_PATHS(self) -> list[str]:  # noqa: N802
        """Full paths to all PDB template files."""
        return [
            str(Path(self.INPUT_DIR) / name) for name in self.PDB_TEMPLATE_FILES_NAMES
        ]

    @computed_field
    @property
    def MAIN_PDB_TEMPLATE_PATH(self) -> str:  # noqa: N802
        """Full path to the main (first) PDB template."""
        return self.PDB_TEMPLATE_FILES_PATHS[0]

    @computed_field
    @property
    def MAIN_ALIGN_CODE_TEMPLATE(self) -> str:  # noqa: N802
        """Alignment code for the main (first) template."""
        return self.PDB_TEMPLATE_FILES_NAMES[0]

    @computed_field
    @property
    def ALIGNMENT_FILE(self) -> str:  # noqa: N802
        """Full path to the auto-generated alignment file."""
        return str(
            Path(self.MODELING_RESULTS_DIR)
            / f"{self.MAIN_ALIGN_CODE_TEMPLATE}_{self.ALIGN_CODE_SEQUENCE}.ali"
        )

    @computed_field
    @property
    def ALIGNMENT_CDE_FILE(self) -> str:  # noqa: N802
        """Full path to the auto-generated CDE alignment file."""
        return str(
            Path(self.MODELING_RESULTS_DIR)
            / f"{self.MAIN_ALIGN_CODE_TEMPLATE}_{self.ALIGN_CODE_SEQUENCE}_cde.ali"
        )

    @computed_field
    @property
    def FINAL_RANKING_CSV(self) -> str:  # noqa: N802
        """Full path to the final model ranking CSV."""
        return str(Path(self.MODELING_RESULTS_DIR) / "final_ranking.csv")

    # ---- Persistence -------------------------------------------

    def save_settings(self, yaml_path: str = "config.yaml") -> None:  # noqa: PLR0912
        """Save settings while PRESERVING comments and structure.

        Line-by-line update to maintain the YAML structure with comments and sections.
        """
        full_yaml_path = Path(self.PROJECT_ROOT) / yaml_path
        computed_names = set(PrismConfig.model_computed_fields)
        data = self.model_dump(exclude=computed_names)
        data["PSIPRED_EMAIL"] = str(data["PSIPRED_EMAIL"])

        if not full_yaml_path.exists():
            # Fallback for new files — use standard dump
            with full_yaml_path.open("w", encoding="utf-8") as fh:
                yaml.dump(data, fh, default_flow_style=False, sort_keys=False)
            return

        lines = full_yaml_path.read_text(encoding="utf-8").splitlines()
        new_lines: list[str] = []
        i = 0
        while i < len(lines):
            line = lines[i]
            # Match top-level keys
            match = re.match(r"^([A-Z_]+):", line)
            if match:
                key = match.group(1)
                if key in data:
                    val = data[key]
                    # Handle complex types (lists and nested dicts)
                    if isinstance(val, (dict, list)) and val:
                        # Preserve comment on the parent key line if it exists
                        parent_comment = ""
                        if "#" in line:
                            parent_comment = "  #" + line.split("#", 1)[1]

                        # Generate the YAML block for this key
                        dump_data = {key: val}
                        dumped_block = yaml.dump(
                            dump_data, default_flow_style=False, sort_keys=False
                        )
                        block_lines = dumped_block.splitlines()
                        if parent_comment and block_lines:
                            block_lines[0] = f"{block_lines[0]}{parent_comment}"

                        new_lines.extend(block_lines)

                        # Skip the original block in the source file
                        i += 1
                        while i < len(lines) and (
                            not lines[i].strip()
                            or lines[i].startswith(" ")
                            or lines[i].startswith("-")
                        ):
                            # Stop if we hit a blank line followed by a new key
                            # or if we hit a line that looks like a new key.
                            if not lines[i].strip() and (
                                i + 1 < len(lines)
                                and re.match(r"^[A-Z_]+:", lines[i + 1].strip())
                            ):
                                break
                            i += 1
                        continue

                    # Handle scalar types
                    if isinstance(val, bool):
                        val_str = "true" if val else "false"
                    elif val is None:
                        val_str = "null"
                    else:
                        val_str = str(val)

                    # Preserve trailing comment if present
                    comment = ""
                    if "#" in line:
                        comment = "  #" + line.split("#", 1)[1]

                    new_lines.append(f"{key}: {val_str}{comment}")
                    i += 1
                    continue

            new_lines.append(line)
            i += 1

        full_yaml_path.write_text("\n".join(new_lines) + "\n", encoding="utf-8")


# ============================================================================
# INSTANTIATION LOGIC
# ============================================================================


def load_settings(yaml_path: str | None = None) -> PrismConfig:
    """Load and validate the YAML configuration file.

    Args:
        yaml_path: Optional path to the YAML file.  When *None*, defaults to
            ``<project_root>/config.yaml``.  Relative paths are resolved
            against the project root.

    Returns:
        A fully validated ``PrismConfig`` instance.

    Raises:
        FileNotFoundError: If the configuration file does not exist.
    """
    project_root = Path(__file__).resolve().parent.parent

    if yaml_path is None:
        resolved = project_root / "config.yaml"
    else:
        resolved = Path(yaml_path)
        if not resolved.is_absolute():
            resolved = project_root / resolved

    if not resolved.exists():
        raise FileNotFoundError(f"Configuration missing: {resolved}")

    logger.info("Loading settings from %s", resolved)

    with resolved.open(encoding="utf-8") as fh:
        raw_data = yaml.safe_load(fh)
    return PrismConfig(**raw_data)


settings: PrismConfig | None = None
_CONFIG_LOCK = threading.Lock()


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Sequence cache (defined here to avoid circular imports from utils).
SEQUENCE_FULL: str | None = None


def read_fasta_sequence(file_path: str) -> str:
    """Read a FASTA file and return the concatenated sequence.

    Args:
        file_path: Path to the FASTA file.

    Returns:
        The protein sequence as a single string (header lines excluded).
    """
    with Path(file_path).open(encoding="utf-8") as fh:
        parts = [line.strip() for line in fh if not line.startswith(">")]
    return "".join(parts)


def get_sequence() -> str:
    """Return the target sequence, loading it from disk on first call.

    Returns:
        The full protein sequence string.
    """
    global SEQUENCE_FULL, settings  # noqa: PLW0603
    if settings is None:
        settings = load_settings()
    if SEQUENCE_FULL is None:
        SEQUENCE_FULL = read_fasta_sequence(settings.FASTA_FILE_PATH)
    return SEQUENCE_FULL


# ============================================================================
# MODULE-LEVEL ATTRIBUTE INJECTION
# ============================================================================
# On first access of any config attribute (e.g. ``config.CHAIN_ID``), the
# YAML file is loaded and every validated field is promoted to a module
# global.  This avoids requiring an explicit ``load_settings()`` call in
# every importing module.


def __getattr__(name: str) -> Any:
    global settings  # noqa: PLW0603
    if settings is None:
        with _CONFIG_LOCK:
            if settings is None:
                settings = load_settings()
                for field in settings.model_fields:
                    globals()[field] = getattr(settings, field)
                for computed_name in PrismConfig.model_computed_fields:
                    globals()[computed_name] = getattr(settings, computed_name)
    if name in globals():
        return globals()[name]
    raise AttributeError(f"Missing attribute: {name}")


__all__ = [
    "SEQUENCE_FULL",
    "PrismConfig",
    "get_sequence",
    "load_settings",
    "read_fasta_sequence",
    "settings",
]

if __name__ == "__main__":
    if settings is None:
        settings = load_settings()
    print(settings.model_dump_json(indent=4))  # noqa: T201
