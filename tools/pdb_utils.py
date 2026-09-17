#!/usr/bin/env python3
"""Shared utilities for PRISM tool scripts.

Provides the ``PDBAtom`` class for parsing and writing PDB-format atom
records, and the ``copy_to_input`` helper used by multiple tools to
stage output files into the pipeline's ``input/`` directory.
"""

from __future__ import annotations

import logging
import shutil
import sys
from pathlib import Path

logger = logging.getLogger("tools.pdb_utils")

# PDB Format Constants
PDB_COORD_LINE_MIN_LEN = 54
PDB_OCCUPANCY_END = 60
PDB_TEMP_FACTOR_END = 66
PDB_ELEMENT_START = 76
PDB_ATOM_NAME_FOUR_CHAR_LEN = 4


# Dynamic load of PRISM.config to obtain PROJECT_ROOT
def _get_project_root() -> Path:
    """Load PRISM.config dynamically and return the project_root Path."""
    base_path = Path(__file__).resolve().parent.parent
    try:
        if str(base_path) not in sys.path:
            sys.path.insert(0, str(base_path))

        from PRISM import config  # noqa: PLC0415

        return Path(config.PROJECT_ROOT)
    except Exception as exc:
        logger.warning("Failed to load PRISM.config: %s", exc)
        return base_path


# Path Constants
PROJECT_ROOT = _get_project_root()
INPUT_DIR = PROJECT_ROOT / "input"
OUTPUT_TOOLS_DIR = PROJECT_ROOT / "output_tools"


class PDBAtom:
    """Parser and writer for a single PDB ATOM/HETATM record.

    Parses fixed-width columns according to the PDB file format
    specification and can regenerate a standards-compliant line via
    :meth:`to_pdb_line`.

    Args:
        line: A raw PDB-format line (at least 54 characters).
    """

    def __init__(self, line: str) -> None:
        """Parse a PDB line into its constituent fields."""
        if len(line) < PDB_COORD_LINE_MIN_LEN:
            raise ValueError(
                f"PDB line too short ({len(line)} chars, "
                f"need ≥{PDB_COORD_LINE_MIN_LEN}): {line!r}"
            )
        self.line = line
        self.record_type: str = line[0:6].strip()
        try:
            self.serial: int = int(line[6:11])
        except ValueError:
            self.serial = 0
        self.name: str = line[12:16]
        self.alt_loc: str = line[16]
        self.res_name: str = line[17:20].strip()
        self.chain_id: str = line[21]
        try:
            self.res_seq: int = int(line[22:26])
        except ValueError:
            self.res_seq = 0
        self.i_code: str = line[26]
        self.x: float = float(line[30:38])
        self.y: float = float(line[38:46])
        self.z: float = float(line[46:54])
        self.occ: float = (
            float(line[PDB_COORD_LINE_MIN_LEN:PDB_OCCUPANCY_END])
            if len(line) > PDB_COORD_LINE_MIN_LEN
            and line[PDB_COORD_LINE_MIN_LEN:PDB_OCCUPANCY_END].strip()
            else 1.00
        )
        self.temp: float = (
            float(line[PDB_OCCUPANCY_END:PDB_TEMP_FACTOR_END])
            if (
                len(line) > PDB_OCCUPANCY_END
                and line[PDB_OCCUPANCY_END:PDB_TEMP_FACTOR_END].strip()
            )
            else 0.00
        )
        self.element: str = (
            line[PDB_ELEMENT_START:78].strip() if len(line) > PDB_ELEMENT_START else ""
        )

    def to_pdb_line(self) -> str:
        """Render this atom as a PDB-format fixed-width line.

        Returns:
            An 80-character PDB record string (no trailing newline).
        """
        if len(self.name) == PDB_ATOM_NAME_FOUR_CHAR_LEN:
            name_str = f"{self.name}"
        elif self.name[0].isdigit():
            name_str = f"{self.name:<4}"
        else:
            name_str = f" {self.name:<3}"
        return (
            f"{self.record_type:<6}{self.serial:>5} {name_str:4}"
            f"{self.alt_loc}{self.res_name:>3} "
            f"{self.chain_id}{self.res_seq:>4}{self.i_code}   "
            f"{self.x:>8.3f}{self.y:>8.3f}{self.z:>8.3f}"
            f"{self.occ:>6.2f}{self.temp:>6.2f}          {self.element:>2}  "
        )


def copy_to_input(filename: str | Path) -> None:
    """Copy a file to the pipeline's ``input/`` directory.

    Searches for the ``input/`` directory relative to the repository root.
    Falls back to a local ``input/`` in the current working directory if
    the project structure is not found.

    Args:
        filename: Path to the file to copy.
    """
    filename = Path(filename)

    if INPUT_DIR.is_dir():
        target = INPUT_DIR
    else:
        local_input = Path.cwd() / "input"
        if local_input.is_dir():
            target = local_input
        else:
            logger.warning(
                "[COPY] No input/ directory found. Skipping copy of %s.", filename.name
            )
            return

    try:
        shutil.copy2(str(filename), str(target / filename.name))
        logger.info("[COPY] Copied %s to %s", filename.name, target)
    except Exception as exc:
        logger.warning("[COPY] Could not copy %s to %s: %s", filename.name, target, exc)


def get_chains(pdb_path: Path) -> list[str]:
    """Extract all unique chain IDs from a PDB file.

    Args:
        pdb_path: Path to the PDB file.

    Returns:
        Sorted list of unique chain IDs.
    """
    chains = set()
    with pdb_path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21].strip()
                if chain:
                    chains.add(chain)
    return sorted(chains)


def get_residues(pdb_path: Path, chain_id: str) -> list[tuple[str, str]]:
    """Extract all unique residues (res_seq, i_code) for a specific chain.

    Args:
        pdb_path: Path to the PDB file.
        chain_id: The chain ID to extract residues for.

    Returns:
        List of unique (res_seq, i_code) tuples.
    """
    residues = set()
    with pdb_path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                atom = PDBAtom(line)
                if atom.chain_id == chain_id:
                    residues.add((str(atom.res_seq), atom.i_code))
    return sorted(residues, key=lambda x: (int(x[0]), x[1]))


def count_residues(pdb_path: Path, chain_id: str) -> int:
    """Count unique residues in a specific chain.

    Args:
        pdb_path: Path to the PDB file.
        chain_id: The chain ID to count residues for.

    Returns:
        Number of unique residues.
    """
    return len(get_residues(pdb_path, chain_id))


def remove_chain(input_pdb: Path, output_pdb: Path, chain_id: str) -> None:
    """Create a new PDB file excluding a specific chain.

    Args:
        input_pdb: Path to the source PDB.
        output_pdb: Path to the destination PDB.
        chain_id: Chain ID to exclude.
    """
    with (
        input_pdb.open(encoding="utf-8") as f_in,
        output_pdb.open("w", encoding="utf-8") as f_out,
    ):
        for line in f_in:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21].strip()
                if chain == chain_id:
                    continue
            f_out.write(line)


def extract_protein(
    input_pdb: Path,
    output_pdb: Path,
    target_chain: str,
    exclude_chains: list[str] | None = None,
) -> None:
    """Extract protein ATOM/HETATM lines, excluding specific ligand chains.

    If exclude_chains is provided (even if empty), extracts all ATOM and HETATM
    lines for all chains EXCEPT those specified in exclude_chains.
    Otherwise, extracts only ATOM lines for the target_chain.

    Args:
        input_pdb: Path to the source PDB.
        output_pdb: Path to the destination PDB.
        target_chain: Target chain ID.
        exclude_chains: Optional list of chain IDs to exclude.
    """
    lines = []
    exclude_set = set(exclude_chains) if exclude_chains else set()
    with input_pdb.open(encoding="utf-8") as fh:
        for line in fh:
            if (
                line.startswith(("ATOM", "HETATM"))
                and len(line) >= PDB_COORD_LINE_MIN_LEN
            ):
                atom = PDBAtom(line)
                if exclude_chains is not None:
                    if atom.chain_id not in exclude_set:
                        lines.append(line)
                elif line.startswith("ATOM") and atom.chain_id == target_chain:
                    lines.append(line)
    with output_pdb.open("w", encoding="utf-8") as fh:
        fh.writelines(lines)


def extract_residue_pdb(
    input_pdb: Path, output_pdb: Path, chain_id: str, res_seq: int, i_code: str
) -> str:
    """Extract all atom/hetatm records for a specific residue.

    Args:
        input_pdb: Path to the source PDB.
        output_pdb: Path to the destination PDB.
        chain_id: Chain ID of the residue.
        res_seq: Residue sequence number.
        i_code: Insertion code.

    Returns:
        The residue name (str).
    """
    lines = []
    res_name = "UNK"
    with input_pdb.open(encoding="utf-8") as fh:
        for line in fh:
            if (
                line.startswith(("ATOM", "HETATM"))
                and len(line) >= PDB_COORD_LINE_MIN_LEN
            ):
                atom = PDBAtom(line)
                if (
                    atom.chain_id == chain_id
                    and atom.res_seq == res_seq
                    and atom.i_code == i_code
                ):
                    lines.append(line)
                    res_name = atom.res_name
    with output_pdb.open("w", encoding="utf-8") as fh:
        fh.writelines(lines)
    return res_name


def extract_chain_pdb(input_pdb: Path, output_pdb: Path, chain_id: str) -> str:
    """Extract all atom/hetatm records for a specific chain.

    Args:
        input_pdb: Path to the source PDB.
        output_pdb: Path to the destination PDB.
        chain_id: Chain ID of the residue.

    Returns:
        The residue name of the first residue in the chain.
    """
    lines = []
    res_name = "UNK"
    first_seen = False
    with input_pdb.open(encoding="utf-8") as fh:
        for line in fh:
            if (
                line.startswith(("ATOM", "HETATM"))
                and len(line) >= PDB_COORD_LINE_MIN_LEN
            ):
                atom = PDBAtom(line)
                if atom.chain_id == chain_id:
                    lines.append(line)
                    if not first_seen:
                        res_name = atom.res_name
                        first_seen = True
    with output_pdb.open("w", encoding="utf-8") as fh:
        fh.writelines(lines)
    return res_name
