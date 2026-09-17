#!/usr/bin/env python3
"""PRISM Tool: Template Overlap Resolution.

Resolves sequence overlaps between multiple templates in a MODELLER
alignment.  Trims protein residues from lower-priority templates while
maintaining a configurable overlap buffer for structural continuity.
BLK/HETATM residues and chain breaks are always preserved.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Any

from modeller import Alignment, Environ

try:
    from pdb_utils import PDBAtom, copy_to_input
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from pdb_utils import PDBAtom, copy_to_input

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("tools.unify_templates")

MIN_TEMPLATES_FOR_UNIFICATION = 2


def renumber_pdb(  # noqa: PLR0912, PLR0915
    pdb_path: str, residues_to_keep: set[int], output_path: str | Path
) -> None:
    """Read a PDB, keep only specified protein residues plus all BLK, and renumber.

    Protein residues are renumbered starting at 1.  BLK/HETATM residues
    continue numbering sequentially from the last protein residue.

    Args:
        pdb_path: Path to the input PDB file.
        residues_to_keep: 1-indexed set of protein residue indices to
            retain.
        output_path: Path for the renumbered output PDB.
    """
    protein_atoms: list[Any] = []
    blk_atoms: list[Any] = []

    with Path(pdb_path).open(encoding="utf-8") as fh:
        prot_res_count = 0
        last_prot_res_id: tuple[str, int, str] | None = None
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                atom = PDBAtom(line)
                res_id = (atom.chain_id, atom.res_seq, atom.i_code)

                if atom.record_type == "HETATM" or atom.res_name == "BLK":
                    blk_atoms.append(atom)
                else:
                    if res_id != last_prot_res_id:
                        prot_res_count += 1
                        last_prot_res_id = res_id
                    if prot_res_count in residues_to_keep:
                        protein_atoms.append(atom)

    new_lines: list[str] = []
    serial = 1
    new_res_seq = 0

    # Renumber protein
    last_old_res_id: tuple[str, int, str] | None = None
    last_chain_id: str | None = None

    for atom in protein_atoms:
        if last_chain_id is not None and atom.chain_id != last_chain_id:
            new_lines.append("TER")
            new_res_seq = 0
            last_old_res_id = None

        old_res_id = (atom.chain_id, atom.res_seq, atom.i_code)
        if old_res_id != last_old_res_id:
            new_res_seq += 1
            last_old_res_id = old_res_id

        atom.serial = serial
        atom.res_seq = new_res_seq
        new_lines.append(atom.to_pdb_line())
        serial += 1
        last_chain_id = atom.chain_id

    if protein_atoms:
        new_lines.append("TER")

    # Renumber BLK
    last_old_res_id = None
    for atom in blk_atoms:
        if last_chain_id is not None and atom.chain_id != last_chain_id:
            new_lines.append("TER")
            last_old_res_id = None

        old_res_id = (atom.chain_id, atom.res_seq, atom.i_code)
        if old_res_id != last_old_res_id:
            new_res_seq += 1
            last_old_res_id = old_res_id

        atom.serial = serial
        atom.res_seq = new_res_seq
        new_lines.append(atom.to_pdb_line())
        serial += 1
        last_chain_id = atom.chain_id

    if blk_atoms:
        new_lines.append("TER")
    new_lines.append("END")

    with Path(output_path).open("w", encoding="utf-8") as fh:
        for line in new_lines:
            fh.write(line + "\n")


def unify_templates(  # noqa: PLR0912, PLR0915
    align_file: str, overlap_limit: int
) -> None:
    """Resolve sequence overlaps between templates in an alignment file.

    Templates appearing earlier in the alignment have higher priority.
    Protein residues from lower-priority templates that overlap with
    already-covered positions are trimmed, retaining only a buffer of
    *overlap_limit* residues at each junction.

    Args:
        align_file: Path to the MODELLER alignment file (``.ali``).
        overlap_limit: Number of overlap residues to preserve at
            junctions.
    """
    env = Environ()
    aln = Alignment(env, file=align_file)

    templates = [s for s in aln if s.prottyp.startswith("structure")]
    if len(templates) < MIN_TEMPLATES_FOR_UNIFICATION:
        logger.info(
            f"Less than {MIN_TEMPLATES_FOR_UNIFICATION} templates found. "
            "Nothing to unify."
        )
        return

    def is_aa(c: str) -> bool:
        return "A" <= c <= "Z"

    # Write temp alignment and read back sequences as strings
    ali_dir = str(Path(align_file).parent)
    tmp_fd, tmp_ali = tempfile.mkstemp(suffix=".ali", dir=ali_dir)
    os.close(tmp_fd)
    aln.write(file=tmp_ali)
    ali_strings: dict[str, str] = {}

    try:
        with Path(tmp_ali).open(encoding="utf-8") as fh:
            curr_code: str | None = None
            for line in fh:
                if line.startswith(">P1;"):
                    curr_code = line[4:].strip()
                    ali_strings[curr_code] = ""
                elif curr_code and not line.startswith(
                    ("structure", "sequence", "C;", " ")
                ):
                    cleaned = re.sub(r"\s+", "", line).rstrip("*")
                    ali_strings[curr_code] += cleaned
    finally:
        if Path(tmp_ali).exists():
            Path(tmp_ali).unlink()

    covered_positions: set[int] = set()
    new_sequences: dict[str, str] = {}

    target_code = next((s.code for s in aln if s.prottyp == "sequence"), None)
    target_seq_str = ali_strings.get(target_code) if target_code else None
    if target_seq_str:
        logger.info("Target sequence: %s. Mismatch-to-gap enabled.", target_code)

    for i, temp in enumerate(templates):
        logger.info("Processing template %d: %s", i + 1, temp.code)

        temp_seq_str = ali_strings.get(temp.code, "")
        if not temp_seq_str:
            logger.warning("  Sequence for %s not found in alignment.", temp.code)
            continue

        orig_pdb_basename = temp.atom_file.strip() if temp.atom_file else ""
        if not orig_pdb_basename:
            logger.error("  Missing atom_file in alignment. Cannot modify structure.")
            continue

        resolved_pdb_path = Path(align_file).parent / orig_pdb_basename
        if not resolved_pdb_path.exists():
            logger.error("  %s not found. Cannot modify structure.", resolved_pdb_path)
            continue

        pdb_residues = []
        seen_res = set()
        with resolved_pdb_path.open(encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    atom = PDBAtom(line)
                    res_key = (atom.chain_id, atom.res_seq, atom.i_code)
                    if res_key not in seen_res:
                        seen_res.add(res_key)
                        is_blk = atom.record_type == "HETATM" or atom.res_name == "BLK"
                        pdb_residues.append({"res_key": res_key, "is_blk": is_blk})

        prot_count = 0
        for r in pdb_residues:
            if not r["is_blk"]:
                prot_count += 1
                r["prot_res_idx"] = prot_count
            else:
                r["prot_res_idx"] = None

        non_gap_idx = 0
        align_to_pdb = {}
        for j, char in enumerate(temp_seq_str):
            if char not in ("-", "/"):
                if non_gap_idx < len(pdb_residues):
                    align_to_pdb[j] = pdb_residues[non_gap_idx]
                    non_gap_idx += 1
                else:
                    align_to_pdb[j] = None

        temp_aa_positions = []
        mismatch_positions: set[int] = set()
        for j, char in enumerate(temp_seq_str):
            res = align_to_pdb.get(j)
            if res is not None and not res["is_blk"]:
                if is_aa(char) and target_seq_str and j < len(target_seq_str):
                    target_res = target_seq_str[j]
                    if is_aa(target_res):
                        if char.upper() != target_res.upper():
                            mismatch_positions.add(j)
                            continue
                    else:
                        mismatch_positions.add(j)
                        continue
                temp_aa_positions.append(j)

        overlapping_aa_pos = [
            pos for pos in temp_aa_positions if pos in covered_positions
        ]

        if not overlapping_aa_pos or len(overlapping_aa_pos) <= overlap_limit:
            keep_aa_pos = set(temp_aa_positions)
            logger.info(
                f"  Template {i + 1} low overlap ({len(overlapping_aa_pos)}). "
                "Keeping available residues."
            )
        else:
            non_overlapping_pos = {
                pos for pos in temp_aa_positions if pos not in covered_positions
            }
            if overlap_limit > 0:
                expanded_nop: set[int] = set()
                for nop in non_overlapping_pos:
                    for offset in range(-overlap_limit, overlap_limit + 1):
                        expanded_nop.add(nop + offset)
                keep_aa_pos = {pos for pos in temp_aa_positions if pos in expanded_nop}
            else:
                keep_aa_pos = non_overlapping_pos

        residues_to_keep = set()
        for j in range(len(temp_seq_str)):
            res = align_to_pdb.get(j)
            if res is not None and not res["is_blk"] and j in keep_aa_pos:
                residues_to_keep.add(res["prot_res_idx"])

        new_seq_list: list[str] = []
        for j, char in enumerate(temp_seq_str):
            res = align_to_pdb.get(j)
            if res is not None and not res["is_blk"]:
                prot_res_idx = res["prot_res_idx"]
                if prot_res_idx in residues_to_keep and j not in mismatch_positions:
                    covered_positions.add(j)
                    new_seq_list.append(char)
                else:
                    new_seq_list.append("-")
            else:
                new_seq_list.append(char)

        # Reorder tail characters (gaps, chain breaks, BLK dots)
        last_aa_pos = -1
        for j in range(len(new_seq_list) - 1, -1, -1):
            if is_aa(new_seq_list[j]):
                last_aa_pos = j
                break

        tail_start = last_aa_pos + 1
        tail_chars = new_seq_list[tail_start:]
        if tail_chars:
            num_dots = tail_chars.count(".")
            has_slash = "/" in tail_chars
            num_gaps = len(tail_chars) - num_dots - (1 if has_slash else 0)

            new_tail: list[str] = ["-"] * num_gaps
            if has_slash:
                new_tail.append("/")
            new_tail.extend(["."] * num_dots)
            new_seq_list[tail_start:] = new_tail

        new_sequences[temp.code] = "".join(new_seq_list)

        # Generate renumbered PDB
        output_pdb_name = f"{temp.code.split('.')[0]}_unified.pdb"
        output_pdb_path = Path(align_file).parent / output_pdb_name

        renumber_pdb(str(resolved_pdb_path), residues_to_keep, output_pdb_path)
        logger.info("  Generated %s", output_pdb_path)
        temp.atom_file = output_pdb_name

        copy_to_input(output_pdb_path)

    # Write unified alignment
    ali_path = Path(align_file)
    output_ali_name = ali_path.with_stem(ali_path.stem + "_unified").name
    output_ali_path = ali_path.parent / output_ali_name

    with (
        ali_path.open(encoding="utf-8") as f_in,
        output_ali_path.open("w", encoding="utf-8") as f_out,
    ):
        current_code: str | None = None
        for line in f_in:
            if line.startswith(">P1;"):
                current_code = line[4:].strip()
                if current_code.endswith(".pdb"):
                    base_name = Path(current_code).stem
                    unified_name = f"{base_name}_unified.pdb"
                    f_out.write(f">P1;{unified_name}\n")
                else:
                    f_out.write(line)
            elif current_code in new_sequences and not line.startswith(
                ("structure", "sequence", "C;", " ")
            ):
                if "*" in line:
                    f_out.write(new_sequences[current_code] + "*\n")
                    current_code = None
            elif current_code in new_sequences and line.startswith(
                ("structure", "sequence")
            ):
                base_name = Path(current_code).stem
                unified_name = f"{base_name}_unified.pdb"
                f_out.write(line.replace(current_code, unified_name))
            else:
                f_out.write(line)

    logger.info("\nSuccess! Unified alignment written to %s", output_ali_path)
    copy_to_input(output_ali_path)


def main() -> None:
    """Entry point for the template unification tool."""
    parser = argparse.ArgumentParser(
        description="Unify templates by resolving sequence overlaps."
    )
    parser.add_argument("alignment", help="Path to the Modeller alignment file (.ali)")
    parser.add_argument(
        "--overlap",
        type=int,
        default=3,
        help="Number of residues allowed to overlap for continuity (default: 3)",
    )

    args = parser.parse_args()

    if not Path(args.alignment).exists():
        logger.error("Alignment file %s not found.", args.alignment)
        sys.exit(1)

    unify_templates(args.alignment, args.overlap)


if __name__ == "__main__":
    main()
