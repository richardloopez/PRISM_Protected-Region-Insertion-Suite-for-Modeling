#!/usr/bin/env python3
"""PRISM PDB Pre/Post-Processor.

Prepares raw PDB files for the PRISM pipeline (``prep`` mode) by splitting
protein and ligand chains, renumbering atoms, and generating a detailed
JSON log.  Restores original metadata after modeling (``retro`` mode) using
the same log.

Supports post-translational modifications (PTMs) as rigid-body attachments
whose relative coordinates are preserved across the modeling cycle.

Usage:
    python3 tools/prep_prism_pdb.py prep  input.pdb A B
    python3 tools/prep_prism_pdb.py retro model.pdb data.json
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np

try:
    from pdb_utils import OUTPUT_TOOLS_DIR, PDBAtom, copy_to_input
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from pdb_utils import OUTPUT_TOOLS_DIR, PDBAtom, copy_to_input

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("tools.prep_prism_pdb")

# Mathematical and parsing constants
NORM_EPSILON = 1e-10
EXPECTED_REMAPPING_PARTS = 2


# ============================================================================
#                                 MATH HELPERS
# ============================================================================


def get_local_frame(
    n_pos: tuple[float, ...],
    ca_pos: tuple[float, ...],
    c_pos: tuple[float, ...],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute an orthonormal basis from backbone N, CA, C coordinates.

    Args:
        n_pos: Coordinates of the backbone nitrogen.
        ca_pos: Coordinates of the C-alpha.
        c_pos: Coordinates of the backbone carbon.

    Returns:
        A 3-tuple of unit vectors ``(v_x, v_y, v_z)`` defining the local
        reference frame as np.ndarray objects.
    """
    n = np.array(n_pos)
    ca = np.array(ca_pos)
    c = np.array(c_pos)

    # Unit vector from CA to C
    diff_xc = c - ca
    norm_xc = np.linalg.norm(diff_xc)
    v_x = diff_xc / norm_xc if norm_xc > NORM_EPSILON else np.zeros(3)

    # Unit vector from CA to N
    diff_nc = n - ca
    norm_nc = np.linalg.norm(diff_nc)
    v_nc = diff_nc / norm_nc if norm_nc > NORM_EPSILON else np.zeros(3)

    # Orthonormal basis
    v_z_raw = np.cross(v_x, v_nc)
    norm_z = np.linalg.norm(v_z_raw)
    v_z = v_z_raw / norm_z if norm_z > NORM_EPSILON else np.zeros(3)

    v_y = np.cross(v_z, v_x)

    return v_x, v_y, v_z


# ============================================================================
#                                 HELPER FUNCTIONS
# ============================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for prep/retro modes.

    Returns:
        Parsed ``Namespace`` with mode-specific attributes.
    """
    parser = argparse.ArgumentParser(description="PRISM PDB Pre/Post-Processor.")
    subparsers = parser.add_subparsers(
        dest="mode", required=True, help="Mode of operation"
    )

    # --- PREP MODE ---
    parser_prep = subparsers.add_parser(
        "prep", help="Prepare PDB for PRISM (Generate Chain A/B and Log)."
    )
    parser_prep.add_argument("input_pdb", help="Original input PDB file.")
    parser_prep.add_argument(
        "protein_chains", help="Chains to merge into Chain A (e.g. 'A,C')."
    )
    parser_prep.add_argument(
        "ligand_chains", help="Chains to convert to BLK Chain B (e.g. 'B,D')."
    )
    parser_prep.add_argument(
        "ptm_args",
        nargs="*",
        help="Optional PTM arguments (e.g. posttranslational=2 E-1,E-2:A-265 ...)",
    )

    # --- RETRO MODE ---
    parser_retro = subparsers.add_parser(
        "retro", help="Restore original ligand info to PRISM output."
    )
    parser_retro.add_argument(
        "prism_output_pdb", help="The output PDB from PRISM (with BLK ligands)."
    )
    parser_retro.add_argument(
        "log_file", help="The .json log file generated during the 'prep' stage."
    )
    parser_retro.add_argument(
        "ptm_args",
        nargs="*",
        help=(
            "Optional remapping arguments "
            "(e.g. posttranslational=2 A-209-new-A-232 ...)"
        ),
    )

    return parser.parse_args()


def format_atom_name_blk(counter: int) -> str:
    """Generate a BLK-style atom name (X1, X2, ...).

    Args:
        counter: Sequential atom index within the BLK residue.

    Returns:
        Atom name string.
    """
    return f"X{counter}"


def parse_ptm_string(
    ptm_str: str,
) -> tuple[list[tuple[str, int]], tuple[str, int]] | None:
    """Parse a PTM specification string.

    Accepts strings like ``'E-1,E-2:A-265'`` and returns structured IDs.

    Args:
        ptm_str: Raw PTM specification.

    Returns:
        A 2-tuple of ``(ptm_residue_ids, attachment_id)`` or *None* on
        parse failure.
    """
    if ":" not in ptm_str:
        return None
    ptm_part, attach_part = ptm_str.split(":")

    def _parse_id(s: str) -> tuple[str, int] | None:
        if "-" not in s:
            return None
        c, seq = s.split("-", 1)
        return (c.strip(), int(seq.strip()))

    ptms = [_parse_id(s) for s in ptm_part.split(",")]
    ptms = [p for p in ptms if p is not None]
    attach = _parse_id(attach_part)
    if not ptms or not attach:
        return None
    return ptms, attach


def parse_remap_string(
    s: str,
) -> tuple[tuple[str, int], tuple[str, int]] | None:
    """Parse a residue remapping string for retro mode.

    Accepts strings like ``'A-209-new-A-232'``.

    Args:
        s: Raw remap specification.

    Returns:
        A 2-tuple of ``(original_id, model_id)`` or *None* on failure.
    """
    if "-new-" not in s:
        return None
    parts = s.split("-new-")
    if len(parts) != EXPECTED_REMAPPING_PARTS:
        return None

    def _parse_id(ss: str) -> tuple[str, int] | None:
        if "-" not in ss:
            return None
        c, seq = ss.split("-", 1)
        return (c.strip(), int(seq.strip()))

    o_id = _parse_id(parts[0])
    m_id = _parse_id(parts[1])
    if not o_id or not m_id:
        return None
    return o_id, m_id


# ============================================================================
#                                   PREP LOGIC
# ============================================================================


def run_prep(  # noqa: PLR0912, PLR0915
    input_path: str,
    prot_chains_str: str,
    lig_chains_str: str,
    ptm_args: list[str] | None = None,
) -> None:
    """Prepare a PDB file for PRISM by splitting chains and generating a log.

    Args:
        input_path: Path to the original PDB file.
        prot_chains_str: Comma-separated protein chain IDs.
        lig_chains_str: Comma-separated ligand chain IDs.
        ptm_args: Optional list of PTM specification strings.
    """
    if ptm_args is None:
        ptm_args = []

    input_file = Path(input_path)
    if not input_file.exists():
        raise FileNotFoundError(f"File {input_path} not found.")

    prot_chains = [c.strip() for c in prot_chains_str.split(",")]
    lig_chains = [c.strip() for c in lig_chains_str.split(",")]

    logger.info("[PREP] Processing %s", input_path)
    logger.info("[PREP] Protein Chains -> A: %s", prot_chains)
    logger.info("[PREP] Ligand Chains  -> B: %s", lig_chains)

    protein_atoms: list[Any] = []
    ligand_atoms: list[Any] = []

    # 1. Read Atoms
    with Path(input_path).open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                atom = PDBAtom(line)
                chain = atom.chain_id
                if chain in prot_chains:
                    protein_atoms.append(atom)
                elif chain in lig_chains:
                    ligand_atoms.append(atom)

    # 2. Handle PTMs
    ptm_specs: list[Any] = []
    for arg in ptm_args:
        if arg.startswith("posttranslational="):
            continue
        spec = parse_ptm_string(arg)
        if spec:
            ptm_specs.append(spec)

    ptm_map: list[dict[str, Any]] = []
    with Path(input_path).open(encoding="utf-8") as fh:
        full_atom_list = [
            PDBAtom(line) for line in fh if line.startswith(("ATOM", "HETATM"))
        ]

    res_dict: dict[tuple[str, int], list[Any]] = {}
    for atom in full_atom_list:
        key = (atom.chain_id, atom.res_seq)
        if key not in res_dict:
            res_dict[key] = []
        res_dict[key].append(atom)

    for ptms, attach in ptm_specs:
        if attach not in res_dict:
            logger.warning(
                "[PREP] Attachment residue %s not found. Skipping PTM.",
                attach,
            )
            continue

        attach_atoms = res_dict[attach]
        n_at = next((a for a in attach_atoms if a.name.strip() == "N"), None)
        ca_at = next((a for a in attach_atoms if a.name.strip() == "CA"), None)
        c_at = next((a for a in attach_atoms if a.name.strip() == "C"), None)

        if not (n_at and ca_at and c_at):
            logger.warning(
                "[PREP] N, CA, or C missing for %s. Cannot compute local frame.",
                attach,
            )
            continue

        n_p = (n_at.x, n_at.y, n_at.z)
        ca_p = (ca_at.x, ca_at.y, ca_at.z)
        c_p = (c_at.x, c_at.y, c_at.z)
        vx, vy, vz = get_local_frame(n_p, ca_p, c_p)

        ptm_residue_data: list[dict[str, Any]] = []
        for ptm_res_id in ptms:
            if ptm_res_id not in res_dict:
                logger.warning("[PREP] PTM residue %s not found.", ptm_res_id)
                continue

            for p_at in res_dict[ptm_res_id]:
                rel_p = np.array((p_at.x, p_at.y, p_at.z)) - ca_p
                rx = np.dot(rel_p, vx)
                ry = np.dot(rel_p, vy)
                rz = np.dot(rel_p, vz)

                ptm_residue_data.append(
                    {
                        "name": p_at.name.strip(),
                        "res_name": p_at.res_name,
                        "res_seq": p_at.res_seq,
                        "chain_id": p_at.chain_id,
                        "element": p_at.element,
                        "record_type": p_at.record_type,
                        "rel_coords": (rx, ry, rz),
                    }
                )

                # Remove from output lists if present
                protein_atoms = [
                    a
                    for a in protein_atoms
                    if not (a.chain_id == p_at.chain_id and a.res_seq == p_at.res_seq)
                ]
                ligand_atoms = [
                    a
                    for a in ligand_atoms
                    if not (a.chain_id == p_at.chain_id and a.res_seq == p_at.res_seq)
                ]

        ptm_map.append(
            {
                "orig_attach_chain": attach[0],
                "orig_attach_res_seq": attach[1],
                "ptm_atoms": ptm_residue_data,
            }
        )
        logger.info(
            "[PREP] PTM: %s attached to %s -> Processed %d atoms.",
            ptms,
            attach,
            len(ptm_residue_data),
        )

    # 3. Prepare Data Structures
    new_chain_a_lines: list[str] = []
    new_chain_b_lines: list[str] = []

    log_data: dict[str, Any] = {
        "original_filename": Path(input_path).name,
        "protein_map": [],
        "ligand_map": {},
        "ptm_map": ptm_map,
    }

    current_serial = 1

    # --- PROCESS PROTEIN (CHAIN A) ---
    current_res_seq = 0
    prev_id: tuple[str, int, str] | None = None

    for atom in protein_atoms:
        curr_id = (atom.chain_id, atom.res_seq, atom.i_code)
        if curr_id != prev_id:
            current_res_seq += 1

        new_atom = PDBAtom(atom.line)
        new_atom.serial = current_serial
        new_atom.chain_id = "A"
        new_atom.res_seq = current_res_seq
        new_atom.record_type = "ATOM  "

        log_data["protein_map"].append(
            {
                "new_chain": "A",
                "new_res_seq": current_res_seq,
                "new_atom_name": new_atom.name.strip(),
                "orig_chain": atom.chain_id,
                "orig_res_name": atom.res_name,
                "orig_res_seq": atom.res_seq,
                "orig_atom_name": atom.name.strip(),
            }
        )

        new_chain_a_lines.append(new_atom.to_pdb_line())
        current_serial += 1
        prev_id = curr_id

    # --- PROCESS LIGAND (CHAIN B) ---
    prev_id = None
    atom_counter = 1
    log_res_seq = 0

    for atom in ligand_atoms:
        curr_id = (atom.chain_id, atom.res_seq, atom.i_code)
        if curr_id != prev_id:
            current_res_seq += 1
            log_res_seq += 1
            atom_counter = 1

        new_name = format_atom_name_blk(atom_counter)

        new_atom = PDBAtom(atom.line)
        new_atom.serial = current_serial
        new_atom.chain_id = "B"
        new_atom.res_name = "BLK"
        new_atom.res_seq = current_res_seq
        new_atom.name = new_name
        new_atom.element = "X"
        new_atom.temp = 99.99
        new_atom.record_type = "HETATM"

        key = f"{log_res_seq}_{new_name}"
        log_data["ligand_map"][key] = {
            "orig_chain": atom.chain_id,
            "orig_res_name": atom.res_name,
            "orig_res_seq": atom.res_seq,
            "orig_atom_name": atom.name.strip(),
            "orig_element": atom.element,
            "orig_record": atom.record_type.strip(),
        }

        new_chain_b_lines.append(new_atom.to_pdb_line())
        current_serial += 1
        atom_counter += 1
        prev_id = curr_id

    # 4. Output
    base = input_file.stem
    out_pdb_path = OUTPUT_TOOLS_DIR / f"{base}_prism_prep.pdb"
    out_log_path = OUTPUT_TOOLS_DIR / f"{base}_prism_data.json"

    with out_pdb_path.open("w", encoding="utf-8") as fh:
        for line in new_chain_a_lines:
            fh.write(line + "\n")
        if new_chain_a_lines:
            fh.write("TER\n")
        for line in new_chain_b_lines:
            fh.write(line + "\n")
        if new_chain_b_lines:
            fh.write("TER\n")
        fh.write("END\n")

    with out_log_path.open("w", encoding="utf-8") as fh:
        json.dump(log_data, fh, indent=4)

    logger.info("[PREP] Success!")
    logger.info(" > Generated PDB: %s", out_pdb_path)
    logger.info(" > Generated Log: %s", out_log_path)

    copy_to_input(out_pdb_path)


# ============================================================================
#                                   RETRO LOGIC
# ============================================================================


def run_retro(  # noqa: PLR0912, PLR0915
    model_path: str,
    log_path: str,
    extra_args: list[str] | None = None,
) -> None:
    """Restore original ligand and protein metadata to a PRISM output PDB.

    Args:
        model_path: Path to the PRISM-generated model PDB.
        log_path: Path to the JSON log from the ``prep`` stage.
        extra_args: Optional list of remapping specification strings.
    """
    if extra_args is None:
        extra_args = []

    if not Path(model_path).exists():
        raise FileNotFoundError(f"Model file {model_path} not found.")
    if not Path(log_path).exists():
        raise FileNotFoundError(f"Log file {log_path} not found.")

    logger.info("[RETRO] Restoring original ligand info...")
    logger.info(" > Model: %s", model_path)
    logger.info(" > Log:   %s", log_path)

    with Path(log_path).open(encoding="utf-8") as fh:
        log_data = json.load(fh)

    ligand_map = log_data.get("ligand_map", {})
    ptm_map_data = log_data.get("ptm_map", [])

    # 1. Parse Extra Args (Remappings)
    model_to_orig_explicit: dict[tuple[str, int], tuple[str, int]] = {}
    for arg in extra_args:
        if arg.startswith("posttranslational="):
            continue
        remap = parse_remap_string(arg)
        if remap:
            o_id, m_id = remap
            model_to_orig_explicit[m_id] = o_id

    # 2. Build Lookups from Log
    prot_id_lookup: dict[int, tuple[str, int, str]] = {}
    for p in log_data.get("protein_map", []):
        seq = p["new_res_seq"]
        if seq not in prot_id_lookup:
            prot_id_lookup[seq] = (
                p["orig_chain"],
                p["orig_res_seq"],
                p["orig_res_name"],
            )

    found_attachments: dict[tuple[str, int], dict[str, Any]] = {}
    target_orig_ids = {
        (p["orig_attach_chain"], p["orig_attach_res_seq"]) for p in ptm_map_data
    }

    restored_lines: list[str] = []
    log_res_seq = 0
    prev_id: tuple[str, int, str] | None = None
    current_serial = 1

    with Path(model_path).open(encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                lstrip = line.strip()
                if lstrip and not lstrip.startswith(("TER", "END", "CONECT")):
                    restored_lines.append(lstrip)
                continue

            atom = PDBAtom(line)
            curr_id = (atom.chain_id, atom.res_seq, atom.i_code)

            orig_id: tuple[str, int] | None = None
            m_id = (atom.chain_id, atom.res_seq)

            if m_id in model_to_orig_explicit:
                orig_id = model_to_orig_explicit[m_id]
            elif atom.chain_id == "A":
                if atom.res_seq in prot_id_lookup:
                    info = prot_id_lookup[atom.res_seq]
                    orig_id = (info[0], info[1])
            elif atom.chain_id == "B":
                if curr_id != prev_id:
                    log_res_seq += 1
                key = f"{log_res_seq}_{atom.name.strip()}"
                if key in ligand_map:
                    info_lig = ligand_map[key]
                    orig_id = (info_lig["orig_chain"], info_lig["orig_res_seq"])
                    atom.chain_id = info_lig["orig_chain"]
                    atom.res_name = info_lig["orig_res_name"]
                    atom.res_seq = info_lig["orig_res_seq"]
                    atom.name = info_lig["orig_atom_name"]
                    atom.element = info_lig["orig_element"]
                    atom.record_type = f"{info_lig['orig_record']:<6}"
                    atom.temp = 0.00

            if orig_id in target_orig_ids:
                if orig_id not in found_attachments:
                    found_attachments[orig_id] = {}
                a_name = atom.name.strip()
                if a_name in ["N", "CA", "C"]:
                    found_attachments[orig_id][a_name] = atom

            atom.serial = current_serial
            restored_lines.append(atom.to_pdb_line())
            current_serial += 1
            prev_id = curr_id

    final_output = [
        line for line in restored_lines if not line.startswith(("ATOM", "HETATM"))
    ]

    all_atoms = [
        PDBAtom(line) for line in restored_lines if line.startswith(("ATOM", "HETATM"))
    ]

    # Process PTMs and add to all_atoms
    ptm_atoms: list[Any] = []
    if ptm_map_data:
        for ptm_spec in ptm_map_data:
            orig_attach_id = (
                ptm_spec["orig_attach_chain"],
                ptm_spec["orig_attach_res_seq"],
            )
            if orig_attach_id not in found_attachments:
                logger.warning(
                    "[RETRO] Could not find attachment %s in model.",
                    orig_attach_id,
                )
                continue

            bits = found_attachments[orig_attach_id]
            if not ("N" in bits and "CA" in bits and "C" in bits):
                continue

            n_at, ca_at, c_at = bits["N"], bits["CA"], bits["C"]
            vx, vy, vz = get_local_frame(
                (n_at.x, n_at.y, n_at.z),
                (ca_at.x, ca_at.y, ca_at.z),
                (c_at.x, c_at.y, c_at.z),
            )

            for p_at_info in ptm_spec["ptm_atoms"]:
                rel_coords = np.array(p_at_info["rel_coords"])
                # Project back from local frame to global coordinates
                new_pos = (
                    np.array((ca_at.x, ca_at.y, ca_at.z))
                    + rel_coords[0] * vx
                    + rel_coords[1] * vy
                    + rel_coords[2] * vz
                )
                new_x, new_y, new_z = new_pos[0], new_pos[1], new_pos[2]

                dummy_line = (
                    f"{'HETATM':<6}{0:>5} {'X':<4} {'PTM':>3} {'C'}{0:>4}"
                    f"    {0.0:>8.3f}{0.0:>8.3f}{0.0:>8.3f}{1.0:>6.2f}{0.0:>6.2f}"
                )
                new_atom = PDBAtom(dummy_line)
                new_atom.record_type = f"{p_at_info['record_type']:<6}"
                new_atom.name = p_at_info["name"]
                new_atom.res_name = p_at_info["res_name"]
                new_atom.chain_id = p_at_info["chain_id"]
                new_atom.res_seq = p_at_info["res_seq"]
                new_atom.x, new_atom.y, new_atom.z = new_x, new_y, new_z
                new_atom.element = p_at_info["element"]
                ptm_atoms.append(new_atom)

    full_atom_list = all_atoms + ptm_atoms

    chains_in_order: list[str] = []
    for a in full_atom_list:
        if a.chain_id not in chains_in_order:
            chains_in_order.append(a.chain_id)

    serial = 1
    final_atom_lines: list[str] = []
    for c in sorted(chains_in_order):
        for a in full_atom_list:
            if a.chain_id == c:
                a.serial = serial
                final_atom_lines.append(a.to_pdb_line())
                serial += 1
        final_atom_lines.append("TER")

    final_output.extend(final_atom_lines)
    final_output.append("END")

    # Output
    base_model = Path(model_path).stem
    out_path = OUTPUT_TOOLS_DIR / f"{base_model}_restored.pdb"

    with out_path.open("w", encoding="utf-8") as fh:
        for line in final_output:
            fh.write(line + "\n")

    logger.info("[RETRO] Success! Restored file saved as: %s", out_path)


# ============================================================================
#                                      MAIN
# ============================================================================

if __name__ == "__main__":
    args = parse_args()

    if args.mode == "prep":
        run_prep(args.input_pdb, args.protein_chains, args.ligand_chains, args.ptm_args)
    elif args.mode == "retro":
        run_retro(args.prism_output_pdb, args.log_file, args.ptm_args)
