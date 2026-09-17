#!/usr/bin/env python3
"""PRISM Tool: Experimental Structure Merging.

Merges multiple experimental structures aligned to a predicted model
(e.g. AlphaFold) using Biopython.  Creates a unified experimental
template with continuous, gapless renumbering and a merged alignment file.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

from Bio.PDB import PDBIO, Chain, Model, PDBParser, Structure, Superimposer

try:
    from pdb_utils import INPUT_DIR, OUTPUT_TOOLS_DIR
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from pdb_utils import INPUT_DIR, OUTPUT_TOOLS_DIR

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("tools.merge_experimental_templates")


def parse_pir_sequences(align_file: str) -> tuple[dict[str, str], list[str]]:
    """Parse a PIR alignment file and extract full alignment strings.

    Args:
        align_file: Path to the PIR alignment file.

    Returns:
        A 2-tuple of ``(sequences_dict, all_codes)`` where *sequences_dict*
        maps alignment codes to their full sequences (including gaps) and
        *all_codes* is an ordered list of code prefixes.
    """
    sequences: dict[str, str] = {}
    all_codes: list[str] = []
    curr_code: str | None = None
    curr_seq: list[str] = []

    with Path(align_file).open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">P1;"):
                if curr_code:
                    sequences[curr_code] = (
                        "".join(curr_seq).replace("\n", "").replace(" ", "")
                    )
                curr_code = line[4:].strip()
                descript_code = curr_code.split("_")[0]
                all_codes.append(descript_code)
                curr_seq = []
            elif curr_code and not line.startswith(("structure", "sequence", " ")):
                curr_seq.append(line.strip().rstrip("*"))
        if curr_code:
            sequences[curr_code] = "".join(curr_seq).replace("\n", "").replace(" ", "")

    return sequences, all_codes


def get_all_residues(structure: Any) -> list[Any]:
    """Extract all residues from the first model of a Biopython structure.

    Args:
        structure: A Biopython ``Structure`` object.

    Returns:
        List of residue objects.
    """
    return list(structure[0].get_residues())


def merge_structures(  # noqa: PLR0912, PLR0915
    align_file: str,
    ref_code: str,
    output_pdb: str | None = None,
) -> None:
    """Merge multiple experimental structures into a unified template PDB.

    Superimposes each experimental template onto the reference model,
    then creates a merged PDB by selecting experimental residues at each
    alignment position.

    Args:
        align_file: Path to the alignment file (``.ali``).
        ref_code: Code of the predictive model used as the reference.
        output_pdb: Optional output filename for the merged PDB.  If
            *None*, a descriptive default is generated.
    """
    pir_seqs, all_codes = parse_pir_sequences(align_file)
    logger.info("PIR sequences found: %s", list(pir_seqs.keys()))

    if not output_pdb:
        ref_prefix = ref_code.split("_", maxsplit=1)[0]
        template_codes = [c for c in all_codes if c not in (ref_prefix, "FullSeq")]
        output_pdb_name = "_".join(template_codes) + "_merged_experimental.pdb"
    else:
        output_pdb_name = Path(output_pdb).name

    output_pdb_path = OUTPUT_TOOLS_DIR / output_pdb_name

    target_code: str | None = None
    templates: list[str] = []

    with Path(align_file).open(encoding="utf-8") as fh:
        lines = fh.readlines()
        for i, line in enumerate(lines):
            if line.startswith(">P1;"):
                code = line[4:].strip()
                type_line = lines[i + 1]
                if type_line.startswith("sequence:"):
                    target_code = code
                elif type_line.startswith("structure") and code != ref_code:
                    templates.append(code)

    if ref_code not in pir_seqs:
        logger.error("Reference model '%s' not found in alignment.", ref_code)
        return
    if not target_code:
        logger.error("Target sequence not found in alignment.")
        return

    logger.info("Reference: %s", ref_code)
    logger.info("Target: %s", target_code)
    logger.info("Experimental Templates: %s", templates)

    target_seq = pir_seqs[target_code]
    ref_seq = pir_seqs[ref_code]

    parser = PDBParser(QUIET=True)

    ref_pdb_path = INPUT_DIR / ref_code
    if not ref_pdb_path.exists():
        logger.error("Could not find PDB for reference %s", ref_pdb_path)
        return
    ref_structure = parser.get_structure(ref_code, str(ref_pdb_path))
    ref_residues = get_all_residues(ref_structure)

    transformed_structures: dict[str, Any] = {}
    superimposer = Superimposer()

    for temp_code in templates:
        logger.info("Superposing %s onto %s...", temp_code, ref_code)
        temp_pdb = INPUT_DIR / temp_code

        if not temp_pdb.exists():
            logger.warning("  Could not find PDB for %s. Skipping.", temp_code)
            continue

        temp_structure = parser.get_structure(temp_code, str(temp_pdb))
        temp_residues = get_all_residues(temp_structure)
        temp_seq = pir_seqs[temp_code]

        ref_atoms: list[Any] = []
        temp_atoms: list[Any] = []
        r_idx, t_idx = 0, 0

        for r_char, t_char in zip(ref_seq, temp_seq, strict=False):
            is_r_valid = r_char != "-"
            is_t_valid = t_char != "-"

            if (
                is_r_valid
                and is_t_valid
                and r_idx < len(ref_residues)
                and t_idx < len(temp_residues)
            ):
                r_res = ref_residues[r_idx]
                t_res = temp_residues[t_idx]
                if "CA" in r_res and "CA" in t_res:
                    ref_atoms.append(r_res["CA"])
                    temp_atoms.append(t_res["CA"])

            if is_r_valid:
                r_idx += 1
            if is_t_valid:
                t_idx += 1

        if len(ref_atoms) == 0:
            logger.warning("  No aligned CA atoms found for %s.", temp_code)
            continue

        superimposer.set_atoms(ref_atoms, temp_atoms)
        superimposer.apply(temp_structure.get_atoms())

        transformed_structures[temp_code] = temp_structure
        logger.info("  RMSD: %.3f Å", superimposer.rms)

    # --- BUILD MERGED STRUCTURE ---
    logger.info("Merging superposed template structures into merged model...")
    merged_structure = Structure.Structure("Merged")
    merged_model = Model.Model(0)
    merged_chain = Chain.Chain("A")
    merged_model.add(merged_chain)
    merged_structure.add(merged_model)

    current_res_seq = 0

    for pos in range(len(target_seq)):
        if target_seq[pos] in ("-", ".", "/"):
            continue

        chosen_temp: str | None = None
        for t_code in templates:
            if t_code in transformed_structures and pir_seqs[t_code][pos] != "-":
                chosen_temp = t_code
                break

        if chosen_temp:
            temp_seq_up_to_pos = pir_seqs[chosen_temp][:pos]
            temp_res_idx = len(temp_seq_up_to_pos.replace("-", ""))

            temp_structure = transformed_structures[chosen_temp]
            temp_residues = get_all_residues(temp_structure)

            if temp_res_idx < len(temp_residues):
                current_res_seq += 1
                res_to_copy = temp_residues[temp_res_idx].copy()
                res_to_copy.id = (" ", current_res_seq, " ")
                merged_chain.add(res_to_copy)

    # --- MERGE NON-A CHAINS (LIGANDS) ---
    if templates and templates[0] in transformed_structures:
        first_temp = transformed_structures[templates[0]]
        merged_chain_b = Chain.Chain("B")
        for chain in first_temp[0]:
            if chain.id != "A":
                for res in chain:
                    current_res_seq += 1
                    res_copy = res.copy()
                    res_copy.id = (res_copy.id[0], current_res_seq, " ")
                    merged_chain_b.add(res_copy)
        if len(merged_chain_b) > 0:
            merged_model.add(merged_chain_b)

    io = PDBIO()
    io.set_structure(merged_structure)
    io.save(str(output_pdb_path))
    logger.info("Merged experimental PDB written to %s", output_pdb_path)

    # --- WRITE MERGED ALIGNMENT ---
    ali_path = Path(align_file)
    output_ali_name = ali_path.with_stem(ali_path.stem + "_merged").name
    output_ali_path = OUTPUT_TOOLS_DIR / output_ali_name

    merged_seq = list("-" * len(target_seq))
    for pos in range(len(target_seq)):
        for t_code in templates:
            if pir_seqs[t_code][pos] != "-":
                merged_seq[pos] = pir_seqs[t_code][pos]
                break

    with output_ali_path.open("w", encoding="utf-8") as fh:
        fh.write(f">P1;{output_pdb_name}\n")
        fh.write(f"structure:{output_pdb_name}:FIRST:@:END:@::::\n")
        fh.write("".join(merged_seq) + "*\n")

        fh.write(f">P1;{ref_code}\n")
        fh.write(f"structure:{ref_code}:FIRST:@:END:@::::\n")
        fh.write(ref_seq + "*\n")

        fh.write(f">P1;{target_code}\n")
        fh.write(f"sequence:{target_code}:FIRST:@:END:@::::\n")
        fh.write(target_seq + "*\n")

    logger.info("Merged alignment written to %s", output_ali_path)


def main() -> None:
    """Entry point for the experimental structure merging tool."""
    parser = argparse.ArgumentParser(
        description=(
            "Merge multiple experimental structures aligned to a "
            "predicted model using Biopython."
        )
    )
    parser.add_argument("alignment", help="Path to the alignment file (.ali)")
    parser.add_argument(
        "reference",
        help="Code of the predictive model in the alignment to use as reference",
    )
    parser.add_argument(
        "--output_pdb", default=None, help="Output filename for the merged PDB"
    )

    args = parser.parse_args()

    if not Path(args.alignment).exists():
        logger.error("Alignment file %s not found.", args.alignment)
        sys.exit(1)

    merge_structures(args.alignment, args.reference, args.output_pdb)


if __name__ == "__main__":
    main()
