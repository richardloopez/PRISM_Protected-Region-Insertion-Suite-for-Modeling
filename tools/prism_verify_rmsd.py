#!/usr/bin/env python3
"""PRISM Coordinate Fidelity Verification.

Verifies that experimental coordinates have been preserved (0.0 Å RMSD)
in a PRISM output model by comparing CA positions between the original
template PDB and the generated model.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from modeller import Environ, Model, log

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("tools.prism_verify_rmsd")

RMSD_IDENTITY_THRESHOLD = 0.01  # Å


class PrismVerify:
    """Verifier for experimental coordinate preservation in PRISM models.

    Args:
        original_pdb: Path to the original experimental PDB.
        model_pdb: Path to the PRISM-generated model PDB.
        align_file: Path to the PIR alignment file.
    """

    def __init__(self, original_pdb: Path, model_pdb: Path, align_file: Path) -> None:
        """Initialize the verifier with paths to experimental and model PDBs."""
        self.original_pdb = original_pdb
        self.model_pdb = model_pdb
        self.align_file = align_file
        self.env = Environ()
        self._setup_modeller()

    def _setup_modeller(self) -> None:
        """Configure a minimal MODELLER environment for structure loading."""
        self.env.libs.topology.read(file="$(LIB)/top.lib")
        self.env.libs.parameters.read(file="$(LIB)/par.lib")
        self.env.io.atom_files_directory = ["."]
        log.none()

    def parse_pir_alignment(self, template_code: str) -> tuple[str, str]:
        """Parse a PIR file and return the template and target sequences.

        Args:
            template_code: Alignment code identifying the template entry.

        Returns:
            A 2-tuple of ``(template_sequence, target_sequence)``.

        Raises:
            SystemExit: If the alignment file is missing or codes not found.
        """
        if not Path(self.align_file).exists():
            raise FileNotFoundError(f"Alignment file {self.align_file} not found.")

        with Path(self.align_file).open(encoding="utf-8") as fh:
            content = fh.read()

        blocks = content.split(">P1;")
        template_seq = ""
        target_seq = ""

        for block in blocks[1:]:
            lines = block.strip().splitlines()
            code = lines[0].strip()
            header = lines[1].strip()
            sequence = "".join(
                line.strip()
                for line in lines[2:]
                if line.strip() and not line.startswith("#")
            )

            if code == template_code:
                template_seq = sequence
            elif header.startswith("sequence:"):
                target_seq = sequence

        if not template_seq or not target_seq:
            raise ValueError(
                f"Could not find template '{template_code}' "
                f"or target sequence in alignment."
            )

        return template_seq, target_seq

    def build_mapping_from_alignment(self, template_code: str) -> dict[int, int]:
        """Build a residue mapping from aligned positions.

        Positions where both template and target have actual residues (not
        gaps) define the mapping.

        Args:
            template_code: Alignment code identifying the template entry.

        Returns:
            Dictionary mapping template residue numbers to target residue
            numbers.
        """
        t_seq, q_seq = self.parse_pir_alignment(template_code)

        mapping: dict[int, int] = {}
        t_res_count = 0
        q_res_count = 0

        for t_char, q_char in zip(t_seq, q_seq, strict=False):
            if t_char not in ("-", "/", "*"):
                t_res_count += 1
            if q_char not in ("-", "/", "*"):
                q_res_count += 1
            if t_char not in ("-", "/", ".") and q_char not in ("-", "/", "."):
                mapping[t_res_count] = q_res_count

        return mapping

    def run_rmsd_check(
        self, mapping: dict[int, int], orig_chain: str, mod_chain: str
    ) -> None:
        """Calculate and report per-residue CA RMSD.

        Args:
            mapping: Template-to-model residue mapping.
            orig_chain: Chain ID in the original PDB.
            mod_chain: Chain ID in the model PDB.
        """
        mdl_orig = Model(self.env)
        mdl_orig.read(file=self.original_pdb)

        mdl_mod = Model(self.env)
        mdl_mod.read(file=self.model_pdb)

        results: list[float] = []
        count = 0

        logger.info("\n%s", "=" * 80)
        logger.info(
            "%s | %s | %s | %s",
            "RESIDUE".ljust(15),
            "ORIG POS".ljust(10),
            "MOD POS".ljust(10),
            "RMSD (Å)".ljust(10),
        )
        logger.info("-" * 80)

        for t_idx, q_idx in sorted(mapping.items()):
            try:
                res_o = mdl_orig.residues[f"{t_idx}:{orig_chain}"]
                res_m = mdl_mod.residues[f"{q_idx}:{mod_chain}"]

                ca_o = res_o.atoms["CA"]
                ca_m = res_m.atoms["CA"]

                dist = (
                    (ca_o.x - ca_m.x) ** 2
                    + (ca_o.y - ca_m.y) ** 2
                    + (ca_o.z - ca_m.z) ** 2
                ) ** 0.5

                status = "✓" if dist < RMSD_IDENTITY_THRESHOLD else "!"
                logger.info(
                    "%s %-13s | %-10s | %-10s | %.6f",
                    status,
                    res_o.pdb_name,
                    t_idx,
                    q_idx,
                    dist,
                )

                results.append(dist)
                count += 1
            except (KeyError, AttributeError) as exc:  # noqa: PERF203
                logger.debug("Skipping residue %d->%d: %s", t_idx, q_idx, exc)
                continue

        if count == 0:
            logger.error("No matching residues found for comparison.")
            return

        avg_rmsd = sum(results) / count
        max_rmsd = max(results)

        logger.info("-" * 80)
        logger.info("SUMMARY for %d residues:", count)
        logger.info(" > Average RMSD: %.6f Å", avg_rmsd)
        logger.info(" > Maximum RMSD: %.6f Å", max_rmsd)

        if max_rmsd < RMSD_IDENTITY_THRESHOLD:
            logger.info("\nSUCCESS: Experimental coordinates are FIXED.")
        else:
            logger.warning("\nWARNING: Coordinate drift detected (> 0.01 Å).")


def parse_manual_segments(segment_str: str) -> dict[int, int]:
    """Parse a manual segment specification string.

    Args:
        segment_str: Format ``'start-end:start-end'`` (e.g. ``'1-191:2-192'``).

    Returns:
        Residue mapping dictionary.

    Raises:
        SystemExit: If the format is invalid.
    """
    mapping: dict[int, int] = {}
    try:
        orig_part, mod_part = segment_str.split(":")
        o_start, o_end = map(int, orig_part.split("-"))
        m_start = int(mod_part.split("-")[0])

        for i in range(o_end - o_start + 1):
            mapping[o_start + i] = m_start + i
    except ValueError:
        sys.exit(
            "Error: Manual segments format must be 'start-end:start-end' "
            "(e.g. 1-191:2-192)"
        )
    return mapping


def main() -> None:
    """Entry point for the RMSD verification tool."""
    parser = argparse.ArgumentParser(
        description="PRISM Coordinate Fidelity Verification"
    )
    parser.add_argument("original_pdb", help="Original experimental PDB file.")
    parser.add_argument("model_pdb", help="Generated PRISM model PDB.")
    parser.add_argument("alignment", help="PIR alignment file used for modeling.")
    parser.add_argument(
        "--manual", help="Manual segments (e.g., '1-191:2-192'). Overrides alignment."
    )
    parser.add_argument(
        "--orig-chain", default="A", help="Chain in original PDB (default: A)."
    )
    parser.add_argument(
        "--mod-chain", default="A", help="Chain in model PDB (default: A)."
    )

    args = parser.parse_args()

    template_code = Path(args.original_pdb).name
    raw_code = (
        template_code.replace(".pdb", "")
        if template_code.endswith(".pdb")
        else template_code
    )

    verifier = PrismVerify(args.original_pdb, args.model_pdb, args.alignment)

    if args.manual:
        logger.info("[VERIFY] Using manual segment definition: %s", args.manual)
        mapping = parse_manual_segments(args.manual)
    else:
        logger.info(
            "[VERIFY] Automatically building mapping from %s...", args.alignment
        )
        mapping = verifier.build_mapping_from_alignment(template_code)
        if not mapping:
            mapping = verifier.build_mapping_from_alignment(raw_code)

    if not mapping:
        sys.exit(
            f"Error: No mapping generated. Check if '{template_code}' "
            f"is in {args.alignment}."
        )

    verifier.run_rmsd_check(mapping, args.orig_chain, args.mod_chain)


if __name__ == "__main__":
    main()
