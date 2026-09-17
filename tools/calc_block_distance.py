#!/usr/bin/env python3
"""PRISM Tool: BLK Distance Calculator.

Computes the distance between protein CA atoms and BLK gravity centres.
This matches the physics used in PRISM's ``add_hetatm_repulsion_shield``
and helps determine an appropriate ``BLOCK_REPULSION_RADIUS``.
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
from pathlib import Path
from typing import Any

from modeller import Environ, Model
from modeller.selection import Selection

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("tools.calc_block_distance")

MIN_SAFE_ATOMIC_DISTANCE = 2.0  # Å


def get_gravity_center(residue: Any) -> tuple[float, float, float]:
    """Compute the gravity centre of all atoms in a residue.

    Args:
        residue: A MODELLER residue object.

    Returns:
        A 3-tuple ``(x, y, z)`` of the centre-of-mass coordinates.
    """
    atoms = residue.atoms
    n = len(atoms)
    x = sum(a.x for a in atoms) / n
    y = sum(a.y for a in atoms) / n
    z = sum(a.z for a in atoms) / n
    return (x, y, z)


def calculate_distances(
    pdb_file: str,
    protein_chain: str,
    blk_chain: str,
) -> list[dict[str, Any]]:
    """Load a PDB and compute CA-to-BLK gravity-centre distances.

    Args:
        pdb_file: Path to the PDB file.
        protein_chain: Chain ID containing the protein.
        blk_chain: Chain ID containing the BLK/HETATM residues.

    Returns:
        List of dictionaries with per-residue distance data, or an empty
        list if no BLK residues or CA atoms are found.
    """
    env = Environ()
    env.io.atom_files_directory = ["."]
    env.io.hetatm = True
    env.libs.topology.read(file="$(LIB)/top_heav.lib")
    env.libs.parameters.read(file="$(LIB)/par.lib")

    mdl = Model(env, file=pdb_file)

    het_residues = [
        r
        for r in mdl.residues
        if r.chain.name == blk_chain
        and (r.hetatm or r.name == "BLK")
        and r.name != "HOH"
    ]
    if not het_residues:
        available_chains = sorted({r.chain.name for r in mdl.residues})
        logger.error("No HETATM/BLK residues found in chain %s.", blk_chain)
        logger.error("Available chains in PDB: %s", ", ".join(available_chains))
        return []

    try:
        prot_chain_obj = mdl.chains[protein_chain]
    except KeyError:
        available_chains = sorted(c.name for c in mdl.chains)
        logger.exception("Chain %s not found in %s.", protein_chain, pdb_file)
        logger.error("Available chains in PDB: %s", ", ".join(available_chains))
        return []

    protein_ca = Selection(prot_chain_obj).only_atom_types("CA")
    if not protein_ca:
        logger.error("No CA atoms found in chain %s.", protein_chain)
        return []

    results: list[dict[str, Any]] = []

    for ca in protein_ca:
        ca_coord = (ca.x, ca.y, ca.z)
        min_dist = float("inf")
        closest_het = None

        for res in het_residues:
            het_coord = get_gravity_center(res)
            dist = math.dist(ca_coord, het_coord)
            if dist < min_dist:
                min_dist = dist
                closest_het = res

        results.append(
            {
                "residue_index": ca.residue.index,
                "residue_num": ca.residue.num,
                "residue_name": ca.residue.name,
                "min_distance": min_dist,
                "closest_het": (
                    f"{closest_het.name}:{closest_het.num}:{closest_het.chain.name}"
                    if closest_het
                    else "N/A"
                ),
            }
        )

    logger.info(
        "[calc_block_distance] Computed distances for %d CA atoms.",
        len(results),
    )
    return results


def main() -> None:
    """Entry point for the BLK distance calculator."""
    parser = argparse.ArgumentParser(
        description=(
            "Calculate distances between protein CA atoms and BLK gravity centres."
        )
    )
    parser.add_argument("pdb_file", help="Path to the PDB file.")
    parser.add_argument(
        "--protein_chain", default="A", help="Chain ID for the protein (default: A)."
    )
    parser.add_argument(
        "--blk_chain", default="B", help="Chain ID for the BLK residues (default: B)."
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=10.0,
        help="Output only residues closer than this threshold (default: 10.0).",
    )

    args = parser.parse_args()

    if not Path(args.pdb_file).exists():
        logger.error("File %s not found.", args.pdb_file)
        sys.exit(1)

    logger.info("Analyzing %s...", args.pdb_file)
    logger.info("Protein Chain: %s, BLK Chain: %s", args.protein_chain, args.blk_chain)

    distances = calculate_distances(args.pdb_file, args.protein_chain, args.blk_chain)

    if not distances:
        logger.info("No distances calculated.")
        return

    distances.sort(key=lambda x: x["min_distance"])

    logger.info("\n" + "=" * 80)
    logger.info(
        "%s | %s | %s",
        "Residue".ljust(15),
        "Min Distance (A)".ljust(18),
        "Closest BLK Group".ljust(20),
    )
    logger.info("-" * 80)

    min_all = distances[0]["min_distance"]

    for d in distances:
        if d["min_distance"] <= args.threshold:
            res_str = f"{d['residue_name']} {d['residue_num']}"
            logger.info(
                "%s | %-18.3f | %s",
                res_str.ljust(15),
                d["min_distance"],
                d["closest_het"],
            )

    logger.info("=" * 80)
    logger.info("\nGLOBAL MINIMUM DISTANCE: %.3f A", min_all)
    logger.info(
        "SUGGESTED BLOCK_REPULSION_RADIUS: %.1f A (Example suggestion)",
        min_all,
    )
    if min_all < MIN_SAFE_ATOMIC_DISTANCE:
        logger.warning(
            "WARNING: The minimum distance is less than "
            "%s A.\n"
            "Make sure that the BLK residues are not too close to "
            "the protein.\n"
            "BLK is going to be fixed during the modeling process, but the "
            "protein could move its relative position to the BLK.",
            MIN_SAFE_ATOMIC_DISTANCE,
        )
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
