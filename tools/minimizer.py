"""OpenMM Protein-Ligand Minimization Pipeline with GAFF2 and Stepwise Restraints.

This script performs a rigorous energy minimization of a target protein in the
presence of ligands or cofactors.
Positional restraints are applied in a two-phase approach to safely resolve
steric clashes introduced by hydrogen addition before performing a deep relaxation.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import openmm
from openmm import app, unit

try:
    from pdb_utils import OUTPUT_TOOLS_DIR
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from pdb_utils import OUTPUT_TOOLS_DIR

LOG_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
logger = logging.getLogger(__name__)


@dataclass
class MinimizerConfig:
    """Configuration parameters for the minimization process.

    Attributes:
        input_pdb (Path): Path to the input PDB file.
        output_pdb (Path): Path to the output minimized PDB file.
        log_file (Path): Path to the PRISM log file containing frozen residues.
        prmtop_path (Path): Path to the Amber prmtop topology file.
        inpcrd_path (Path): Path to the Amber inpcrd coordinate file.
        leap_pdb_path (Path): Path to the Amber-generated PDB file.
        target_chain_id (str): Chain ID for the target protein.
        k_phase1_all (unit.Quantity): Stiffness for initial hydrogen/clash relaxation.
        k_phase2_target (unit.Quantity): Final stiffness for the frozen backbone.
        k_phase2_env (unit.Quantity): Final stiffness for the environment.
        target_atoms (frozenset[str]): Atoms to restrain in the target chain.
    """

    input_pdb: Path
    output_pdb: Path
    log_file: Path
    prmtop_path: Path
    inpcrd_path: Path
    leap_pdb_path: Path
    target_chain_id: str = "A"

    k_phase1_all: unit.Quantity = field(
        default_factory=lambda: 1000.0 * unit.kilocalories_per_mole / unit.angstrom**2
    )
    k_phase2_target: unit.Quantity = field(
        default_factory=lambda: 100.0 * unit.kilocalories_per_mole / unit.angstrom**2
    )
    k_phase2_env: unit.Quantity = field(
        default_factory=lambda: 1000.0 * unit.kilocalories_per_mole / unit.angstrom**2
    )

    target_atoms: frozenset[str] = frozenset({"N", "CA", "C", "CB"})


def extract_target_residues(log_file: Path) -> set[int]:
    """Parse the PRISM log file to extract the list of frozen target residues.

    Args:
        log_file (Path): The log file containing the 'Freezing residues' list.

    Returns:
        set[int]: A set of integer residue indices that should be restrained.
    """
    if not log_file.exists():
        logger.error("Log file not found: %s", log_file)
        sys.exit(1)

    content = log_file.read_text(encoding="utf-8")
    match = re.search(r"Freezing residues:\s*\{([^}]+)\}", content)

    if not match:
        logger.error("No 'Freezing residues: {...}' pattern found in %s", log_file)
        sys.exit(1)

    residues_str = match.group(1)
    if not residues_str.strip():
        logger.error("The frozen residues list is empty.")
        sys.exit(1)

    try:
        target_residues = {int(x.strip()) for x in residues_str.split(",") if x.strip()}
    except ValueError:
        logger.exception("Failed to parse target residues from log file.")
        sys.exit(1)

    if not target_residues:
        logger.error("No valid residues parsed from the log file.")
        sys.exit(1)

    logger.info("Extracted %d target residues to restrain.", len(target_residues))
    return target_residues


class PDBResidueInfo:
    """Stores information about a residue parsed from a PDB file."""

    def __init__(self, chain_id: str, res_seq: int, i_code: str, res_name: str) -> None:
        """Initialize PDBResidueInfo.

        Args:
            chain_id: The chain ID of the residue.
            res_seq: The residue sequence number.
            i_code: The insertion code of the residue.
            res_name: The name of the residue.
        """
        self.chain_id = chain_id
        self.res_seq = res_seq
        self.i_code = i_code
        self.res_name = res_name


AMBER_MAP = {
    "HIE": "HIS",
    "HID": "HIS",
    "HIP": "HIS",
    "CYX": "CYS",
    "CYS": "CYS",
    "CYM": "CYS",
    "ASH": "ASP",
    "GLH": "GLU",
    "LYN": "LYS",
    "ARN": "ARG",
    "NMET": "MET",
    "CMET": "MET",
    "WAT": "HOH",
    "HOH": "HOH",
    "TIP3": "HOH",
    "SOL": "HOH",
}

AMBER_NAME_LEN = 4
PDB_COORD_LINE_MIN_LEN = 54


def normalize_res_name(name: str) -> str:
    """Normalize Amber residue names to standard residue names."""
    name = name.strip().upper()

    if len(name) > 1 and name[-1] in ("3", "5"):
        base_candidates = (
            "DA",
            "DT",
            "DG",
            "DC",
            "A",
            "U",
            "G",
            "C",
            "RA",
            "RU",
            "RG",
            "RC",
        )
        if name[:-1] in base_candidates:
            name = name[:-1]

    if name in AMBER_MAP:
        return AMBER_MAP[name]
    if len(name) == AMBER_NAME_LEN and name.startswith(("N", "C")):
        sub = name[1:]
        if sub in AMBER_MAP:
            return AMBER_MAP[sub]
        return sub
    return name


def parse_pdb_residues(pdb_path: Path) -> list[PDBResidueInfo]:
    """Parse residues sequentially from a PDB file."""
    residues = []
    seen = set()
    with pdb_path.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                if len(line) < PDB_COORD_LINE_MIN_LEN:
                    continue
                res_name = line[17:20].strip()
                chain_id = line[21]
                try:
                    res_seq = int(line[22:26])
                except ValueError:
                    res_seq = 0
                i_code = line[26]

                key = (chain_id, res_seq, i_code)
                if key not in seen:
                    seen.add(key)
                    residues.append(PDBResidueInfo(chain_id, res_seq, i_code, res_name))
    return residues


def align_residues(
    input_residues: list[PDBResidueInfo],
    leap_residues: list[PDBResidueInfo],
) -> dict[int, PDBResidueInfo]:
    """Align residues in the leap PDB to the original input PDB residues."""
    mapping = {}
    input_idx = 0
    num_input = len(input_residues)

    for leap_idx, leap_res in enumerate(leap_residues):
        leap_name = normalize_res_name(leap_res.res_name)
        matched = False
        for offset in range(30):
            curr_idx = input_idx + offset
            if curr_idx >= num_input:
                break
            input_res = input_residues[curr_idx]
            input_name = normalize_res_name(input_res.res_name)
            if leap_name == input_name:
                mapping[leap_idx] = input_res
                input_idx = curr_idx + 1
                matched = True
                break
        if not matched:
            logger.warning(
                "Could not match leap residue %s (index %d) to any input residue "
                "near index %d.",
                leap_res.res_name,
                leap_idx,
                input_idx,
            )
    return mapping


def setup_simulation(
    config: MinimizerConfig, target_residues: set[int]
) -> tuple[app.Simulation, app.Topology]:
    """Load pre-parametrized Amber system and apply restraints."""
    logger.info("Loading pre-parametrized Amber system...")
    prmtop = app.AmberPrmtopFile(str(config.prmtop_path))
    inpcrd = app.AmberInpcrdFile(str(config.inpcrd_path))

    pdb = app.PDBFile(str(config.leap_pdb_path))

    logger.info("Creating system from Amber topology...")
    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, constraints=None)

    logger.info("Aligning Leap topology residues with original PDB residues...")
    input_residues = parse_pdb_residues(config.input_pdb)
    leap_residues = parse_pdb_residues(config.leap_pdb_path)
    residue_mapping = align_residues(input_residues, leap_residues)

    openmm_force_unit = unit.kilojoules_per_mole / unit.nanometer**2
    k_phase1_val = config.k_phase1_all.value_in_unit(openmm_force_unit)

    target_force = openmm.CustomExternalForce("k_target*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    target_force.addGlobalParameter("k_target", k_phase1_val)
    target_force.addPerParticleParameter("x0")
    target_force.addPerParticleParameter("y0")
    target_force.addPerParticleParameter("z0")

    env_force = openmm.CustomExternalForce("k_env*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    env_force.addGlobalParameter("k_env", k_phase1_val)
    env_force.addPerParticleParameter("x0")
    env_force.addPerParticleParameter("y0")
    env_force.addPerParticleParameter("z0")

    logger.info("Assigning atoms to restraint potentials...")
    target_count, env_count = 0, 0

    positions = inpcrd.positions.value_in_unit(unit.nanometers)

    for atom in pdb.topology.atoms():
        if atom.element.symbol == "H":
            continue

        pos = positions[atom.index]

        original_res = residue_mapping.get(atom.residue.index)
        if original_res is not None:
            if original_res.chain_id == config.target_chain_id:
                if (
                    original_res.res_seq in target_residues
                    and atom.name in config.target_atoms
                ):
                    target_force.addParticle(atom.index, pos)
                    target_count += 1
            else:
                env_force.addParticle(atom.index, pos)
                env_count += 1
        else:
            env_force.addParticle(atom.index, pos)
            env_count += 1

    system.addForce(target_force)
    system.addForce(env_force)
    logger.info(
        "Restrained %d target atoms and %d environment atoms.",
        target_count,
        env_count,
    )

    integrator = openmm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)

    return simulation, pdb.topology


def main() -> None:
    """Execute the main minimization pipeline."""
    parser = argparse.ArgumentParser(
        description=("Protein minimization with dynamic GAFF2 ligand parameterization.")
    )
    parser.add_argument(
        "input_pdb",
        type=Path,
        help="Path to the input PDB file (e.g., protein_ligand.pdb).",
    )
    parser.add_argument(
        "--prmtop",
        type=Path,
        help=(
            "Path to the Amber prmtop topology file "
            "(default: output_tools/<input_pdb_stem>_parm.prmtop)."
        ),
    )
    parser.add_argument(
        "--inpcrd",
        type=Path,
        help=(
            "Path to the Amber inpcrd coordinate file "
            "(default: output_tools/<input_pdb_stem>_parm.inpcrd)."
        ),
    )
    parser.add_argument(
        "--leap-pdb",
        type=Path,
        help=(
            "Path to the Amber-generated PDB file "
            "(default: output_tools/<input_pdb_stem>_parm_leap.pdb)."
        ),
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("../modeling_results/logs/S2_automodel_1.log"),
        help="Path to the PRISM modeling log file specifying frozen residues.",
    )
    parser.add_argument(
        "--target-chain",
        type=str,
        default="A",
        help="Chain ID for the target protein.",
    )
    args = parser.parse_args()

    input_pdb = args.input_pdb
    output_pdb = input_pdb.with_name(f"{input_pdb.stem}_minimized{input_pdb.suffix}")

    prmtop_path = (
        args.prmtop
        if args.prmtop
        else OUTPUT_TOOLS_DIR / f"{input_pdb.stem}_parm.prmtop"
    )
    inpcrd_path = (
        args.inpcrd
        if args.inpcrd
        else OUTPUT_TOOLS_DIR / f"{input_pdb.stem}_parm.inpcrd"
    )
    leap_pdb_path = (
        args.leap_pdb
        if args.leap_pdb
        else OUTPUT_TOOLS_DIR / f"{input_pdb.stem}_parm_leap.pdb"
    )

    min_cfg = MinimizerConfig(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        log_file=args.log_file,
        prmtop_path=prmtop_path,
        inpcrd_path=inpcrd_path,
        leap_pdb_path=leap_pdb_path,
        target_chain_id=args.target_chain,
    )

    try:
        target_residues = extract_target_residues(min_cfg.log_file)
        simulation, final_topology = setup_simulation(min_cfg, target_residues)

        unit_str = unit.kilocalories_per_mole / unit.angstrom**2
        k_phase1_val = min_cfg.k_phase1_all.value_in_unit(unit_str)
        logger.info(
            "Phase 1: Local relaxation with K=%f for 500 steps...",
            k_phase1_val,
        )
        simulation.minimizeEnergy(maxIterations=500)

        k_target_val = min_cfg.k_phase2_target.value_in_unit(unit_str)
        k_env_val = min_cfg.k_phase2_env.value_in_unit(unit_str)

        logger.info(
            "Phase 2: Deep relaxation (Target K=%f, Env K=%f)...",
            k_target_val,
            k_env_val,
        )

        openmm_force_unit = unit.kilojoules_per_mole / unit.nanometer**2
        simulation.context.setParameter(
            "k_target", min_cfg.k_phase2_target.value_in_unit(openmm_force_unit)
        )
        simulation.context.setParameter(
            "k_env", min_cfg.k_phase2_env.value_in_unit(openmm_force_unit)
        )

        simulation.minimizeEnergy(
            tolerance=0.1 * unit.kilojoules_per_mole / unit.nanometer,
            maxIterations=0,
        )

        logger.info("Saving minimized coordinates...")
        min_pos = simulation.context.getState(getPositions=True).getPositions()

        with min_cfg.output_pdb.open("w", encoding="utf-8") as f:
            app.PDBFile.writeFile(final_topology, min_pos, f)

        logger.info("Process completed successfully. Saved to: %s", min_cfg.output_pdb)

    except Exception:
        logger.exception("Fatal error during execution.")
        sys.exit(1)


if __name__ == "__main__":
    main()
