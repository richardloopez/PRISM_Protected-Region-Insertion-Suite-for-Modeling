#!/usr/bin/env python3
"""Automated Protein-Ligand Parameterization Utility.

Extracts the target protein chain and specified ligand chains/residues,
runs reduce, antechamber, parmchk2, and tleap to generate Amber parameter
and coordinate files dynamically named after the input PDB.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

try:
    from pdb_utils import (
        OUTPUT_TOOLS_DIR,
        extract_chain_pdb,
        extract_protein,
        extract_residue_pdb,
        get_residues,
    )
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from pdb_utils import (
        OUTPUT_TOOLS_DIR,
        extract_chain_pdb,
        extract_protein,
        extract_residue_pdb,
        get_residues,
    )

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("parameterizer")


def parse_charges(charges_str: str) -> dict[str, int]:
    """Parse ligand charges mapping from string like ATP:-4,XXX:3.

    Args:
        charges_str: A string mapping residue names to integer charges.

    Returns:
        A dictionary mapping residue names (in uppercase) or 'DEFAULT'
        to their integer charges.

    Raises:
        ValueError: If the charge format is invalid or cannot be parsed.
    """
    charges = {}
    if not charges_str:
        return charges
    for raw_part in charges_str.split(","):
        part = raw_part.strip()
        if not part:
            continue
        if ":" not in part:
            try:
                val = int(part)
                charges["DEFAULT"] = val
            except ValueError as err:
                raise ValueError(
                    f"Invalid charge format (expected RES:CHARGE): {part}"
                ) from err
        else:
            res_name, val_str = part.split(":", 1)
            charges[res_name.strip().upper()] = int(val_str.strip())
    return charges


def _resolve_executable(name: str) -> str:
    """Resolve the absolute path of an executable.

    Args:
        name: The name of the executable to search for.

    Returns:
        The absolute path of the executable.

    Raises:
        FileNotFoundError: If the executable is not found in the system PATH.
    """
    path = shutil.which(name)
    if not path:
        raise FileNotFoundError(f"Required executable '{name}' not found in PATH.")
    return path


def run_cmd(
    cmd: list[str], cwd: Path | None = None
) -> subprocess.CompletedProcess[str]:
    """Run system command and handle failure.

    Args:
        cmd: The command to execute as a list of strings.
        cwd: The working directory for the command.

    Returns:
        The completed process object.

    Raises:
        subprocess.CalledProcessError: If the command returns a non-zero exit code.
    """
    logger.info("Executing: %s", " ".join(cmd))
    res = subprocess.run(  # noqa: S603
        cmd,
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )
    if res.returncode != 0:
        logger.error("Command failed: %s", " ".join(cmd))
        logger.error("STDOUT:\n%s", res.stdout)
        logger.error("STDERR:\n%s", res.stderr)
        raise subprocess.CalledProcessError(
            res.returncode, cmd, output=res.stdout, stderr=res.stderr
        )
    return res


def _parameterize_single_ligand(
    res_name: str,
    work_dir: Path,
    charges: dict[str, int],
    var_name: str,
    ligand_base: str,
) -> dict[str, Any]:
    """Parameterize a single ligand structure using reduce, antechamber, and parmchk2.

    Args:
        res_name: Residue name of the ligand.
        work_dir: Path to the working directory.
        charges: Mapping of residue names to charges.
        var_name: Variable name to use in tleap.
        ligand_base: Base filename for output files.

    Returns:
        A dictionary containing ligand records for tleap.
    """
    ligand_pdb = work_dir / f"{ligand_base}.pdb"

    # Count atoms in the residue/chain to check if it's a single atom (e.g. ion)
    with ligand_pdb.open(encoding="utf-8") as f:
        atom_lines = [line for line in f if line.startswith(("ATOM", "HETATM"))]

    is_single_atom = len(atom_lines) == 1

    if is_single_atom:
        logger.info(
            "Residue/Chain %s (var %s) has only 1 atom. "
            "Skipping reduce, antechamber, and parmchk2.",
            res_name,
            var_name,
        )
        return {
            "var_name": var_name,
            "is_single_atom": True,
            "pdb": ligand_pdb.name,
        }

    ligand_reduced_pdb = work_dir / f"{ligand_base}_reduced.pdb"
    ligand_mol2 = work_dir / f"{ligand_base}.mol2"
    ligand_frcmod = work_dir / f"{ligand_base}.frcmod"

    # Run reduce to add hydrogens
    try:
        reduce_bin = _resolve_executable("reduce")
        with ligand_reduced_pdb.open("w", encoding="utf-8") as out_f:
            res = subprocess.run(  # noqa: S603
                [reduce_bin, "-NOFLIP", str(ligand_pdb)],
                stdout=out_f,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
        if not ligand_reduced_pdb.exists() or ligand_reduced_pdb.stat().st_size == 0:
            logger.error("reduce failed to produce any output. stderr:\n%s", res.stderr)
            sys.exit(1)
        logger.info("Reduced ligand using reduce: %s", ligand_reduced_pdb.name)
    except Exception:
        logger.exception("Failed to run reduce")
        raise

    charge = charges.get(res_name.upper(), charges.get("DEFAULT", 0))
    logger.info("Using charge %d for residue %s", charge, res_name)

    # Run antechamber
    antechamber_bin = _resolve_executable("antechamber")
    run_cmd(
        [
            antechamber_bin,
            "-i",
            ligand_reduced_pdb.name,
            "-fi",
            "pdb",
            "-o",
            ligand_mol2.name,
            "-fo",
            "mol2",
            "-c",
            "bcc",
            "-at",
            "gaff2",
            "-nc",
            str(charge),
            "-rn",
            res_name,
        ],
        cwd=work_dir,
    )

    # Run parmchk2
    parmchk2_bin = _resolve_executable("parmchk2")
    run_cmd(
        [
            parmchk2_bin,
            "-i",
            ligand_mol2.name,
            "-f",
            "mol2",
            "-o",
            ligand_frcmod.name,
        ],
        cwd=work_dir,
    )

    return {
        "var_name": var_name,
        "is_single_atom": False,
        "mol2": ligand_mol2.name,
        "frcmod": ligand_frcmod.name,
    }


def _extract_and_parameterize_ligands(
    input_pdb: Path,
    work_dir: Path,
    ligand_chains_list: list[str],
    keep_intact_chains: list[str],
    charges: dict[str, int],
) -> list[dict[str, Any]]:
    """Extract and parameterize all ligand chains from the input PDB.

    Args:
        input_pdb: Path to the input PDB file.
        work_dir: Path to the working directory.
        ligand_chains_list: List of chain IDs containing ligands.
        keep_intact_chains: List of chain IDs to keep intact.
        charges: Mapping of residue names to charges.

    Returns:
        A list of ligand record dictionaries.
    """
    ligand_records = []
    for chain in ligand_chains_list:
        if chain in keep_intact_chains:
            # Process the entire chain intact as a single ligand/molecule
            temp_ligand_pdb = work_dir / f"temp_ligand_{chain}.pdb"
            res_name = extract_chain_pdb(input_pdb, temp_ligand_pdb, chain)
            if not temp_ligand_pdb.exists() or temp_ligand_pdb.stat().st_size == 0:
                logger.warning("No residues found in ligand chain %s", chain)
                if temp_ligand_pdb.exists():
                    temp_ligand_pdb.unlink()
                continue

            ligand_base = f"ligand_{chain}_{res_name}"
            ligand_pdb = work_dir / f"{ligand_base}.pdb"

            if ligand_pdb.exists():
                ligand_pdb.unlink()
            temp_ligand_pdb.rename(ligand_pdb)
            logger.info(
                "Extracted entire chain %s (representative res_name %s) to %s",
                chain,
                res_name,
                ligand_pdb.name,
            )

            record = _parameterize_single_ligand(
                res_name=res_name,
                work_dir=work_dir,
                charges=charges,
                var_name=f"LIG_{chain}",
                ligand_base=ligand_base,
            )
            ligand_records.append(record)
        else:
            # Default behavior: split chain into individual residues
            unique_res = get_residues(input_pdb, chain)
            if not unique_res:
                logger.warning("No residues found in ligand chain %s", chain)
                continue

            for res_seq_str, i_code in unique_res:
                res_seq = int(res_seq_str)

                temp_ligand_pdb = work_dir / f"temp_ligand_{chain}_{res_seq}.pdb"
                res_name = extract_residue_pdb(
                    input_pdb,
                    temp_ligand_pdb,
                    chain,
                    res_seq,
                    i_code,
                )

                suffix = f"_{i_code.strip()}" if i_code.strip() else ""
                ligand_base = f"ligand_{chain}_{res_seq}{suffix}_{res_name}"
                ligand_pdb = work_dir / f"{ligand_base}.pdb"

                if ligand_pdb.exists():
                    ligand_pdb.unlink()
                temp_ligand_pdb.rename(ligand_pdb)
                logger.info(
                    "Extracted residue %s (Chain %s, seq %s%s) to %s",
                    res_name,
                    chain,
                    res_seq,
                    suffix,
                    ligand_pdb.name,
                )

                record = _parameterize_single_ligand(
                    res_name=res_name,
                    work_dir=work_dir,
                    charges=charges,
                    var_name=f"LIG_{chain}_{res_seq}{suffix}",
                    ligand_base=ligand_base,
                )
                ligand_records.append(record)

    return ligand_records


def _run_tleap(  # noqa: PLR0913, PLR0912, PLR0915
    stem: str,
    work_dir: Path,
    protein_pdb: Path,
    ligand_records: list[dict[str, Any]],
    extra_preps: list[Path] | None = None,
    extra_frcmods: list[Path] | None = None,
    extra_offs: list[Path] | None = None,
) -> None:
    """Run tleap to generate Amber parameter and coordinate files.

    Args:
        stem: Naming stem for output files.
        work_dir: Path to the working directory.
        protein_pdb: Path to the protein PDB file.
        ligand_records: List of records for parameterized ligands.
        extra_preps: List of paths to extra Amber prep files to load.
        extra_frcmods: List of paths to extra frcmod files to load.
        extra_offs: List of paths to extra off files to load.
    """
    tleap_in = work_dir / f"{stem}_parm.in"
    tleap_lines = [
        "source leaprc.protein.ff19SB",
        "source leaprc.DNA.OL21",
        "source leaprc.RNA.OL3",
        "source leaprc.gaff2",
        "source leaprc.water.tip3p",
    ]

    if extra_preps:
        for prep in extra_preps:
            prep_path = prep.resolve()
            if prep_path.exists():
                shutil.copy2(prep_path, work_dir)
                tleap_lines.append(f"loadamberprep {prep_path.name}")
            else:
                logger.warning("Extra prep file not found: %s", prep_path)

    if extra_frcmods:
        for frcmod in extra_frcmods:
            frcmod_path = frcmod.resolve()
            if frcmod_path.exists():
                shutil.copy2(frcmod_path, work_dir)
                tleap_lines.append(f"loadamberparams {frcmod_path.name}")
            else:
                logger.warning("Extra frcmod file not found: %s", frcmod_path)

    if extra_offs:
        for off in extra_offs:
            off_path = off.resolve()
            if off_path.exists():
                shutil.copy2(off_path, work_dir)
                tleap_lines.append(f"loadoff {off_path.name}")
            else:
                logger.warning("Extra off file not found: %s", off_path)

    for lig in ligand_records:
        if lig["is_single_atom"]:
            tleap_lines.append(f"{lig['var_name']} = loadpdb {lig['pdb']}")
        else:
            tleap_lines.append(f"loadamberparams {lig['frcmod']}")
            tleap_lines.append(f"{lig['var_name']} = loadmol2 {lig['mol2']}")

    tleap_lines.append(f"PROT = loadpdb {protein_pdb.name}")

    if ligand_records:
        combine_vars = ["PROT"] + [lig["var_name"] for lig in ligand_records]
        tleap_lines.append(f"COMPLEX = combine {{{' '.join(combine_vars)}}}")
    else:
        tleap_lines.append("COMPLEX = PROT")

    tleap_lines.append(f"saveamberparm COMPLEX {stem}_parm.prmtop {stem}_parm.inpcrd")
    tleap_lines.append(f"savepdb COMPLEX {stem}_parm_leap.pdb")
    tleap_lines.append("quit")

    with tleap_in.open("w", encoding="utf-8") as f:
        f.write("\n".join(tleap_lines) + "\n")

    tleap_bin = _resolve_executable("tleap")
    run_cmd([tleap_bin, "-f", tleap_in.name], cwd=work_dir)

    expected_outputs = [
        f"{stem}_parm.prmtop",
        f"{stem}_parm.inpcrd",
        f"{stem}_parm_leap.pdb",
        f"{stem}_parm.in",
    ]
    for out_file in expected_outputs:
        src_path = work_dir / out_file
        if src_path.exists():
            dest_path = OUTPUT_TOOLS_DIR / out_file
            shutil.copy2(src_path, dest_path)
            logger.info("Copied output to: %s", dest_path)
        else:
            logger.error("Expected output file not generated: %s", out_file)
            sys.exit(1)


def main() -> None:
    """Automate Amber parameterization for protein and small molecules."""
    parser = argparse.ArgumentParser(
        description="Automate Amber parameterization for protein and ligand structures."
    )
    parser.add_argument(
        "input_pdb",
        type=Path,
        help="Path to the input PDB file (e.g., input.pdb).",
    )
    parser.add_argument(
        "--target-chain",
        type=str,
        default="A",
        help="Chain ID for the target protein (default: A).",
    )
    parser.add_argument(
        "--ligand-chains",
        type=str,
        default="",
        help=(
            "Comma-separated list of chain IDs containing ligands "
            "(e.g. --ligand-chains B,C)."
        ),
    )
    parser.add_argument(
        "--ligand-charges",
        type=str,
        default="",
        help="Comma-separated residue name to charge mapping (e.g., ATP:-4,XXX:3).",
    )
    parser.add_argument(
        "--keep-intact-chains",
        type=str,
        default="",
        help=(
            "Comma-separated list of ligand chain IDs that should not "
            "be split into individual residues (e.g. --keep-intact-chains B)."
        ),
    )
    parser.add_argument(
        "--extra-preps",
        type=str,
        default="",
        help=(
            "Comma-separated list of extra Amber prep files "
            "(e.g. extra1.prep,extra2.prep)."
        ),
    )
    parser.add_argument(
        "--extra-frcmods",
        type=str,
        default="",
        help=(
            "Comma-separated list of extra frcmod files "
            "(e.g. extra1.frcmod,extra2.frcmod)."
        ),
    )
    parser.add_argument(
        "--extra-offs",
        type=str,
        default="",
        help=("Comma-separated list of extra off files (e.g. extra1.off,extra2.off)."),
    )

    args = parser.parse_args()

    if not args.input_pdb.exists():
        logger.error("Input PDB file not found: %s", args.input_pdb)
        sys.exit(1)

    charges = parse_charges(args.ligand_charges)
    logger.info("Loaded ligand charges: %s", charges)

    stem = args.input_pdb.stem
    logger.info("Using naming stem: %s_parm", stem)

    OUTPUT_TOOLS_DIR.mkdir(parents=True, exist_ok=True)
    work_dir = OUTPUT_TOOLS_DIR / "parameterizer"
    work_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Intermediate work directory: %s", work_dir)

    ligand_chains_list = [c.strip() for c in args.ligand_chains.split(",") if c.strip()]
    keep_intact_chains = [
        c.strip() for c in args.keep_intact_chains.split(",") if c.strip()
    ]

    # Ensure any chain in keep_intact_chains is also treated as a ligand chain
    for chain in keep_intact_chains:
        if chain not in ligand_chains_list:
            ligand_chains_list.append(chain)

    # 1. Extract protein (and ligands with dedicated forcefields)
    protein_pdb = work_dir / "protein.pdb"
    extract_protein(
        args.input_pdb,
        protein_pdb,
        args.target_chain,
        exclude_chains=ligand_chains_list,
    )

    # 2. Extract and parameterize ligands
    ligand_records = _extract_and_parameterize_ligands(
        input_pdb=args.input_pdb,
        work_dir=work_dir,
        ligand_chains_list=ligand_chains_list,
        keep_intact_chains=keep_intact_chains,
        charges=charges,
    )

    extra_preps = [Path(p.strip()) for p in args.extra_preps.split(",") if p.strip()]
    extra_frcmods = [
        Path(f.strip()) for f in args.extra_frcmods.split(",") if f.strip()
    ]
    extra_offs = [Path(o.strip()) for o in args.extra_offs.split(",") if o.strip()]

    # 3. tleap execution
    _run_tleap(
        stem,
        work_dir,
        protein_pdb,
        ligand_records,
        extra_preps=extra_preps,
        extra_frcmods=extra_frcmods,
        extra_offs=extra_offs,
    )

    logger.info("Parameterization completed successfully!")


if __name__ == "__main__":
    main()
