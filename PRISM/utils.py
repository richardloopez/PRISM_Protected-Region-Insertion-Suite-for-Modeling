#!/usr/bin/env python3
"""Utility functions for the PRISM pipeline.

Handles environment setup, alignment processing, secondary structure
parsing, loop detection, and model ranking.
"""

from __future__ import annotations

import csv
import logging
import re
from pathlib import Path
from typing import Any

from modeller import Alignment, Environ, Model, log
from modeller.parallel import Job, LocalWorker
from modeller.scripts import complete_pdb
from modeller.selection import Selection

from . import config

logger = logging.getLogger("PRISM.utils")

# SS2 file format constants
MIN_SS2_COLUMNS = 3
SS2_CHAR_COL_INDEX = 2


# ============================================================================
#                        ENVIRONMENT SETUP
# ============================================================================


def setup_environment() -> tuple[Environ, Job]:
    """Initialise the MODELLER environment and parallel worker pool.

    Returns:
        A 2-tuple of ``(Environ, Job)`` ready for modeling.
    """
    logger.info("Current Working Directory: %s", Path.cwd())
    logger.info("\n" + "=" * 80)
    logger.info("PRISM: PROTECTED-REGION INSERTION SUITE FOR MODELING")
    logger.info("=" * 80 + "\n")

    env = Environ()

    logger.info("Loading MODELLER topology and parameter libraries...")
    env.libs.topology.read(file="$(LIB)/top_heav.lib")
    env.libs.parameters.read(file="$(LIB)/par.lib")
    logger.info(
        "[ENVIRONMENT] MODELLER topology and parameter libraries loaded successfully."
    )

    env.io.atom_files_directory = [
        ".",
        config.INPUT_DIR,
        "../atom_files",
        config.MODELING_RESULTS_DIR,
    ]
    env.io.hetatm = True

    log.verbose()

    env.jobs = config.MODELLER_CORES
    job = Job()
    logger.info("[ENVIRONMENT] Using %d local workers.", env.jobs)
    for _ in range(env.jobs):
        job.append(LocalWorker())
    job.start()

    return env, job


# ============================================================================
#                        ALIGNMENT GENERATION
# ============================================================================


def run_prereq_cde(env: Environ) -> None:
    """Generate the CDE-annotated alignment file.

    Depending on ``USE_MANUAL_ALIGNMENT``, either maps secondary structure
    onto a user-provided alignment or generates one automatically with
    MODELLER's ``salign()``.

    Args:
        env: MODELLER environment.

    Raises:
        FileNotFoundError: If the required alignment file is missing.
    """
    logger.info("[ENVIRONMENT] Running PREREQ CDE...")

    if config.USE_MANUAL_ALIGNMENT:
        logger.info("[ENVIRONMENT] Using manual alignment.")
        clean_ali = config.MANUAL_ALIGNMENT_FILE
        if not Path(clean_ali).exists():
            raise FileNotFoundError(f"Clean alignment file not found: {clean_ali}")

        logger.info("  > Input: %s", clean_ali)
        logger.info("  > Output: %s", config.MANUAL_ALIGNMENT_CDE_FILE)

        add_cde_line_to_pir(
            clean_ali_path=clean_ali,
            cde_ali_path=config.MANUAL_ALIGNMENT_CDE_FILE,
            ss2_path=config.SS2_FILE_PATH,
            target_seq=config.get_sequence(),
            align_code=config.ALIGN_CODE_SEQUENCE,
        )
    else:
        logger.info("[ENVIRONMENT] Using automatic alignment.")
        generate_pir_files_auto(
            env=env,
            ali_file_clean=config.ALIGNMENT_FILE,
            ali_file_cde=config.ALIGNMENT_CDE_FILE,
        )

    logger.info("[ENVIRONMENT] Alignment files generated successfully.")


def run_prerequisites(  # noqa: PLR0912, PLR0915
    env: Environ,  # noqa: ARG001
) -> tuple[str, list[tuple[int, int]], set[int], set[int]]:
    """Run all prerequisites for the automodeling stage.

    Determines which residues are experimental, identifies loops for
    refinement, and selects fixed residues based on the active mode
    (automatic or manual).

    Args:
        env: MODELLER environment.

    Returns:
        A 4-tuple of:
            - ``ali_file_automodel``: Path to the alignment file for AutoModel.
            - ``loop_ranges``: List of ``(start, end)`` residue ranges for
              loop refinement.
            - ``truly_fixed_residues``: Set of residue numbers that are
              unconditionally frozen.
            - ``automodel_fixed``: Set of residue numbers frozen during
              AutoModel (may include or exclude flanks).

    Raises:
        FileNotFoundError: If alignment files are missing.
    """
    logger.info("\n" + "=" * 80)
    logger.info("STEP 1: Alignment Preparation (Read-Only)")
    logger.info("=" * 80 + "\n")

    if config.USE_MANUAL_ALIGNMENT:
        ali_file_automodel = config.MANUAL_ALIGNMENT_FILE
        ali_file_analysis = config.MANUAL_ALIGNMENT_CDE_FILE
    else:
        ali_file_automodel = config.ALIGNMENT_FILE
        ali_file_analysis = config.ALIGNMENT_CDE_FILE

    if not Path(ali_file_automodel).exists():
        raise FileNotFoundError(f"Clean alignment file not found: {ali_file_automodel}")
    if not Path(ali_file_analysis).exists():
        raise FileNotFoundError(f"CDE alignment file not found: {ali_file_analysis}")

    seq_map = read_sequences_from_ali(ali_file_analysis)
    main_tmpl_seq = seq_map[config.MAIN_ALIGN_CODE_TEMPLATE]
    target_seq = seq_map[config.ALIGN_CODE_SEQUENCE]

    logger.info("[ENVIRONMENT] Experimental residue identification")
    experimental_residues = identify_experimental_residues(main_tmpl_seq, target_seq)

    logger.info("[ENVIRONMENT] Loop detection")

    # MANUAL OVERRIDE MODE
    if config.USE_MANUAL_OPTIMIZATION_SELECTION:
        logger.info("[MANUAL_MODE] Flag USE_MANUAL_OPTIMIZATION_SELECTION is True.")
        logger.info(
            "[MANUAL_MODE] Ignoring automatic flank detection and SS2 coil detection."
        )

        manual_list = set(config.MANUAL_OPTIMIZATION_RESIDUES)
        truly_fixed_residues = experimental_residues - manual_list
        loop_ranges = group_ranges(list(manual_list))
        automodel_fixed = truly_fixed_residues

        logger.info(
            "[MANUAL_MODE] User selected %d residues to optimize.",
            len(manual_list),
        )
        logger.info(
            "[MANUAL_FIXATION] Total fixed residues (Experimental + Manual): %d",
            len(experimental_residues),
        )

    # AUTOMATIC DETECTION MODE
    else:
        logger.info("[AUTO_MODE] Running automatic detection with connectivity checks.")

        raw_seq_str = config.get_sequence()
        full_seq_str = re.sub(r"[^A-Z]", "", raw_seq_str.upper())
        max_len = len(full_seq_str)

        all_coil = get_coil_residues(config.SS2_FILE_PATH, full_seq_str)
        experimental_ranges = group_ranges(list(experimental_residues))
        flank_residues: set[int] = set()
        n_flanks = config.MOBILE_FLANK_RESIDUES

        for start, end in experimental_ranges:
            if start > 1:
                for i in range(n_flanks):
                    if start + i <= end:
                        flank_residues.add(start + i)

            if end < max_len:
                for i in range(n_flanks):
                    if end - i >= start:
                        flank_residues.add(end - i)

        truly_fixed_residues = experimental_residues - flank_residues
        refinable_residues = all_coil - truly_fixed_residues
        loop_ranges = group_ranges(list(refinable_residues))

        if config.REFINE_FLANKS_DURING_AUTOMODEL:
            automodel_fixed = truly_fixed_residues
            logger.info("  > Fixed residues (Auto): %d", len(automodel_fixed))
        else:
            automodel_fixed = experimental_residues
            logger.info("  > Fixed residues (Exp Only): %d", len(automodel_fixed))

    # OUTPUT SUMMARY
    ranges_str = [f"[{s}-{e}]" for s, e in loop_ranges]
    if not ranges_str:
        logger.info("\n[ENVIRONMENT] No regions selected for optimization.")
    else:
        logger.info(
            "\n[ENVIRONMENT] Regions selected for optimization: %s",
            ", ".join(ranges_str),
        )

    return ali_file_automodel, loop_ranges, truly_fixed_residues, automodel_fixed


# ============================================================================
#                               HELPER FUNCTIONS
# ============================================================================


def flatten_ali_file(filepath: str | Path) -> None:
    """Flatten an alignment file by joining multi-line sequences.

    Sequence lines between header entries are concatenated into a single
    line, terminated by ``*``.

    Args:
        filepath: Path to the alignment file to flatten.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    lines = path.read_text(encoding="utf-8").splitlines()
    new_lines: list[str] = []
    seq_buffer: list[str] = []

    for line in lines:
        stripped = line.strip()
        if stripped.startswith((">P1;", "structure", "sequence", "#")):
            if seq_buffer:
                full_seq = "".join(seq_buffer).replace(" ", "")
                if not full_seq.endswith("*"):
                    full_seq += "*"
                new_lines.append(full_seq)
                seq_buffer = []
            new_lines.append(stripped)
        elif stripped:
            seq_buffer.append(stripped)

    if seq_buffer:
        full_seq = "".join(seq_buffer).replace(" ", "")
        if not full_seq.endswith("*"):
            full_seq += "*"
        new_lines.append(full_seq)

    path.write_text("\n".join(new_lines) + "\n")
    logger.info("[ENVIRONMENT] Ali file flattened: %s", filepath)


def read_ss2_file(ss2_path: str | Path, seq_full: str) -> str:
    """Read an SS2 file and return the secondary structure string.

    Args:
        ss2_path: Path to the SS2 file.
        seq_full: Full target sequence string.

    Returns:
        A string of the same length as *seq_full*, with one SS character
        per residue (``H``, ``E``, or ``C``).

    Raises:
        ValueError: If the SS2 string is longer than the sequence, or if
            the missing portion contains actual residues.
    """
    logger.info("[ENVIRONMENT] Reading ss2 file: %s", ss2_path)
    lines = Path(ss2_path).read_text(encoding="utf-8").splitlines()
    ss_parts = [
        line.split()[SS2_CHAR_COL_INDEX]
        for line in lines
        if (
            line.strip()
            and not line.startswith("#")
            and len(line.split()) >= MIN_SS2_COLUMNS
        )
    ]
    ss_string = "".join(ss_parts)

    len_ss, len_seq = len(ss_string), len(seq_full)

    if len_ss > len_seq:
        raise ValueError(f"SS2 length exceeds sequence: {len_ss} > {len_seq}")

    if len_ss < len_seq:
        missing = len_seq - len_ss
        if seq_full[len_ss:].strip("./"):
            raise ValueError("SS2 too short for target sequence.")
        logger.warning(
            "[ENVIRONMENT] Padding short SS2 with %d '.' characters.",
            missing,
        )
        ss_string += "." * missing

    return ss_string[:len_seq]


def read_sequences_from_ali(ali_file: str | Path) -> dict[str, str]:
    """Read a PIR alignment file and return a mapping of code → sequence.

    Args:
        ali_file: Path to the alignment file.

    Returns:
        Dictionary mapping sequence codes to their aligned sequences
        (gaps included).

    Raises:
        ValueError: If aligned sequences have different lengths.
    """
    sequences: dict[str, str] = {}
    current_code: str | None = None
    current_seq_parts: list[str] = []

    lines = Path(ali_file).read_text().splitlines()
    for raw_line in lines:
        stripped = raw_line.strip()
        if stripped.startswith(">P1;"):
            if current_code:
                sequences[current_code] = "".join(current_seq_parts)
            current_code = stripped.split(";")[1].strip()
            current_seq_parts = []
        elif current_code and not stripped.startswith(
            ("structure", "sequence", "CDE", "#")
        ):
            current_seq_parts.append(stripped)

    if current_code:
        sequences[current_code] = "".join(current_seq_parts)

    clean_map = {
        k: re.sub(r"[^A-Z\-\./]", "", v.upper().rstrip("*"))
        for k, v in sequences.items()
    }

    lengths = {len(s) for s in clean_map.values()}
    if len(lengths) > 1:
        raise ValueError(f"Inconsistent sequence lengths: {lengths}")

    return clean_map


def add_cde_line_to_pir(
    clean_ali_path: str | Path,
    cde_ali_path: str | Path,
    ss2_path: str | Path,
    target_seq: str,
    align_code: str,
) -> None:
    """Add a CDE (secondary structure annotation) line to a PIR file.

    The CDE line maps PSIPRED predictions onto the aligned target
    sequence, inserting ``'.'`` at gap positions.

    Args:
        clean_ali_path: Path to the clean alignment file (input).
        cde_ali_path: Path to the CDE alignment file (output).
        ss2_path: Path to the SS2 file.
        target_seq: Target sequence string.
        align_code: Alignment code identifying the target entry.

    Raises:
        ValueError: If *align_code* is not found in the alignment.
        IndexError: If the SS2 string is shorter than expected.
    """
    logger.info("[ENVIRONMENT] Adding CDE line to PIR file: %s", clean_ali_path)

    ss2_string = read_ss2_file(ss2_path, target_seq)
    aligned_seqs = read_sequences_from_ali(clean_ali_path)

    if align_code not in aligned_seqs:
        raise ValueError(
            f"Alignment code {align_code} not found in alignment file {clean_ali_path}"
        )

    align_target = aligned_seqs[align_code]
    cde_chars: list[str] = []
    ss_idx = 0
    gap_count = 0

    for char in align_target:
        if char in ("-", "/", "."):
            cde_chars.append(".")
            gap_count += 1
        elif ss_idx < len(ss2_string):
            cde_chars.append(ss2_string[ss_idx])
            ss_idx += 1
        else:
            raise IndexError("SS2 string is shorter than aligned target sequence")

    cde_line = "# CDE:" + "".join(cde_chars)

    input_lines = Path(clean_ali_path).read_text().splitlines()
    output_lines: list[str] = []
    sequences_counter = 0
    add_cde_sequences_counter = 0
    found_target = False

    for line in input_lines:
        output_lines.append(line)
        if line.strip().startswith(">P1;"):
            sequences_counter += 1
        if line.strip().startswith(f">P1;{align_code}"):
            found_target = True
        elif found_target and line.strip().startswith("sequence:"):
            output_lines.append(cde_line)
            add_cde_sequences_counter += 1
            found_target = False

    Path(cde_ali_path).write_text("\n".join(output_lines) + "\n")
    flatten_ali_file(cde_ali_path)
    logger.info("[ENVIRONMENT] CDE line added to PIR file: %s", cde_ali_path)
    logger.info(
        "[ENVIRONMENT] %d '.' inserted into CDE line "
        "(gap positions in aligned target sequence)",
        gap_count,
    )
    logger.info(
        "[ENVIRONMENT] %d sequences in PIR file: %s",
        sequences_counter,
        clean_ali_path,
    )
    logger.info(
        "[ENVIRONMENT] %d CDE sequences added to PIR file: %s | "
        "The sequence code is %s",
        add_cde_sequences_counter,
        cde_ali_path,
        align_code,
    )


def generate_pir_files_auto(
    env: Environ,
    ali_file_clean: str,
    ali_file_cde: str,
) -> None:
    """Generate PIR alignment files automatically via MODELLER salign.

    Args:
        env: MODELLER environment.
        ali_file_clean: Output path for the clean alignment file.
        ali_file_cde: Output path for the CDE-annotated alignment file.
    """
    logger.info(
        "[ENVIRONMENT] Generating PIR files for %s and %s",
        ali_file_clean,
        ali_file_cde,
    )
    aln = Alignment(env)

    for pdb_path, code in zip(
        config.PDB_TEMPLATE_FILES_PATHS, config.PDB_TEMPLATE_FILES_NAMES, strict=True
    ):
        logger.info(
            "[ENVIRONMENT] Adding template: %s (%s)",
            Path(pdb_path).name,
            code,
        )
        aln.append_model(
            mdl=Model(env, file=pdb_path),
            align_codes=code,
            atom_files=Path(pdb_path).name,
        )

    aln.append_sequence(config.get_sequence())
    aln[len(config.PDB_TEMPLATE_FILES_PATHS)].code = config.ALIGN_CODE_SEQUENCE

    aln.salign()
    aln.write(file=ali_file_clean, alignment_format="PIR")
    flatten_ali_file(ali_file_clean)

    add_cde_line_to_pir(
        ali_file_clean,
        ali_file_cde,
        config.SS2_FILE_PATH,
        config.get_sequence(),
        config.ALIGN_CODE_SEQUENCE,
    )


def prepare_prism_power_files(phase: str) -> tuple[list[str], str]:
    """Generate a virtual power alignment for the prism-power paradigm.

    Creates replica entries pointing to the same underlying PDB file,
    biasing MODELLER's coordinate averaging toward the experimental
    template.

    Args:
        phase: Either ``'precalculation'`` or ``'precomputed'``.

    Returns:
        A 2-tuple of:
            - ``all_expanded_knowns``: List of virtual template codes (the
              target sequence is excluded).
            - ``power_ali_path``: Path to the generated power alignment file.
    """
    power_map = getattr(config.PRISM_POWER_SETTINGS, phase.upper())
    input_dir = Path(config.INPUT_DIR)
    original_ali_path = input_dir / config.MANUAL_ALIGNMENT_BASENAME
    power_ali_path = input_dir / f"prism_power_{phase}.ali"

    all_expanded_knowns: list[str] = []

    with Path(original_ali_path).open(encoding="utf-8") as fh:
        content = fh.read()

    blocks = content.split(">P1;")
    header = blocks[0]
    ali_map: dict[str, list[str]] = {}
    for b in blocks[1:]:
        if not b.strip():
            continue
        lines = b.splitlines()
        code = lines[0].strip()
        ali_map[code] = lines[1:]

    new_ali_content = header

    for base_pdb, count in power_map.items():
        base_code: str | None = None
        if base_pdb in ali_map:
            base_code = base_pdb
        elif base_pdb.replace(".pdb", "") in ali_map:
            base_code = base_pdb.replace(".pdb", "")

        if not base_code:
            logger.error(
                "[ERROR] Template %s not found in alignment IDs: %s",
                base_pdb,
                list(ali_map.keys()),
            )
            continue

        original_lines = ali_map[base_code]

        for i in range(count):
            replica_id = f"{base_code}_{i:03d}" if i > 0 else base_code
            all_expanded_knowns.append(replica_id)

            block_lines = list(original_lines)
            if block_lines and block_lines[0].startswith("structure"):
                fields = block_lines[0].split(":")
                fields[1] = base_pdb
                block_lines[0] = ":".join(fields)

            new_ali_content += f">P1;{replica_id}\n" + "\n".join(block_lines) + "\n"

    target_code = config.ALIGN_CODE_SEQUENCE
    if target_code in ali_map:
        new_ali_content += (
            f">P1;{target_code}\n" + "\n".join(ali_map[target_code]) + "\n"
        )
    else:
        logger.error(
            "[ERROR] Target code '%s' not found in original alignment.", target_code
        )

    with Path(power_ali_path).open("w", encoding="utf-8") as fh:
        fh.write(new_ali_content)

    logger.info("[ENVIRONMENT] Virtual power alignment generated: %s", power_ali_path)
    logger.info(
        "%d virtual templates defined (Target excluded).",
        len(all_expanded_knowns),
    )
    logger.info("[ENVIRONMENT] Virtual templates defined: %s", all_expanded_knowns)

    return all_expanded_knowns, str(power_ali_path)


# ============================================================================
#                        RESIDUE ANALYSIS
# ============================================================================


def identify_experimental_residues(
    aligned_template_seq: str,
    aligned_target_seq: str,
) -> set[int]:
    """Identify target residues that map to experimental regions in the template.

    Positions where both template and target have actual residues (not gaps)
    indicate experimental coordinates that must be preserved.

    Args:
        aligned_template_seq: Aligned sequence from the main template.
        aligned_target_seq: Aligned sequence from the target protein.

    Returns:
        Set of 1-indexed residue numbers in the target that correspond to
        experimental positions in the template.

    Raises:
        ValueError: If aligned sequences have different lengths.
    """
    experimental_residues: set[int] = set()
    target_res_num = 0

    if len(aligned_template_seq) != len(aligned_target_seq):
        raise ValueError("Aligned sequence mismatch.")

    for template_res, target_res in zip(
        aligned_template_seq, aligned_target_seq, strict=False
    ):
        if target_res not in ("-", ".", "/"):
            target_res_num += 1
        if template_res not in ("-", ".", "/") and target_res not in ("-", ".", "/"):
            experimental_residues.add(target_res_num)

    if config.USE_MANUAL_FIXATION_SELECTION:
        logger.info("[MANUAL_FIXATION] Flag USE_MANUAL_FIXATION_SELECTION is True.")
        manual_fix = set(config.MANUAL_FIXATION_RESIDUES)
        if manual_fix:
            logger.info(
                "[MANUAL_FIXATION] Adding %d user-specified residues "
                "to experimental set.",
                len(manual_fix),
            )
            experimental_residues.update(manual_fix)
        else:
            logger.warning(
                "[MANUAL_FIXATION] Flag is True but MANUAL_FIXATION_RESIDUES list is "
                "empty."
            )

    logger.info(
        "[FIXED_REGION] Identified %d experimental residues mapped from template",
        len(experimental_residues),
    )
    logger.info(
        "[FIXED_REGION] These residues will NOT be optimized or refined "
        "(unless flank size takes them)"
    )

    return experimental_residues


def group_ranges(residues: list[int]) -> list[tuple[int, int]]:
    """Group a sorted list of residue numbers into contiguous ranges.

    The input list does not need to be pre-sorted or deduplicated;
    duplicates are removed and values are sorted internally before
    grouping.

    Args:
        residues: List of residue numbers (may be unsorted/duplicated).

    Returns:
        List of ``(start, end)`` tuples representing contiguous ranges.
    """
    if not residues:
        return []

    residues = sorted(set(residues))
    ranges: list[tuple[int, int]] = []
    start = end = residues[0]

    for res in residues[1:]:
        if res == end + 1:
            end = res
        else:
            ranges.append((start, end))
            start = end = res
    ranges.append((start, end))
    return ranges


def get_coil_residues(ss2_path: str, seq_full: str) -> set[int]:
    """Extract residue indices predicted as coil ('C') from the SS2 file.

    Args:
        ss2_path: Path to the SS2 file.
        seq_full: Full target sequence string.

    Returns:
        Set of 1-indexed residue numbers predicted as coil.

    Raises:
        ValueError: If the sequence and SS2 string have different lengths.
    """
    ss2_string = read_ss2_file(ss2_path, seq_full)

    if len(seq_full) != len(ss2_string):
        raise ValueError("SS2/Sequence length mismatch.")

    coil_residues: list[int] = []
    for i, char in enumerate(ss2_string):
        if char == "C":
            coil_residues.append(i + 1)
    return set(coil_residues)


# ============================================================================
#                            RANKING & EVALUATION
# ============================================================================


def run_rank_automodel_models(env: Environ) -> None:
    """Rank initial homology models by DOPE-HR score and rename them.

    Models matching ``<ALIGN_CODE_SEQUENCE>.B*.pdb`` are evaluated,
    sorted by score, and renamed to ``AUTO_1.pdb``, ``AUTO_2.pdb``, etc.
    The top ``TOP_MODELS_FOR_REFINEMENT`` are selected for loop refinement.

    Args:
        env: MODELLER environment.
    """
    logger.info("[AUTOMODEL_RANKING] Starting automodel ranking")

    cwd = Path.cwd()
    pattern = re.compile(
        rf"^{re.escape(config.ALIGN_CODE_SEQUENCE)}\.B[0-9]{{5,}}\.pdb$"
    )

    raw_models = [f for f in cwd.iterdir() if f.is_file() and pattern.match(f.name)]

    if not raw_models:
        logger.error(
            "[ERROR] No models found in %s (Pattern: %s.B*.pdb).",
            cwd,
            config.ALIGN_CODE_SEQUENCE,
        )
        return
    logger.info(
        "[AUTOMODEL_RANKING] Found %d models. Beginning ranking...",
        len(raw_models),
    )

    results: list[dict[str, Any]] = []
    for model_path in raw_models:
        mdl = complete_pdb(env, str(model_path))
        score = Selection(mdl.chains[config.CHAIN_ID]).assess_dopehr()
        results.append({"path": model_path, "score": score})
        logger.info(
            "[AUTOMODEL_RANKING] Evaluated -> %s: | DOPEHR score: %s",
            model_path,
            score,
        )

    results.sort(key=lambda x: x["score"])
    logger.info("[AUTOMODEL_RANKING] Ranked %d models... renaming", len(results))

    selected_for_refinement: list[str] = []

    for rank, data in enumerate(results, 1):
        old_path = data["path"]
        new_name = f"AUTO_{rank}.pdb"
        new_path = cwd / new_name

        old_path.rename(new_path)
        if rank <= config.TOP_MODELS_FOR_REFINEMENT:
            selected_for_refinement.append(new_name)

    logger.info("Ranking complete")
    logger.info("Selected %d models for refinement", len(selected_for_refinement))
    for name in selected_for_refinement:
        logger.info("Selected -> %s", name)


def final_evaluation_and_ranking(
    env: Environ,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Perform final evaluation and ranking of all models.

    Assesses all ``AUTO_*.pdb`` and ``*LOOP*.pdb`` files using DOPE-HR
    and normalised DOPE-HR Z-scores, then writes the results to
    ``final_ranking.csv``.

    Args:
        env: MODELLER environment.

    Returns:
        A 2-tuple of:
            - Full sorted results list.
            - Dictionary for the best model (or empty dict if none).
    """
    logger.info("[FINAL_EVALUATION] Starting final evaluation and ranking")

    template_names = set(config.PDB_TEMPLATE_FILES_NAMES)
    pdbs = [
        p
        for p in Path().glob("*.pdb")
        if p.name not in template_names and ("AUTO_" in p.name or "LOOP" in p.name)
    ]
    if not pdbs:
        logger.warning("No valid PDB files found for final evaluation.")
        return [], {}

    results: list[dict[str, Any]] = []
    for pdb_path in pdbs:
        mdl = complete_pdb(env, str(pdb_path))
        score = Selection(mdl.chains[config.CHAIN_ID]).assess_dopehr()
        zscore = mdl.assess_normalized_dopehr()
        results.append(
            {
                "name": pdb_path,
                "DOPEHR_score": score,
                "DOPEHR_zscore": zscore,
            }
        )
        logger.info(
            "[FINAL_EVALUATION] %s: DOPEHR score: %s, DOPEHR zscore: %s",
            pdb_path,
            score,
            zscore,
        )

    results.sort(key=lambda x: x["DOPEHR_score"])
    val = config.NUM_BEST_FINAL_MODELS
    limit = None if val == float("inf") else int(val)
    best_models = results[:limit]

    with Path(config.FINAL_RANKING_CSV).open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["Rank", "Model_Name", "DOPEHR_score", "DOPEHR_zscore"]
        )
        writer.writeheader()
        for rank, data in enumerate(best_models, start=1):
            writer.writerow(
                {
                    "Rank": rank,
                    "Model_Name": data["name"],
                    "DOPEHR_score": data["DOPEHR_score"],
                    "DOPEHR_zscore": data["DOPEHR_zscore"],
                }
            )
    logger.info(
        "[FINAL_EVALUATION] Best models ranking saved to %s",
        config.FINAL_RANKING_CSV,
    )

    return results, (best_models[0] if best_models else {})


__all__ = [
    "add_cde_line_to_pir",
    "final_evaluation_and_ranking",
    "flatten_ali_file",
    "generate_pir_files_auto",
    "get_coil_residues",
    "group_ranges",
    "identify_experimental_residues",
    "prepare_prism_power_files",
    "read_sequences_from_ali",
    "read_ss2_file",
    "run_prereq_cde",
    "run_prerequisites",
    "run_rank_automodel_models",
    "setup_environment",
]
