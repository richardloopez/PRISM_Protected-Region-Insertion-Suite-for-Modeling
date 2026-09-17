#!/usr/bin/env python3
"""Core modeling engine for PRISM.

Provides custom MODELLER classes and execution logic for both homology
modeling (``FixedRegionAutoModel``) and loop refinement
(``FixedRegionLoopModel``).  The key innovation is per-residue coordinate
fixation combined with HETATM repulsion shields.
"""

from __future__ import annotations

import logging
import re
import shutil
from pathlib import Path
from typing import Any

from modeller import Alignment, Environ, features, forms, physical, pseudo_atom
from modeller.automodel import (
    AutoModel,
    DOPEHRLoopModel,
    assess,
    autosched,
    refine,
)
from modeller.parallel import Job
from modeller.selection import Selection

from . import config

logger = logging.getLogger("PRISM.modeling_engine")

# ============================================================================
#                        MODULE-LEVEL CONSTANTS
# ============================================================================

MIN_LOOP_LENGTH: int = 1
"""Minimum residue count for a valid loop range."""

MAX_LOOP_LENGTH: int = 1000
"""Maximum residue count for a valid loop range."""


def _get_residue_list(selection: Selection) -> list[dict[str, str]]:
    """Extract unique residue numbers from a MODELLER selection.

    Args:
        selection: A MODELLER ``Selection`` object.

    Returns:
        List of dictionaries with ``'res_num'`` keys for each unique
        residue in the selection.
    """
    residues: list[dict[str, str]] = []
    seen: set[str] = set()
    for sel_atom in selection:
        atom_str = str(sel_atom)
        try:
            _, _, res_num, _ = re.split(r"[: ]+", atom_str)
        except ValueError:
            logger.debug("Skipping unparseable atom string: %s", atom_str)
            continue
        if res_num not in seen:
            seen.add(res_num)
            residues.append({"res_num": res_num})
    return residues


# ============================================================================
#                           CUSTOM MODEL CLASSES
# ============================================================================


class FixedRegionAutoModel(AutoModel):
    """AutoModel subclass that freezes experimental residues during optimisation.

    Residues in ``experimental_residues`` and the entire BLK chain (if
    present) are excluded from the MODELLER objective function.  Per-residue
    HETATM repulsion shields are applied via ``nonstd_restraints``.
    """

    def __init__(
        self,
        env: Environ,
        experimental_residues: set[int],
        chain_id: str | None = None,
        blk_chain_id: str | None = None,
        **kwargs: Any,
    ) -> None:
        """Initialise the FixedRegionAutoModel.

        Args:
            env: MODELLER environment.
            experimental_residues: Set of 1-indexed residue numbers to freeze.
            chain_id: Protein chain identifier.  Defaults to
                ``config.CHAIN_ID`` when *None*.
            blk_chain_id: HETATM/ligand chain identifier.  Defaults to
                ``config.BLK_CHAIN_ID`` when *None*.  Pass ``"None"`` to
                skip BLK fixation explicitly.
            **kwargs: Forwarded to ``AutoModel.__init__``.
        """
        super().__init__(env, **kwargs)
        self.experimental_residues = experimental_residues
        self.chain_id = chain_id if chain_id is not None else config.CHAIN_ID
        self.blk_chain_id = (
            blk_chain_id if blk_chain_id is not None else config.BLK_CHAIN_ID
        )

    def select_atoms(self) -> Selection:
        """Select atoms for optimisation, **excluding** experimental residues and BLK.

        Returns:
            A ``Selection`` containing only the atoms that should be
            optimised by MODELLER.
        """
        all_atoms = Selection(self)
        fixed_selection_protein = Selection()
        fixed_selection_blk = Selection()

        if self.experimental_residues:
            for res_num in sorted(self.experimental_residues):
                curr_res = self.residue_range(
                    f"{res_num}:{self.chain_id}", f"{res_num}:{self.chain_id}"
                )
                fixed_selection_protein.add(curr_res)

        if self.blk_chain_id and str(self.blk_chain_id).lower() not in ("none", "null"):
            try:
                blk_chain_selection = Selection(self.chains[self.blk_chain_id])
                fixed_selection_blk.add(blk_chain_selection)
            except KeyError:
                logger.warning(
                    "[FixedRegionAutoModel] BLK chain '%s' "
                    "not found in model. Skipping BLK fixation.",
                    self.blk_chain_id,
                )

        optimizable = all_atoms - fixed_selection_protein - fixed_selection_blk
        logger.info(
            "[FixedRegionAutoModel] Optimizing %d atoms "
            "(Fixed Protein: %d atoms, BLK: %d atoms)",
            len(optimizable),
            len(fixed_selection_protein),
            len(fixed_selection_blk),
        )

        fixed_prot_ids = _get_residue_list(fixed_selection_protein)
        fixed_blk_ids = _get_residue_list(fixed_selection_blk)
        optimizable_ids = _get_residue_list(optimizable)

        logger.info("\n" + "=" * 80)
        logger.info("PRISM OPTIMIZATION SELECTION REPORT ([FixedRegionAutoModel])")
        logger.info("=" * 80)
        logger.info(
            "Frozen Protein Residues (Chain: %s) (Total: %d): %s",
            self.chain_id,
            len(fixed_prot_ids),
            ", ".join(item["res_num"] for item in fixed_prot_ids),
        )
        if fixed_blk_ids:
            logger.info(
                "Frozen BLK/HETATM Residues (Chain: %s) (Total: %d): %s",
                self.blk_chain_id,
                len(fixed_blk_ids),
                ", ".join(item["res_num"] for item in fixed_blk_ids),
            )
        else:
            logger.info("Frozen BLK/HETATM Residues: None")
        logger.info(
            "Mobile Optimizable Residues (Chain: %s) (Total: %d): %s",
            self.chain_id,
            len(optimizable_ids),
            ", ".join(item["res_num"] for item in optimizable_ids),
        )

        return optimizable

    def nonstd_restraints(self, aln: Alignment) -> None:
        """Add HETATM repulsion shield restraints.

        Args:
            aln: MODELLER alignment object (passed by the framework).
        """
        super().nonstd_restraints(aln)
        add_hetatm_repulsion_shield(self, config.BLOCK_REPULSION_RADIUS)


class FixedRegionLoopModel(DOPEHRLoopModel):
    """LoopModel subclass for loop refinement with frozen experimental residues.

    Validates that the loop range does not overlap with experimental
    coordinates, and applies per-residue repulsion shields to avoid clashes
    with bound HETATM groups.
    """

    def __init__(  # noqa: PLR0913
        self,
        env: Environ,
        inimodel: str,
        sequence: str,
        loop_start: int,
        loop_end: int,
        chain_id: str,
        experimental_residues: set[int],
        **kwargs: Any,
    ) -> None:
        """Initialise the FixedRegionLoopModel.

        Args:
            env: MODELLER environment.
            inimodel: Path to the initial model PDB.
            sequence: Target sequence code.
            loop_start: First residue of the loop (1-indexed).
            loop_end: Last residue of the loop (1-indexed).
            chain_id: Protein chain identifier.
            experimental_residues: Residues that must remain frozen.
            **kwargs: Forwarded to ``DOPEHRLoopModel.__init__``.

        Raises:
            ValueError: If the loop range overlaps with fixed residues.
        """
        super().__init__(env, inimodel=inimodel, sequence=sequence, **kwargs)
        self.loop_start = loop_start
        self.loop_end = loop_end
        self.chain_id = chain_id
        self.experimental_residues = experimental_residues
        logger.info(
            "[FixedRegionLoopModel] Initialized with experimental FIXED residues: %s",
            self.experimental_residues,
        )

        loop_res = set(range(loop_start, loop_end + 1))
        overlap = loop_res.intersection(experimental_residues)
        if overlap:
            raise ValueError(f"Loop overlap: {sorted(overlap)}")
        logger.info(
            "[FixedRegionLoopModel] Loop [%d-%d] does not overlap with fixed residues.",
            loop_start,
            loop_end,
        )

    def select_loop_atoms(self) -> Selection:
        """Select loop atoms for optimisation.

        Returns:
            A ``Selection`` spanning the configured loop range.
        """
        rng_start = f"{self.loop_start}:{self.chain_id}"
        rng_end = f"{self.loop_end}:{self.chain_id}"
        return Selection(self.residue_range(rng_start, rng_end))

    def nonstd_restraints(self, aln: Alignment) -> None:
        """Add HETATM repulsion shield restraints for loop atoms.

        Args:
            aln: MODELLER alignment object (passed by the framework).
        """
        super().nonstd_restraints(aln)
        add_hetatm_repulsion_shield(
            self, config.BLOCK_REPULSION_RADIUS, only_loop_atoms=True
        )


# ============================================================================
#                        HETATM REPULSION SHIELD
# ============================================================================


def add_hetatm_repulsion_shield(
    model: Any,
    min_dist: float,
    only_loop_atoms: bool = False,
) -> None:
    """Add lower-bound distance restraints between CA atoms and HETATM centres.

    Shared logic between ``FixedRegionAutoModel`` and
    ``FixedRegionLoopModel``.

    Args:
        model: A MODELLER model instance (AutoModel or LoopModel subclass).
        min_dist: Minimum allowed distance in Ångströms between a CA atom
            and a HETATM gravity centre.
        only_loop_atoms: If *True*, only apply restraints to the loop
            selection; otherwise apply to all optimisable atoms.
    """
    rsr = model.restraints

    het_residues = [r for r in model.residues if r.hetatm and r.name != "HOH"]
    if not het_residues:
        return

    target_sel = model.select_loop_atoms() if only_loop_atoms else model.select_atoms()

    target_ca = target_sel.only_atom_types("CA")
    if len(target_ca) == 0:
        return

    logger.info(
        "[add_hetatm_repulsion_shield] Adding repulsion: %d CA atoms vs %d HET groups.",
        len(target_ca),
        len(het_residues),
    )

    het_centers: list[Any] = []
    for res in het_residues:
        center = pseudo_atom.GravityCenter(Selection(res))
        rsr.pseudo_atoms.append(center)
        het_centers.append(center)

    count = 0
    for ca in target_ca:
        for center in het_centers:
            rsr.add(
                forms.LowerBound(
                    group=physical.xy_distance,
                    feature=features.Distance(ca, center),
                    mean=min_dist,
                    stdev=1.0,
                )
            )
            count += 1

    logger.info(
        "[add_hetatm_repulsion_shield] Added %d repulsion restraints (Min Dist: %s Å).",
        count,
        min_dist,
    )


# ============================================================================
#                        HOMOLOGY MODELING EXECUTION
# ============================================================================


def run_automodel(  # noqa: PLR0913
    env: Environ,
    align_file: str,
    job: Job,
    residues_to_freeze: set[int],
    start_model: int,
    end_model: int,
    knowns: list[str],
    input_mode: str,
) -> list[dict[str, Any]]:
    """Run AutoModel for homology modeling.

    Args:
        env: MODELLER environment.
        align_file: Path to the alignment file.
        job: MODELLER parallel Job object.
        residues_to_freeze: Residues to exclude from optimisation.
        start_model: First model index in this job slice.
        end_model: Last model index in this job slice.
        knowns: Template PDB codes.
        input_mode: One of ``'precalculation'``, ``'precomputed'``, or
            ``'normal'``.

    Returns:
        List of model output dictionaries from MODELLER, or an empty list
        if no models were generated.
    """
    if start_model > end_model:
        logger.warning(
            "[run_automodel] Start model %d > end model %d. Skipping.",
            start_model,
            end_model,
        )
        return []

    num_models = (end_model - start_model) + 1

    logger.info("\n" + "=" * 80)
    logger.info(
        "[run_automodel] Running AutoModel (Job range: %d-%d | Count: %d models)",
        start_model,
        end_model,
        num_models,
    )
    logger.info("[run_automodel] Templates: %s", knowns)
    logger.info("[run_automodel] Freezing residues: %s", residues_to_freeze)
    logger.info("=" * 80)

    if input_mode == "precalculation":
        logger.info("[run_automodel] PRECALCULATION MODE: Generating fresh inputs.")
        extra_inputs: dict[str, str] = {}
    elif input_mode == "precomputed":
        logger.info(
            "[run_automodel] Using precomputed files: %s, %s",
            config.CUSTOM_INIFILE_PATH,
            config.CUSTOM_RSRFILE_PATH,
        )
        extra_inputs = {
            "inifile": config.CUSTOM_INIFILE_PATH,
            "csrfile": config.CUSTOM_RSRFILE_PATH,
        }
    elif input_mode == "normal":
        logger.info("[run_automodel] Normal mode: Generating inputs as usual.")
        extra_inputs = {}
    else:
        logger.error("[run_automodel] Unknown input_mode: '%s'.", input_mode)
        return []

    a = FixedRegionAutoModel(
        env,
        experimental_residues=residues_to_freeze,
        chain_id=config.CHAIN_ID,
        alnfile=align_file,
        knowns=knowns,
        sequence=config.ALIGN_CODE_SEQUENCE,
        assess_methods=(assess.DOPEHR, assess.GA341),
        **extra_inputs,
    )
    a.use_parallel_job(job)
    a.starting_model = start_model
    a.ending_model = end_model

    a.library_schedule = autosched.slow
    a.max_var_iterations = 1000

    if input_mode == "precalculation":
        logger.info("[run_automodel] Running in PRECALCULATION mode (exit_stage=1).")
        a.make(exit_stage=1)

        generated_ini = f"{config.ALIGN_CODE_SEQUENCE}.ini"
        generated_rsr = f"{config.ALIGN_CODE_SEQUENCE}.rsr"

        ini_path = Path(generated_ini)
        rsr_path = Path(generated_rsr)
        if ini_path.exists() and rsr_path.exists():
            shutil.move(str(ini_path), config.CUSTOM_INIFILE_PATH)
            shutil.move(str(rsr_path), config.CUSTOM_RSRFILE_PATH)
            logger.info(
                "[run_automodel] Precalculation complete. "
                "Restraints (.rsr) and Initial (.ini) files generated."
            )
            return []
        logger.error(
            "[run_automodel] Precalculation failed. "
            "Restraints (.rsr) and Initial (.ini) files not generated."
        )
        return []
    a.make()

    if not a.outputs:
        logger.error("[run_automodel] No models were generated.")
        return []

    logger.info("[run_automodel] Generated %d models.", len(a.outputs))
    return a.outputs


# ============================================================================
#                        LOOP REFINEMENT EXECUTION
# ============================================================================


def run_loop_model(  # noqa: PLR0912
    env: Environ,
    job: Job,
    initial_models_names: list[str],
    loop_ranges: list[tuple[int, int]],
    experimental_residues: set[int],
) -> None:
    """Run LoopModel for loop refinement on selected models.

    Each model is refined sequentially across all valid loop ranges.  The
    best refinement output for each loop becomes the input for the next.

    Args:
        env: MODELLER environment.
        job: MODELLER parallel Job object.
        initial_models_names: PDB filenames of models to refine.
        loop_ranges: List of ``(start, end)`` residue ranges.
        experimental_residues: Residues that must remain frozen.
    """
    if not loop_ranges:
        logger.info(
            "[run_loop_model] No loop ranges provided. Skipping loop refinement."
        )
        return

    logger.info(
        "[run_loop_model] Validating loop candidates. "
        "Min loop length: %d | Max loop length: %d",
        MIN_LOOP_LENGTH,
        MAX_LOOP_LENGTH,
    )

    valid_loops: list[tuple[int, int]] = []
    invalid_loops: list[tuple[int, int]] = []
    for start, end in loop_ranges:
        length = end - start + 1
        if MIN_LOOP_LENGTH <= length <= MAX_LOOP_LENGTH:
            valid_loops.append((start, end))
        else:
            invalid_loops.append((start, end))

    if invalid_loops:
        invalid_str = ", ".join(f"[{s}-{e}]" for s, e in invalid_loops)
        logger.warning(
            "[run_loop_model] Skipping invalid loops (%d): %s",
            len(invalid_loops),
            invalid_str,
        )

    if valid_loops:
        valid_str = ", ".join(f"[{s}-{e}]" for s, e in valid_loops)
        logger.info(
            "[run_loop_model] Proceeding with valid loop ranges (%d): %s",
            len(valid_loops),
            valid_str,
        )

    if not valid_loops:
        logger.info(
            "[run_loop_model] No valid loop ranges remain. Aborting loop refinement."
        )
        return

    logger.info(
        "[run_loop_model] Loop refinement on %d models",
        len(initial_models_names),
    )
    logger.info(
        "[run_loop_model] Loop refinement for %d loops per model.",
        len(valid_loops),
    )

    for model_idx, pdb_file in enumerate(initial_models_names):
        base_name = Path(pdb_file).stem
        current_best_pdb = pdb_file

        logger.info(
            "[run_loop_model] Loop refinement for model %d: %s",
            model_idx + 1,
            pdb_file,
        )

        for j, (start, end) in enumerate(valid_loops):
            logger.info(
                "[run_loop_model] Loop refinement for loop %d: Residues %d-%d",
                j + 1,
                start,
                end,
            )

            ml = FixedRegionLoopModel(
                env,
                inimodel=current_best_pdb,
                sequence=config.ALIGN_CODE_SEQUENCE,
                loop_start=start,
                loop_end=end,
                chain_id=config.CHAIN_ID,
                experimental_residues=experimental_residues,
            )

            ml.use_parallel_job(job)
            ml.loop.starting_model = 1
            ml.loop.ending_model = config.LOOP_MODELS_PER_TARGET
            ml.loop.md_level = refine.slow_large

            ml.md_level = None

            ml.loop.assess_methods = (assess.DOPEHR, assess.GA341)
            ml.max_var_iterations = 1000

            ml.make()

            results = ml.loop.outputs
            if results:
                results.sort(key=lambda x: x.get("DOPE-HR score", 9e9))
                for m, info in enumerate(results):
                    old = info["name"]
                    new = f"{base_name}_LOOP{j + 1}_R{m + 1}.pdb"
                    old_path = Path(old)
                    if old_path.exists():
                        old_path.rename(new)
                    else:
                        logger.error(
                            "[run_loop_model] Loop model %s failed to generate. "
                            "Skipping.",
                            old,
                        )

                current_best_pdb = f"{base_name}_LOOP{j + 1}_R1.pdb"
            else:
                logger.error(
                    "[run_loop_model] No models generated for loop %d. "
                    "Keeping previous model.",
                    j + 1,
                )

            logger.info(
                "[run_loop_model] Loop refinement for model %d completed.",
                model_idx + 1,
            )


__all__ = [
    "FixedRegionAutoModel",
    "FixedRegionLoopModel",
    "add_hetatm_repulsion_shield",
    "run_automodel",
    "run_loop_model",
]
