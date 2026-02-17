#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM (modeling engine)

'''
Core modeling engine for PRISM.

This module combines custom Modeller classes and execution logic for both
homology modeling (AutoModel) and loop refinement (DOPEHRLoopModel).
'''

import os 
import shutil
import logging
from typing import List, Tuple, Set, Dict, Any, Union
import re

from modeller import *
from modeller.automodel import *
from modeller.selection import Selection
from modeller.parallel import Job

from . import config


# Setup logging
logger = logging.getLogger("PRISM.modeling_engine")

# ============================================================================
#                           CUSTOM MODEL CLASSES
# ============================================================================

class FixedRegionAutoModel(AutoModel):
    '''
    Custom AutoModel class for homology modeling with fixed experimental residues.
    Applies per-residue repulsion shields to avoid clashes with heteroatoms.
    '''
    DEFAULT_CHAIN = config.CHAIN_ID
    def __init__(self, env: Environ, experimental_residues: Set[int], chain_id: str = DEFAULT_CHAIN, blk_chain_id: str = config.BLK_CHAIN_ID, **kwargs) -> None:
        '''
        Initialize the FixedRegionAutoModel class.
        '''
        super().__init__(env, **kwargs)
        self.experimental_residues = experimental_residues
        self.chain_id = chain_id
        self.blk_chain_id = blk_chain_id
    

    def select_atoms(self) -> Selection:
        '''
        Select atoms for optimization, EXCLUDING experimental residues and BLK.
        '''
        all_atoms = Selection(self)
        fixed_selection_protein = Selection()
        fixed_selection_blk = Selection()

        if self.experimental_residues:
            for res_num in sorted(self.experimental_residues):
                curr_res = self.residue_range(f'{res_num}:{self.chain_id}', f'{res_num}:{self.chain_id}')
                fixed_selection_protein.add(curr_res)
        if self.blk_chain_id is not None:
            blk_chain_selection = Selection(self.chains[self.blk_chain_id])
            fixed_selection_blk.add(blk_chain_selection)

        optimizable = all_atoms - fixed_selection_protein - fixed_selection_blk
        logger.info(f"[ENVIRONMENT][FixedRegionAutoModel] Optimizing {len(optimizable)} atoms (Fixed Protein: {len(fixed_selection_protein)} atoms, BLK: {len(fixed_selection_blk)} atoms)")
        
    
        def get_residue_list(selection):
            residues = []
            seen = set()
            for atom in selection:
                atom_str = str(atom)
                atom, atom_name, res_num, chain_id = re.split(r'[: ]+', atom_str)
                if res_num not in seen:
                    seen.add(res_num)
                    residues.append({
                        "res_num": res_num
                    })
            return residues
        fixed_prot_ids = get_residue_list(fixed_selection_protein)
        fixed_blk_ids = get_residue_list(fixed_selection_blk)
        optimizable_ids = get_residue_list(optimizable)
        
        logger.info("\n" + "="*80)
        logger.info("PRISM OPTIMIZATION SELECTION REPORT ([ENVIRONMENT][FixedRegionAutoModel])")
        ("\n" + "="*80)
        logger.info(f"Frozen Protein Residues (Chain: {self.chain_id}) (Total: {len(fixed_prot_ids)}): {', '.join(item['res_num'] for item in fixed_prot_ids)}")
        logger.info(f"Frozen BLK/HETATM Residues (Chain: {self.blk_chain_id}) (Total: {len(fixed_blk_ids)}): {', '.join(item['res_num'] for item in fixed_blk_ids)}")
        logger.info(f"Mobile Optimizable Residues (Chain: {self.chain_id}) (Total: {len(optimizable_ids)}): {', '.join(item['res_num'] for item in optimizable_ids)}")

        return optimizable


    def nonstd_restraints(self, aln: Alignment) -> None:
        '''
        Add HETATM repulsion shield restraints.
        '''
        super().nonstd_restraints(aln)
        add_hetatm_repulsion_shield(self, config.BLOCK_REPULSION_RADIUS)


class FixedRegionLoopModel(DOPEHRLoopModel):
    '''
    Custom LoopModel class for loop refinement with fixed experimental residues.
    Applies per-residue repulsion shields to avoid clashes with heteroatoms.
    '''
    def __init__(self, env: Environ, inimodel: str, sequence: str, loop_start: int, loop_end: int, chain_id: str, experimental_residues: Set[int], **kwargs) -> None:
        '''
        Initialize the LoopModel with fixed experimental residues.
        '''
        super().__init__(env, inimodel=inimodel, sequence=sequence, **kwargs)
        self.loop_start = loop_start
        self.loop_end = loop_end
        self.chain_id = chain_id
        self.experimental_residues = experimental_residues
        logger.info(f"[ENVIRONMENT][FixedRegionLoopModel] Initialized LoopModel with Experimental FIXED residues: {self.experimental_residues}")

        loop_res = set(range(loop_start, loop_end + 1))
        overlap = loop_res.intersection(experimental_residues)
        if overlap:
            raise ValueError(f"[ENVIRONMENT][FixedRegionLoopModel] Loop [{loop_start}-{loop_end}] overlaps with fixed residues {sorted(overlap)}.")
        logger.info(f"[ENVIRONMENT][FixedRegionLoopModel] Loop [{loop_start}-{loop_end}] does not overlap with fixed residues.")
        

    def select_loop_atoms(self) -> Selection:
        '''
        Select loop atoms for optimization.
        '''
        rng_start = f'{self.loop_start}:{self.chain_id}'
        rng_end = f'{self.loop_end}:{self.chain_id}'
        return Selection(self.residue_range(rng_start, rng_end))
        

    def nonstd_restraints(self, aln: Alignment) -> None:
        '''
        Add HETATM repulsion shield restraints.
        '''
        super().nonstd_restraints(aln)
        add_hetatm_repulsion_shield(self, config.BLOCK_REPULSION_RADIUS, only_loop_atoms=True)


def add_hetatm_repulsion_shield(model: Any, min_dist: float, only_loop_atoms: bool = False) -> None:
    '''
    Helper to add repulsion restraints between CA atoms and HETATM centers.
    Shared logic between AutoModel and LoopModel.
    '''
    rsr = model.restraints

    het_residues = [r for r in model.residues if r.hetatm and r.name != 'HOH']
    if not het_residues:
        return
    
    if only_loop_atoms:
        target_sel = model.select_loop_atoms()
    else:
        target_sel = model.select_atoms()
    
    target_ca = target_sel.only_atom_types('CA')
    if len(target_ca) == 0:
        return

    print(f"\n[ENVIRONMENT][add_hetatm_repulsion_shield] Adding repulsion: {len(target_ca)} CA atoms vs {len(het_residues)} HET groups.")

    het_centers = []
    for res in het_residues:
        center = pseudo_atom.GravityCenter(Selection(res))
        rsr.pseudo_atoms.append(center)
        het_centers.append(center)
    
    count = 0
    for ca in target_ca:
        for center in het_centers:
            rsr.add(forms.LowerBound(
                group=physical.xy_distance,
                feature=features.Distance(ca, center),
                mean=min_dist,
                stdev=1.0
            ))
            count += 1
        
    print(f"\n[ENVIRONMENT][add_hetatm_repulsion_shield] Added {count} repulsion restraints (Min Dist: {min_dist}A).")
    

# ============================================================================
#                        HOMOLOGY MODELING EXECUTION
# ============================================================================

def run_automodel(env: Environ, align_file: str, job: Job,
                  residues_to_freeze: Set[int],
                  start_model: int, end_model: int, knowns: List[str], input_mode: str) -> List[Dict[str, Any]]:
    '''
    Run AutoModel for homology modeling. Returns a list of models.
    '''
    if start_model > end_model:
        logger.warning(f"[ENVIRONMENT][run_automodel] Start model {start_model} is greater than end model {end_model}. Skipping.")
        return []
    num_models = (end_model - start_model) + 1

    logger.info("\n" + "="*80)
    logger.info(f"[ENVIRONMENT][run_automodel] Running AutoModel (Job range: {start_model}-{end_model} | Count: {num_models} models)")
    logger.info(f"[ENVIRONMENT][run_automodel] Templates: {knowns}")
    logger.info(f"[ENVIRONMENT][run_automodel] Freezing residues: {residues_to_freeze}")
    logger.info("\n" + "="*80)

    if input_mode == 'precalculation':
        logger.info("[ENVIRONMENT][run_automodel] PRECALCULATION MODE: Ignoring precomputed inputs to generate fresh ones.")
        extra_inputs = {}

    elif input_mode == 'precomputed':
        logger.info(f"[ENVIRONMENT][run_automodel] Using precomputed files: {config.CUSTOM_INIFILE_PATH}, {config.CUSTOM_RSRFILE_PATH}")
        extra_inputs = {
            'inifile': config.CUSTOM_INIFILE_PATH,
            'csrfile': config.CUSTOM_RSRFILE_PATH
        }
    elif input_mode == 'normal':
        logger.info("[ENVIRONMENT][run_automodel] Normal mode: Generating inputs as usual.")
        extra_inputs = {}
        
    a = FixedRegionAutoModel(env,
                            experimental_residues=residues_to_freeze,
                            chain_id=config.CHAIN_ID,
                            alnfile=align_file,
                            knowns=knowns,
                            sequence=config.ALIGN_CODE_SEQUENCE,
                            assess_methods=(assess.DOPEHR, assess.GA341),
                            **extra_inputs
                            )
    a.use_parallel_job(job)
    a.starting_model = start_model
    a.ending_model = end_model

    a.library_schedule = autosched.slow
    a.max_var_iterations = 1000
    
    if input_mode == "precalculation": 
        logger.info("[ENVIRONMENT][run_automodel] running in PRECALCULATION mode (exit_stage=1).")
        a.make(exit_stage=1)

        generated_ini = f"{config.ALIGN_CODE_SEQUENCE}.ini"
        generated_rsr = f"{config.ALIGN_CODE_SEQUENCE}.rsr"

        if os.path.exists(generated_ini) and os.path.exists(generated_rsr):
            shutil.move(generated_ini, config.CUSTOM_INIFILE_PATH)
            shutil.move(generated_rsr, config.CUSTOM_RSRFILE_PATH)
            logger.info("[ENVIRONMENT][run_automodel] Precalculation complete. Restraints (.rsr) and Initial (.ini) files generated.")
            return []
        else:
            logger.error("[ENVIRONMENT][run_automodel] Precalculation failed. Restraints (.rsr) and Initial (.ini) files not generated.")
            return []
    else:
        a.make()
    
    if not a.outputs:
        logger.error("[ENVIRONMENT][run_automodel] No models were generated.")
        return []

    logger.info(f"[ENVIRONMENT][run_automodel] Generated {len(a.outputs)} models.")
    return a.outputs


# ============================================================================
#                        LOOP REFINEMENT EXECUTION
# ============================================================================


def run_loop_model(env: Environ, job: Job, initial_models_names: List[str],
                        loop_ranges : List[Tuple[int, int]],
                        experimental_residues: Set[int]) -> None:
    min_loop_length = 1
    max_loop_length = 1000
    '''
    Run LoopModel for loop refinement. Returns a list of models.
    '''
    if not loop_ranges:
        logger.info("[ENVIRONMENT][run_loop_model] No loop ranges provided. Skipping loop refinement.")
        return

    logger.info(f"[ENVIRONMENT][run_loop_model] Validating loop candidates. Min loop length: {min_loop_length} | Max loop length: {max_loop_length}")

    valid_loops = []
    invalid_loops = []
    for start, end in loop_ranges:
        length = end - start + 1
        if min_loop_length <= length <= max_loop_length:
            valid_loops.append((start, end))
        else:
            invalid_loops.append((start, end))
    
    if invalid_loops:
        invalid_str = ', '.join([f"[{s}-{e}]" for s, e in invalid_loops])
        logger.warning(f"[ENVIRONMENT][run_loop_model] Skipping invalid loops ({len(invalid_loops)}): {invalid_str}")
    
    if valid_loops:
        valid_str = ', '.join([f"[{s}-{e}]" for s, e in valid_loops])
        logger.info(f"[ENVIRONMENT][run_loop_model] Proceeding with valid loop ranges ({len(valid_loops)}): {valid_str}")

    if not valid_loops:
        logger.info("[ENVIRONMENT][run_loop_model] No valid loop ranges remain. Aborting loop refinement stage.")
        return
    
    logger.info(f"[ENVIRONMENT][run_loop_model] Loop refinement on {len(initial_models_names)} models")
    logger.info(f"[ENVIRONMENT][run_loop_model] Loop refinement for {len(valid_loops)} loops per model.")

    for model_idx, pdb_file in enumerate(initial_models_names):
        base_name = pdb_file.replace('.pdb', '')
        current_best_pdb = pdb_file

        logger.info(f"[ENVIRONMENT][run_loop_model] Loop refinement for model {model_idx+1}: {pdb_file}")

        for j, (start, end) in enumerate(valid_loops):
            logger.info(f"[ENVIRONMENT][run_loop_model] Loop refinement for loop {j+1}: Residues {start}-{end}")

            ml = FixedRegionLoopModel(
                env,
                inimodel=current_best_pdb,
                sequence=config.ALIGN_CODE_SEQUENCE,
                loop_start=start,
                loop_end=end,
                chain_id=config.CHAIN_ID,
                experimental_residues=experimental_residues
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
                results.sort(key=lambda x: x.get('DOPE-HR score', 9e9))
                for m, info in enumerate(results):
                    old = info['name']
                    new = f'{base_name}_LOOP{j+1}_R{m+1}.pdb'
                    if os.path.exists(old):
                        os.rename(old,new)
                    else:
                        logger.error(f"[ENVIRONMENT][run_loop_model] Loop model {old} failed to generate. Skipping.")
                        
                current_best_pdb = f'{base_name}_LOOP{j+1}_R1.pdb'
            else:
                logger.error(f"[ENVIRONMENT][run_loop_model] No models were generated for loop {j+1}. Keeping previous model.")
            
            logger.info(f"[ENVIRONMENT][run_loop_model] Loop refinement for model {model_idx+1} completed.")












