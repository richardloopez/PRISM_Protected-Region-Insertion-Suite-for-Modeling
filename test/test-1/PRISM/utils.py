#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM (utils)

'''
Utility functions for PRISM pipeline.

Handles environment setup, alignment processing, secondary structure parsing, loop detection,
and model ranking
'''

import re
import csv
import logging
from pathlib import Path
from typing import List, Tuple, Set, Dict, Any, Union

from modeller import *
from modeller.automodel import *
from modeller.parallel import Job, LocalWorker
from modeller.scripts import complete_pdb
from modeller.selection import Selection

from . import config

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("PRISM.utils")

def setup_environment() -> Tuple[Environ, Job]:
    '''
    Setup the environment for MODELLER.

    Returns:
        Tuple of Environ and Job objects
    '''
    logger.info(f"Current Working Directory: {Path.cwd()}")
    logger.info("\n" + '='*80)  
    logger.info("PRISM: PROTECTED-REGION INSERTION SUITE FOR MODELING")
    logger.info("="*80 + "\n")
    
    env = Environ()

    logger.info("Loading MODELLER topology and parameter libraries...")
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    logger.info("[ENVIRONMENT] MODELLER topology and parameter libraries loaded successfully.")
    
    env.io.atom_files_directory = ['.', config.INPUT_DIR, '../atom_files', config.MODELING_RESULTS_DIR]
    env.io.hetatm = True
    
    log.verbose()

    env.jobs = config.MODELLER_CORES
    job = Job()
    logger.info(f"[ENVIRONMENT] Using {env.jobs} local workers.")
    for _ in range(env.jobs):
        job.append(LocalWorker())
    job.start()

    return env, job


def run_prereq_cde(env: Environ) -> None:
    '''
    Run PREREQ CDE.
    '''
    logger.info("[ENVIRONMENT] Running PREREQ CDE...")

    if config.USE_MANUAL_ALIGNMENT:
        logger.info("[ENVIRONMENT] Using manual alignment.")
        clean_ali = config.MANUAL_ALIGNMENT_FILE
        if not Path(clean_ali).exists():
            raise FileNotFoundError(f"Clean alignment file not found: {clean_ali}")
        
        logger.info(f"  > Input: {clean_ali}")
        logger.info(f"  > Output: {config.MANUAL_ALIGNMENT_CDE_FILE}")

        add_cde_line_to_pir(
            clean_ali_path = clean_ali,
            cde_ali_path = config.MANUAL_ALIGNMENT_CDE_FILE,
            ss2_path = config.SS2_FILE_PATH,
            target_seq = config.get_sequence(),
            align_code = config.ALIGN_CODE_SEQUENCE
            )

    else:
        logger.info("[ENVIRONMENT] Using automatic alignment.")
        generate_pir_files_AUTO(
        env = env,
        ali_file_clean = config.ALIGNMENT_FILE, 
        ali_file_cde = config.ALIGNMENT_CDE_FILE
        )
        
    logger.info("[ENVIRONMENT] Alignment files generated successfully.")


def run_prerequisites(env: Environ) -> Tuple[str, List[Tuple[int, int]], Set[int], Set[int]]:
    '''
    Run prerequisites for automodeling.
    '''
    logger.info("\n" + "="*80)  
    logger.info("STEP 1: Alignment Preparation (Read-Only)")
    logger.info("="*80 + "\n")

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
        logger.info("[MANUAL_MODE] Ignoring automatic flank detection and SS2 coil detection.")
            
        manual_list = set(config.MANUAL_OPTIMIZATION_RESIDUES)
        truly_fixed_residues = experimental_residues - manual_list
        loop_ranges = group_ranges(list(manual_list))
        automodel_fixed = truly_fixed_residues
            
        logger.info(f"[MANUAL_MODE] User selected {len(manual_list)} residues to optimize.")
        logger.info(f"[MANUAL_MODE] Resulting Fixed Residues: {len(truly_fixed_residues)}")

    # AUTOMATIC DETECTION MODE
    else:
        logger.info("[AUTO_MODE] Running automatic detection with connectivity checks.")
            
        raw_seq_str = config.get_sequence()
        full_seq_str = re.sub(r'[^A-Z]', '', raw_seq_str.upper())
        MAX_LEN = len(full_seq_str)

        if experimental_residues:
            last_experimental_res = max(experimental_residues)
        else:
            last_experimental_res = MAX_LEN

        all_coil = get_coil_residues(config.SS2_FILE_PATH, full_seq_str)
        experimental_ranges = group_ranges(list(experimental_residues))
        flank_residues = set()
        N = config.MOBILE_FLANK_RESIDUES

        for start, end in experimental_ranges:
            if start > 1:
                for i in range(N):
                    if start + i <= end: flank_residues.add(start + i)
                
            if end < MAX_LEN:
                for i in range(N):
                    if end - i >= start: flank_residues.add(end - i)

        truly_fixed_residues = experimental_residues - flank_residues
        refinable_residues = all_coil - truly_fixed_residues
        loop_ranges = group_ranges(list(refinable_residues))

        if config.REFINE_FLANKS_DURING_AUTOMODEL:
            automodel_fixed = truly_fixed_residues
            logger.info(f"  > Fixed residues (Auto): {len(automodel_fixed)}")
        else:
            automodel_fixed = experimental_residues
            logger.info(f"  > Fixed residues (Exp Only): {len(automodel_fixed)}")

        # OUTPUT SUMMARY
    ranges_str = [f"[{s}-{e}]" for s, e in loop_ranges]
    if not ranges_str:
        logger.info("\n[ENVIRONMENT] No regions selected for optimization.")
    else:
        logger.info(f"\n[ENVIRONMENT] Regions selected for optimization: {', '.join(ranges_str)}")

    return ali_file_automodel, loop_ranges, truly_fixed_residues, automodel_fixed


# ============================================================================
#                               HELPER FUNCTIONS
# ============================================================================


def flatten_ali_file(filepath: Union[str, Path]) -> None:
    '''
    Flatten an alignment file by removing all lines between sequences.

    Args:
        filepath: Path to the alignment file to flatten
    '''
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    lines = path.read_text().splitlines()
    new_lines, seq_buffer = [], []

    for line in lines:
        stripped = line.strip()
        if stripped.startswith(('>P1;', 'structure', 'sequence', '#')):
            if seq_buffer:
                full_seq = ''.join(seq_buffer).replace(' ', '')
                if not full_seq.endswith('*'):
                    full_seq += '*'
                new_lines.append(full_seq)
                seq_buffer = []
            new_lines.append(stripped)
        elif stripped:
            seq_buffer.append(stripped)

    if seq_buffer:
        full_seq = ''.join(seq_buffer).replace(' ', '')
        if not full_seq.endswith('*'):
            full_seq += '*'
        new_lines.append(full_seq)

    path.write_text('\n'.join(new_lines) + '\n') 
    logger.info(f"[ENVIRONMENT] Ali file flattened: {filepath}")

        
def read_ss2_file(ss2_path: Union[str, Path], seq_full: str) -> str:
    '''
    Read an ss2 file and return the secondary structure string.

    Args:
        ss2_path: Path to the ss2 file
        seq_full: Full sequence string

    Returns:
        Secondary structure string
    '''
    logger.info(f"[ENVIRONMENT] Reading ss2 file: {ss2_path}")
    lines = Path(ss2_path).read_text().splitlines()
    ss_parts = []
    for line in lines:
        if line.strip() and not line.startswith('#') and len(line.split()) >= 3:
            ss_parts.append(line.split()[2])
    ss_string = ''.join(ss_parts)
    
    len_ss, len_seq = len(ss_string), len(seq_full)

    if len_ss > len_seq:
        raise ValueError(f"SS2 and sequence lengths do not match. SS2: {len_ss}, Sequence: {len_seq}")

    if len_ss < len_seq:
        missing = len_seq - len_ss
        if seq_full[len_ss:].strip('./'):
            raise ValueError(f"SS2 too short (SS2: {len_ss}, Sequence   : {len_seq}) and missing part contains residues.")

        logger.warning(f"[ENVIRONMENT] Padding short SS2 with {missing} '.' characters.")
        ss_string += '.' * missing
        
    return ss_string[:len_seq]
    

def read_sequences_from_ali(ali_file: Union[str, Path]) -> Dict[str, str]:
    '''
    Read an alignment file and return a dictionary of sequences.

    Args:
        ali_file: Path to the alignment file

    Returns:
        Dictionary of sequences
    '''
    sequences = {}
    current_code = None
    current_seq_parts = []

    lines = Path(ali_file).read_text().splitlines()
    for line in lines:
        line = line.strip()
        if line.startswith('>P1;'):
            if current_code:
                sequences[current_code] = ''.join(current_seq_parts)
            current_code = line.split(';')[1].strip()
            current_seq_parts = []
        elif current_code and not line.startswith(('structure', 'sequence', 'CDE', '#')):
            current_seq_parts.append(line)

    if current_code:
        sequences[current_code] = ''.join(current_seq_parts)

    clean_map = {k: re.sub(r'[^A-Z\-\.\/]', '', v.upper()) for k, v in sequences.items()}

    lengths = {len(s) for s in clean_map.values()}
    if len(lengths) > 1:
        raise ValueError(f"All sequences must have the same length. Found lengths: {lengths}")
        
    return clean_map


def add_cde_line_to_pir(clean_ali_path: Union[str, Path], cde_ali_path: Union[str, Path], 
                        ss2_path: Union[str, Path], target_seq: str, align_code: str) -> None:
    '''
    Add a CDE line to a PIR file.

    Args:
        clean_ali_path: Path to the clean alignment file
        cde_ali_path: Path to the CDE alignment file
        ss2_path: Path to the SS2 file
        target_seq: Target sequence string
        align_code: Alignment code
    '''
    logger.info(f"[ENVIRONMENT] Adding CDE line to PIR file: {clean_ali_path}")

    ss2_string = read_ss2_file(ss2_path, target_seq)
    aligned_seqs = read_sequences_from_ali(clean_ali_path)

    if align_code not in aligned_seqs:
        raise ValueError(f"Alignment code {align_code} not found in alignment file {clean_ali_path}")
    
    align_target = aligned_seqs[align_code]
    cde_chars = []
    ss_idx = 0
    ss_blk = 0

    for char in align_target:
        if char in ('-', '/', '.'):
            cde_chars.append('.')
            ss_blk += 1
        else:
            if ss_idx < len(ss2_string):
                cde_chars.append(ss2_string[ss_idx])
                ss_idx += 1
            else:
                raise IndexError(f"SS2 string is shorter than aligned target sequence")
    
    cde_line = "# CDE:" + "".join(cde_chars)

    input_lines = Path(clean_ali_path).read_text().splitlines()
    output_lines = []
    sequences_counter = 0
    add_cde_sequences_counter = 0
    found_target = False

    for line in input_lines:
        output_lines.append(line)
        if line.strip().startswith('>P1;'):
            sequences_counter += 1
        if line.strip().startswith(f'>P1;{align_code}'):
            found_target = True
        elif found_target and line.strip().startswith('sequence:'):
            output_lines.append(cde_line)
            add_cde_sequences_counter += 1
            found_target = False
    
    Path(cde_ali_path).write_text('\n'.join(output_lines) + '\n')
    flatten_ali_file(cde_ali_path)
    logger.info(f"[ENVIRONMENT] CDE line added to PIR file: {cde_ali_path}")
    logger.info(f"[ENVIRONMENT] {ss_blk} '.' added to CDE line (BLK residues in aligned target sequence)")
    logger.info(f"[ENVIRONMENT] {sequences_counter} sequences in PIR file: {clean_ali_path}")
    logger.info(f"[ENVIRONMENT] {add_cde_sequences_counter} CDE sequences added to PIR file: {cde_ali_path} | The sequence code is {align_code}")
    

def generate_pir_files_AUTO(env: Environ, ali_file_clean: str, ali_file_cde: str) -> None:
    '''
    Generate PIR files for the target and template sequences.

    Args:
        env: Modeller environment
        ali_file_clean: Path to the clean alignment file
        ali_file_cde: Path to the CDE alignment file
    '''
    logger.info(f"[ENVIRONMENT] Generating PIR files for {ali_file_clean} and {ali_file_cde}")
    aln = Alignment(env)

    for pdb_path, code in zip(config.PDB_TEMPLATE_FILES_PATHS, config.PDB_TEMPLATE_FILES_NAMES):
        logger.info(f"[ENVIRONMENT] Adding template: {Path(pdb_path).name} ({code})")
        aln.append_model(Model(env, file=pdb_path), align_codes=code, atom_files=Path(pdb_path).name)

    aln.append_sequence(config.get_sequence())
    aln[len(config.PDB_TEMPLATE_FILES_PATHS)].code = config.ALIGN_CODE_SEQUENCE

    aln.salign()
    aln.write(file=ali_file_clean, alignment_format='PIR')
    flatten_ali_file(ali_file_clean)

    add_cde_line_to_pir(ali_file_clean, ali_file_cde, config.SS2_FILE_PATH, config.get_sequence(), config.ALIGN_CODE_SEQUENCE)
    

def prepare_prism_power_files(phase: str) -> Tuple[List[str], str]:
    '''
    Optimized Power Alignment:
    1. Points all replicas to the original PDB file (no disk copying).
    2. Ensures PIR structure headers are correctly formatted.
    3. Excludes the target sequence from the 'knowns' list to avoid MODELLER errors.
    '''
    power_map = getattr(config.PRISM_POWER_SETTINGS, phase)
    input_dir = Path(config.INPUT_DIR)
    original_ali_path = input_dir / config.MANUAL_ALIGNMENT_BASENAME
    power_ali_path = input_dir / f"prism_power_{phase}.ali"
    
    all_expanded_knowns = []

    with open(original_ali_path, 'r') as f:
        content = f.read()
    
    blocks = content.split('>P1;')
    header = blocks[0] 
    ali_map = {}
    for b in blocks[1:]:
        if not b.strip(): continue
        lines = b.splitlines()
        code = lines[0].strip()
        ali_map[code] = lines[1:] 

    new_ali_content = header
    
    for base_pdb, count in power_map.items():
        base_code = None
        if base_pdb in ali_map:
            base_code = base_pdb
        elif base_pdb.replace('.pdb', '') in ali_map:
            base_code = base_pdb.replace('.pdb', '')
        
        if not base_code:
            logger.error(f"[ERROR] Template {base_pdb} not found in alignment IDs: {list(ali_map.keys())}")
            continue

        original_lines = ali_map[base_code]
        
        for i in range(count):
            replica_id = f"{base_code}_{i:03d}" if i > 0 else base_code
            all_expanded_knowns.append(replica_id)
            
            block_lines = list(original_lines)
            if block_lines and block_lines[0].startswith('structure'):
                fields = block_lines[0].split(':')
                fields[1] = base_pdb 
                block_lines[0] = ':'.join(fields)
            
            new_ali_content += f">P1;{replica_id}\n" + '\n'.join(block_lines) + "\n"

    target_code = config.ALIGN_CODE_SEQUENCE
    if target_code in ali_map:
        new_ali_content += f">P1;{target_code}\n" + '\n'.join(ali_map[target_code]) + "\n"
    else:
        logger.error(f"[ERROR] Target code '{target_code}' not found in original alignment.")

    with open(power_ali_path, 'w') as f:
        f.write(new_ali_content)
    
    logger.info(f"[ENVIRONMENT] Virtual power alignment generated: {power_ali_path}")
    logger.info(f"{len(all_expanded_knowns)} virtual templates defined (Target excluded).")
    logger.info(f"[ENVIRONMENT] Virtual templates defined: {all_expanded_knowns}")

    return all_expanded_knowns, str(power_ali_path) 


def identify_experimental_residues(aligned_template_seq: str, aligned_target_seq: str) -> Set[int]:
    '''
    Identify target residues that map to experimental regions in the main template.

    These residues correspond to positions where both the template and target
    have actual residues (not gaps), indicating experimental coordinates that
    must be preserved during modeling.

    Args:
        aligned_template_seq: Aligned sequence from the main template (MAIN_PDB)
        aligned_target_seq: Aligned sequence from the target protein

    Returns:
        Set of residue numbers (1-indexed) in the target sequence that map
        to experimental positions in the template (not gaps in either sequence)

    Raises:
        ValueError: If aligned sequences have different lengths
    '''
    experimental_residues = set()
    target_res_num = 0

    if len(aligned_template_seq) != len(aligned_target_seq):
        logger.error("[ERROR] Aligned sequence lengths do not match in identify_experimental_residues.")
        raise ValueError("Aligned template and target sequences must have equal lengths.")
    
    for template_res, target_res in zip(aligned_template_seq, aligned_target_seq):
        if template_res == '.' or template_res == '/':
            continue
        if target_res != '-':
            target_res_num += 1
        if template_res != '-' and target_res != '-':
            experimental_residues.add(target_res_num)
    
    logger.info(f"[FIXED_REGION] Identified {len(experimental_residues)} experimental residues mapped from template")
    logger.info("[FIXED_REGION] These residues will NOT be optimized or refined (unless flank size takes them)")

    return experimental_residues


def group_ranges(residues: List[int]) -> List[Tuple[int, int]]:
    '''
    Group a list of residues into ranges.
    
    Args:
        residues: List of residues
    
    Returns:
        List of ranges
    '''
    if not residues:
        return []
    
    residues = sorted(set(residues))
    ranges = []
    start = end = residues[0]

    for res in residues[1:]:
        if res == end + 1:
            end = res
        else:
            ranges.append((start, end))
            start = end = res
    ranges.append((start, end))
    return ranges


def get_coil_residues(ss2_path: str, seq_full: str) -> Set[int]:
    '''
    Get the coil residues from the SS2 file.
    
    Args:
        ss2_path: Path to the SS2 file
        seq_full: Full sequence string
    
    Returns:
        Set of coil residues
    '''
    ss2_string = read_ss2_file(ss2_path, seq_full)

    if len(seq_full) != len(ss2_string):
        raise ValueError("Sequence and SS2 string must have the same length")
    
    coil_residues = []
    for i, char in enumerate(ss2_string):
        if char == 'C':
            coil_residues.append(i+1)
    return set(coil_residues)


# ============================================================================
#                            RANKING & EVALUATION
# ============================================================================

def run_rank_automodel_models(env: Environ) -> None:
    '''
    Run automodel ranking.
    
    Args:
        env: Modeller environment
    '''
    logger.info("[AUTOMODEL_RANKING] Starting automodel ranking")

    cwd = Path.cwd()
    pattern = re.compile(rf"^{re.escape(config.ALIGN_CODE_SEQUENCE)}\.B[0-9]{{5,}}\.pdb$")

    raw_models = [f for f in cwd.iterdir() if f.is_file() and pattern.match(f.name)]

    if not raw_models:
        logger.error(f"[ERROR] No models found in {cwd} (Pattern: {config.ALIGN_CODE_SEQUENCE}.B*.pdb).")
        return
    logger.info(f"[AUTOMODEL_RANKING] Found {len(raw_models)} models. Beginning ranking...")

    results = []
    for model_path in raw_models:
        mdl = complete_pdb(env, str(model_path))
        score = Selection(mdl.chains[config.CHAIN_ID]).assess_dopehr()
        results.append({'path': model_path, 'score': score})
        logger.info(f"[AUTOMODEL_RANKING] Evaluated -> {model_path}: | DOPEHR score: {score}")
    
    results.sort(key=lambda x: x['score'])
    logger.info(f"[AUTOMODEL_RANKING] Ranked {len(results)} models... renaming")

    selected_for_refinement = []

    for rank, data in enumerate(results, 1):
        old_path = data['path']
        new_name = f"AUTO_{rank}.pdb"
        new_path = cwd / new_name

        old_path.rename(new_path)
        if rank <= config.TOP_MODELS_FOR_REFINEMENT:
            selected_for_refinement.append(new_name)
    
    logger.info("Ranking complete")
    logger.info(f"Selected {len(selected_for_refinement)} models for refinement")
    for name in selected_for_refinement:
        logger.info(f"Selected -> {name}")


def final_evaluation_and_ranking(env) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    '''
    Final evaluation and ranking of the models.
    
    Args:
        env: environment
    
    Returns:
        Tuple of List of results and best model
    '''

    logger.info("[FINAL_EVALUATION] Starting final evaluation and ranking")

    template_names = set(config.PDB_TEMPLATE_FILES_NAMES)
    pdbs = [
        p for p in Path('.').glob('*.pdb')
        if p.name not in template_names and ('AUTO_' in p.name or 'LOOP_' in p.name)
    ]
    if not pdbs:
        logger.warning('No valid PDB files found for final evaluation.')
        return [], {}
    
    results = []
    for pdb_path in pdbs:
        mdl = complete_pdb(env, str(pdb_path))
        score = Selection(mdl.chains[config.CHAIN_ID]).assess_dopehr()
        zscore = mdl.assess_normalized_dopehr()
        results.append({
            'name': pdb_path,
            'DOPEHR_score': score,
            'DOPEHR_zscore': zscore
        })
        logger.info(f"[FINAL_EVALUATION] {pdb_path}: DOPEHR score: {score}, DOPEHR zscore: {zscore}")
    results.sort(key=lambda x: x['DOPEHR_score'])
    val = config.NUM_BEST_FINAL_MODELS
    limit = None if (val == 'inf' or val == float('inf') or val is None) else int(val)
    best_models = results[:limit]

    with open (config.FINAL_RANKING_CSV, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['Rank', 'Model_Name', 'DOPEHR_score', 'DOPEHR_zscore'])
        writer.writeheader()
        for rank, data in enumerate(best_models, start=1):
            writer.writerow({
                'Rank': rank,
                'Model_Name': data['name'],
                'DOPEHR_score': data['DOPEHR_score'],
                'DOPEHR_zscore': data['DOPEHR_zscore']
            })
    logger.info(f"[FINAL_EVALUATION] Best models ranking saved to {config.FINAL_RANKING_CSV}")
    
    return results, (best_models[0] if best_models else {})
    
    return results, (best_models[0] if best_models else {})


    


