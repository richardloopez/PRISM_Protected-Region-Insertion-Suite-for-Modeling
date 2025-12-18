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

import sys
import re
import csv
from pathlib import Path
from typing import List, Tuple, Set, Dict, Any, Union

from modeller import *
from modeller.automodel import *
from modeller.parallel import Job, LocalWorker
from modeller.scripts import complete_pdb
from modeller.selection import Selection

from . import config

# ============================================================================
#                               SETUP & STAGE FUNCTIONS
# ============================================================================


def setup_environment() -> Tuple[Environ, Job]:
    '''
    Setup the environment for MODELLER.

    Returns:
        Tuple of Environ and Job objects
    '''
    print(f"\n[ENVIRONMENT] Current Working Directory: {Path.cwd()}")
    print("\n" + '='*80)
    print("PRISM: PROTECTED-REGION INSERTION SUITE FOR MODELING")
    print("="*80 + "\n")
    
    env = Environ()

    print(f"\n[ENVIRONMENT] Loading MODELLER topology and parameter libraries...")
    try: 
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        print(f"\n[ENVIRONMENT] MODELLER topology and parameter libraries loaded successfully.")
    except Exception as e:
        print(f"\n[ENVIRONMENT] Failed to load MODELLER topology and parameter libraries: {e}")
        sys.exit(1)
    
    env.io.atom_files_directory = ['.', config.INPUT_DIR, '../atom_files', config.MODELING_RESULTS_DIR]
    env.io.hetatm = True
    
    log.verbose()

    env.jobs = config.NUM_PROCESSORS
    job = Job()
    print(f"\n[ENVIRONMENT] Using {env.jobs} local workers.")
    for _ in range(env.jobs):
        job.append(LocalWorker())
    job.start()

    return env, job


def run_prereq_cde(env: Environ):
    '''
    Run PREREQ CDE.
    '''
    print(f"\n[ENVIRONMENT] Running PREREQ CDE...")

    try: 
        if config.USE_MANUAL_ALIGNMENT:
            print(f"\n[ENVIRONMENT] Using manual alignment.")
            clean_ali = config.MANUAL_ALIGNMENT_FILE
            cde_ali_name = Path(config.MANUAL_ALIGNMENT_CDE_FILE).name
            cde_ali_path = Path(config.MODELING_RESULTS_DIR) / cde_ali_name
            if not Path(clean_ali).exists():
                raise FileNotFoundError(f"Clean alignment file not found: {clean_ali}")
        
            print(f"  > Input: {clean_ali}")
            print(f"  > Output: {cde_ali_path}")

            add_cde_line_to_pir(
                clean_ali_path = clean_ali,
                cde_ali_path = str(cde_ali_path),
                ss2_path = config.SS2_FILE_PATH,
                target_seq = config.get_sequence(),
                align_code = config.ALIGN_CODE_SEQUENCE
                )

        else:
            print(f"\n[ENVIRONMENT] Using automatic alignment.")
            generate_pir_files_AUTO(
            env = env,
            ali_file_clean = config.ALIGNMENT_FILE, 
            ali_file_cde = config.ALIGNMENT_CDE_FILE
            )
        
        print(f"\n[ENVIRONMENT] Alignment files generated successfully.")

    except Exception as e:
        print(f"\n[ENVIRONMENT] Failed to prepare CDE alignment files: {e}")
        sys.exit(1)


def run_prerequisites(env: Environ) -> Tuple[str, List[Tuple[int, int]], Set[int], Set[int]]:
    '''
    Run prerequisites for automodeling.
    '''
    print(f"\n{'='*80}\n[STEP 1] Alignment Preparation (Read-Only)\n{'='*80}\n")

    if config.USE_MANUAL_ALIGNMENT:
        ali_file_automodel = config.MANUAL_ALIGNMENT_FILE
        cde_name = Path(config.MANUAL_ALIGNMENT_CDE_FILE).name
        ali_file_analysis = str(Path(config.MODELING_RESULTS_DIR) / cde_name)
    else:
        ali_file_automodel = config.ALIGNMENT_FILE
        ali_file_analysis = config.ALIGNMENT_CDE_FILE
    
    if not Path(ali_file_automodel).exists():
        raise FileNotFoundError(f"Clean alignment file not found: {ali_file_automodel}")
    if not Path(ali_file_analysis).exists():
        raise FileNotFoundError(f"CDE alignment file not found: {ali_file_analysis}")
    
    try:
        seq_map = read_sequences_from_ali(ali_file_analysis)
        main_tmpl_seq = seq_map[config.MAIN_ALIGN_CODE_TEMPLATE]
        target_seq = seq_map[config.ALIGN_CODE_SEQUENCE]

        print(f"\n[ENVIRONMENT] Experimental residue identification\n{'='*80}\n")
        experimental_residues = identify_experimental_residues(main_tmpl_seq, target_seq)

    except Exception as e:
        print(f"\n[ENVIRONMENT] Failed to process alignment files: {e}")
        sys.exit(1)
    
    print(f"\n[ENVIRONMENT] Loop detection\n{'='*80}\n")
    try:
        all_coil = get_coil_residues(config.SS2_FILE_PATH, config.get_sequence())
        experimental_ranges = group_ranges(list(experimental_residues))
        flank_residues = set()
        N = config.EXPERIMENTAL_FLANK_SIZE

        for start, end in experimental_ranges:
            for i in range(N):
                if start + i <= end: flank_residues.add(start + i)
            for i in range(N):
                if end - i >= start: flank_residues.add(end - i)
        truly_fixed_residues = experimental_residues - flank_residues
        refinable_residues = all_coil - truly_fixed_residues
        loop_ranges = group_ranges(list(refinable_residues))

        if config.REFINE_FLANKS_DURING_AUTOMODEL:
            automodel_fixed = truly_fixed_residues
            print(f"\n[ENVIRONMENT] Refining flanks during automodeling\n{'='*80}\n")
            print(f"  > Fixed residues: {automodel_fixed}")
        else:
            automodel_fixed = experimental_residues
            print(f"\n[ENVIRONMENT] Not refining flanks during automodeling\n{'='*80}\n")
            print(f"  > Fixed residues: {automodel_fixed}")
        
        ranges_str = [f"[{s}-{e}]" for s, e in loop_ranges]
        print(f"\n[ENVIRONMENT] Detected coil regions from SS2 (Candidates for refinement): {', '.join(ranges_str)}")

        return ali_file_automodel, loop_ranges, truly_fixed_residues, automodel_fixed
    
    except Exception as e:
        print(f"\n[ENVIRONMENT] Loop detection failed: {e}")
        sys.exit(1)


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
    print(f"\n[ALIGNMENT] Ali file flattened: {filepath}")   

        
def read_ss2_file(ss2_path: Union[str, Path], seq_full: str) -> str:
    '''
    Read an ss2 file and return the secondary structure string.

    Args:
        ss2_path: Path to the ss2 file
        seq_full: Full sequence string

    Returns:
        Secondary structure string
    '''
    print(f"\n[ALIGNMENT] Reading ss2 file: {ss2_path}")
    try:
        lines = Path(ss2_path).read_text().splitlines()
        ss_parts = []
        for line in lines:
            if line.strip() and not line.startswith('#') and len(line.split()) >= 3:
                ss_parts.append(line.split()[2])
        ss_string = ''.join(ss_parts)
    except Exception as e:
        raise IOError(f"Failed to read ss2 file: {ss2_path}") from e
    
    len_ss, len_seq = len(ss_string), len(seq_full)

    if len_ss > len_seq:
        raise ValueError(f"SS2 and sequence lengths do not match. SS2: {len_ss}, Sequence: {len_seq}")

    if len_ss < len_seq:
        missing = len_seq - len_ss
        if seq_full[len_ss:].strip('.'):
            raise ValueError(f"SS2 too short (SS2: {len_ss}, Sequence   : {len_seq}) and missing part contains residues.")

        print(f"[SS2] Padding short SS2 with {missing} '.' characters.")
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

    try: 
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

    except Exception as e:
        raise IOError(f"Failed to read alignment file: {ali_file}") from e


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
    print(f"\n[ALIGNMENT] Adding CDE line to PIR file: {clean_ali_path}")

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
    print(f"[ALIGNMENT] CDE line added to PIR file: {cde_ali_path}")
    print(f"[ALIGNMENT] {ss_blk} '.' added to CDE line (BLK residues in aligned target sequence)")
    print(f"[ALIGNMENT] {sequences_counter} sequences in PIR file: {clean_ali_path}")
    print(f"[ALIGNMENT] {add_cde_sequences_counter} CDE sequences added to PIR file: {cde_ali_path} | The sequence code is {align_code}")
    

def generate_pir_files_AUTO(env, ali_file_clean: str, ali_file_cde: str):
    '''
    Generate PIR files for the target and template sequences.

    Args:
        env: environment
        ali_file_clean: Path to the clean alignment file
        ali_file_cde: Path to the CDE alignment file
    '''
    print(f"\n[ALIGNMENT] Generating PIR files for {ali_file_clean} and {ali_file_cde}")
    aln = Alignment(env)

    for pdb_path, code in zip(config.PDB_TEMPLATE_FILES_PATHS, config.PDB_TEMPLATE_FILES_NAMES):
        print(f"\n[ALIGNMENT] Adding template: {Path(pdb_path).name} ({code})")
        aln.append_model(Model(env, file=pdb_path), align_codes=code, atom_files=Path(pdb_path).name)

    aln.append_sequence(config.get_sequence())
    aln[len(config.PDB_TEMPLATE_FILES_PATHS)].code = config.ALIGN_CODE_SEQUENCE

    aln.salign()
    aln.write(file=ali_file_clean, alignment_format='PIR')
    flatten_ali_file(ali_file_clean)

    try:
        add_cde_line_to_pir(ali_file_clean, ali_file_cde, config.SS2_FILE_PATH, config.get_sequence(), config.ALIGN_CODE_SEQUENCE)
    except Exception as e:
        raise RuntimeError(f"Failed to add CDE line to PIR file: {e}")
    

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
        print("[ERROR] Aligned sequence lengths do not match in identify_experimental_residues.")
        raise ValueError("Aligned template and target sequences must have equal lengths.")
    
    for template_res, target_res in zip(aligned_template_seq, aligned_target_seq):
        if template_res == '.' or template_res == '/':
            continue
        if target_res != '-':
            target_res_num += 1
        if template_res != '-' and target_res != '-':
            experimental_residues.add(target_res_num)
    
    print(f"\n[FIXED_REGION] Identified {len(experimental_residues)} experimental residues mapped from template")
    print(f"[FIXED_REGION] These residues will NOT be optimized or refined (unless flank size takes them, which will be specified later)")

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

def run_rank_automodel_models(env):
    '''
    Run automodel ranking.
    
    Args:
        env: environment
    '''
    print("\n[AUTOMODEL_RANKING] Starting automodel ranking")

    cwd = Path.cwd()
    pattern = re.compile(rf"^{re.escape(config.ALIGN_CODE_SEQUENCE)}\.B[0-9]{{5,}}\.pdb$")

    raw_models = [f for f in cwd.iterdir() if f.is_file() and pattern.match(f.name)]

    if not raw_models:
        print(f"[AUTOMODEL_RANKING] No models found in {cwd} (Pattern: {config.ALIGN_CODE_SEQUENCE}.B*.pdb).")
        sys.exit(1)
    print(f"[AUTOMODEL_RANKING] Found {len(raw_models)} models. Beginning ranking...")

    results = []
    for model_path in raw_models:
        try:
            mdl = complete_pdb(env, str(model_path))
            score = Selection(mdl.chains[config.CHAIN_ID]).assess_dopehr()
            results.append({'path': model_path, 'score': score})
            print(f"[AUTOMODEL_RANKING] Evaluated -> {model_path}: | DOPEHR score: {score}")
        except Exception as e:
            print(f"[AUTOMODEL_RANKING] Failed to evaluate {model_path}: {e}")
            results.append({'path': model_path, 'score': float('inf')})
    
    results.sort(key=lambda x: x['score'])
    print(f"[AUTOMODEL_RANKING] Ranked {len(results)} models... renaming")

    selected_for_refinement = []

    for rank, data in enumerate(results, 1):
        old_path = data['path']
        new_name = f"AUTO_{rank}.pdb"
        new_path = cwd / new_name

        try:
            old_path.rename(new_path)
            if rank <= config.NUM_MODELS_TO_REFINE:
                selected_for_refinement.append(new_name)
        except OSError as e:
            print(f"[AUTOMODEL_RANKING] Failed to rename {old_path} to {new_path}: {e}")
    
    print(f"[AUTOMODEL_RANKING] Ranking complete")
    print(f"[AUTOMODEL_RANKING] Selected {len(selected_for_refinement)} models for refinement")
    for name in selected_for_refinement:
        print(f"[AUTOMODEL_RANKING] Selected -> {name}")


def final_evaluation_and_ranking(env) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    '''
    Final evaluation and ranking of the models.
    
    Args:
        env: environment
    
    Returns:
        Tuple of List of results and best model
    '''

    print("\n[FINAL_EVALUATION] Starting final evaluation and ranking")

    template_names = set(config.PDB_TEMPLATE_FILES_NAMES)
    pdbs = [
        p for p in Path('.').glob('*.pdb')
        if p.name not in template_names and ('AUTO_' in p.name or 'LOOP_' in p.name)
    ]
    if not pdbs:
        raise ValueError('No valid PDB files found in the current directory')
    
    results = []
    for pdb_path in pdbs:
        try:
            mdl = complete_pdb(env, str(pdb_path))
            score = Selection(mdl.chains[config.CHAIN_ID]).assess_dopehr()
            zscore = mdl.assess_normalized_dopehr()
            results.append({
                'name': pdb_path,
                'DOPEHR_score': score,
                'DOPEHR_zscore': zscore
            })
            print(f"[FINAL_EVALUATION] {pdb_path}: DOPEHR score: {score}, DOPEHR zscore: {zscore}")
        except Exception as e:
            print(f"[FINAL_EVALUATION] {pdb_path}: Failed to evaluate: {e}")
            results.append({
                'name': pdb_path,
                'DOPEHR_score': float('inf'),
                'DOPEHR_zscore': float('inf')
            })
    results.sort(key=lambda x: x['DOPEHR_score'])
    best_models = results[:config.NUM_BEST_FINAL_MODELS]

    try: 
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
            print(f"[FINAL_EVALUATION] Best models ranking saved to {config.FINAL_RANKING_CSV}")
    except Exception as e:
        print(f"[FINAL_EVALUATION] Failed to save best models ranking: {e}")
    
    return results, (best_models[0] if best_models else {})


    

