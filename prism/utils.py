#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM (Protected-Region Insertion Suite for Modeling)
#
# Citation:
# If you use this software in your research, please cite:
# Lopez-Corbalan, R.

"""
Utility functions for PRISM pipeline.

This module provides essential utilities for alignment processing, secondary
structure handling, loop detection, and model evaluation.
"""

import os
import re
import csv
from typing import List, Tuple, Dict, Any, Set

from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from modeller.selection import Selection

from . import config
from .config import (
    SS2_FILE_PATH, sequence_full,
    PDB_TEMPLATE_FILES_PATHS, PDB_TEMPLATE_FILES_NAMES, ALIGN_CODES_TEMPLATES,
    MAIN_ALIGN_CODE_TEMPLATE, ALIGN_CODE_SEQUENCE, CHAIN_ID,
    FINAL_RANKING_CSV, NUM_BEST_FINAL_MODELS
)

def _flatten_ali_file(filepath: str):
    """
    Reads an ALI file (filepath) and rewrites it in place,
    with all sequence blocks flattened to a single line.

    Args:
        filepath (str): Path to the .ali file to modify.
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        new_lines = []
        seq_buffer = []
        for line in lines:
            line_stripped = line.strip()
            
            # Check for header or description lines
            if line_stripped.startswith(('>P1;', 'structureX:', 'structure:', 'sequence:', '# CDE:')):
                # If we have a sequence in the buffer, write it first
                if seq_buffer:
                    # Join, remove all newlines/spaces, and ensure it ends with *
                    full_seq = "".join(seq_buffer).replace('\n', '').replace(' ', '')
                    if not full_seq.endswith('*'):
                        full_seq += '*'
                    new_lines.append(full_seq + '\n')
                
                seq_buffer = [] # Reset buffer
                
                # We strip() the original line and add our own \n
                # This ensures consistent newlines on all header lines.
                new_lines.append(line_stripped + '\n') 
                # ---------------------------------
            
            elif not line_stripped:
                continue # Skip empty lines
            
            else:
                # This must be a sequence line
                seq_buffer.append(line_stripped)
        
        # Write the last sequence block remaining in the buffer
        if seq_buffer:
            full_seq = "".join(seq_buffer).replace('\n', '').replace(' ', '')
            if not full_seq.endswith('*'):
                full_seq += '*'
            new_lines.append(full_seq + '\n')

        # Rewrite the file with flattened sequences
        with open(filepath, 'w') as f:
            f.writelines(new_lines)
        
        print(f"  > Flattened sequence lines in: {os.path.basename(filepath)}")

    except Exception as e:
        print(f"  > [Warning] Failed to flatten {filepath}. Skipping. Error: {e}")
        
def _read_ss2_file(ss2_file_path: str, target_seq_full: str) -> str:
    """
    Read and validate a PSIPRED SS2 file against the full target sequence.

    Validation logic:
    1. If len(ss) == len(seq): Perfect match
    2. If len(ss) > len(seq): ERROR
    3. If len(ss) < len(seq):
       a. If extra sequence part contains ONLY '.' (BLK), accept and pad with '.'
       b. If extra sequence part contains ANY amino acid, ERROR

    Args:
        ss2_file_path: Path to PSIPRED SS2 secondary structure file
        target_seq_full: Complete target protein sequence

    Returns:
        Secondary structure string validated and padded if necessary

    Raises:
        FileNotFoundError: If SS2 file doesn't exist
        ValueError: If SS2 length validation fails
    """
    print(f"\n[SS2] Reading secondary structure file: {ss2_file_path}")
    ss_string = ""
    try:
        with open(ss2_file_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 3:
                    ss_string += parts[2]
    except FileNotFoundError:
        print(f"  [ERROR] PSIPRED SS2 file '{ss2_file_path}' not found.")
        raise
    except Exception as e:
        print(f"  [ERROR] Could not read SS2 file. Error: {e}")
        raise

    len_ss = len(ss_string)
    len_seq = len(target_seq_full)

    if len_ss == len_seq:
        print(f"[SS2] Secondary structure length ({len_ss}) matches sequence ({len_seq}).")
        return ss_string

    elif len_ss > len_seq:
        print(f"  [ERROR] SS2 file is LONGER ({len_ss}) than sequence ({len_seq}).")
        raise ValueError(f"Secondary structure (len {len_ss}) is longer than sequence (len {len_seq}).")

    else:
        missing_len = len_seq - len_ss
        extra_seq_part = target_seq_full[len_ss:]

        if extra_seq_part.strip('.') == "":
            print(f"  [INFO] SS length ({len_ss}) is shorter than sequence ({len_seq}).")
            print(f"         Missing part of sequence ({missing_len} res) contains only BLK characters ('.').")
            print("         Padding missing SS with '.' (BLK).")
            ss_string += '.' * missing_len
            return ss_string
        else:
            print(f"  [ERROR] SS length ({len_ss}) is shorter than sequence ({len_seq}).")
            print(f"          Missing part of sequence contains amino acids.")
            print(f"          Please regenerate the .ss2 file to cover the complete sequence.")
            raise ValueError(f"SS2 file (len {len_ss}) is shorter than sequence (len {len_seq}) and missing part is not all BLK.")


def extract_ss_from_ss2(ss2_file_path: str, seq_full: str) -> str:
    """
    Extract validated secondary structure string from SS2 file.

    Args:
        ss2_file_path: Path to PSIPRED SS2 file
        seq_full: Complete target sequence

    Returns:
        Secondary structure string (H, E, C, or .) for each residue
    """
    ss_string_full = _read_ss2_file(ss2_file_path, seq_full)

    if not ss_string_full:
        return ""

    ss_string_sliced = ss_string_full[:len(seq_full)]
    print(f"[CDE] Secondary structure extracted. Length used: {len(ss_string_sliced)} for CDE line.")
    return ss_string_sliced


def read_sequences_from_ali_temp(ali_file: str) -> Dict[str, str]:
    """
    Read all aligned sequences (including gaps) from a PIR/ALI file robustly.

    This function correctly handles multi-line sequences, chain breaks ('/'),
    and the final '*' terminator, ensuring all sequences in the map have
    the same length.

    Args:
        ali_file (str): Path to PIR alignment file.

    Returns:
        Dict[str, str]: Dictionary mapping align_code to its full aligned sequence string.

    Raises:
        IOError: If file cannot be read or parsed.
        ValueError: If no sequences are found or lengths mismatch.
    """
    sequences_map: Dict[str, str] = {}
    current_code = ""
    current_raw_sequence = ""
    
    try:
        with open(ali_file, 'r') as f:
            lines = f.readlines()

        # Pass 1: Read all sequence blocks into a raw map
        for line in lines:
            line_stripped = line.strip()
            
            if line_stripped.startswith('>P1;'):
                # 1. Save the previous sequence block (if one exists)
                if current_code and current_raw_sequence:
                    sequences_map[current_code] = current_raw_sequence
                
                # 2. Start the new sequence
                current_code = line_stripped.split(';')[1].strip()
                current_raw_sequence = "" # Reset sequence buffer
            
            elif current_code and line_stripped:
                # 3. Append sequence lines, ignoring headers/comments
                if not line_stripped.startswith(('structureX', 'structure', 'sequence', 'CDE:', '#')):
                    current_raw_sequence += line_stripped
        
        # 4. Save the very last sequence in the file (after loop finishes)
        if current_code and current_raw_sequence:
            sequences_map[current_code] = current_raw_sequence

        if not sequences_map:
            raise ValueError(f"No sequences found in ALI file: {ali_file}")

        # Pass 2: Clean all sequences and validate length
        cleaned_map = {}
        
        # Regex: Allow A-Z, gaps (-), BLKs (.), and CHAIN BREAKS (/)
        allowed_chars_re = re.compile(r'[A-Z\-\.\/]') 
        
        len_check = -1
        first_code = ""

        for code, raw_seq in sequences_map.items():
            # Join all matched characters
            cleaned_seq = "".join(allowed_chars_re.findall(raw_seq.upper()))
            
            if not cleaned_seq:
                print(f"[Warning] No sequence data found for code '{code}' in {ali_file}")
                continue

            # Check for consistent alignment lengths
            if len_check == -1:
                len_check = len(cleaned_seq)
                first_code = code
            elif len_check != len(cleaned_seq):
                # This is the error detected in the log
                print(f"    [FATAL ALIGNMENT ERROR] Sequence length mismatch found in {ali_file}:")
                print(f"    Code '{first_code}' has length: {len_check}")
                print(f"    Code '{code}' has length: {len(cleaned_seq)}")
                print(f"    This is often due to special characters (like '/') being missed by the parser.")
                # This will be caught by the check in fixed_region_utils.py,
                # but we raise it here to be explicit.
                raise ValueError(f"Alignment length mismatch: {code} ({len(cleaned_seq)}) vs {first_code} ({len_check})")
            
            cleaned_map[code] = cleaned_seq

        if not cleaned_map:
            raise ValueError(f"Failed to parse any valid sequence data from {ali_file}")

        return cleaned_map
        
    except Exception as e:
        raise IOError(f"Error reading PIR alignment file '{ali_file}'. Check format. Error: {e}")


def extract_cde_line(ali_file_path: str) -> str:
    """
    Extract the CDE line from an alignment file if present.
    (No change)
    """
    try:
        with open(ali_file_path, 'r') as f:
            for line in f:
                if line.startswith("# CDE:"):
                    return line.strip()
    except FileNotFoundError:
        print(f"[CDE] Warning: File not found at {ali_file_path}")
    except Exception as e:
        print(f"[CDE] Warning: Could not read {ali_file_path}. Error: {e}")
    return ""


def add_cde_line_to_pir(clean_ali_path: str, cde_ali_path: str,
                        ss2_file_path: str, target_seq_full: str,
                        align_code_sequence: str):
    """
    Generate a PIR file with CDE line from a clean PIR file.

    The CDE (Computed Disulfide Environment) line contains secondary
    structure information aligned to the target sequence.

    Args:
        clean_ali_path (str): Path to clean PIR file (input).
        cde_ali_path (str): Path to output PIR file with CDE line.
        ss2_file_path (str): Path to PSIPRED SS2 file.
        target_seq_full (str): Complete target sequence.
        align_code_sequence (str): Alignment code for the target sequence.

    Raises:
        RuntimeError: If CDE line generation fails.
        ValueError: If target sequence not found in alignment.
        IndexError: If alignment/SS2 length mismatch occurs.
    """
    print(f"\n[CDE] Generating CDE line for {os.path.basename(clean_ali_path)}...")

    ss_string_full = extract_ss_from_ss2(ss2_file_path, target_seq_full)
    if not ss_string_full:
        raise RuntimeError("Could not generate CDE line from SS2 file.")

    all_aligned_seqs_map = read_sequences_from_ali_temp(clean_ali_path)

    if align_code_sequence not in all_aligned_seqs_map:
            raise ValueError(f"Target '{align_code_sequence}' not found in {clean_ali_path}")

    aligned_target_seq_temp = all_aligned_seqs_map[align_code_sequence]

    cde_content = ""
    ss_index = 0

    for char_in_alignment in aligned_target_seq_temp:
        if char_in_alignment == '-':
            cde_content += '.'
        elif char_in_alignment == '/': # Handle chain break as a gap in SS
            cde_content += '.'
        else: # Is an AA or a BLK
            if ss_index < len(ss_string_full):
                cde_content += ss_string_full[ss_index] # Append H, E, C, or .
                ss_index += 1
            else:
                raise IndexError(f"Error generating CDE: Alignment {align_code_sequence} (len {len(aligned_target_seq_temp)}) "
                                 f"is longer than the SS2 sequence (len {len(ss_string_full)}).")

    cde_line_full = "CDE:" + cde_content
    cde_line_commented = "# " + cde_line_full

    # Write the new CDE file
    with open(clean_ali_path, 'r') as f:
        pir_content_modeller = f.readlines()

    pir_content_cde = []
    found_target_header = False
    for line in pir_content_modeller:
        pir_content_cde.append(line)
        
        # Look for the >P1; header of the target sequence
        if line.strip() == f">P1;{align_code_sequence}":
            found_target_header = True
            continue # Go to the next line (which should be the 'sequence:' line)
        
        # If we found the header, this must be the 'sequence:' line
        if found_target_header and line.startswith("sequence:"):
            pir_content_cde.append(cde_line_commented + "\n")
            found_target_header = False # Reset flag

    with open(cde_ali_path, 'w') as f:
        f.writelines(pir_content_cde)
    
    _flatten_ali_file(cde_ali_path)

    print(f"[CDE] PIR file (with CDE) generated: {os.path.basename(cde_ali_path)}")

def generate_pir_files_AUTO(env: Environ, align_file_modeller: str, align_file_cde: str):
    """
    Generate final PIR files (with and without CDE) in AUTOMATIC mode.

    Uses MODELLER's salign() to automatically align multiple templates
    with the target sequence.

    Args:
        env: MODELLER environment
        align_file_modeller: Output path for clean PIR file
        align_file_cde: Output path for PIR file with CDE line

    Raises:
        ValueError: If template lists have mismatched lengths
        Exception: If salign() or CDE generation fails
    """
    print("\n[STEP 2] Generating automatic alignment with MODELLER salign()")

    num_templates = len(PDB_TEMPLATE_FILES_PATHS)

    aln = Alignment(env)

    if len(PDB_TEMPLATE_FILES_PATHS) != len(ALIGN_CODES_TEMPLATES):
        raise ValueError("PDB_TEMPLATE_FILES_PATHS and ALIGN_CODES_TEMPLATES in config.py must have the same length.")

    print("[STEP 2.1] Aligning multiple templates:")
    for pdb_path, align_code in zip(PDB_TEMPLATE_FILES_PATHS, ALIGN_CODES_TEMPLATES):
        print(f"  -> Adding template: {pdb_path} (Code: {align_code})")
        mdl = Model(env, file=pdb_path)
        aln.append_model(mdl, align_codes=align_code, atom_files=pdb_path)

    aln.append_sequence(sequence_full)
    aln[num_templates].code = ALIGN_CODE_SEQUENCE

    print("[STEP 2.1] Running salign() to align all templates and target...")
    aln.salign()

    aln.write(file=align_file_modeller, alignment_format='PIR')
    
    _flatten_ali_file(align_file_modeller)

    try:
        add_cde_line_to_pir(
            clean_ali_path=align_file_modeller,
            cde_ali_path=align_file_cde,
            ss2_file_path=SS2_FILE_PATH,
            target_seq_full=sequence_full,
            align_code_sequence=ALIGN_CODE_SEQUENCE
        )
    except Exception as e:
        print(f"[ERROR] CDE file generation failed in AUTO mode. Error: {e}")
        raise

    print(f"\n[STEP 2] Alignment PIR file (CLEAN) generated: {align_file_modeller}")


def group_ranges(residue_list: List[int]) -> List[Tuple[int, int]]:
    """
    Group consecutive residue numbers into continuous ranges.

    Args:
        residue_list: List of residue numbers

    Returns:
        List of (start, end) tuples representing continuous ranges

    Example:
        [1, 2, 3, 7, 8, 10] -> [(1, 3), (7, 8), (10, 10)]
    """
    if not residue_list:
        return []
    residue_list = sorted(residue_list)
    ranges = []
    start = residue_list[0]
    end = residue_list[0]
    for n in residue_list[1:]:
        if n == end + 1:
            end = n
        else:
            ranges.append((start, end))
            start = n
            end = n
    ranges.append((start, end))
    return ranges


def get_coil_residues(ss2_file_path: str, seq_full: str) -> Set[int]:
    """
    Identify all residues in coil ('C') conformation from PSIPRED prediction.

    Uses validated SS2 reading that pads with '.' for missing regions,
    automatically excluding non-protein residues (BLK).

    Args:
        ss2_file_path: Path to PSIPRED SS2 file
        seq_full: Complete target sequence

    Returns:
        Set of residue numbers (1-indexed) predicted to be in coil

    Raises:
        ValueError: If SS2 and sequence lengths don't match after processing
    """
    coil_residues = set()

    ss_string = _read_ss2_file(ss2_file_path, seq_full)
    seq_length = len(seq_full)

    if len(ss_string) != seq_length:
        raise ValueError(f"Post-processed SS length ({len(ss_string)}) doesn't match sequence ({seq_length}).")

    for i in range(seq_length):
        if ss_string[i] == 'C':
            coil_residues.add(i + 1)

    print(f"[SS2] Detected {len(coil_residues)} 'C' (Coil) residues for refinement.")
    return coil_residues


def final_evaluation_and_ranking(env: Environ) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """
    Evaluate all generated models and create a ranking by DOPE-HR score.

    Assesses all PDB models generated during the pipeline (excluding templates)
    and ranks them by their DOPE-HR score (lower/more negative is better).

    Args:
        env: MODELLER environment

    Returns:
        Tuple of (full_ranking_list, best_model_dict)

    Note:
        Results are saved to CSV file specified in config.FINAL_RANKING_CSV
    """
    print(f"\n{'='*75}\n[STEP 7] INITIATING FINAL EVALUATION OF ALL PDB MODELS\n")
    template_files_set = set(PDB_TEMPLATE_FILES_NAMES)
    pdbs_to_evaluate = [
        f for f in os.listdir()
        if f.endswith(".pdb")
        and f not in template_files_set
        and (f.startswith("AUTO_") or "_LOOP" in f)
    ]
    if not pdbs_to_evaluate:
        print("[FINAL] No generated PDB files found for evaluation.")
        return [], {}

    final_results: List[Dict[str, Any]] = []
    for filename in pdbs_to_evaluate:
        try:
            mdl = complete_pdb(env, filename)
            atmsel = Selection(mdl.chains[CHAIN_ID])
            dopeHR_score = atmsel.assess_dopehr()
            normalized_dopeHR_zscore = mdl.assess_normalized_dopehr()
            final_results.append({
                'name': filename,
                'DOPEHR score': dopeHR_score,
                'DOPEHR Z-score': normalized_dopeHR_zscore
            })
            print(f"  -> Evaluated {filename:<40} | DOPEHR: {dopeHR_score:.3f} | Z-score: {normalized_dopeHR_zscore:.3f}")
        except Exception as e:
            print(f"  [ERROR] Failed to evaluate {filename}. Error: {e}")
            final_results.append({
                'name': filename,
                'DOPEHR score': float('inf'),
                'DOPEHR Z-score': float('inf')
            })
            continue

    final_ranking = sorted(final_results, key=lambda x: x['DOPEHR score'], reverse=False)
    best_final_models = final_ranking[:NUM_BEST_FINAL_MODELS]

    print(f"\n{'='*75}")
    print(f"FINAL RANKING - Top {NUM_BEST_FINAL_MODELS} Models by DOPEHR Score")
    print(f"{'='*75}\n")
    print(f"{'Rank':<6} {'Model Name':<45} {'DOPEHR':<12} {'Z-score':<12}")
    print(f"{'-'*75}")
    for rank, model_data in enumerate(best_final_models, start=1):
        print(f"{rank:<6} {model_data['name']:<45} {model_data['DOPEHR score']:<12.3f} {model_data['DOPEHR Z-score']:<12.3f}")
    print(f"\n{'='*75}\n")

    csv_filename = FINAL_RANKING_CSV
    try:
        with open(csv_filename, 'w', newline='') as csvfile:
            fieldnames = ['Rank', 'Model Name', 'DOPEHR Score', 'DOPEHR Z-score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for rank, model_data in enumerate(best_final_models, start=1):
                writer.writerow({
                    'Rank': rank,
                    'Model Name': model_data['name'],
                    'DOPEHR Score': f"{model_data['DOPEHR score']:.3f}",
                    'DOPEHR Z-score': f"{model_data['DOPEHR Z-score']:.3f}"
                })
        print(f"[FINAL] Ranking exported to: {csv_filename}\n")
    except Exception as e:
        print(f"[WARNING] Could not write ranking CSV file '{csv_filename}'. Error: {e}")

    best_model = best_final_models[0] if best_final_models else {}
    return final_ranking, best_model
