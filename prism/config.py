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
Configuration module for PRISM pipeline.

This module contains all configurable parameters for the protein modeling pipeline,
including file paths, modeling parameters, and execution settings.
"""

import os
from typing import List, Tuple, Dict, Any


def read_fasta_sequence(fasta_file: str) -> str:
    """
    Read a FASTA file and return the complete sequence as a single string.

    Ignores header lines (starting with '>') and joins all sequence lines.

    Args:
        fasta_file: Path to the FASTA file

    Returns:
        Complete protein sequence as a string

    Raises:
        FileNotFoundError: If the FASTA file doesn't exist
        ValueError: If no sequence is found in the file
    """
    sequence_lines = []
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    continue
                sequence_lines.append(line.strip())

        full_sequence = "".join(sequence_lines)
        if not full_sequence:
            raise ValueError(f"No sequence found in file: {fasta_file}")
        return full_sequence

    except FileNotFoundError:
        print(f"ERROR: Could not find FASTA file at path: {fasta_file}")
        raise
    except Exception as e:
        print(f"ERROR: Failed to read FASTA file: {e}")
        raise


try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
except NameError:
    PROJECT_ROOT = os.getcwd()
    print("WARNING: __file__ not defined. Assuming current working directory is project root.")


INPUT_DIR = os.path.join(PROJECT_ROOT, 'input')
MODELING_RESULTS_DIR = os.path.join(PROJECT_ROOT, 'modeling_results')
PSIPRED_RESULTS_DIR = os.path.join(PROJECT_ROOT, 'psipred_results')

FASTA_FILE_PATH = os.path.join(INPUT_DIR, 'sequence_full.fasta')

SS2_FILE_BASENAME = 'P1_NCL_secondary_structure.ss2'
SS2_FILE_PATH = os.path.join(INPUT_DIR, SS2_FILE_BASENAME)

MANUAL_ALIGNMENT_FILE = os.path.join(INPUT_DIR, 'manual_template_FullSeq.ali')
MANUAL_ALIGNMENT_CDE_FILE = os.path.join(INPUT_DIR, 'manual_template_FullSeq_cde.ali')

ALIGN_CODE_SEQUENCE = 'FullSeq'
CHAIN_ID = "A"

sequence_full = read_fasta_sequence(FASTA_FILE_PATH)

NUM_PROCESSORS = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
NUM_MODELS_AUTO = 32
NUM_MODELS_TO_REFINE = 2
NUM_MODELS_LOOP = 16
NUM_BEST_FINAL_MODELS = 1000

MIN_DIST_FROM_NUCLEOTIDE_COM = 5.0

EXPERIMENTAL_FLANK_SIZE: int = 5
REFINE_FLANKS_DURING_AUTOMODEL: bool = True

USE_MANUAL_ALIGNMENT = True

PERFORM_PSIPRED_PREDICTION = False
PSIPRED_EMAIL = "richard.lopezc@uah.es"
PSIPRED_POLL_INTERVAL = 60

PDB_TEMPLATE_FILES_NAMES: List[str] = [
    '9CB5_DS_renum_HETATM-B.pdb',
    "2fc9_DS_renum.pdb",
    "2fc8_DS_renum.pdb",
    "AF-P19338-F1-model_v6_DS_renum.pdb",
    "9CB5_DS_renum_HETATM-B-1.pdb",
    "9CB5_DS_renum_HETATM-B-2.pdb",
    "9CB5_DS_renum_HETATM-B-3.pdb",
    "9CB5_DS_renum_HETATM-B-4.pdb",
    "9CB5_DS_renum_HETATM-B-5.pdb",
    "9CB5_DS_renum_HETATM-B-6.pdb",
    "9CB5_DS_renum_HETATM-B-7.pdb",
    "9CB5_DS_renum_HETATM-B-8.pdb",
    "9CB5_DS_renum_HETATM-B-9.pdb",
    "9CB5_DS_renum_HETATM-B-10.pdb"
]

ALIGN_CODES_TEMPLATES: List[str] = PDB_TEMPLATE_FILES_NAMES

PDB_TEMPLATE_FILES_PATHS: List[str] = [os.path.join(INPUT_DIR, name) for name in PDB_TEMPLATE_FILES_NAMES]

MAIN_PDB_TEMPLATE_PATH: str = PDB_TEMPLATE_FILES_PATHS[0]
MAIN_ALIGN_CODE_TEMPLATE: str = ALIGN_CODES_TEMPLATES[0]

ALIGNMENT_FILE = os.path.join(MODELING_RESULTS_DIR, f'{MAIN_ALIGN_CODE_TEMPLATE}_{ALIGN_CODE_SEQUENCE}.ali')
ALIGNMENT_CDE_FILE = os.path.join(MODELING_RESULTS_DIR, f'{MAIN_ALIGN_CODE_TEMPLATE}_{ALIGN_CODE_SEQUENCE}_cde.ali')
FINAL_RANKING_CSV = os.path.join(MODELING_RESULTS_DIR, 'final_models_ranking.csv')

__all__ = [
    'PROJECT_ROOT', 'INPUT_DIR', 'MODELING_RESULTS_DIR', 'PSIPRED_RESULTS_DIR',
    'FASTA_FILE_PATH', 'SS2_FILE_PATH', 'PDB_TEMPLATE_FILES_PATHS',
    'MANUAL_ALIGNMENT_FILE', 'MANUAL_ALIGNMENT_CDE_FILE', 'MAIN_PDB_TEMPLATE_PATH',
    'SS2_FILE_BASENAME', 'PDB_TEMPLATE_FILES_NAMES', 'ALIGN_CODES_TEMPLATES',
    'MAIN_ALIGN_CODE_TEMPLATE', 'ALIGN_CODE_SEQUENCE', 'CHAIN_ID',
    'sequence_full',
    'NUM_PROCESSORS', 'NUM_MODELS_AUTO', 'NUM_MODELS_TO_REFINE',
    'NUM_MODELS_LOOP', 'NUM_BEST_FINAL_MODELS',
    'MIN_DIST_FROM_NUCLEOTIDE_COM',
    'USE_MANUAL_ALIGNMENT', 'PERFORM_PSIPRED_PREDICTION',
    'PSIPRED_EMAIL', 'PSIPRED_POLL_INTERVAL',
    'ALIGNMENT_FILE', 'ALIGNMENT_CDE_FILE', 'FINAL_RANKING_CSV',
    'EXPERIMENTAL_FLANK_SIZE', 'REFINE_FLANKS_DURING_AUTOMODEL'
]
