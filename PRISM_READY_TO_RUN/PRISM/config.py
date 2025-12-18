#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
# PRISM (Protected-Region Insertion Suite for Modeling) (config)
#

'''
PRISM: Protected-Region Insertion Suite for Modeling

A MODELLER pipeline for high-fidelity homology modeling and loop
refinement while maintaining experimental core coordinates completely fixed.

This package prevents coordinate drift of experimental regions during optimization
and uses HETATM repulsion shields to model loops in complex environments.
'''

# Citation:
# If you use this software in your research, please cite:
# [1] R. Lopez-Corbalan, "PRISM: Protected-Region Insertion Suite for Modeling", GitHub, 2025. [Online]. Available: https://github.com/richardloopez/PRISM_Protected-Region-Insertion-Suite-for-Modeling. [Accessed: Date].

'''
Configuration module for PRISM pipeline.

This module contains all configurable parameters for PRISM pipeline,
including file paths, model parameters, and execution settings.

'''

import os
from typing import List

# File paths
ALIGN_CODE_SEQUENCE = 'FullSeq'
CHAIN_ID = 'A'
CUSTOM_INIFILE_BASENAME = 'precomputed_ini.pdb'
CUSTOM_RSRFILE_BASENAME = 'precomputed_rsr.rsr'
FASTA_FILE_BASENAME = 'sequence_full.fasta'
SS2_FILE_BASENAME = 'P1_NCL_secondary_structure.ss2'
MANUAL_ALIGNMENT_BASENAME = 'manual_template_FullSeq.ali'
MANUAL_ALIGNMENT_CDE_BASENAME = 'manual_template_FullSeq_cde.ali'

# Parallelization and number of models
NUM_PROCESSORS = int(10)
NUM_SPLITTED_PROCESSES = 10
NUM_MODELS_AUTO = 100
NUM_MODELS_TO_REFINE = 10
NUM_MODELS_LOOP = 10
NUM_BEST_FINAL_MODELS = 1000000000000000000000000000000000000000000000000000000000000000000000000000

# Model parameters and execution settings
PERFORM_PSIPRED_PREDICTION = False
PSIPRED_EMAIL = "richard.lopezc@uah.es"
PSIPRED_POLL_INTERVAL = 60
USE_MANUAL_ALIGNMENT = True
MIN_DIST_FROM_NUCLEOTIDE_COM = 5.0
REFINE_FLANKS_DURING_AUTOMODEL = True
EXPERIMENTAL_FLANK_SIZE = 5
USE_PRECOMPUTED_FILES = True

RSR_INI_PRECALCULATION = False

# Templates: First one is the main template
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
    "9CB5_DS_renum_HETATM-B-10.pdb",
    "2fc9_DS_renum-1.pdb",
    "2fc9_DS_renum-2.pdb",
    "2fc9_DS_renum-3.pdb",
    "2fc9_DS_renum-4.pdb",
    "2fc9_DS_renum-5.pdb",
    "2fc8_DS_renum-1.pdb",
    "2fc8_DS_renum-2.pdb",
    "2fc8_DS_renum-3.pdb",
    "2fc8_DS_renum-4.pdb",
    "2fc8_DS_renum-5.pdb"
]

# FILE PATHS

# main directory
try:
    _current_dir = os.path.dirname(os.path.abspath(__file__))
    PROJECT_ROOT = os.path.dirname(_current_dir)
except NameError:
    PROJECT_ROOT = os.getcwd()
    print(f'Warning: PROJECT_ROOT not found, using current working directory: {PROJECT_ROOT}')


# input directory
INPUT_DIR = os.path.join(PROJECT_ROOT, 'input')
MODELING_RESULTS_DIR = os.path.join(PROJECT_ROOT, 'modeling_results')
CUSTOM_INIFILE_PATH = os.path.join(INPUT_DIR, CUSTOM_INIFILE_BASENAME)
CUSTOM_RSRFILE_PATH = os.path.join(INPUT_DIR, CUSTOM_RSRFILE_BASENAME)
FASTA_FILE_PATH = os.path.join(INPUT_DIR, FASTA_FILE_BASENAME)
PSIPRED_RESULTS_DIR = os.path.join(PROJECT_ROOT, 'psipred_results')
SS2_FILE_PATH = os.path.join(INPUT_DIR, SS2_FILE_BASENAME)
MANUAL_ALIGNMENT_FILE = os.path.join(INPUT_DIR, MANUAL_ALIGNMENT_BASENAME)
MANUAL_ALIGNMENT_CDE_FILE = os.path.join(INPUT_DIR, MANUAL_ALIGNMENT_CDE_BASENAME)

# sequence (def here to avoid circular imports from utils)
sequence_full = None

def read_fasta_sequence(file_path: str) -> str:
    """
    Read a FASTA file and return the sequence.
    """
    sequence = []
    try: 
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    continue
                else:
                    sequence.append(line.strip())
        return ''.join(sequence)
    except FileNotFoundError:
        print(f"Error: File not found: {file_path}")
        return None
    except Exception as e:
        print(f"Error reading file: {file_path}") 
        print(str(e))
        return None


def get_sequence():
    """
    Get the sequence from the FASTA file.
    """
    global sequence_full
    if sequence_full is None:
        sequence_full = read_fasta_sequence(FASTA_FILE_PATH)
    return sequence_full


# templates
PDB_TEMPLATE_FILES_PATHS: List[str] = [os.path.join(INPUT_DIR, name) for name in PDB_TEMPLATE_FILES_NAMES]
MAIN_PDB_TEMPLATE_PATH: str = PDB_TEMPLATE_FILES_PATHS[0]
MAIN_ALIGN_CODE_TEMPLATE: str = PDB_TEMPLATE_FILES_NAMES[0]

# alignment files
ALIGNMENT_FILE = os.path.join(MODELING_RESULTS_DIR, f'{MAIN_ALIGN_CODE_TEMPLATE}_{ALIGN_CODE_SEQUENCE}.ali')
ALIGNMENT_CDE_FILE = os.path.join(MODELING_RESULTS_DIR, f'{MAIN_ALIGN_CODE_TEMPLATE}_{ALIGN_CODE_SEQUENCE}_cde.ali')
FINAL_RANKING_CSV = os.path.join(MODELING_RESULTS_DIR, 'final_models_ranking.csv')

__all__ = [
    'ALIGN_CODE_SEQUENCE', 'CHAIN_ID', 'CUSTOM_INIFILE_BASENAME', 'CUSTOM_RSRFILE_BASENAME',
    'FASTA_FILE_BASENAME', 'SS2_FILE_BASENAME', 'MANUAL_ALIGNMENT_BASENAME',
    'MANUAL_ALIGNMENT_CDE_BASENAME',

    'NUM_PROCESSORS', 'NUM_SPLITTED_PROCESSES', 'NUM_MODELS_AUTO', 'NUM_MODELS_TO_REFINE',
    'NUM_MODELS_LOOP', 'NUM_BEST_FINAL_MODELS',

    'PERFORM_PSIPRED_PREDICTION', 'PSIPRED_EMAIL', 'PSIPRED_POLL_INTERVAL', 'USE_MANUAL_ALIGNMENT',
    'MIN_DIST_FROM_NUCLEOTIDE_COM', 'REFINE_FLANKS_DURING_AUTOMODEL', 'EXPERIMENTAL_FLANK_SIZE', 'USE_PRECOMPUTED_FILES',
    'RSR_INI_PRECALCULATION',

    'PDB_TEMPLATE_FILES_NAMES',

    'PROJECT_ROOT',

    'INPUT_DIR', 'MODELING_RESULTS_DIR', 'CUSTOM_INIFILE_PATH', 'CUSTOM_RSRFILE_PATH', 'FASTA_FILE_PATH',
    'PSIPRED_RESULTS_DIR', 'SS2_FILE_PATH', 'MANUAL_ALIGNMENT_FILE', 'MANUAL_ALIGNMENT_CDE_FILE',

    'sequence_full',

    'PDB_TEMPLATE_FILES_PATHS', 'MAIN_PDB_TEMPLATE_PATH', 'MAIN_ALIGN_CODE_TEMPLATE',

    'ALIGNMENT_FILE', 'ALIGNMENT_CDE_FILE', 'FINAL_RANKING_CSV'

]