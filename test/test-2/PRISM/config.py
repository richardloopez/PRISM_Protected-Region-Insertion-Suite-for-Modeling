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
import yaml
from typing import List, Union, Literal, Dict, Optional, Tuple
from pydantic import BaseModel, EmailStr, Field, computed_field, model_validator

class PrismPowerConfig(BaseModel):
    precalculation: Dict[str, int]
    precomputed: Dict[str, int]

class PrismConfig(BaseModel):

    # Filenames (only basenames)
    ALIGN_CODE_SEQUENCE: str
    CHAIN_ID: str
    BLK_CHAIN_ID: str
    CUSTOM_INIFILE_BASENAME: str
    CUSTOM_RSRFILE_BASENAME: str
    FASTA_FILE_BASENAME: str
    SS2_FILE_BASENAME: str
    MANUAL_ALIGNMENT_BASENAME: str
    MANUAL_ALIGNMENT_CDE_BASENAME: str

    # Parallelization and number of models
    MODELLER_CORES: int = Field(gt=0, alias="MODELLER_CORES")
    TOTAL_PARALLEL_JOBS: int = Field(gt=0, alias="TOTAL_PARALLEL_JOBS")
    TOTAL_HOMOLOGY_MODELS: int = Field(gt=0, alias="TOTAL_HOMOLOGY_MODELS")
    TOP_MODELS_FOR_REFINEMENT: int = Field(ge=0, alias="TOP_MODELS_FOR_REFINEMENT") 
    LOOP_MODELS_PER_TARGET: int = Field(gt=0, alias="LOOP_MODELS_PER_TARGET")
    NUM_BEST_FINAL_MODELS: Union[int, str]

    # Model parameters and execution settings
    PERFORM_PSIPRED_PREDICTION: bool
    PSIPRED_EMAIL: EmailStr
    PSIPRED_POLL_INTERVAL: int
    USE_MANUAL_ALIGNMENT: bool
    BLOCK_REPULSION_RADIUS: float
    REFINE_FLANKS_DURING_AUTOMODEL: bool
    MOBILE_FLANK_RESIDUES: int
    EXECUTION_PARADIGM: Literal["precalculation", "precomputed", "normal", "prism-power"]
    USE_MANUAL_OPTIMIZATION_SELECTION: bool
    MANUAL_OPTIMIZATION_RESIDUES: List[int]
    

    # Templates: First one is the main template
    PDB_TEMPLATE_FILES_NAMES: List[str]
    PRISM_POWER_SETTINGS: Optional[PrismPowerConfig] = None

    # Fixed Directories (Internal Constants)
    INPUT_DIR_NAME: str = "input"
    MODELING_RESULTS_DIR_NAME: str = "modeling_results"
    PSIPRED_RESULTS_DIR_NAME: str = "psipred_results"

    # Validators
    @model_validator(mode="after")
    def check_num_models_to_refine(self):
        if self.TOP_MODELS_FOR_REFINEMENT > self.TOTAL_HOMOLOGY_MODELS:
            raise ValueError("TOP_MODELS_FOR_REFINEMENT cannot be larger than TOTAL_HOMOLOGY_MODELS")
        return self

    @model_validator(mode="after")
    def check_num_best_final_models(self):
        if (isinstance(self.NUM_BEST_FINAL_MODELS, int) and self.NUM_BEST_FINAL_MODELS < 1) or (isinstance(self.NUM_BEST_FINAL_MODELS, str) and self.NUM_BEST_FINAL_MODELS.lower() != "inf"):
            raise ValueError("num_best_final_models must be a positive integer or 'inf'")
        return self

    # Calculated Fields

    # Program dirs 
    @computed_field
    @property
    def PROJECT_ROOT(self) -> str:
        return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    @computed_field
    @property
    def INPUT_DIR(self) -> str:
        return os.path.join(self.PROJECT_ROOT, self.INPUT_DIR_NAME)
    @computed_field
    @property
    def MODELING_RESULTS_DIR(self) -> str:
        return os.path.join(self.PROJECT_ROOT, self.MODELING_RESULTS_DIR_NAME)
    @computed_field
    @property
    def PSIPRED_RESULTS_DIR(self) -> str:
        return os.path.join(self.PROJECT_ROOT, self.PSIPRED_RESULTS_DIR_NAME)
    
    # Files dirs
    @computed_field
    @property
    def CUSTOM_INIFILE_PATH(self) -> str:
        return os.path.join(self.INPUT_DIR, self.CUSTOM_INIFILE_BASENAME)
    @computed_field
    @property
    def CUSTOM_RSRFILE_PATH(self) -> str:
        return os.path.join(self.INPUT_DIR, self.CUSTOM_RSRFILE_BASENAME)
    @computed_field
    @property
    def FASTA_FILE_PATH(self) -> str:
        return os.path.join(self.INPUT_DIR, self.FASTA_FILE_BASENAME)
    @computed_field
    @property
    def SS2_FILE_PATH(self) -> str:
        return os.path.join(self.INPUT_DIR, self.SS2_FILE_BASENAME)
    @computed_field
    @property
    def MANUAL_ALIGNMENT_FILE(self) -> str:
        return os.path.join(self.INPUT_DIR, self.MANUAL_ALIGNMENT_BASENAME)
    @computed_field
    @property
    def MANUAL_ALIGNMENT_CDE_FILE(self) -> str:
        return os.path.join(self.INPUT_DIR, self.MANUAL_ALIGNMENT_CDE_BASENAME)
    @computed_field
    @property
    def PDB_TEMPLATE_FILES_PATHS(self) -> List[str]:
        return [os.path.join(self.INPUT_DIR, pdb_file) for pdb_file in self.PDB_TEMPLATE_FILES_NAMES]
    @computed_field
    @property
    def MAIN_PDB_TEMPLATE_PATH(self) -> str:
        return self.PDB_TEMPLATE_FILES_PATHS[0]
    @computed_field
    @property
    def MAIN_ALIGN_CODE_TEMPLATE(self) -> str:
        return self.PDB_TEMPLATE_FILES_NAMES[0]
    @computed_field
    @property
    def ALIGNMENT_FILE(self) -> str:
        return os.path.join(self.MODELING_RESULTS_DIR, f'{self.MAIN_ALIGN_CODE_TEMPLATE}_{self.ALIGN_CODE_SEQUENCE}.ali')
    @computed_field
    @property
    def ALIGNMENT_CDE_FILE(self) -> str:
        return os.path.join(self.MODELING_RESULTS_DIR, f'{self.MAIN_ALIGN_CODE_TEMPLATE}_{self.ALIGN_CODE_SEQUENCE}_cde.ali')
    @computed_field
    @property
    def FINAL_RANKING_CSV(self) -> str:
        return os.path.join(self.MODELING_RESULTS_DIR, f'final_ranking.csv')

# ============================================================================
# INSTANTIATION LOGIC
# ============================================================================

    def save_settings(self, yaml_path: str = "config.yaml"):
        '''
        Saves current settings back to a YAML file.
        '''
        data = self.model_dump(exclude={"PROJECT_ROOT", "INPUT_DIR", "MODELING_RESULTS_DIR", "PSIPRED_RESULTS_DIR",
                                        "CUSTOM_INIFILE_PATH", "CUSTOM_RSRFILE_PATH", "FASTA_FILE_PATH", "SS2_FILE_PATH",
                                        "MANUAL_ALIGNMENT_FILE", "MANUAL_ALIGNMENT_CDE_FILE", "PDB_TEMPLATE_FILES_PATHS",
                                        "MAIN_PDB_TEMPLATE_PATH", "MAIN_ALIGN_CODE_TEMPLATE", "ALIGNMENT_FILE",
                                        "ALIGNMENT_CDE_FILE", "FINAL_RANKING_CSV"})
        
        data['PSIPRED_EMAIL'] = str(data['PSIPRED_EMAIL'])
        data = {k: v for k, v in data.items() if v is not None}

        full_yaml_path = os.path.join(self.PROJECT_ROOT, yaml_path)
        with open(full_yaml_path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

def load_settings(yaml_path: Optional[str] = None) -> PrismConfig:
    '''
    Loads the YAML file and returns a validated PrismConfig object.
    '''
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    
    if yaml_path is None:
        yaml_path = os.path.join(project_root, "config.yaml")
    elif not os.path.isabs(yaml_path):
        yaml_path = os.path.join(project_root, yaml_path)

    if not os.path.exists(yaml_path):
        raise FileNotFoundError(f"Config file not found at {yaml_path}")
        
    with open(yaml_path, "r") as f:
        raw_data = yaml.safe_load(f)
    return PrismConfig(**raw_data)

settings = load_settings()


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# sequence (def here to avoid circular imports from utils)
sequence_full = None

def read_fasta_sequence(file_path: str) -> str:
    '''
    Read a FASTA file and return the sequence.
    '''
    sequence = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            else:
                sequence.append(line.strip())
    return ''.join(sequence)

def get_sequence():
    '''
    Get the sequence from the FASTA file.
    '''
    global sequence_full
    if sequence_full is None:
        sequence_full = read_fasta_sequence(settings.FASTA_FILE_PATH)
    return sequence_full

# ============================================================================
# INJECTION
# ============================================================================

for field_name in PrismConfig.model_fields:
    globals()[field_name] = getattr(settings, field_name)

for computed_name in PrismConfig.model_computed_fields:
    try:
        globals()[computed_name] = getattr(settings, computed_name)
    except AttributeError:
        continue

__all__ = ['settings', 'sequence_full'] + list(PrismConfig.model_fields.keys()) + list(PrismConfig.model_computed_fields.keys())

if __name__ == "__main__":
    print(settings.model_dump_json())