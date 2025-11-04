# PRISM Project Memory

## Project Overview
PRISM (Protected-Region Insertion Suite for Modeling) is a professional Python package for MODELLER-based protein modeling with protected experimental regions and loop refinement.

## Architecture
- **Package Structure**: Modular Python package in `prism/` directory
- **Main Entry Point**: `prism/controller.py`
- **Core Modules**:
  - `config.py`: Centralized configuration
  - `utils.py`: Alignment, SS2, and evaluation utilities
  - `fixed_region_utils.py`: Experimental residue identification
  - `custom_models.py`: FixedRegionAutoModel and FixedRegionLoopRefiner classes
  - `homology_modeling.py`: Parallel homology modeling
  - `loop_refinement.py`: Parallel loop refinement

## Key Features
1. Keeps experimental coordinates completely fixed during optimization
2. Per-nucleotide HETATM repulsion shields
3. Parallel processing support via MODELLER
4. Traceable model nomenclature (AUTO_1_LOOP_2_R1.pdb)
5. Comprehensive DOPE-HR evaluation and ranking

## Recent Changes
- Complete English translation from Spanish proof-of-concept
- Professional docstrings in Google/reST style
- PEP 8 compliance
- Logical step numbering in pipeline output
- Comprehensive README.md documentation

## Dependencies
- MODELLER (10.5+) - must be installed separately
- Python 3.7+
- No additional pip packages required

## Input Requirements
- `input/sequence_full.fasta`: Target sequence
- `input/*.pdb`: Template structures
- `input/P1_NCL_secondary_structure.ss2`: PSIPRED predictions
- Optional: `input/manual_template_FullSeq.ali`: Manual alignment

## Output
All results in `modeling_results/`:
- PDB models with traceable names
- `final_models_ranking.csv`: Complete ranking
- Alignment files (.ali)

## Author
Richard Lopez Corbalan
GitHub: github.com/richardloopez
