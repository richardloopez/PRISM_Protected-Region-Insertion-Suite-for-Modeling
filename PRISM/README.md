# PRISM: Protected-Region Insertion Suite for Modeling

**A professional MODELLER pipeline for high-fidelity protein homology modeling and loop refinement with guaranteed preservation of experimental coordinates.**

Author: Richard Lopez Corbalan  
GitHub: [github.com/richardloopez](https://github.com/richardloopez)

---

## Overview

PRISM is a specialized protein modeling pipeline built on MODELLER that addresses a critical challenge in structural biology: **preventing coordinate drift of experimental regions during computational optimization**. 

The pipeline enables high-fidelity modeling of loop insertions and flexible regions while maintaining the structural integrity of experimentally-determined core coordinates with 100% fidelity.

### Key Features

- **🔒 Absolute Core Protection**: Experimental coordinates remain completely fixed during all optimization steps
- **🧬 HETATM Repulsion Shields**: Per-nucleotide repulsion prevents loop clashes with bound molecules (DNA/RNA/ligands)
- **⚡ Parallel Processing**: Full support for multi-processor execution via MODELLER's parallel framework
- **📊 Traceable Nomenclature**: Models are systematically named (e.g., `AUTO_1_LOOP_2_R1.pdb`) for complete traceability
- **🎯 Smart Loop Detection**: Combines PSIPRED secondary structure predictions with experimental boundary analysis
- **📈 Comprehensive Evaluation**: Automatic ranking using DOPE-HR scores with detailed CSV output

---

## The Problem PRISM Solves

### 1. Experimental Core Drift

During standard homology modeling and loop refinement, optimization algorithms can inadvertently modify experimentally-determined coordinates. Even small deviations accumulate and compromise the reliability of the structural core. PRISM prevents this by:

- Identifying all residues that map to experimental coordinates in the template
- Completely excluding these residues from optimization (`select_atoms()`)
- Maintaining coordinate integrity across all refinement cycles

### 2. Loop-Ligand Clashes

When modeling loops in the presence of bound molecules (DNA, RNA, ligands), standard approaches often produce steric clashes. PRISM addresses this through:

- **Per-nucleotide gravity center calculation**: Each HETATM residue gets a pseudo-atom at its center of mass
- **Lower-bound distance restraints**: C-alpha atoms in flexible regions maintain minimum separation from bound molecules
- **Configurable repulsion strength**: Adjustable minimum distance (default: 5.0 Å) balances accuracy with sampling

### 3. Experimental Flank Refinement

The boundaries between experimental and modeled regions often require special treatment. PRISM enables:

- **Configurable flank size**: Define how many edge residues can be refined if predicted as coil
- **Secondary structure awareness**: Only refine flank residues predicted to be flexible ('C')
- **Deep core preservation**: Inner experimental residues remain absolutely fixed

---

## Installation

### Prerequisites

1. **MODELLER** (version 10.5 or higher)
   - Download from: https://salilab.org/modeller/
   - Academic license required (free): https://salilab.org/modeller/registration.html
   - Ensure MODELLER is in your Python path

2. **Python** (3.7 or higher)

3. **PSIPRED** (optional, for secondary structure prediction)
   - If using pre-computed `.ss2` files, PSIPRED is not required
   - Set `PERFORM_PSIPRED_PREDICTION = False` in `config.py`

### Setup

```bash
# Clone or download PRISM
git clone https://github.com/richardloopez/PRISM.git
cd PRISM

# Verify MODELLER installation
python3 -c "import modeller; print(modeller.__version__)"

# No additional Python packages required (uses MODELLER's environment)
```

---

## Quick Start

### 1. Prepare Input Files

Place the following files in the `input/` directory:

- **`sequence_full.fasta`**: Your target protein sequence in FASTA format
- **Template PDB files**: Experimental structures to use as templates
- **`P1_NCL_secondary_structure.ss2`**: PSIPRED secondary structure prediction
- **`manual_template_FullSeq.ali`** (optional): Manual alignment in PIR format

### 2. Configure Parameters

Edit `prism/config.py` to set:

```python
# Template files (must be in input/ directory)
PDB_TEMPLATE_FILES_NAMES = [
    '9CB5_DS_renum_HETATM-B.pdb',  # Main template (first in list)
    '2fc9_DS_renum.pdb',
    # ... additional templates
]

# Modeling parameters
NUM_MODELS_AUTO = 16          # Number of initial homology models
NUM_MODELS_TO_REFINE = 2      # Top models to refine
NUM_MODELS_LOOP = 2           # Models per loop refinement
EXPERIMENTAL_FLANK_SIZE = 5   # Refinable edge residues (0 to disable)

# Processing
NUM_PROCESSORS = 4            # Parallel workers (adjust for your system)

# Alignment mode
USE_MANUAL_ALIGNMENT = True   # True = use manual, False = auto salign()
```

### 3. Run the Pipeline

```bash
# Run as a module (recommended for proper package imports)
python3 -m prism.controller

# Alternative: Run from within the prism directory
cd prism
python3 -c "from controller import main_workflow; main_workflow()"
```

### 4. Check Results

All outputs are written to `modeling_results/`:

- **`AUTO_*.pdb`**: Initial homology models (ranked by quality)
- **`AUTO_*_LOOP*_R*.pdb`**: Loop-refined models with traceable nomenclature
- **`final_models_ranking.csv`**: Complete ranking with DOPE-HR scores

---

## Configuration Reference

### Core Parameters

| Parameter | Description | Default | Notes |
|-----------|-------------|---------|-------|
| `NUM_MODELS_AUTO` | Initial homology models to generate | 16 | Higher = better sampling, longer runtime |
| `NUM_MODELS_TO_REFINE` | Top models selected for loop refinement | 2 | Typically 1-5 models |
| `NUM_MODELS_LOOP` | Models generated per loop refinement | 2 | Should match `NUM_PROCESSORS` for efficiency |
| `NUM_PROCESSORS` | Parallel workers | 1 | Set to CPU count for maximum speed |
| `NUM_BEST_FINAL_MODELS` | Models included in final ranking | 1000 | Limits CSV output size |

### Experimental Region Control

| Parameter | Description | Default | Impact |
|-----------|-------------|---------|--------|
| `EXPERIMENTAL_FLANK_SIZE` | Edge residues to consider for refinement | 5 | 0 = no flank refinement, 10+ = aggressive |
| `MIN_DIST_FROM_NUCLEOTIDE_COM` | Minimum C-alpha to HETATM distance (Å) | 5.0 | Lower = tighter packing, higher = safer |

### Alignment Mode

| Parameter | Description | Options |
|-----------|-------------|---------|
| `USE_MANUAL_ALIGNMENT` | Use manual or automatic alignment | `True` (manual) / `False` (auto salign) |

**Manual Mode**: Provide `input/manual_template_FullSeq.ali` in PIR format  
**Automatic Mode**: MODELLER generates alignment using `salign()` with multiple templates

### PSIPRED Integration

| Parameter | Description | Default |
|-----------|-------------|---------|
| `PERFORM_PSIPRED_PREDICTION` | Query PSIPRED web server | `False` |
| `PSIPRED_EMAIL` | Email for PSIPRED submission | Required if `True` |
| `SS2_FILE_BASENAME` | Secondary structure filename | `P1_NCL_secondary_structure.ss2` |

**Recommended**: Run PSIPRED offline and place `.ss2` file in `input/` directory

---

## Pipeline Workflow

```
┌─────────────────────────────────────────────────────────┐
│ STEP 1: Alignment Preparation                          │
│  • Manual mode: Use provided PIR alignment             │
│  • Auto mode: Generate with MODELLER salign()          │
│  • Add CDE line with secondary structure                │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│ STEP 2: Experimental Residue Identification            │
│  • Map template coordinates to target sequence          │
│  • Identify core regions (both seqs have residues)      │
│  • Calculate experimental flanks                        │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│ STEP 3: Loop Detection                                 │
│  • Extract 'C' (coil) residues from PSIPRED            │
│  • Exclude deep core (inner experimental residues)     │
│  • Include flank 'C' residues if enabled               │
│  • Group into continuous loop ranges                    │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│ STEP 4: Homology Modeling (FixedRegionAutoModel)       │
│  • Use all templates for gap filling                    │
│  • Freeze experimental core via select_atoms()          │
│  • Apply HETATM repulsion shields                       │
│  • Generate and rank N models → AUTO_1.pdb, AUTO_2.pdb │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│ STEP 5: Loop Refinement (FixedRegionLoopRefiner)       │
│  • Sequential refinement of each loop region            │
│  • Maintain deep core as absolutely fixed               │
│  • Apply HETATM shields per loop                        │
│  • Traceable naming: AUTO_1_LOOP1_R1.pdb, etc.         │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│ STEP 6: Final Evaluation                               │
│  • Assess all models with DOPE-HR                       │
│  • Rank by score (more negative = better)              │
│  • Export to final_models_ranking.csv                   │
└─────────────────────────────────────────────────────────┘
```

---

## Model Nomenclature

PRISM uses systematic, traceable naming:

| Filename Pattern | Meaning | Example |
|-----------------|---------|---------|
| `AUTO_N.pdb` | Initial homology model (rank N by DOPE-HR) | `AUTO_1.pdb` |
| `AUTO_N_LOOPJ_RK.pdb` | Loop-refined model: base N, loop J, refinement K | `AUTO_1_LOOP2_R1.pdb` |

**Reading `AUTO_1_LOOP2_R3.pdb`**:
- Started from `AUTO_1.pdb` (best initial model)
- Refined loop #2 sequentially
- 3rd ranked refinement for that loop

**Best Model Strategy**: The final ranking includes all intermediate models, allowing you to select the globally best model regardless of refinement stage.

---

## Output Files

### In `modeling_results/`

| File | Description |
|------|-------------|
| `AUTO_*.pdb` | Initial homology models |
| `AUTO_*_LOOP*_R*.pdb` | Loop-refined models |
| `final_models_ranking.csv` | Complete model ranking with scores |
| `*.ali` | Alignment files (clean and with CDE line) |

### Ranking CSV Format

```csv
Rank,Model Name,DOPEHR Score,DOPEHR Z-score
1,AUTO_1_LOOP3_R1.pdb,-52847.234,-1.234
2,AUTO_2_LOOP3_R1.pdb,-52523.456,-1.156
...
```

---

## Advanced Usage

### Custom Template Selection

Edit `PDB_TEMPLATE_FILES_NAMES` in `config.py`:

```python
PDB_TEMPLATE_FILES_NAMES = [
    'main_template.pdb',      # FIRST = main template (defines experimental core)
    'template2.pdb',          # Additional templates for better coverage
    'template3.pdb',
]
```

**Important**: The FIRST template defines the experimental core that will be frozen.

### Adjusting Repulsion Strength

```python
# In config.py
MIN_DIST_FROM_NUCLEOTIDE_COM = 5.0  # Ångströms

# Tighter packing (risk of clashes):
MIN_DIST_FROM_NUCLEOTIDE_COM = 4.0

# Conservative (more breathing room):
MIN_DIST_FROM_NUCLEOTIDE_COM = 6.0
```

### Disabling Flank Refinement

```python
# In config.py
EXPERIMENTAL_FLANK_SIZE = 0  # No flank residues will be refined
```

### HPC/Cluster Usage

```python
# In config.py
NUM_PROCESSORS = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
```

Use with SLURM:
```bash
#SBATCH --cpus-per-task=16
python3 controller.py
```

---

## Troubleshooting

### Issue: "No sequence found in file"

**Solution**: Ensure `sequence_full.fasta` is properly formatted:
```
>YourProteinName
MAKTQVLVL...
```

### Issue: "SS2 file is longer/shorter than sequence"

**Solution**: Regenerate PSIPRED prediction to match your exact target sequence length. Ensure no extra residues or missing regions.

### Issue: "Loop contains experimental residues"

**Solution**: This is intentional protection. The loop overlaps with the experimental core and should not be refined. Check your alignment and secondary structure predictions.

### Issue: Models have steric clashes with DNA/ligands

**Solution**: Increase `MIN_DIST_FROM_NUCLEOTIDE_COM` in `config.py`:
```python
MIN_DIST_FROM_NUCLEOTIDE_COM = 6.0  # More conservative
```

### Issue: Long runtime

**Solutions**:
- Reduce `NUM_MODELS_AUTO` (try 8 instead of 16)
- Reduce `NUM_MODELS_TO_REFINE` (try 1 instead of 2)
- Increase `NUM_PROCESSORS` to use more CPU cores
- Use faster refinement: `refine.fast` instead of `refine.slow_large` in `loop_refinement.py`

---

## Citation

If you use PRISM in your research, please cite:

```
Lopez-Corbalan, R. (2025). PRISM: Protected-Region Insertion Suite for Modeling.
GitHub: https://github.com/richardloopez/PRISM
```

---

## License

PRISM is distributed under the terms specified by the author. For commercial use or redistribution, please contact Richard Lopez Corbalan.

**Note**: PRISM requires MODELLER, which has its own licensing terms. Academic users can obtain free licenses from https://salilab.org/modeller/registration.html.

---

## Technical Details

### Architecture

PRISM is organized into modular components:

- **`controller.py`**: Main pipeline orchestrator
- **`config.py`**: Centralized configuration management
- **`utils.py`**: Alignment, secondary structure, and evaluation utilities
- **`fixed_region_utils.py`**: Experimental residue identification
- **`custom_models.py`**: Specialized AutoModel and LoopModel classes
- **`homology_modeling.py`**: Parallel homology modeling execution
- **`loop_refinement.py`**: Parallel loop refinement execution

### Core Innovation: Fixed-Region Selection

```python
def select_atoms(self):
    """Exclude experimental residues from optimization"""
    all_residues = Selection(self).only_std_residues()
    fixed_selection = Selection()
    
    for res_num in self.experimental_residues:
        res_range = self.residue_range(f'{res_num}:{self.chain_id}', ...)
        fixed_selection = fixed_selection | Selection(res_range)
    
    return all_residues - fixed_selection  # Only optimize non-experimental
```

### HETATM Repulsion Implementation

```python
def nonstd_restraints(self, aln):
    """Apply per-nucleotide repulsion shield"""
    for het_res in hetatm_residues:
        pseudo = pseudo_atom.GravityCenter(Selection(het_res))
        rsr.pseudo_atoms.append(pseudo)
        
        for ca_atom in optimizable_ca_atoms:
            rsr.add(forms.LowerBound(
                group=physical.xy_distance,
                feature=features.Distance(ca_atom, pseudo),
                mean=MIN_DIST_FROM_NUCLEOTIDE_COM,
                stdev=1.0
            ))
```

---

## Support

For questions, issues, or contributions, please contact:

- **Author**: Richard Lopez Corbalan
- **GitHub**: https://github.com/richardloopez
- **Email**: richardlopezcorbalan@gmail.com

---

## Acknowledgments

Built on the powerful [MODELLER](https://salilab.org/modeller/) comparative modeling framework developed by the Sali Lab at UCSF.

PSIPRED secondary structure predictions: Jones, D.T. (1999) J. Mol. Biol. 292: 195-202.
