# PRISM: Protected-Region Insertion Suite for Modeling

**A professional MODELLER pipeline for high-fidelity protein homology modeling and loop refinement with guaranteed preservation of experimental coordinates.**

Author: Richard Lopez-Corbalan  
GitHub: [github.com/richardloopez](https://github.com/richardloopez)

---

## 0. Overview

PRISM is a specialized protein modeling pipeline built on MODELLER that addresses a critical challenge in structural biology: **preventing coordinate drift of experimental regions during computational optimization**. 

The pipeline enables high-fidelity modeling of loop insertions and flexible regions while maintaining the structural integrity of experimentally-determined core coordinates with 100% fidelity.

### Key Features

- **🔒 Absolute Core Protection**: Experimental coordinates remain completely fixed during all optimization steps
- **🧬 HETATM Repulsion Shields**: Per-nucleotide repulsion prevents loop clashes with bound molecules (DNA/RNA/ligands)
- **⚡ Parallel Processing**: Full support for multi-processor execution via MODELLER's parallel framework
- **🔮 Automatic PSIPRED Prediction**: Integrates a client (psipred_client.py) to automatically submit, poll, and download secondary structure predictions from the PSIPRED web server.
- **📊 Traceable Nomenclature**: Models are systematically named (e.g., `AUTO_1_LOOP_2_R1.pdb`) for complete traceability
- **🚀 HPC-Ready**: Includes a run_prism.sh script for easy submission to SLURM queuing systems.
- **🔬 Advanced Flank Control**: Granular user control over the behavior of experimental flanks (the junction between fixed and modeled regions) using EXPERIMENTAL_FLANK_SIZE and REFINE_FLANKS_DURING_AUTOMODEL flags.
- **🎯 Smart Loop Detection**: Combines PSIPRED secondary structure predictions with experimental boundary analysis
- **📈 Comprehensive Evaluation**: Automatic ranking using DOPE-HR scores with detailed CSV output

---

## 1. The Problem PRISM Solves

### 1.1. Experimental Core Drift

During standard homology modeling and loop refinement, optimization algorithms can inadvertently modify experimentally-determined coordinates. Even small deviations accumulate and compromise the reliability of the structural core. PRISM prevents this by:

- Identifying all residues that map to experimental coordinates in the template
- Completely excluding these residues from optimization (`select_atoms()`)
- Maintaining coordinate integrity across all refinement cycles

### 1.2. Loop-Ligand Clashes

When modeling loops in the presence of bound molecules (DNA, RNA, ligands), standard approaches often produce steric clashes. PRISM addresses this through:

- **Per-nucleotide gravity center calculation**: Each HETATM residue gets a pseudo-atom at its center of mass
- **Lower-bound distance restraints**: C-alpha atoms in flexible regions maintain minimum separation from bound molecules
- **Configurable repulsion strength**: Adjustable minimum distance (default: 5.0 Å) balances accuracy with sampling


### 1.3. The Template Averaging Problem

MODELLER generates its initial structure by averaging the coordinates of all provided templates. This creates a critical problem:

- **The Problem**: When mixing an experimental PDB (e.g., residues 306-472) with a full-length AlphaFold model (e.g., 1-710), the 306-472 region will be averaged. Its coordinates will move (RMSD > 0), even if marked as "fixed" for later optimization.
- **The Solution**: PRISM uses a "replication trick" to give overwhelming weight to the main template. By providing 100 or 1000 copies of the main PDB, its contribution to the average dominates, ensuring the experimental core's RMSD remains at 0.0. The tools/replicate_template.py script automates this process.

### 1.4. Experimental Flank Refinement

The junction between the 0.0 RMSD experimental core and the new (e.g., AlphaFold) model can be highly strained. PRISM provides granular control over this junction (the "flank").

- **EXPERIMENTAL_FLANK_SIZE**: Defines a buffer zone (e.g., 5 residues) at the edge of the experimental core.
- **REFINE_FLANKS_DURING_AUTOMODEL**: This boolean flag provides full control:
    - False (Default / Safe Mode): AutoModel (Step 4) freezes the entire experimental region (core + flanks). This provides a maximally stable scaffold. Then, LoopModel (Step 5) intelligently refines only the flank residues that are also predicted as 'C' (coil) by PSIPRED, providing a data-driven tension release.
    - True (Advanced / Aggressive Mode): AutoModel (Step 4) freezes only the "deep core" (core - flanks). The flanks are optimized immediately (regardless of secondary structure) to resolve tension from the very beginning. This may be useful for severe clashes but risks destroying stable secondary structures in the flank.


---

## 2. Best Practices

### 2.1. PDB Preparation & Cleaning

The code is designed to receive **clean inputs**.

- **Cleaning**: Remove purification tags, non-target proteins, and unwanted solvent molecules. Tools like Discovery Studio Visualizer are useful for this.
    - **Main Template** (MAIN_PDB): This is the PDB whose coordinates (and optionally ligand) will be preserved.
        - The protein MUST be in Chain A.
        - The desired ligand should be placed in a separate chain (e.g., Chain B).
    - **Other Templates**:
        - The protein MUST be in Chain A.
        - They MUST NOT contain ligands.
- **Utilities**:
    - tools/test_hetatm.py: Use this to analyze your PDBs. It reports which chains are ATOM vs. HETATM and how many BLK (.) characters .ali file needs.
    - tools/renumber_pdb.py: Use this to renumber residues sequentially and, critically, to convert ATOM records to HETATM for specific chains (e.g., convert Chain B ligand to HETATM).

### 2.2. The Averaging Problem & The Replication Trick

As mentioned in Section 1.4, MODELLER averages templates. To achieve a 0.0 RMSD on your experimental core, the replication trick should be used.

- Use **tools/replicate_template.py** to generate many copies (e.g., 1000) of your main PDB.
    - This ensures the experimental coordinates are not distorted by other templates (like AlphaFold models) during the initial AutoModel averaging.
    - If RMSD > 0.0 when using tools/verify_rmsd.py, the cause is likely an insufficient number of replicas.

### 2.3. Alignment Verification

The code can run in a fully automatic mode (USE_MANUAL_ALIGNMENT = False) by generating an alignment with salign().

- **WARNING!** The auto-generated alignment is not always accurate.
    - It is strongly recommended to **first** run the pipeline with **USE_MANUAL_ALIGNMENT = False** and PERFORM_PSIPRED_PREDICTION = False. The script will quickly generate an alignment file in modeling_results/.
    - **Manually inspect** this .ali file. If it's incorrect, edit it, save it as input/manual_template_FullSeq.ali, and switch to USE_MANUAL_ALIGNMENT = True in config.py.

---

## 3. Installation

### 3.1. Prerequisites

1. **MODELLER** (version 10.5 or higher)
   - Download from: https://salilab.org/modeller/
   - Academic license required (free): https://salilab.org/modeller/registration.html
   - Ensure MODELLER is in your Python path

2. **Python** (3.7 or higher)

3. **PSIPRED** (optional, for secondary structure prediction)
   - If using pre-computed `.ss2` files, PSIPRED is not required
   - Set `PERFORM_PSIPRED_PREDICTION = False` in `config.py`

### 3.2. Setup

```bash
# Clone or download PRISM
git clone https://github.com/richardloopez/PRISM.git
cd PRISM

# Verify MODELLER installation
python3 -c "import modeller; print(modeller.__version__)"

# No additional Python packages required (uses MODELLER's environment)
```
### 3.3. Directory Structure

PRISM/
├── input/                 # Input files (PDBs, FASTA, .ss2, .ali)
├── prism/                 # The main package source code (.py)
├── tools/                 # Utility scripts for preparation & analysis
├── modeling_results/      # Default output directory
├── psipred_results/       # Default output from psipred_client
└── run_prism.sh           # SLURM execution script

---

## 4. Quick Start

### 4.1. Prepare Input Files

Place the following files in the `input/` directory:

- **`sequence_full.fasta`**: Your target protein sequence in FASTA format
- **Template PDB files**: Experimental structures to use as templates
- **`P1_NCL_secondary_structure.ss2`**(Optional): PSIPRED secondary structure prediction. Required if PERFORM_PSIPRED_PREDICTION = False.
- **`manual_template_FullSeq.ali`** (optional): Manual alignment in PIR format

### 4.2. Configure Parameters

Edit `prism/config.py` to set:

```python
# Use the manual .ali from 'input/'? (True) or generate one automatically? (False)
USE_MANUAL_ALIGNMENT = True

# Request a new .ss2 prediction from PSIPRED? (True) or use the .ss2 file from 'input/'? (False)
PERFORM_PSIPRED_PREDICTION = False
PSIPRED_EMAIL = "your_email@domain.com" # Required if True

# --- Template Files ---
PDB_TEMPLATE_FILES_NAMES = [
    '9CB5_DS_renum_HETATM-B.pdb',  # Main Template (first in list)
    '2fc9_DS_renum.pdb',
    # ... additional templates
    # ... 1000 replicas of the main PDB (if using replicate_template.py)
]

# --- Modeling Parameters ---
NUM_MODELS_AUTO = 16          # Number of initial homology models
NUM_MODELS_TO_REFINE = 2      # Top N models to refine
NUM_MODELS_LOOP = 2           # Models per loop refinement

# --- Flank Control ---
EXPERIMENTAL_FLANK_SIZE = 5   # Refinable edge residues (0 to disable)
REFINE_FLANKS_DURING_AUTOMODEL = False # Optimize flanks in AutoModel? (Default: False)

# --- Processing ---
# Auto-reads SLURM CPUs, or defaults to 1
NUM_PROCESSORS = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
```

### 4.3. Run the Pipeline

- Option A: **Locally**:

```bash
# Run as a module (recommended for correct package imports)
python3 -m prism.controller
```

- Option B: With **SLURM** (HPC Cluster)

```bash
# Adjust --cpus-per-task in run_prism.sh if needed
sbatch run_prism.sh
```

### 4.4. Check Results
All outputs are written to `modeling_results/`:

- **`AUTO_*.pdb`**: Initial homology models (ranked by quality)
- **`AUTO_*_LOOP*_R*.pdb`**: Loop-refined models with traceable nomenclature
- **`final_models_ranking.csv`**: Complete ranking with DOPE-HR scores

---

## 5. Configuration Reference

### 5.1. Core Parameters

| Parameter | Description | Recommended | Notes |
|-----------|-------------|---------|-------|
| `NUM_MODELS_AUTO` | Initial homology models to generate | 1000 - 10000 | Higher = better sampling, longer runtime |
| `NUM_MODELS_TO_REFINE` | Top models selected for loop refinement | 20 | Higher = better sampling, longer runtime |
| `NUM_MODELS_LOOP` | Models generated per loop refinement | 48 | Should match `NUM_PROCESSORS` for efficiency |
| `NUM_PROCESSORS` | Parallel workers | X | Set to CPU count for maximum speed |
| `NUM_BEST_FINAL_MODELS` | Models included in final ranking | 1000 | Limits CSV output size |

### 5.2. Experimental Region Control

| Parameter | Description | Default | Impact |
|-----------|-------------|---------|--------|
| `EXPERIMENTAL_FLANK_SIZE` | Edge residues to consider for refinement | 5 | 0 = no flank refinement, 10+ = aggressive |
| `MIN_DIST_FROM_NUCLEOTIDE_COM` | Minimum C-alpha to HETATM distance (Å) | 5.0 | Lower = tighter packing, higher = safer |

### 5.3. Alignment Mode

| Parameter | Description | Options |
|-----------|-------------|---------|
| `USE_MANUAL_ALIGNMENT` | Use manual or automatic alignment | `True` (manual) / `False` (auto salign) |

**Manual Mode**: Provide `input/manual_template_FullSeq.ali` in PIR format  
**Automatic Mode**: MODELLER generates alignment using `salign()` with multiple templates

### 5.4. PSIPRED Integration

| Parameter | Description | Default |
|-----------|-------------|---------|
| `PERFORM_PSIPRED_PREDICTION` | Query PSIPRED web server | `False` |
| `PSIPRED_EMAIL` | Email for PSIPRED submission | Required if `True` |
| `SS2_FILE_BASENAME` | Secondary structure filename | `P1_NCL_secondary_structure.ss2` |

---

## 6. Pipeline Workflow

```
┌─────────────────────────────────────────────────────────┐
│ STEP 0.5: PSIPRED Prediction (Optional)                 │
│  • Submits FASTA to PSIPRED server                │
│  • Polls for results and downloads .ss2 file       │
└─────────────────────────────────────────────────────────┘
                          ↓
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

## 7. Model Nomenclature

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

## 8. Output Files

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

## 9. Test Case Walkthrough

This repository is provided with the necessary files to run a test case (see input/). The default config.py is set up for this test:

- **PERFORM_PSIPRED_PREDICTION = False**: The pipeline will use the provided P1_NCL_secondary_structure.ss2 file to run immediately.
    - To test the client: Delete the .ss2 file from input/, set this to True, and add your email to PSIPRED_EMAIL. The pipeline will pause until PSIPRED is finished (can take hours).

- **USE_MANUAL_ALIGNMENT = True**: The pipeline will use the provided manual_template_FullSeq.ali, which has been manually vetted.
    - To test auto-alignment: Delete the .ali file from input/ and set this to False. The pipeline will generate a new .ali in modeling_results/. Check this file to see if it matches the manual one.

- **Replicas & RMSD**: The test case includes only 10 replicas of the main PDB (9CB5) to keep the repository lightweight.
    - Consequence: If you run tools/verify_rmsd.py on the results, you will see an RMSD > 0.0! This is expected.

- **Publication-Ready Solution**: For a real project, use tools/replicate_template.py to generate 1000 replicas (indicative). This will ensure a low RMSD on the deep core.

- **Template Selection**: All experimental PDBs for this protein were used, except 2KRR.
    - Reason: 2KRR almost completely overlaps with 9CB5 (the main PDB). Since 9CB5 is ligand-bound and may be slightly deformed, we do not want the apo 2KRR structure to distort its coordinates during the initial averaging step.

---

## 10. Citation

If you use PRISM in your research, please cite:

```
Lopez-Corbalan, R. (2025). PRISM: Protected-Region Insertion Suite for Modeling.
GitHub: https://github.com/richardloopez/PRISM
```

## 11. Technical Details

### 11.1. Architecture

#### 11.1.1. Core package (prism/)

This directory contains all the main pipeline logic, organized as a Python package.

- **`__init__.py`**: Defines the prism package, making its modules importable.
- **`controller.py`**: Main pipeline orchestrator
- **`psipred_client.py`**: Handles communication with the PSIPRED web server.
- **`config.py`**: Centralized configuration management
- **`utils.py`**: Alignment, secondary structure, and evaluation utilities
- **`fixed_region_utils.py`**: Experimental residue identification
- **`custom_models.py`**: Specialized FixedRegionAutoModel and FixedRegionLoopRefiner classes
- **`homology_modeling.py`**: Parallel homology modeling execution using FixedRegionAutoModel
- **`loop_refinement.py`**: Parallel loop refinement execution using FixedRegionLoopRefiner

#### 11.1.2. Preparation and analysis toolkit (tools/)

This directory contains standalone helper scripts that support the main pipeline.

- **`replicate_template.py`**: Create N copies of the main PDB and returns the update of config.py and .ali files to prevent the "template averaging" problem.
- **`verify_rmsd.py`**: Superimpose a final model onto the original template and verify that the experimental core's RMSD remains low.
- **`test_hetatm.py`**: Analyze PDB files for ATOM vs. HETATM content, helpful when preparing PDBs with ligands or DNA/RNA.
- **`renuamber_pdb.py`**: Renumber residues sequentially and convert specific chains (e.g., ligands) from ATOM to HETATM.
- **`extract_results.py`**: Gather final LOOP models and .csv reports from multiple parallel job directories into a single {_processed} folder.


### 11.2. Core Innovation: Fixed-Region Selection

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

### 11.3. HETATM Repulsion Implementation

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

### 11.4. Parallelism and Performance

PRISM is built on MODELLER's parallel job architecture to accelerate model generation. When **NUM_PROCESSORS** is set to a value greater than 1, the main controller spawns that number of **controller.worker** processes, each responsible for building a single structure. This parallel design has two important implications for performance tuning.

#### 11.4.1. NUM_PROCESSORS and HPC (SLURM) Configuration

The NUM_PROCESSORS parameter in config.py dictates how many parallel workers MODELLER will launch. The provided run_prism.sh script correctly reads the #SBATCH --cpus-per-task value and passes it to the environment, where config.py can read it.

However, on modern HPC systems that use hyperthreading (e.g., 2 threads per physical core), the --cpus-per-task value can be ambiguous—it may refer to physical cores or logical cores (threads).

   - **The Challenge**: A single MODELLER worker (controller.worker) is a single-threaded process. If 32 physical cores are requested (--cpus-per-task=32) on a 2-way hyperthreaded node (64 total threads) but set NUM_PROCESSORS = 32, 32 workers will be launched. The operating system will likely schedule these 32 workers on the 32 physical cores, leaving the other 32 "hyper-threads" idle and "wasting" half of the node's processing potential.

   - **The Recommendation**: For maximum resource utilization, NUM_PROCESSORS should be set to match the total number of logical cores (threads) available.

      - **Example**: SLURM script requests --cpus-per-task=32 on a node with 2-way hyperthreading.
      - SLURM Allocation: 32 Physical Cores (64 Logical Cores/Threads).
      - Optimal config.py setting: NUM_PROCESSORS = 64.
      - This will launch 64 workers, allowing the HPC's scheduler to utilize all 64 threads, maximizing throughput and ensuring you are not competing with other jobs on the node for those idle threads.

#### 11.4.2. NUM_MODELS_LOOP vs. NUM_MODELS_AUTO

The optimal strategy for setting the number of models differs between the initial homology modeling and the loop refinement phases.

   - **NUM_MODELS_AUTO** (Homology Modeling): This number is typically large (e.g., 1000+). Because the time to generate each .pdb can vary significantly (depending on the model, thread load, etc.), and the total number of jobs is high, a precise match between NUM_MODELS_AUTO and NUM_PROCESSORS is not critical for efficiency. A large number of models will naturally balance the load across the workers.

   - **NUM_MODELS_LOOP** (Loop Refinement): This phase is a significant performance bottleneck. The number of models is small, and each one takes a long time to generate. Crucially, the PRISM controller processes loops sequentially: it must wait for all NUM_MODELS_LOOP to finish for a given base model (e.g., AUTO_1) before it can begin the next loop refinement.

   - **The Recommendation**: For maximum efficiency, **NUM_MODELS_LOOP** should be set to an exact multiple (or equal to) **NUM_PROCESSORS**.
      - Justification: This setting prevents "idle worker" scenarios.
      - **Bad Example**: NUM_PROCESSORS = 32 and NUM_MODELS_LOOP = 33.
      - Result: The pipeline will run 32 models in parallel (e.g., 1 hour). It will then force 31 workers to sit idle for another hour while a single worker processes the 33rd model. This effectively doubles the runtime.
      - **Good Example**: NUM_PROCESSORS = 32 and NUM_MODELS_LOOP = 64.   
      - Result: The pipeline runs 32 models in parallel (1 hour), and then immediately runs the next 32 models in parallel (1 hour). The total time is 2 hours, but with zero worker idle time and maximum resource utilization.

---

## 12. Support

For questions, issues, or contributions, please contact:

- **Author**: Richard Lopez-Corbalan
- **GitHub**: https://github.com/richardloopez
- **Email**: richard.lopezc@uah.es

---

## 13. Acknowledgments

Built on the powerful [MODELLER](https://salilab.org/modeller/) comparative modeling framework developed by the Sali Lab at UCSF.

PSIPRED secondary structure predictions: Jones, D.T. (1999) J. Mol. Biol. 292: 195-202.

---

## 14. License

**Note**: PRISM requires MODELLER, which has its own licensing terms. Academic users can obtain free licenses from https://salilab.org/modeller/registration.html.

---