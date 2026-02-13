# 🛡️ PRISM: Protected-Region Insertion Suite for Modeling

**A pipeline for high-fidelity protein homology modeling and loop refinement with guaranteed preservation of experimental coordinates.**

Author: Richard Lopez-Corbalan  
GitHub: [github.com/richardloopez](https://github.com/richardloopez)

---

## 0. Overview

PRISM is a specialized protein modeling pipeline built on MODELLER that addresses a critical challenge in structural biology: **preventing coordinate drift of experimental regions during computational optimization**.

The pipeline enables high-fidelity modeling of loop insertions and flexible regions while maintaining the structural integrity of experimentally-determined core coordinates with 100% fidelity (0.0 Å RMSD).

### Key Features

- **🔒 Absolute Core Protection**: Experimental coordinates remain completely fixed during all optimization steps
- **🧬 HETATM Repulsion Shields**: Per-residue repulsion prevents loop clashes with bound molecules (DNA/RNA/ligands)
- **⚡ Parallel Processing**: Full support for multi-processor execution via MODELLER's parallel framework
- **🔮 Automatic PSIPRED Prediction**: Integrated client to automatically submit, poll, and download secondary structure predictions from the PSIPRED web server
- **📊 Traceable Nomenclature**: Models are systematically named (e.g., `AUTO_1_LOOP2_R1.pdb`) for complete traceability
- **🚀 HPC-Ready**: Nextflow-based orchestration with native SLURM support
- **🔬 Advanced Flank Control**: Granular user control over the behavior of experimental flanks (the junction between fixed and modeled regions)
- **🎯 Smart Loop Detection**: Combines PSIPRED secondary structure predictions with experimental boundary analysis
- **📈 Comprehensive Evaluation**: Automatic ranking using DOPE-HR scores with detailed CSV output
- **🖥️ Interactive Dashboard**: A Streamlit-based GUI (`PRISM_Dashboard.py`) for configuration, 3D visualization, and real-time monitoring

---

## 1. The Problem PRISM Solves

### 1.1. Experimental Core Drift

During standard homology modeling and loop refinement, optimization algorithms can inadvertently modify experimentally-determined coordinates. Even small deviations accumulate and compromise the reliability of the structural core. PRISM prevents this by:

- Identifying all residues that map to experimental coordinates in the template
- Completely excluding these residues from optimization (`select_atoms()`)
- Maintaining coordinate integrity across all refinement cycles

### 1.2. Loop-Ligand Clashes

When modeling loops in the presence of bound molecules (DNA, RNA, ligands), standard approaches often produce steric clashes. PRISM addresses this through:

- **Per-residue gravity center calculation**: Each HETATM residue gets a pseudo-atom at its center of mass
- **Lower-bound distance restraints**: C-alpha atoms in flexible regions maintain minimum separation from bound molecules
- **Configurable repulsion strength**: Adjustable minimum distance (`BLOCK_REPULSION_RADIUS`, default: 5.0 Å) balances accuracy with sampling

### 1.3. The Template Averaging Problem

MODELLER generates its initial structure by averaging the coordinates of all provided templates. This creates a critical problem:

- **The Problem**: When mixing an experimental PDB (e.g., residues 306-472) with a full-length AlphaFold model (e.g., 1-710), the 306-472 region will be averaged. Its coordinates will move (RMSD > 0), even if marked as "fixed" for later optimization.
- **The Solution**: PRISM employs a "template weighting" strategy via `PRISM_POWER_SETTINGS`. By utilizing numerous replicas of the main template and a two-phase execution paradigm (`prism-power`), the starting `.ini` and `.rsr` files are heavily biased toward the experimental coordinates.

### 1.4. Experimental Flank Refinement

The junction between the 0.0 RMSD experimental core and the new (e.g., AlphaFold) model can be highly strained. PRISM provides granular control over this junction (the "flank").

- **`MOBILE_FLANK_RESIDUES`**: Defines a buffer zone (e.g., 5 residues) at the edge of the experimental core.
- **`REFINE_FLANKS_DURING_AUTOMODEL`**: This boolean flag provides full control:
    - `false` (Default / Safe Mode): AutoModel (Stage 2) freezes the entire experimental region (core + flanks). LoopModel (Stage 4) then intelligently refines only flank residues predicted as 'C' (coil) by PSIPRED, providing a data-driven tension release.
    - `true` (Advanced / Aggressive Mode): AutoModel freezes only the "deep core" (core - flanks). Flanks are optimized immediately to resolve tension from the beginning. Useful for severe clashes, but risks destroying stable secondary structures in the flank.

---

## 2. Best Practices

### 2.1. PDB Preparation & Cleaning

The code is designed to receive **clean inputs**. Use the included `tools/prep_prism_pdb.py` utility (see [Section 10](#10-utility-tools-tools)) for automated preparation.

- **Main Template** (first entry in `PDB_TEMPLATE_FILES_NAMES`): This is the PDB whose coordinates (and optionally ligand) will be preserved.
    - The protein MUST be in Chain A.
    - The desired ligand/HETATM group should be placed in a separate chain (e.g., Chain B, controlled by `BLK_CHAIN_ID`).
- **Other Templates**:
    - The protein MUST be in Chain A.
    - They MUST NOT contain ligands.

### 2.2. The Averaging Problem & The Replication Trick

As mentioned in Section 1.3, MODELLER averages templates. To achieve a 0.0 RMSD on your experimental core, the `prism-power` execution paradigm should be used.

- **Use `EXECUTION_PARADIGM: prism-power`**:
    - Configure `PRISM_POWER_SETTINGS` in `config.yaml` with a `precalculation` phase (high replicas, e.g., 1000:1) and a `precomputed` phase (lower replicas, e.g., 100:1).
    - The `precalculation` phase generates the initial `.ini` and `.rsr` files with extreme bias toward the experimental PDB.
    - The `precomputed` phase then uses those files for the actual model generation with a more balanced replica count.
    - This ensures the experimental coordinates are not distorted by other templates (like AlphaFold models) during the initial AutoModel averaging.
    - If RMSD > 0.0 when using `tools/prism_verify_rmsd.py`, the cause is likely an insufficient number of replicas.

### 2.3. Alignment Verification

The code can run in a fully automatic mode (`USE_MANUAL_ALIGNMENT: false`) by generating an alignment with `salign()`.

- ⚠️ **WARNING!** The auto-generated alignment is not always accurate.
    - It is strongly recommended to **first** run the pipeline with `USE_MANUAL_ALIGNMENT: false` and `PERFORM_PSIPRED_PREDICTION: false`. The script will quickly generate an alignment file in `input/`.
    - **Manually inspect** this `.ali` file. If it's incorrect, edit it, save it as `input/manual_template_FullSeq.ali`, and switch to `USE_MANUAL_ALIGNMENT: true` in `config.yaml`.

---

## 3. Installation

### 3.1. Prerequisites

**Pixi** (environment manager) — this is the only prerequisite. Pixi handles everything else (Python, MODELLER, Nextflow, Streamlit, etc.) automatically, with zero root privileges required.

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### 3.2. Setup

```bash
# Clone PRISM
git clone https://github.com/richardloopez/PRISM.git
cd PRISM

# Apply MODELLER license and verify environment
pixi run setup
```

> **Note**: `pixi run setup` automatically applies the MODELLER academic license key (free for academic use) to the environment. Place your license key in `pixi.toml` under `[activation.env]` before running setup.

### 3.3. Directory Structure

```
PRISM/                           # Repository root
├── PRISM/                       # Core Python package
│   ├── config.py                # Centralized Pydantic configuration (loads config.yaml)
│   ├── controller.py            # Pipeline stage controller (dispatches stages)
│   ├── modeling_engine.py       # Custom MODELLER classes (FixedRegionAutoModel, FixedRegionLoopModel)
│   ├── psipred_client.py        # PSIPRED web API client
│   ├── utils.py                 # Alignment, secondary structure, ranking utilities
│   ├── ui_utils.py              # Streamlit helper functions (visualization, file management)
│   └── logo.png                 # PRISM branding logo
├── tools/                       # Utility scripts for PDB preparation & analysis
│   ├── prep_prism_pdb.py        # Prepares raw PDB files for PRISM (chain splitting, renaming)
│   ├── prism_verify_rmsd.py     # Verifies experimental coordinate fidelity (CA RMSD)
│   ├── calc_block_distance.py   # Calculates CA-to-BLK distances for BLOCK_REPULSION_RADIUS
│   └── unify_templates.py       # Resolves template overlaps and generates unified PDBs
├── input/                       # Input files (PDBs, FASTA, .ss2, .ali)
├── modeling_results/            # Output directory (models, rankings, logs)
├── test/                        # Self-contained test case (ready to run)
│   ├── config.yaml              # Pre-configured test settings
│   ├── input/                   # Test PDBs, FASTA, alignment
│   └── ...
├── config.yaml                  # Main pipeline configuration
├── pixi.toml                    # Environment & dependency manager (Pixi)
├── orchestrator.nf              # Nextflow DSL2 workflow (SLURM orchestration)
├── nextflow.config              # Nextflow resource allocation (CPUs, memory)
├── PRISM_Dashboard.py           # Streamlit GUI (interactive dashboard)
└── README.md                    # This file
```

---

## 4. Quick Start

### 4.1. Prepare Input Files

Place the following files in the `input/` directory:

| File | Required | Description |
|------|----------|-------------|
| `sequence_full.fasta` | ✅ Yes | Target protein sequence in FASTA format |
| Template PDB files (`.pdb`) | ✅ Yes | One or more experimental structures as templates |
| `P1_NCL_secondary_structure.ss2` | Conditional | PSIPRED prediction. Required if `PERFORM_PSIPRED_PREDICTION: false` |
| `manual_template_FullSeq.ali` | Conditional | Manual alignment in PIR format. Required if `USE_MANUAL_ALIGNMENT: true` |
| `precomputed_ini.pdb` | Optional | Precomputed average structure (for `precomputed` paradigm) |
| `precomputed_rsr.rsr` | Optional | Precomputed restraint file (for `precomputed` paradigm) |

### 4.2. Configure Parameters

Edit `config.yaml` to match your project. See [Section 7](#7-configuration-reference) for a full parameter reference.

```yaml
# --- Template Selection ---
PDB_TEMPLATE_FILES_NAMES:
  - 9CB5_renum_HETATM.pdb     # Main template (first = experimental core)
  - AF_renum.pdb               # AlphaFold or additional template

# --- Modeling Scale ---
MODELLER_CORES: 10              # CPU threads for MODELLER
TOTAL_PARALLEL_JOBS: 5          # Nextflow parallel job splits
TOTAL_HOMOLOGY_MODELS: 10000    # Total AutoModel structures to generate
TOP_MODELS_FOR_REFINEMENT: 20   # Top N models selected for loop refinement
LOOP_MODELS_PER_TARGET: 10      # Loop models generated per top model

# --- Execution Mode ---
EXECUTION_PARADIGM: prism-power # Recommended: 'prism-power', 'precalculation', 'precomputed', or 'normal'

# --- PRISM Power Settings (only for 'prism-power' paradigm) ---
PRISM_POWER_SETTINGS:
  precalculation:               # Phase 1: Generate .ini/.rsr with heavy experimental bias
    9CB5_renum_HETATM.pdb: 1000
    AF_renum.pdb: 1
  precomputed:                  # Phase 2: Use .ini/.rsr for actual model generation
    9CB5_renum_HETATM.pdb: 100
    AF_renum.pdb: 1
```

### 4.3. Run the Pipeline

**Option A: Interactive Dashboard**
```bash
pixi run gui-prism          # Standard launch (attached to terminal)
pixi run gui-prism-persist  # Persistent launch (safe to close terminal)
```
This launches the Streamlit Dashboard where you can:
- Edit all configuration parameters in a visual interface
- Upload/delete input files
- Launch and monitor the pipeline in real-time
- Visualize results in interactive 3D

**Option B: Terminal / HPC Cluster**
```bash
pixi run terminal-prism          # Standard launch
pixi run terminal-prism-persist  # Persistent launch
```
This runs the full Nextflow pipeline directly. On a SLURM cluster, Nextflow automatically submits each stage as a separate job via the configuration in `nextflow.config`.


### 4.4. Check Results

All outputs are written to `modeling_results/`:

| File | Description |
|------|-------------|
| `AUTO_*.pdb` | Initial homology models (ranked by DOPE-HR quality) |
| `AUTO_*_LOOP*_R*.pdb` | Loop-refined models with traceable nomenclature |
| `final_ranking.csv` | Complete ranking with DOPE-HR scores and Z-scores |
| `modeling_results/logs/` | Per-stage log files for debugging |

---

## 5. The Streamlit Dashboard

The `PRISM_Dashboard.py` provides a full-featured graphical interface:

| Tab | Features |
|-----|----------|
| **⚙️ Config** | Edit all `config.yaml` parameters. Includes sections for Modeling & Refinement, Templates, PSIPRED, Execution Paradigm, and Manual Overrides. Save/reload configuration in real-time. |
| **📁 Files** | Upload PDB templates, FASTA sequences, alignment files, and any additional files. View a file inventory with one-click deletion. |
| **🔬 Visualization** | Interactive 3D visualization of PDB files using py3Dmol. Multiple rendering styles: ribbon, VDW, surface, line, spacefill. |
| **⚡ Execution** | Launch the Nextflow pipeline and monitor real-time progress with automatic log updates. |
| **📊 Results** | View ranked models, DOPE-HR score distributions, Z-score plots, and download the ranking CSV. |
| **🔧 Tools** | Run `prep_prism_pdb.py`, `prism_verify_rmsd.py`, `calc_block_distance.py`, and `unify_templates.py` directly from the GUI. |

**Remote Access (HPC):**  
If running on a remote cluster, use SSH port forwarding to access the dashboard locally:
```bash
ssh -L 8501:localhost:8501 your_hpc_address
```
Then open `http://localhost:8501` in your browser.

---

## 6. Pipeline Workflow

```
┌──────────────────────────────────────────────────────────────────┐
│ STAGE 0.5: PSIPRED Prediction (Optional)                         │
│  • Submits FASTA to UCL PSIPRED API via psipred_client.py        │
│  • Polls for completion and downloads .ss2 file                  │
│  • Skipped if PERFORM_PSIPRED_PREDICTION: false                  │
└──────────────────────────────────────────────────────────────────┘
                             ↓
┌──────────────────────────────────────────────────────────────────┐
│ STAGE 1: Prereq-CDE (Alignment & Secondary Structure)            │
│  • Manual Mode: Maps .ss2 to provided PIR alignment              │
│  • Auto Mode: Performs MODELLER salign() across templates        │
│  • Generates _cde.ali with secondary structure annotation (CDE)  │
└──────────────────────────────────────────────────────────────────┘
                             ↓
┌──────────────────────────────────────────────────────────────────┐
│ STAGE 2: AutoModel (Parallelized via Nextflow)                   │
│  • Nextflow splits TOTAL_HOMOLOGY_MODELS into TOTAL_PARALLEL_JOBS│
│  • FixedRegionAutoModel freezes experimental core residues       │
│  • Applies HETATM repulsion shields to CA atoms                 │
│  • prism-power: Runs precalculation phase first, then precomp.  │
└──────────────────────────────────────────────────────────────────┘
                             ↓
┌──────────────────────────────────────────────────────────────────┐
│ STAGE 3: Rank-AutoModel                                          │
│  • Assesses initial models using DOPE-HR (assess_dopehr())       │
│  • Renames top candidates to AUTO_1.pdb, AUTO_2.pdb, etc.       │
│  • Selects TOP_MODELS_FOR_REFINEMENT for the next stage          │
└──────────────────────────────────────────────────────────────────┘
                             ↓
┌──────────────────────────────────────────────────────────────────┐
│ STAGE 4: Loop Refinement (Parallelized via Nextflow)             │
│  • Refines detected coil regions on best AutoModels              │
│  • FixedRegionLoopModel maintains core/template stability        │
│  • Sequential refinement per loop per model                      │
│  • Nomenclature: AUTO_1_LOOP1_R1.pdb, AUTO_1_LOOP2_R1.pdb, etc. │
└──────────────────────────────────────────────────────────────────┘
                             ↓
┌──────────────────────────────────────────────────────────────────┐
│ STAGE 5: Final Ranking                                           │
│  • Final DOPE-HR assessment of all Auto and Loop models          │
│  • Calculates Normalized DOPE-HR Z-scores                        │
│  • Exports results to final_ranking.csv                          │
└──────────────────────────────────────────────────────────────────┘
```

---

## 7. Configuration Reference

All parameters are set in `config.yaml` (or via the Dashboard).

### 7.1. Modeling Parameters

| Parameter | Description | Recommended | Notes |
|-----------|-------------|-------------|-------|
| `TOTAL_HOMOLOGY_MODELS` | Initial homology models to generate | 1,000 – 10,000 | Higher = better sampling, longer runtime |
| `TOP_MODELS_FOR_REFINEMENT` | Top models selected for loop refinement | 20 | Based on DOPE-HR ranking |
| `LOOP_MODELS_PER_TARGET` | Models generated per loop refinement | 10–48 | Should be a multiple of `MODELLER_CORES` for efficiency |
| `NUM_BEST_FINAL_MODELS` | Models included in final ranking | `inf` | Set to `inf` to rank all models |

### 7.2. Execution & Parallelization

| Parameter | Description | Notes |
|-----------|-------------|-------|
| `MODELLER_CORES` | MODELLER parallel workers per Nextflow task | Set to match available CPU threads |
| `TOTAL_PARALLEL_JOBS` | Number of Nextflow job splits for Stage 2 | Determines how many concurrent SLURM tasks |
| `EXECUTION_PARADIGM` | Workflow mode | `prism-power` (recommended), `precalculation`, `precomputed`, or `normal` |

### 7.3. Experimental Region Control

| Parameter | Description | Default | Impact |
|-----------|-------------|---------|--------|
| `MOBILE_FLANK_RESIDUES` | Edge residues considered for refinement | `5` | `0` = no flank refinement |
| `REFINE_FLANKS_DURING_AUTOMODEL` | Optimize flanks in AutoModel? | `false` | `true` = aggressive (see Section 1.4) |
| `BLOCK_REPULSION_RADIUS` | Min C-alpha to HETATM distance (Å) | `5.0` | Lower = tighter packing |
| `USE_MANUAL_OPTIMIZATION_SELECTION` | Manually specify residues to optimize | `false` | Overrides automatic detection |
| `MANUAL_OPTIMIZATION_RESIDUES` | List of residue indices to optimize | `[]` | Only used when above is `true` |

### 7.4. Chain & File Configuration

| Parameter | Description | Default |
|-----------|-------------|---------|
| `CHAIN_ID` | Protein chain identifier | `A` |
| `BLK_CHAIN_ID` | Ligand/HETATM chain identifier | `B` |
| `FASTA_FILE_BASENAME` | Target sequence file | `sequence_full.fasta` |
| `SS2_FILE_BASENAME` | PSIPRED prediction file | `P1_NCL_secondary_structure.ss2` |
| `MANUAL_ALIGNMENT_BASENAME` | Manual alignment file | `manual_template_FullSeq.ali` |
| `CUSTOM_INIFILE_BASENAME` | Precomputed initial structure | `precomputed_ini.pdb` |
| `CUSTOM_RSRFILE_BASENAME` | Precomputed restraint file | `precomputed_rsr.rsr` |

### 7.5. PSIPRED Integration

| Parameter | Description | Default |
|-----------|-------------|---------|
| `PERFORM_PSIPRED_PREDICTION` | Query PSIPRED web server | `false` |
| `PSIPRED_EMAIL` | Email for PSIPRED submission | Required if `true` |
| `PSIPRED_POLL_INTERVAL` | Seconds between polling attempts | `60` |

### 7.6. Template Configuration

| Parameter | Description |
|-----------|-------------|
| `PDB_TEMPLATE_FILES_NAMES` | Ordered list of PDB template filenames. The **first entry** is always the main experimental template. |
| `PRISM_POWER_SETTINGS` | Per-template replica weights for the `prism-power` paradigm. Contains a `precalculation` and a `precomputed` dictionary mapping template filenames to replica counts. |

---

## 8. Execution Paradigms

PRISM supports four execution paradigms, controlled by `EXECUTION_PARADIGM`:

| Paradigm | Description | When to Use |
|----------|-------------|-------------|
| `normal` | Standard MODELLER run. Templates used as-is from `PDB_TEMPLATE_FILES_NAMES`. | Simple homology modeling without coordinate protection needs. |
| `precalculation` | Runs Stage 2 only to generate `.ini` and `.rsr` files with heavy template weighting. Stops before model generation. | When you need to generate initial files for later `precomputed` runs. |
| `precomputed` | Uses pre-existing `.ini` and `.rsr` files (from `input/`) for Stage 2. Skips initial averaging. | When `.ini`/`.rsr` files have already been generated (e.g., by a previous `precalculation` run). |
| `prism-power` ⭐ | **Two-phase automatic execution.** First runs a `precalculation` phase to generate `.ini`/`.rsr` with extreme experimental bias, then automatically runs `precomputed` model generation with balanced replicas. | **Recommended for all production runs.** Provides the best coordinate fidelity without manual intervention. |

---

## 9. Model Nomenclature

PRISM uses systematic, traceable naming:

| Filename Pattern | Meaning | Example |
|-----------------|---------|---------|
| `AUTO_N.pdb` | Initial homology model (rank N by DOPE-HR) | `AUTO_1.pdb` |
| `AUTO_N_LOOPJ_RK.pdb` | Loop-refined model: base N, loop J, refinement K | `AUTO_1_LOOP2_R1.pdb` |

**Reading `AUTO_1_LOOP2_R3.pdb`:**
- Started from `AUTO_1.pdb` (best initial model)
- Refined loop #2 sequentially
- 3rd ranked refinement for that loop

**Best Model Strategy:** The final ranking includes all intermediate models, allowing you to select the globally best model regardless of refinement stage.

### Ranking CSV Format

```csv
Rank,Model Name,DOPEHR Score,DOPEHR Z-score
1,AUTO_1_LOOP3_R1.pdb,-52847.234,-1.234
2,AUTO_2_LOOP3_R1.pdb,-52523.456,-1.156
...
```

---

## 10. Utility Tools (`tools/`)

### 10.1. `prep_prism_pdb.py` — PDB Preparation

Prepares raw PDB files for PRISM by splitting chains, renumbering atoms, and generating a preparation log.

**Mode: `prep` (Preparation)**
```bash
python3 tools/prep_prism_pdb.py prep \
    --input-pdb raw_structure.pdb \
    --protein-chains A \
    --ligand-chains B
```
- Splits the PDB into protein (Chain A) and ligand (Chain B) sections
- Renumbers atoms sequentially
- Generates `*_renum_HETATM.pdb` (main template) and a JSON log

**Mode: `retro` (Restoration)**
```bash
python3 tools/prep_prism_pdb.py retro \
    --prism-output-pdb modeling_results/AUTO_1.pdb \
    --original-pdb raw_structure.pdb \
    --log-file prep_log.json
```
- Restores original ligand atoms and naming to a PRISM output model
- Uses the JSON log from the `prep` step to reverse-map atom names

### 10.2. `prism_verify_rmsd.py` — Coordinate Fidelity Verification

Verifies that experimental coordinates have been preserved (0.0 Å RMSD) in a PRISM output model.

```bash
python3 tools/prism_verify_rmsd.py \
    input/9CB5_renum_HETATM.pdb \
    modeling_results/AUTO_1.pdb \
    input/manual_template_FullSeq.ali
```

**Optional flags:**
- `--manual 1-191:2-192`: Manually specify residue mapping instead of using the alignment file
- `--orig-chain A`: Chain in the original PDB (default: `A`)
- `--mod-chain A`: Chain in the model PDB (default: `A`)

**Output:**
```
RESIDUE         | ORIG POS   | MOD POS    | RMSD (Å)
--------------------------------------------------------------------------------
✓ ALA           | 306        | 306        | 0.000000
✓ GLY           | 307        | 307        | 0.000000
...
SUMMARY for 167 residues:
 > Average RMSD: 0.000000 Å
 > Maximum RMSD: 0.000000 Å

✅ SUCCESS: Experimental coordinates are FIXED.
```

### 10.3. `calc_block_distance.py` — BLK Distance Calculator

Calculates the minimum distance between each protein CA atom and the gravity centers of BLK/HETATM residue groups. This tool mirrors the internal physics of PRISM's `add_hetatm_repulsion_shield` and helps determine an appropriate `BLOCK_REPULSION_RADIUS`.

```bash
python3 tools/calc_block_distance.py input/9CB5_renum_HETATM.pdb \
    --protein_chain A \
    --blk_chain B \
    --threshold 10.0
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `pdb_file` | *(required)* | Path to the PDB file to analyze |
| `--protein_chain` | `A` | Chain ID for the protein |
| `--blk_chain` | `B` | Chain ID for the BLK/HETATM residues |
| `--threshold` | `10.0` | Only display residues closer than this distance (Å) |

**Output:**
```
Residue         | Min Distance (A)   | Closest BLK Group
--------------------------------------------------------------------------------
PHE 1           | 4.523              | BLK:28:B
ASN 2           | 5.112              | BLK:28:B
...
GLOBAL MINIMUM DISTANCE: 4.523 A
SUGGESTED BLOCK_REPULSION_RADIUS: 4.5 A
```

> [!TIP]
> If the global minimum distance is less than 2.0 Å, the tool will warn you — this could indicate a clash that needs manual inspection before modeling.

### 10.4. `unify_templates.py` — Template Overlap Resolution

Resolves sequence overlaps between multiple templates in a MODELLER alignment. When templates cover the same regions of the target protein, MODELLER averages their coordinates, which can **destroy experimental accuracy**. This tool trims redundant residues from lower-priority templates while maintaining a configurable overlap buffer for structural continuity.

```bash
python3 tools/unify_templates.py input/manual_template_FullSeq.ali --overlap 10
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `alignment` | *(required)* | Path to the Modeller alignment file (`.ali`) |
| `--overlap` | `10` | Number of overlap residues to keep at junctions for structural continuity |

**How it works:**

1. Templates appearing **earlier** in the `.ali` file have higher priority
2. When a later template overlaps with a higher-priority one, its overlapping protein residues are trimmed (replaced with `-` gaps)
3. A buffer of `--overlap` residues is preserved around all junction points between covered and uncovered regions
4. BLK/HETATM residues are **never trimmed** — they are always preserved entirely
5. Residues are renumbered sequentially: protein chains start at 1, BLK residues continue from the last protein residue number
6. `TER` records are inserted at chain ends in the output PDB files

**Output files:**

| File | Description |
|------|-------------|
| `*_unified.ali` | Updated alignment file with trimmed sequences and corrected PDB references |
| `*_unified.pdb` | One per template — renumbered PDB with only the kept residues |

**Example:** If Template A (priority 1) covers residues 1-100 and Template B (priority 2) covers 80-200 with `--overlap 10`:
- Template A is kept as-is (all 100 residues)
- Template B is trimmed to keep residues ~91-200 (10-residue buffer into the overlap)
- Template B's PDB is renumbered starting from 1, with BLK continuing sequentially
- The alignment is updated so Template B has gaps at positions 80-90

> [!IMPORTANT]
> After running this tool, you must **rename the generated `*_unified.pdb` files** (or update the references in the alignment file and `config.yaml`) before running the PRISM pipeline.

---

## 11. HPC Configuration (`nextflow.config`)

The `nextflow.config` file controls SLURM resource allocation:

```groovy
process {
    executor = 'slurm'
    queue    = 'all'    // Your cluster partition name

    withName: 'AUTOMODEL|AUTOMODEL_PRECALC' {
        cpus = config_data.MODELLER_CORES  // Reads from config.yaml
        memory = '80 GB'
    }

    withName: 'LOOP_MODEL' {
        cpus = config_data.MODELLER_CORES
        memory = '80 GB'
    }
    // PREREQ_CDE, RANK_AUTOMODEL, FINAL_RANKING: 1 CPU, 8 GB
}
```

### Performance Tuning

- **`MODELLER_CORES`**: Set this to the total number of logical cores (threads) available on your node, not physical cores. On a 32-core node with 2-way hyperthreading, use `MODELLER_CORES: 64`.
- **`LOOP_MODELS_PER_TARGET`**: For maximum efficiency, set to an exact **multiple** of `MODELLER_CORES`. If `MODELLER_CORES: 32` and `LOOP_MODELS_PER_TARGET: 33`, 31 workers will idle while 1 finishes the 33rd model.
- **`TOTAL_PARALLEL_JOBS`**: Controls how many Stage 2 tasks run *concurrently* as separate SLURM jobs.

---

## 12. Test Case

The `test/` directory contains a self-contained, ready-to-run demo:

```bash
# Copy test files to root (or run from the test/ directory)
cp test/config.yaml config.yaml
cp test/input/* input/

# Run the demo (from terminal)
pixi run terminal-prism

# Run the demo (from Streamlit)
pixi run streamlit-prism
```

**What's included:**
- 4 template PDBs (9CB5 experimental + 2fc8 + 2fc9 + AlphaFold)
- Pre-vetted manual alignment (`manual_template_FullSeq.ali`)
- Target FASTA sequence
- Pre-configured `config.yaml` with `prism-power` paradigm

**Note:** The test case uses a low replica count for speed. Running `tools/prism_verify_rmsd.py` on the results may show RMSD > 0.0 — this is expected with few replicas.

---

## 13. Persistent Sessions & Multi-Working

PRISM includes built-in automation for persistent execution. This is critical for long-running protein modeling jobs on HPC nodes where terminal disconnections are common.

### 13.1. "Fire and Forget" Mode
To launch the pipeline in the background so it survives terminal closure:

- **GUI Mode**: `pixi run gui-prism-persist`
- **CLI Mode**: `pixi run terminal-prism-persist`

### 13.2. Monitoring Progress (Attach)
You can "re-attach" to see the live logs at any time:

- **GUI Logs**: `pixi run attach-gui`
- **CLI Logs**: `pixi run attach-terminal`

### 13.3. Stopping a Session
To safely terminate a background session:

- `pixi run stop-gui`
- `pixi run stop-terminal`

### 13.4. Running Multiple Instances
You can run multiple independent PRISM instances concurrently:
1.  **Isolated Folders**: Use separate folders for each experiment. Each folder will maintain its own hidden `.prism-*.pid` and `.prism-*.log` files.
2.  **Independent Results**: Each run keeps its own `modeling_results/`, `config.yaml`, and Nextflow `work/` directory.


---

## 14. Portability

PRISM is fully portable without requiring root privileges:


- **Pixi** installs everything (Python, MODELLER, Nextflow, Streamlit) in a local environment
- **Supported platforms**: Linux (all distributions), macOS (Intel & Apple Silicon), Windows via WSL2
- **Zero system dependencies**: Copy the repository to a new machine, install Pixi, run `pixi run setup` — done

---

## 14. Technical Details

### 14.1. Core Innovation: Fixed-Region Selection

```python
def select_atoms(self):
    '''
    Select atoms for optimization, EXCLUDING experimental residues.
    '''
    if not self.experimental_residues:
        return Selection(self).only_std_residues()
    all_std = Selection(self).only_std_residues()
    fixed_sel = Selection()

    for res_num in sorted(self.experimental_residues):
        try:
            curr_res = self.residue_range(f'{res_num}:{self.chain_id}', f'{res_num}:{self.chain_id}')
            fixed_sel.add(curr_res)
        except Exception as e:
            continue

    optimizable = all_std - fixed_sel
    return optimizable
```

### 14.2. HETATM Repulsion Implementation

```python
def add_hetatm_repulsion_shield(model, min_dist: float, only_loop_atoms: bool = False):
    '''
    Adds lower-bound distance restraints between CA atoms and HETATM centers.
    '''
    het_residues = [r for r in model.residues if r.hetatm and r.name != 'HOH']
    target_ca = (model.select_loop_atoms() if only_loop_atoms else model.select_atoms()).only_atom_types('CA')

    het_centers = []
    for res in het_residues:
        center = pseudo_atom.GravityCenter(Selection(res))
        model.restraints.pseudo_atoms.append(center)
        het_centers.append(center)

    for ca in target_ca:
        for center in het_centers:
            model.restraints.add(forms.LowerBound(
                group=physical.xy_distance,
                feature=features.Distance(ca, center),
                mean=min_dist, stdev=1.0
            ))
```

---

## 15. Citation

If you use PRISM in your research, please cite:

```
Lopez-Corbalan, R. (2025). PRISM: Protected-Region Insertion Suite for Modeling.
GitHub: https://github.com/richardloopez/PRISM
```

---

## 16. Support

For questions, issues, or contributions:

- **Author**: Richard Lopez-Corbalan
- **GitHub**: [github.com/richardloopez](https://github.com/richardloopez)
- **Email**: richard.lopezc@uah.es

---

## 17. Acknowledgments

Built on the powerful [MODELLER](https://salilab.org/modeller/) comparative modeling framework developed by the Sali Lab at UCSF.

PSIPRED secondary structure predictions: Jones, D.T. (1999) *J. Mol. Biol.* 292: 195-202.

---

## 18. License

**Note**: PRISM requires MODELLER, which has its own licensing terms. Academic users can obtain free licenses from https://salilab.org/modeller/registration.html.
