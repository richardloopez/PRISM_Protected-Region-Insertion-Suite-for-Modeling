#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PRISM Weighted Pre-calculator (.ini y .rsr)

Este script genera DOS archivos esenciales para el pipeline de PRISM:
1.  Un 'weighted_average.pdb' (.ini) ponderado.
2.  Un 'weighted_restraints.rsr' (.rsr) ponderado.

Ambos se basan en la lista de pesos (WEIGHTS_LIST) definida abajo.

¿Cómo funciona?
1.  Crea un alineamiento masivo en memoria donde las plantillas
    se repiten según sus pesos (ej. 1000 copias de 9CB5, 500 de 2fc9...).
2.  Llama a `mdl.transfer_xyz()` con este alineamiento masivo,
    lo que (según el manual, 6.6.26) crea automáticamente
    el promedio ponderado para el .ini.
3.  Llama a `mdl.restraints.make()` con el mismo alineamiento masivo,
    lo que crea restricciones MultiGaussian (Manual, A.3.2) donde
    el peso de cada plantilla es proporcional a su número de copias.
4.  Guarda ambos archivos.

El pipeline principal (`homology_modeling.py`) puede entonces cargar
estos dos archivos usando 'inifile' y 'csrfile', saltándose
todo el pre-cálculo y evitando el "cuelgue" en los jobs paralelos.
"""

import sys
import os
from typing import Set

try:
    from modeller import Environ, Alignment, Model
    from modeller.automodel import assess
    from modeller.selection import Selection
except ImportError:
    print("FATAL ERROR: MODELLER library not found.")
    sys.exit(1)

# --- Configuración del Proyecto ---
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)
sys.path.append(project_root)

try:
    from prism import config
    from prism.config import (
        PDB_TEMPLATE_FILES_NAMES, ALIGN_CODES_TEMPLATES, ALIGN_CODE_SEQUENCE,
        USE_MANUAL_ALIGNMENT, MANUAL_ALIGNMENT_FILE, ALIGNMENT_FILE,
        INPUT_DIR, MODELING_RESULTS_DIR, CUSTOM_INIFILE_PATH, CUSTOM_RSRFILE_PATH,
        MAIN_ALIGN_CODE_TEMPLATE, REFINE_FLANKS_DURING_AUTOMODEL,
        EXPERIMENTAL_FLANK_SIZE, CHAIN_ID
    )
    from prism.custom_models import FixedRegionAutoModel
    from prism.fixed_region_utils import identify_experimental_residues
    from prism.utils import read_sequences_from_ali_temp, group_ranges
    print(f"Loading configuration from: {config.PROJECT_ROOT}/prism/config.py\n")
except ImportError as e:
    print("FATAL ERROR: Could not import 'prism.config' o módulos.")
    print(f"  -> {e}")
    sys.exit(1)

# --- 1. Pesos Definidos por el Usuario ---
# ¡IMPORTANTE! Esta lista DEBE estar en el mismo orden que
# PDB_TEMPLATE_FILES_NAMES en config.py
WEIGHTS_LIST = [
    10.0,  # <-- '9CB5_DS_renum_HETATM-B.pdb'
    5.0,   # <-- "2fc9_DS_renum.pdb"
    5.0,   # <-- "2fc8_DS_renum.pdb"
    1.0      # <-- "AF-P19338-F1-model_v6_DS_renum.pdb"
]

# --- 2. Rutas de Archivos ---
TARGET_CODE = config.ALIGN_CODE_SEQUENCE
TEMPLATE_CODES = list(config.ALIGN_CODES_TEMPLATES)

if USE_MANUAL_ALIGNMENT:
    ALIGNMENT_FILE_PATH = MANUAL_ALIGNMENT_FILE
    print(f"Using MANUAL alignment: {ALIGNMENT_FILE_PATH}")
else:
    ALIGNMENT_FILE_PATH = ALIGNMENT_FILE
    print(f"Using AUTO alignment: {ALIGNMENT_FILE_PATH}")

OUTPUT_INI_PATH = CUSTOM_INIFILE_PATH
OUTPUT_RSR_PATH = CUSTOM_RSRFILE_PATH

os.makedirs(INPUT_DIR, exist_ok=True)
os.makedirs(MODELING_RESULTS_DIR, exist_ok=True)


def get_residues_to_freeze(aln_map) -> Set[int]:
    """Obtiene los residuos a congelar (lógica de controller.py)."""
    print("[INFO] Calculando residuos a congelar...")
    main_template_seq = aln_map[MAIN_ALIGN_CODE_TEMPLATE]
    target_seq = aln_map[TARGET_CODE]
    
    experimental_residues = identify_experimental_residues(main_template_seq, target_seq)
    experimental_ranges = group_ranges(list(experimental_residues))
    flank_residues = set()
    N = EXPERIMENTAL_FLANK_SIZE
    for start, end in experimental_ranges:
        for i in range(N):
            if start + i <= end: flank_residues.add(start + i)
        for i in range(N):
            if end - i >= start: flank_residues.add(end - i)

    truly_fixed_residues = experimental_residues - flank_residues

    if REFINE_FLANKS_DURING_AUTOMODEL:
        fixed_residues = truly_fixed_residues
    else:
        fixed_residues = experimental_residues
    
    print(f"[INFO] Se congelarán {len(fixed_residues)} residuos (modo: {REFINE_FLANKS_DURING_AUTOMODEL}).")
    return fixed_residues


def generate_weighted_precalc():
    """Genera el .ini y .rsr ponderados."""
    
    # --- 3. Validación ---
    if len(TEMPLATE_CODES) != len(WEIGHTS_LIST):
        print("FATAL ERROR: Mismatch in config/script.")
        print("Number of templates does not match number of weights.")
        sys.exit(1)

    print("\nPesos a aplicar (Plantillas:Copias):")
    total_copies = 0
    for code, weight in zip(TEMPLATE_CODES, WEIGHTS_LIST):
        print(f"  - {code}: {weight}")
        total_copies += weight
    print(f"Total 'knowns' en alineamiento en memoria: {total_copies}")

    if not os.path.exists(ALIGNMENT_FILE_PATH):
        print(f"FATAL ERROR: Alignment file not found at: {ALIGNMENT_FILE_PATH}")
        sys.exit(1)

    # --- 4. Setup MODELLER ---
    print("\n[PASO 1] Inicializando MODELLER...")
    env = Environ()
    env.io.atom_files_directory = [config.INPUT_DIR, '.']
    env.io.hetatm = True
    env.io.water = True
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    # --- 5. Cargar Alineamiento Base ---
    print(f"\n[PASO 2] Cargando alineamiento base: {os.path.basename(ALIGNMENT_FILE_PATH)}")
    aln_master_map = {}
    aln_reader = Alignment(env)
    all_codes_to_load = TEMPLATE_CODES + [TARGET_CODE]
    aln_reader.append(file=ALIGNMENT_FILE_PATH, align_codes=all_codes_to_load)
    for seq in aln_reader:
        aln_master_map[seq.code] = seq
    
    # --- 6. Crear Alineamiento Ponderado EN MEMORIA ---
    print("\n[PASO 3] Creando alineamiento ponderado en memoria...")
    aln_weighted = Alignment(env)
    
    # Helper function para reconstruir el string de secuencia (con gaps)
    def get_gapped_string(aln_obj: Alignment, seq_obj: Model) -> str:
        seq_str = []
        for pos in aln_obj.positions:
            res = pos.get_residue(seq_obj)
            if res:
                seq_str.append(res.code)
            else:
                seq_str.append('-')
        return "".join(seq_str)

    base_seq_strings = {}
    for code in TEMPLATE_CODES:
        base_seq_strings[code] = get_gapped_string(aln_reader, aln_master_map[code])
    target_seq_string = get_gapped_string(aln_reader, aln_master_map[TARGET_CODE])
    
    weighted_knowns_codes = [] # Lista de todos los códigos (ej. 2001)
    
    for code, weight in zip(TEMPLATE_CODES, WEIGHTS_LIST):
        base_seq_data = aln_master_map[code]
        base_seq_str = base_seq_strings[code]
            
        for i in range(int(weight)):
            aln_weighted.append_sequence(base_seq_str)
            new_seq = aln_weighted[-1]
            new_seq_code = f"{code}_copy{i+1}" # Código único
            new_seq.code = new_seq_code
            new_seq.atom_file = base_seq_data.atom_file
            new_seq.range = base_seq_data.range
            new_seq.prottyp = base_seq_data.prottyp
            weighted_knowns_codes.append(new_seq_code) # Añadir a la lista
            
    # Añadir el target UNA VEZ al final
    aln_weighted.append_sequence(target_seq_string)
    target_seq_weighted = aln_weighted[-1]
    target_seq_weighted.code = TARGET_CODE
    target_seq_weighted.prottyp = 'sequence'
    
    print(f"  ✓ Alineamiento en memoria creado con {len(aln_weighted)} secuencias.")
    
    # --- 7. Generar .ini y .rsr PONDERADOS usando AutoModel ---
    print(f"\n[PASO 4] Instanciando FixedRegionAutoModel con alineamiento ponderado...")
    
    # Obtenemos los residuos a congelar del alineamiento ORIGINAL
    aln_map_for_freeze = read_sequences_from_ali_temp(ALIGNMENT_FILE_PATH)
    residues_to_freeze = get_residues_to_freeze(aln_map_for_freeze)

    a = FixedRegionAutoModel(env,
                             experimental_residues=residues_to_freeze,
                             chain_id=CHAIN_ID,
                             # --- ¡AQUÍ ESTÁ LA MAGIA! ---
                             alnfile=aln_weighted,             # Pasa el objeto Alignment
                             knowns=weighted_knowns_codes,     # Pasa los 2001 códigos
                             sequence=TARGET_CODE,
                             assess_methods=(assess.DOPEHR))
    
    # Fuerza la cadena en blanco para que coincida con el .ali manual
    # (Manual 4.1.20)
    a.blank_single_chain = True
    
    print("\n[PASO 5] Ejecutando a.make(exit_stage=1) para pre-calcular .ini y .rsr...")
    # Manual 4.1.27: exit_stage=1 genera .ini y .rsr, y para.
    a.make(exit_stage=1)
    print("  ✓ Cálculo de .ini y .rsr completado.")

    # --- 8. Mover y Renombrar Archivos ---
    print(f"\n[PASO 6] Moviendo y renombrando archivos...")
    
    # AutoModel guarda esto en CWD con el nombre de la *secuencia*
    generated_ini_file = f"{TARGET_CODE}.ini"
    generated_rsr_file = f"{TARGET_CODE}.rsr"
    
    print(f"  > Moviendo {generated_ini_file} -> {OUTPUT_INI_PATH}")
    os.rename(generated_ini_file, OUTPUT_INI_PATH)
    
    print(f"  > Moviendo {generated_rsr_file} -> {OUTPUT_RSR_PATH}")
    os.rename(generated_rsr_file, OUTPUT_RSR_PATH)

    # Limpiar .sch
    temp_file = f"{TARGET_CODE}.sch"
    if os.path.exists(temp_file):
        os.remove(temp_file)

    print("\n" + "="*70)
    print(f"¡ÉXITO! Se han pre-calculado:")
    print(f"  1. {os.path.basename(OUTPUT_INI_PATH)} (ponderado)")
    print(f"  2. {os.path.basename(OUTPUT_RSR_PATH)} (ponderado)")
    print("="*70 + "\n")


if __name__ == "__main__":
    # Añade esto a tu prism/config.py si no lo tienes
    if not hasattr(config, 'MAX_CA_CA_DISTANCE'):
        print("Añadiendo MAX_CA_CA_DISTANCE=14.0 a config")
        config.MAX_CA_CA_DISTANCE = 14.0
    if not hasattr(config, 'MAX_N_O_DISTANCE'):
        print("Añadiendo MAX_N_O_DISTANCE=11.0 a config")
        config.MAX_N_O_DISTANCE = 11.0

    generate_weighted_precalc()