#!/usr/bin/env python3
"""
Script to verify that experimental coordinates did NOT move.

Compares the original MAIN_PDB with a generated model.
FIRST, it aligns the model onto the original (using a C-alpha mapping)
and THEN calculates the RMSD of those same residues.
"""

import sys
import os

try:
    from modeller import *
    # We don't need complete_pdb for this script
except ImportError:
    print("ERROR: MODELLER is not installed or not in PYTHONPATH")
    sys.exit(1)

def verify_fixed_coordinates(original_pdb_path: str, model_pdb_path: str, 
                             residue_mapping: dict, 
                             original_chain_id: str = 'A', 
                             model_chain_id: str = 'A'):
    """
    Verifies that experimental residues have not moved,
    using an explicit (dictionary) mapping.
    """
    env = Environ()
    
    # Load standard topology and parameter libraries
    env.libs.topology.read(file='$(LIB)/top.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    
    env.io.atom_files_directory = ['.']
    
    # --- Load FULL models ---
    print(f"\nLoading original PDB: {original_pdb_path} (Chain: {original_chain_id})")
    mdl_original = Model(env)
    mdl_original.read(file=original_pdb_path)
    
    print(f"Loading generated model: {model_pdb_path} (Chain: {model_chain_id})")
    mdl_model = Model(env)
    mdl_model.read(file=model_pdb_path)

    # --- Create 1-to-1 alignment using ONLY the FRAGMENTS ---
    print("\nCreating 1-to-1 alignment for superposition...")
    aln = Alignment(env)
    
    # Original Fragment
    orig_res_min = min(residue_mapping.keys())
    orig_res_max = max(residue_mapping.keys())
    mdl_frag_original = Model(env)
    mdl_frag_original.read(file=original_pdb_path, 
                           model_segment=(f'{orig_res_min}:{original_chain_id}', f'{orig_res_max}:{original_chain_id}'))
    aln.append_model(mdl_frag_original, align_codes='ORIG_FRAG', atom_files=original_pdb_path)

    # Model Fragment
    model_res_min = min(residue_mapping.values())
    model_res_max = max(residue_mapping.values())
    mdl_frag_model = Model(env)
    mdl_frag_model.read(file=model_pdb_path, 
                        model_segment=(f'{model_res_min}:{model_chain_id}', f'{model_res_max}:{model_chain_id}'))
    aln.append_model(mdl_frag_model, align_codes='MOD_FRAG', atom_files=model_pdb_path)
    
    # --- Selection of C-alphas from the REFERENCE FRAGMENT ---
    print(f"\nCreating selection of {len(mdl_frag_original.residues)} C-alphas from the original fragment...")
    sel_frag_original = Selection(mdl_frag_original).only_atom_types('CA')
    
    if not sel_frag_original:
        print("❌ ERROR: Selection for alignment is empty.")
        return False

    print(f"Selected {len(sel_frag_original)} C-alpha atoms from the original PDB.")
    
    # --- Superpose FRAGMENT vs FRAGMENT ---
    print("\nPerforming structural superposition (fragment vs fragment)...")
    r = sel_frag_original.superpose(mdl_frag_model, aln)
    global_rmsd = r.rms
    
    print(f"✓ Fragment superposition OK (RMSD={global_rmsd:.6f} Å).")
    
    # --- Apply that same transformation to the FULL MODEL ---
    print(f"Applying transformation to the full {model_pdb_path}...")
    sel_full_model = Selection(mdl_model)
    sel_full_model.transform(r.rotation)
    sel_full_model.translate(r.translation)
    
    # --- CALCULATION LOOP (NO PRINTS) ---
    # First, calculate everything and store in a list
    
    max_rmsd = 0.0
    total_rmsd = 0.0
    results_list = [] # Store all results here
    calculation_errors_found = False
    
    for res_orig_id, res_model_id in sorted(residue_mapping.items()):
        try:
            res_orig = mdl_original.residues[f'{res_orig_id}:{original_chain_id}']
            res_model = mdl_model.residues[f'{res_model_id}:{model_chain_id}']
            
            aa_orig = res_orig.pdb_name
            aa_mod = res_model.pdb_name
            
            ca_orig = res_orig.atoms['CA']
            ca_model = res_model.atoms['CA']
            
            rmsd = ((ca_orig.x - ca_model.x)**2 + 
                    (ca_orig.y - ca_model.y)**2 + 
                    (ca_orig.z - ca_model.z)**2)**0.5
            
            total_rmsd += rmsd
            max_rmsd = max(max_rmsd, rmsd)
            
            results_list.append((res_orig_id, res_model_id, aa_orig, aa_mod, rmsd))
                    
        except Exception as e:
            print(f"❌  Error calculating {res_orig_id} (Orig) -> {res_model_id} (Mod): {e}")
            results_list.append((res_orig_id, res_model_id, '???', '???', None))
            calculation_errors_found = True
    
    if not residue_mapping:
        print("\n❌ RESULT: The mapping dictionary is empty.")
        return False

    # --- SECTION 1: ORDERED BY POSITION ---
    print(f"\n--------------------------------------------------------------------------------")
    print(f"RESIDUES ORDERED BY POSITION")
    print(f"--------------------------------------------------------------------------------")
    
    for res_o, res_m, aa_o, aa_m, rmsd in results_list: # Already sorted by position
        if rmsd is not None:
            # Add a simple visual indicator
            prefix = "✓  " if rmsd <= 0.01 else "   "
            print(f"{prefix} Res {aa_o:3s} {res_o:4d} (Orig) -> {aa_m:3s} {res_m:4d} (Mod): {rmsd:.6f} Å")
        else:
            print(f"❌   Res {aa_o:3s} {res_o:4d} (Orig) -> {aa_m:3s} {res_m:4d} (Mod): Comparison Error")

    # --- SECTION 2: ORDERED BY RMSD ---
    print(f"\n--------------------------------------------------------------------------------")
    print(f"RESIDUES ORDERED BY RMSD (Top Deviations)")
    print(f"--------------------------------------------------------------------------------")
    
    # Sort by RMSD (descending), Nones at the end
    sorted_results = sorted(results_list, 
                             key=lambda x: x[4] if x[4] is not None else -float('inf'), 
                             reverse=True)
    
    for res_o, res_m, aa_o, aa_m, rmsd in sorted_results:
        if rmsd is not None:
            prefix = "✓  " if rmsd <= 0.01 else "   "
            print(f"{prefix} Res {aa_o:3s} {res_o:4d} (Orig) -> {aa_m:3s} {res_m:4d} (Mod): {rmsd:.6f} Å")
        else:
            print(f"❌   Res {aa_o:3s} {res_o:4d} (Orig) -> {aa_m:3s} {res_m:4d} (Mod): Comparison Error")

    # --- SECTION 3: RMSD SUMMARY ---
    print(f"\n--------------------------------------------------------------------------------")
    print(f"RMSD SUMMARY")
    print(f"--------------------------------------------------------------------------------")
    print(f"  - Average RMSD (individual): {total_rmsd / len(residue_mapping):.6f} Å")
    print(f"  - Maximum RMSD (individual):   {max_rmsd:.6f} Å")
    print(f"  - Global RMSD (from 'superpose'): {global_rmsd:.6f} Å")
    
    # Decide success or failure for the exit code
    rmsd_errors_found = any(r[4] > 0.01 for r in results_list if r[4] is not None)

    if not rmsd_errors_found and not calculation_errors_found:
        print(f"\n✓  SUCCESS: All mapped experimental residues remained fixed (RMSD < 0.01 Å)")
        return True
    else:
        if calculation_errors_found:
            print("\nWARNING: Errors occurred while comparing some residues.")
        if rmsd_errors_found:
             print("\nWARNING: Some residues moved (RMSD > 0.01 Å).")
        return False

if __name__ == '__main__':
    print("="*80)
    print("COORDINATE VERIFICATION (WITH ALIGNMENT)")
    print("="*80)
    
    # --- USER CONFIGURATION ---
    ORIGINAL_PDB = '9CB5_DS_renum_HETATM-B.pdb'
    MODEL_PDB = 'pdb-from-prism.pdb'
    
    # --- CONFIGURE YOUR MAPPING HERE ---
    # Map 1-166 (Orig) to 307-472 (Mod)
    
    orig_res = list(range(1, 167))      # 1 to 166
    mod_res = list(range(307, 473))   # 307 to 472
    
    if len(orig_res) != len(mod_res):
        print("ERROR in script setup: mapping ranges do not match.")
        print(f"Original: {len(orig_res)} residues, Model: {len(mod_res)} residues")
        sys.exit(1)
        
    RESIDUE_MAPPING = dict(zip(orig_res, mod_res))
    
    # -------------------------------------------------------------------
    
    if not os.path.exists(ORIGINAL_PDB):
        print(f"\n❌  ERROR: File not found: {ORIGINAL_PDB}")
        sys.exit(1)
    
    if not os.path.exists(MODEL_PDB):
        print(f"\n❌  ERROR: File not found: {MODEL_PDB}")
        sys.exit(1)
    
    print(f"\nOriginal PDB:  {ORIGINAL_PDB}")
    print(f"Model PDB:    {MODEL_PDB}")
    print(f"Residues to map and align: {len(RESIDUE_MAPPING)}")
    
    success = verify_fixed_coordinates(
        ORIGINAL_PDB, 
        MODEL_PDB, 
        RESIDUE_MAPPING,
        original_chain_id='A', # Chain to use in Original PDB
        model_chain_id='A'   # Chain to use in Model PDB
    )
    
    sys.exit(0 if success else 1)