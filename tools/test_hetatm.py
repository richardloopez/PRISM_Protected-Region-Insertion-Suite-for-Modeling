#!/usr/bin/env python3
"""
Improved script for HETATM/BLK detection.
Automatically detects all chains and analyzes them individually.
Does not require Modeller, only reads the PDB file directly.
"""

import sys
import os
from typing import List, Dict, Any, Set

# --- EXTRACTION FUNCTIONS (Modified to remove try/except) ---

def extract_hetatm_residues_simple(pdb_file: str, chain_id: str) -> List[Dict[str, Any]]:
    """
    Extracts unique HETATM residue info for a specific CHAIN.
    """
    hetatm_residues = []
    last_atom_resnum = 0
    seen_hetatm = set()
    
    with open(pdb_file, 'r') as f:
        for line in f:
            res_chain = line[21:22].strip()
            
            # Only process lines for the requested chain
            if res_chain != chain_id:
                continue

            if line.startswith('ATOM'):
                try:
                    resnum = int(line[22:26].strip())
                    last_atom_resnum = max(last_atom_resnum, resnum)
                except ValueError:
                    continue
            
            elif line.startswith('HETATM'):
                resname = line[17:20].strip()
                try:
                    resnum = int(line[22:26].strip())
                except ValueError:
                    continue
                
                res_key = (resname, resnum, res_chain)
                if res_key not in seen_hetatm:
                    seen_hetatm.add(res_key)
                    hetatm_residues.append({
                        'resname': resname,
                        'resnum': resnum,
                        'chain': res_chain,
                        'position_after_atom_resnum': last_atom_resnum
                    })
    
    hetatm_residues.sort(key=lambda x: x['resnum'])
    return hetatm_residues

def get_all_chains(pdb_file: str) -> Dict[str, Dict[str, int]]:
    """Gets information about all chains in the PDB."""
    chains_info = {}
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain = line[21:22].strip()
                if not chain: # Skip unidentified chains
                    continue
                record_type = 'ATOM' if line.startswith('ATOM') else 'HETATM'
                
                if chain not in chains_info:
                    chains_info[chain] = {'ATOM': 0, 'HETATM': 0}
                chains_info[chain][record_type] += 1
    
    return chains_info

# --- HELPER PRINTING FUNCTIONS ---

def print_hetatm_summary_tables(hetatm_residues_list: List[Dict[str, Any]]):
    """Prints summary tables of HETATM types and details."""
    
    # Group by residue type
    residue_types = {}
    for het in hetatm_residues_list:
        resname = het['resname']
        if resname not in residue_types:
            residue_types[resname] = 0
        residue_types[resname] += 1
    
    print(f"\n   HETATM Residue Types:")
    for resname, count in sorted(residue_types.items()):
        print(f"     - {resname}: {count} residue(s)")
    
    print(f"\n   First 10 HETATM residues detected:")
    print(f"     {'Type':<8} {'ResNum':<8} {'Chain':<8} {'After ATOM #':<20}")
    print(f"     {'-'*50}")
    for het in hetatm_residues_list[:10]:
        print(f"     {het['resname']:<8} {het['resnum']:<8} {het['chain']:<8} {het['position_after_atom_resnum']:<20}")
    
    if len(hetatm_residues_list) > 10:
        print(f"     ... and {len(hetatm_residues_list) - 10} more")

def print_modeller_explanation(hetatm_count: int, chain_id: str):
    """Prints modeling instructions for a specific chain."""
    print(f"\n   --- 💡 Modeling Instructions (Chain {chain_id}) ---")
    print(f"    a) `env.io.hetatm = True` is set ✓")
    print(f"    b) Modeller will read {hetatm_count} HETATM residues automatically for this chain.")
    print(f"    c) In the PIR alignment for CHAIN {chain_id}, YOU MUST ADD:")
    print(f"       - Template: {hetatm_count} '.' (BLK) characters at the end of the chain {chain_id} sequence")
    print(f"       - Target:   {hetatm_count} '.' (BLK) characters at the end of the chain {chain_id} sequence")
    print(f"    d) HETATM residues will be treated as obstacles during modeling.")
    print("   -----------------------------------------------------")


# --- MAIN FUNCTION ---

def main():
    """
    Runs the HETATM detection analysis on all chains of a PDB file.
    """
    
    # --- 1. Setup and File Loading ---
    if len(sys.argv) > 1:
        PDB_FILE = sys.argv[1]
    else:
        print("Usage: python test_hetatm.py <file.pdb>")
        PDB_FILE = '8vx1_DS_renum_HETATM.pdb' # Default test file
        print(f"Using default test file: '{PDB_FILE}'\n")

    if not os.path.exists(PDB_FILE):
        print(f"[ERROR] File not found: {PDB_FILE}")
        return 1

    # Assume 'A' is always the main peptidic chain
    PEPTIDIC_CHAIN_ID = 'A'
    modeled_hetatm_chains = {} # For the final summary

    print("="*70)
    print("HETATM/BLK DETECTION TEST (v2.0 Automatic)")
    print("="*70)
    print(f"\n1. Template PDB File: {PDB_FILE}")

    # --- 2. Get Chain Summary ---
    try:
        chains_info = get_all_chains(PDB_FILE)
    except Exception as e:
        print(f"[ERROR] Could not read PDB file '{PDB_FILE}'. Error: {e}")
        return 1
        
    if not chains_info:
        print("[ERROR] Empty PDB file or no chains could be read.")
        return 1

    print(f"\n2. PDB Chain Summary:")
    print(f"   {'Chain':<10} {'ATOM (lines)':<15} {'HETATM (lines)':<15}")
    print(f"   {'-'*42}")
    for chain, counts in sorted(chains_info.items()):
        print(f"   {chain:<10} {counts['ATOM']:<15} {counts['HETATM']:<15}")

    # --- 3. Detailed Chain Analysis ---
    print("\n" + "="*70)
    print("3. DETAILED CHAIN ANALYSIS")
    print("="*70)

    for chain_id in sorted(chains_info.keys()):
        print(f"\n### 🧬 Analyzing Chain: {chain_id} ###")
        
        try:
            hetatm_residues = extract_hetatm_residues_simple(PDB_FILE, chain_id)
        except Exception as e:
            print(f"  [ERROR] An error occurred while analyzing chain {chain_id}: {e}")
            continue

        num_het_residues = len(hetatm_residues)
        atom_count = chains_info[chain_id]['ATOM']

        # Logic for the main peptidic chain
        if chain_id == PEPTIDIC_CHAIN_ID:
            if num_het_residues > 0:
                print(f"  ❌ WARNING! Detected {num_het_residues} HETATM residues in the main peptide chain '{chain_id}'.")
                print("     This chain should ideally only contain ATOM records. Please review your PDB file.")
                print_hetatm_summary_tables(hetatm_residues)
            else:
                print(f"  ✓ OK: Peptide chain '{chain_id}' contains no HETATM residues.")
                print(f"     Contains {atom_count} ATOM records.")
        
        # Logic for other chains (ligands, DNA, other proteins)
        else:
            if num_het_residues > 0:
                print(f"  ⓘ DETECTED: {num_het_residues} HETATM residues (ligands, DNA, etc.).")
                print(f"     This chain also contains {atom_count} ATOM records.")
                modeled_hetatm_chains[chain_id] = num_het_residues
                
                # Show tables and modeling explanation
                print_hetatm_summary_tables(hetatm_residues)
                print_modeller_explanation(num_het_residues, chain_id)
            else:
                print(f"  ⓘ No HETATM residues detected in this chain.")
                if atom_count > 0:
                    print(f"     Contains {atom_count} ATOM records (considered a protein chain).")
                else:
                    print(f"     This chain contains no protein residues (0 ATOM).")

    # --- 4. Final Summary ---
    print("\n" + "="*70)
    print("✓ TEST SUMMARY COMPLETE")
    print("="*70)
    
    if not modeled_hetatm_chains:
        print("\nThe system is ready. No HETATM residues requiring BLK modeling")
        print("were detected (aside from chain 'A').")
    else:
        # Calculate total by summing the values (counts) from the dictionary
        total_hetatm_count = sum(modeled_hetatm_chains.values())

        print("\nThe system is ready to model with HETATM/BLK.")
        print("Summary of chains that will require BLK alignment:")
            
        for chain_id, count in modeled_hetatm_chains.items():
            print(f"    - Chain {chain_id}: {count} HETATM residues (must add {count} '.' characters to ALI file for this chain)")
            
        print("---") # Separator for clarity
        print(f"TOTAL HETATM (BLK): {total_hetatm_count}")
        print(f"(A total of {total_hetatm_count} '.' characters must be added to the ALI file)")

    print("\n" + "="*70)
    return 0

if __name__ == '__main__':
    sys.exit(main())