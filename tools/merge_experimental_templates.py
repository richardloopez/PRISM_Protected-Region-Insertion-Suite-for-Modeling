#!/usr/bin/env python3
# PRISM Tool: merge_experimental_to_predicted (Biopython version)
# Author: Antigravity (Adapted for Biopython - FIXED)

import sys
import argparse
import os
from typing import Dict

from Bio.PDB import PDBParser, Superimposer, PDBIO, Structure, Model, Chain

def parse_pir_sequences(align_file: str):
    '''Manual parsing of PIR file to get full alignment strings including gaps.'''
    sequences = {}
    all_codes = []
    curr_code = None
    curr_seq = []
    with open(align_file, 'r') as f:
        for line in f:
            if line.startswith('>P1;'):
                if curr_code:
                    sequences[curr_code] = "".join(curr_seq).replace('\n', '').replace(' ', '')
                curr_code = line[4:].strip()
                descript_code = curr_code.split('_')[0]
                all_codes.append(descript_code)
                curr_seq = []
            elif curr_code and not line.startswith(('structure', 'sequence', ' ')):
                curr_seq.append(line.strip().rstrip('*'))
        if curr_code:
            sequences[curr_code] = "".join(curr_seq).replace('\n', '').replace(' ', '')
    return sequences, all_codes

def get_all_residues(structure):
    '''
    Returns all the residues from the first model.
    '''
    return list(structure[0].get_residues())

def merge_structures(align_file: str, ref_code: str, output_pdb: str = None):
    pir_seqs, all_codes = parse_pir_sequences(align_file)
    print(f"PIR sequences found: {list(pir_seqs.keys())}")
    
    if not output_pdb:
        ref_prefix = ref_code.split('_')[0]
        template_codes = [c for c in all_codes if c != ref_prefix and c != "FullSeq"]
        output_pdb = "_".join(template_codes) + "_merged_experimental.pdb"
        print(f"Defaulting output PDB name to: {output_pdb}")

    target_code = None
    templates = []
    
    with open(align_file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith('>P1;'):
                code = line[4:].strip()
                type_line = lines[i+1]
                if type_line.startswith('sequence:'):
                    target_code = code
                elif type_line.startswith('structure') and code != ref_code:
                    templates.append(code)

    if ref_code not in pir_seqs:
        print(f"Error: Reference model '{ref_code}' not found in alignment.")
        return
    if not target_code:
        print(f"Error: Target sequence {target_code} not found in alignment.")
        return

    print(f"Reference: {ref_code}")
    print(f"Target: {target_code}")
    print(f"Experimental Templates: {templates}")

    target_seq = pir_seqs[target_code]
    ref_seq = pir_seqs[ref_code]

    parser = PDBParser(QUIET=True)
    
    ref_pdb_path = os.path.join(os.getcwd(), f"{ref_code}")
    if not os.path.exists(ref_pdb_path):
        print(f"Error: Could not find PDB for reference {ref_pdb_path}")
        return
    ref_structure = parser.get_structure(ref_code, ref_pdb_path)
    ref_residues = get_all_residues(ref_structure)

    transformed_structures = {}
    superimposer = Superimposer()
    
    for temp_code in templates:
        print(f"Superposing {temp_code} onto {ref_code}...")
        temp_pdb = os.path.join(os.getcwd(), f"{temp_code}")
        
        if not os.path.exists(temp_pdb):
            print(f"  Warning: Could not find PDB for {temp_code}. Skipping.")
            continue
            
        temp_structure = parser.get_structure(temp_code, temp_pdb)
        temp_residues = get_all_residues(temp_structure)
        temp_seq = pir_seqs[temp_code]
            
        ref_atoms = []
        temp_atoms = []
        r_idx, t_idx = 0, 0
            
        for r_char, t_char in zip(ref_seq, temp_seq):
            is_r_valid = r_char != '-'
            is_t_valid = t_char != '-'
                
            if is_r_valid and is_t_valid:
                if r_idx < len(ref_residues) and t_idx < len(temp_residues):
                    r_res = ref_residues[r_idx]
                    t_res = temp_residues[t_idx]
                    if 'CA' in r_res and 'CA' in t_res:
                        ref_atoms.append(r_res['CA'])
                        temp_atoms.append(t_res['CA'])
                
            if is_r_valid: r_idx += 1
            if is_t_valid: t_idx += 1
            
        if len(ref_atoms) == 0:
            print(f"  Warning: No aligned CA atoms found for {temp_code}.")
            continue
                
        superimposer.set_atoms(ref_atoms, temp_atoms)
        superimposer.apply(temp_structure.get_atoms())
            
        transformed_structures[temp_code] = temp_structure
        print(f"  RMSD: {superimposer.rms:.3f} Å")
            

    print("Merging superposed template structures into merged model...")
    merged_structure = Structure.Structure("Merged")
    merged_model = Model.Model(0)
    merged_chain = Chain.Chain('A')
    merged_model.add(merged_chain)
    merged_structure.add(merged_model)
    
    current_res_seq = 0
    
    for pos in range(len(target_seq)):
        if target_seq[pos] in ('-', '.', '/'): 
            continue
        
        chosen_temp = None
        for t_code in templates:
            if t_code in transformed_structures and pir_seqs[t_code][pos] != '-':
                chosen_temp = t_code
                break
                
        if chosen_temp:
            temp_seq_up_to_pos = pir_seqs[chosen_temp][:pos]
            temp_res_idx = len(temp_seq_up_to_pos.replace('-', ''))
            
            temp_structure = transformed_structures[chosen_temp]
            temp_residues = get_all_residues(temp_structure)
            
            if temp_res_idx < len(temp_residues):
                current_res_seq += 1
                res_to_copy = temp_residues[temp_res_idx].copy()
                res_to_copy.id = (' ', current_res_seq, ' ')
                merged_chain.add(res_to_copy)

    if templates and templates[0] in transformed_structures:
        first_temp = transformed_structures[templates[0]]
        merged_chain_b = Chain.Chain('B')
        for chain in first_temp[0]:
            if chain.id != 'A':
                for res in chain:
                    current_res_seq += 1
                    res_copy = res.copy()
                    res_copy.id = (res_copy.id[0], current_res_seq, ' ')
                    merged_chain_b.add(res_copy)
        if len(merged_chain_b) > 0:
            merged_model.add(merged_chain_b)

    io = PDBIO()
    io.set_structure(merged_structure)
    io.save(output_pdb)
    print(f"Merged experimental PDB written to {output_pdb}")
    
    output_ali = align_file.replace('.ali', '_merged.ali')
    merged_seq = list('-' * len(target_seq))
    for pos in range(len(target_seq)):
        for t_code in templates:
            if pir_seqs[t_code][pos] != '-':
                merged_seq[pos] = pir_seqs[t_code][pos]
                break
                
    with open(output_ali, 'w') as f:
        f.write(f">P1;{output_pdb}\n")
        f.write(f"structure:{output_pdb}:FIRST:@:END:@::::\n")
        f.write("".join(merged_seq) + "*\n")
        
        f.write(f">P1;{ref_code}\n")
        atom_file = ref_code
        f.write(f"structure:{atom_file}:FIRST:@:END:@::::\n")
        f.write(ref_seq + "*\n")
        
        f.write(f">P1;{target_code}\n")
        f.write(f"sequence:{target_code}:FIRST:@:END:@::::\n")
        f.write(target_seq + "*\n")
        
    print(f"Merged alignment written to {output_ali}")

def main():
    parser = argparse.ArgumentParser(description="Merge multiple experimental structures aligned to a predicted model using Biopython.")
    parser.add_argument("alignment", help="Path to the alignment file (.ali)")
    parser.add_argument("reference", help="Code of the predictive model in the alignment to use as reference")
    parser.add_argument("--output_pdb", default=None, help="Output filename for the merged PDB")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.alignment):
        print(f"Error: Alignment file {args.alignment} not found.")
        sys.exit(1)
        
    merge_structures(args.alignment, args.reference, args.output_pdb)
    

if __name__ == "__main__":
    main()