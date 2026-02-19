#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM Tool: unify_templates
#
'''
Tool to resolve sequence overlaps between multiple templates in a Modeller alignment.
Trims protein residues from lower-priority templates while maintaining a buffer.
Preserves BLK residues and chain breaks (/) in their original alignment positions.
'''

import sys
import argparse
import os
import re
from typing import List, Dict, Any, Tuple, Set

from modeller import *

class PDBAtom:
    def __init__(self, line):
        self.line = line
        self.record_type = line[0:6].strip()
        try:
            self.serial = int(line[6:11])
        except ValueError: self.serial = 0
        self.name = line[12:16] 
        self.alt_loc = line[16]
        self.res_name = line[17:20].strip()
        self.chain_id = line[21]
        try:
            self.res_seq = int(line[22:26])
        except ValueError: self.res_seq = 0
        self.i_code = line[26]
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occ = float(line[54:60]) if len(line) > 54 and line[54:60].strip() else 1.00
        self.temp = float(line[60:66]) if len(line) > 60 and line[60:66].strip() else 0.00
        self.element = line[76:78].strip() if len(line) > 76 else ""

    def to_pdb_line(self):
        if len(self.name) == 4:
            name_str = f"{self.name}"
        else:
            if self.name[0].isdigit(): name_str = f"{self.name:<4}"
            else: name_str = f" {self.name:<3}"
        return (f"{self.record_type:<6}{self.serial:>5} {name_str:4}{self.alt_loc}{self.res_name:>3} {self.chain_id}{self.res_seq:>4}{self.i_code}   "
                f"{self.x:>8.3f}{self.y:>8.3f}{self.z:>8.3f}{self.occ:>6.2f}{self.temp:>6.2f}          {self.element:>2}")

def renumber_pdb(pdb_path: str, residues_to_keep: Set[int], output_path: str):
    '''
    Reads PDB, keeps only protein residues in residues_to_keep SET + all BLK/HETATM,
    and renumbers BOTH protein and BLK/HETATM residues.
    Protein chains reset to 1. BLK/HETATM chains continue numbering from last residue.
    '''
    protein_atoms = []
    blk_atoms = []
    
    with open(pdb_path, 'r') as f:
        prot_res_count = 0
        last_prot_res_id = None
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom = PDBAtom(line)
                res_id = (atom.chain_id, atom.res_seq, atom.i_code)
                
                if atom.record_type == 'HETATM' or atom.res_name == 'BLK':
                    blk_atoms.append(atom)
                else:
                    if res_id != last_prot_res_id:
                        prot_res_count += 1
                        last_prot_res_id = res_id
                    
                    if prot_res_count in residues_to_keep:
                        protein_atoms.append(atom)
    
    new_lines = []
    serial = 1
    new_res_seq = 0
    
    # 1. Renumber Protein
    last_old_res_id = None
    last_chain_id = None
    for atom in protein_atoms:
        if last_chain_id is not None and atom.chain_id != last_chain_id:
            new_lines.append("TER")
            new_res_seq = 0
            last_old_res_id = None
        
        old_res_id = (atom.chain_id, atom.res_seq, atom.i_code)
        if old_res_id != last_old_res_id:
            new_res_seq += 1
            last_old_res_id = old_res_id
        
        atom.serial = serial
        atom.res_seq = new_res_seq
        new_lines.append(atom.to_pdb_line())
        serial += 1
        last_chain_id = atom.chain_id
    
    if protein_atoms:
        new_lines.append("TER")
    
    # 2. Renumber BLK
    last_old_res_id = None
    for atom in blk_atoms:
        if last_chain_id is not None and atom.chain_id != last_chain_id:
            new_lines.append("TER")
            last_old_res_id = None
            
        old_res_id = (atom.chain_id, atom.res_seq, atom.i_code)
        if old_res_id != last_old_res_id:
            new_res_seq += 1
            last_old_res_id = old_res_id
        
        atom.serial = serial
        atom.res_seq = new_res_seq
        new_lines.append(atom.to_pdb_line())
        serial += 1
        last_chain_id = atom.chain_id
        
    if blk_atoms:
        new_lines.append("TER")
    
    new_lines.append("END")
    
    with open(output_path, 'w') as f:
        for l in new_lines:
            f.write(l + "\n")

def unify_templates(align_file: str, overlap_limit: int):
    env = Environ()
    aln = Alignment(env, file=align_file)
    
    templates = [s for s in aln if s.prottyp.startswith('structure')]
    if len(templates) < 2:
        print("Less than 2 templates found. Nothing to unify.")
        return

    def is_aa(c):
        return 'A' <= c <= 'Z'

    tmp_ali = "_tmp_unify.ali"
    aln.write(file=tmp_ali)
    ali_strings = {}
    with open(tmp_ali, 'r') as f:
        curr_code = None
        for line in f:
            if line.startswith('>P1;'):
                curr_code = line[4:].strip()
                ali_strings[curr_code] = ""
            elif curr_code and not line.startswith(('structure', 'sequence', 'C;', ' ')):
                cleaned = re.sub(r'\s+', '', line).rstrip('*')
                ali_strings[curr_code] += cleaned
    if os.path.exists(tmp_ali):
        os.remove(tmp_ali)

    covered_positions = set()
    new_sequences = {} 

    for i, temp in enumerate(templates):
        print(f"Processing template {i+1}: {temp.code}")
        
        temp_seq_str = ali_strings.get(temp.code, "")
        if not temp_seq_str:
            print(f"  Warning: Sequence for {temp.code} not found in alignment.")
            continue

        temp_aa_positions = [j for j in range(len(temp_seq_str)) if is_aa(temp_seq_str[j])]
        overlapping_aa_pos = [pos for pos in temp_aa_positions if pos in covered_positions]
        
        residues_to_keep = set()
        
        if not overlapping_aa_pos or len(overlapping_aa_pos) <= overlap_limit:
            residues_to_keep = set(range(1, len(temp.residues) + 1))
            print(f"  No excessive overlap. Keeping all {len(residues_to_keep)} protein residues.")
        else:
            non_overlapping_pos = set(pos for pos in temp_aa_positions if pos not in covered_positions)
            
            if not non_overlapping_pos:
                print(f"  WARNING: Completely covered! Keeping only buffer.")
                keep_aa_pos = set(temp_aa_positions[:overlap_limit])
            else:
                expanded_nop = set()
                for nop in non_overlapping_pos:
                    for offset in range(-overlap_limit, overlap_limit + 1):
                        expanded_nop.add(nop + offset)
                
                keep_aa_pos = expanded_nop.intersection(set(temp_aa_positions))

            res_idx = 0
            for j in range(len(temp_seq_str)):
                if is_aa(temp_seq_str[j]):
                    res_idx += 1
                    if j in keep_aa_pos:
                        residues_to_keep.add(res_idx)
            
            print(f"  Overlap resolved. Keeping {len(residues_to_keep)} protein residues.")

        # Update covered positions and create new alignment sequence string
        res_idx = 0
        new_seq_list = []
        for j in range(len(temp_seq_str)):
            char = temp_seq_str[j]
            if is_aa(char):
                res_idx += 1
                if res_idx in residues_to_keep:
                    covered_positions.add(j)
                    new_seq_list.append(char)
                else:
                    new_seq_list.append('-')
            else:
                new_seq_list.append(char)
        
        last_aa_pos = -1
        for j in range(len(new_seq_list)-1, -1, -1):
            if is_aa(new_seq_list[j]):
                last_aa_pos = j
                break
        
        tail_start = last_aa_pos + 1
        tail_chars = new_seq_list[tail_start:]
        if tail_chars:
            num_dots = tail_chars.count('.')
            has_slash = '/' in tail_chars
            num_gaps = len(tail_chars) - num_dots - (1 if has_slash else 0)
            
            new_tail = (['-'] * num_gaps)
            if has_slash:
                new_tail.append('/')
            new_tail.extend(['.'] * num_dots)
            new_seq_list[tail_start:] = new_tail
            
        new_sequences[temp.code] = "".join(new_seq_list)

        orig_pdb = temp.atom_file.strip() if temp.atom_file else ""
        if not orig_pdb or not os.path.exists(orig_pdb):
            orig_pdb = f"{temp.code.strip()}.pdb"
            
        if os.path.exists(orig_pdb):
            output_pdb = f"{temp.code.split('.')[0]}_unified.pdb"
            
            renumber_pdb(orig_pdb, residues_to_keep, output_pdb)
            print(f"  Generated {output_pdb}")
            temp.atom_file = output_pdb
        else:
            print(f"  Error: {orig_pdb} not found. Cannot modify structure.")

    output_ali = align_file.replace('.ali', '_unified.ali')
    with open(align_file, 'r') as f_in, open(output_ali, 'w') as f_out:
        current_code = None
        seq_buffer = ""
        for line in f_in:
            if line.startswith('>P1;'):
                current_code = line[4:].strip()
                f_out.write(line)
            elif current_code in new_sequences and not line.startswith(('structure', 'sequence', 'C;', ' ')):
                if '*' in line:
                    f_out.write(new_sequences[current_code] + "*\n")
                    current_code = None
                else:
                    pass
            elif current_code in new_sequences and line.startswith(('structure', 'sequence')):
                f_out.write(line)
            else:
                f_out.write(line)

    print(f"\nSuccess! Unified alignment written to {output_ali}")
    print('''
            REMEMBER: RENAME THE GENERATED .PDBs (preferred)
            OR 
            UPDATE THE NAMES IN:
                - ALIGNMENT FILE
                - CONFIG.YAML (or GUI)
            
            (Names are not overwritten automatically so the user can verify the files are correct.)
    ''')

def main():
    parser = argparse.ArgumentParser(description="Unify templates by resolving sequence overlaps.")
    parser.add_argument("alignment", help="Path to the Modeller alignment file (.ali)")
    parser.add_argument("--overlap", type=int, default=10, help="Number of residues allowed to overlap for continuity (default: 10)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.alignment):
        print(f"Error: Alignment file {args.alignment} not found.")
        sys.exit(1)
        
    unify_templates(args.alignment, args.overlap)

if __name__ == "__main__":
    main()
