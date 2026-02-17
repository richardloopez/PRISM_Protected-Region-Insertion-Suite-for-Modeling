#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# PRISM (Verification Suite)

import os
import sys
import argparse
import re
from typing import Dict, List, Tuple, Set
from modeller import *

class PrismVerify:
    def __init__(self, original_pdb: str, model_pdb: str, align_file: str):
        self.original_pdb = original_pdb
        self.model_pdb = model_pdb
        self.align_file = align_file
        self.env = Environ()
        self._setup_modeller()

    def _setup_modeller(self):
        '''
        Standard Modeller environment setup.
        '''
        self.env.libs.topology.read(file='$(LIB)/top.lib')
        self.env.libs.parameters.read(file='$(LIB)/par.lib')
        self.env.io.atom_files_directory = ['.']
        log.none()

    def parse_pir_alignment(self, template_code: str) -> Tuple[str, str]:
        '''
        Parses PIR file to find the template and target sequences.
        Assumes the target is the entry with the 'sequence' header.
        '''
        if not os.path.exists(self.align_file):
            sys.exit(f"❌ Error: Alignment file {self.align_file} not found.")

        with open(self.align_file, 'r') as f:
            content = f.read()

        blocks = content.split('>P1;')
        template_seq = ""
        target_seq = ""

        for block in blocks[1:]:
            lines = block.strip().splitlines()
            code = lines[0].strip()
            header = lines[1].strip()
            sequence = "".join([l.strip() for l in lines[2:] if l.strip() and not l.startswith('#')])
            
            if code == template_code:
                template_seq = sequence
            elif header.startswith('sequence:'):
                target_seq = sequence

        if not template_seq or not target_seq:
            sys.exit(f"❌ Error: Could not find template '{template_code}' or target sequence in alignment.")

        return template_seq, target_seq

    def build_mapping_from_alignment(self, template_code: str) -> Dict[int, int]:
        '''
        Calculates residue mapping where both sequences have a match (no gaps).
        '''
        t_seq, q_seq = self.parse_pir_alignment(template_code)
        
        mapping = {}
        t_res_count = 0
        q_res_count = 0

        for t_char, q_char in zip(t_seq, q_seq):
            if t_char not in ('-', '/', '*'): t_res_count += 1
            if q_char not in ('-', '/', '*'): q_res_count += 1

            if t_char not in ('-', '/', '.') and q_char not in ('-', '/', '.'):
                mapping[t_res_count] = q_res_count

        return mapping

    def run_rmsd_check(self, mapping: Dict[int, int], orig_chain: str, mod_chain: str):
        '''
        Calculates CA RMSD using Modeller models.
        '''
        mdl_orig = Model(self.env)
        mdl_orig.read(file=self.original_pdb)
        
        mdl_mod = Model(self.env)
        mdl_mod.read(file=self.model_pdb)

        results = []
        total_sq_diff = 0
        count = 0

        print(f"\n{'='*80}")
        print(f"{'RESIDUE':<15} | {'ORIG POS':<10} | {'MOD POS':<10} | {'RMSD (Å)':<10}")
        print("-" * 80)

        for t_idx, q_idx in sorted(mapping.items()):
            try:
                res_o = mdl_orig.residues[f'{t_idx}:{orig_chain}']
                res_m = mdl_mod.residues[f'{q_idx}:{mod_chain}']
                
                ca_o = res_o.atoms['CA']
                ca_m = res_m.atoms['CA']
                
                dist = ((ca_o.x - ca_m.x)**2 + (ca_o.y - ca_m.y)**2 + (ca_o.z - ca_m.z)**2)**0.5
                
                status = "✓" if dist < 0.01 else "!"
                print(f"{status} {res_o.pdb_name:<13} | {t_idx:<10} | {q_idx:<10} | {dist:.6f}")
                
                results.append(dist)
                total_sq_diff += dist**2
                count += 1
            except Exception:
                # Gaps or missing atoms in PDB
                continue

        if count == 0:
            print("❌ No matching residues found for comparison.")
            return

        avg_rmsd = sum(results) / count
        max_rmsd = max(results)
        
        print("-" * 80)
        print(f"SUMMARY for {count} residues:")
        print(f" > Average RMSD: {avg_rmsd:.6f} Å")
        print(f" > Maximum RMSD: {max_rmsd:.6f} Å")
        
        if max_rmsd < 0.01:
            print(f"\n✅ SUCCESS: Experimental coordinates are FIXED.")
        else:
            print(f"\n⚠️  WARNING: Coordinate drift detected (> 0.01 Å).")

def parse_manual_segments(segment_str: str) -> Dict[int, int]:
    '''
    Helper to parse manual flag: '1-100:5-105'
    '''
    mapping = {}
    try:
        orig_part, mod_part = segment_str.split(':')
        o_start, o_end = map(int, orig_part.split('-'))
        m_start = int(mod_part.split('-')[0])
        
        for i in range(o_end - o_start + 1):
            mapping[o_start + i] = m_start + i
    except ValueError:
        sys.exit("❌ Error: Manual segments format must be 'start-end:start-end' (e.g. 1-191:2-192)")
    return mapping

def main():
    parser = argparse.ArgumentParser(description="PRISM Coordinate Fidelity Verification")
    parser.add_argument("original_pdb", help="Original experimental PDB file.")
    parser.add_argument("model_pdb", help="Generated PRISM model PDB.")
    parser.add_argument("alignment", help="PIR alignment file used for modeling.")
    parser.add_argument("--manual", help="Manual segments (e.g., '1-191:2-192'). Overrides alignment.")
    parser.add_argument("--orig-chain", default="A", help="Chain in original PDB (default: A).")
    parser.add_argument("--mod-chain", default="A", help="Chain in model PDB (default: A).")

    args = parser.parse_args()

    template_code = os.path.basename(args.original_pdb)
    if template_code.endswith('.pdb'):
        raw_code = template_code.replace('.pdb', '')
    else:
        raw_code = template_code

    verifier = PrismVerify(args.original_pdb, args.model_pdb, args.alignment)

    if args.manual:
        print(f"[VERIFY] Using manual segment definition: {args.manual}")
        mapping = parse_manual_segments(args.manual)
    else:
        print(f"[VERIFY] Automatically building mapping from {args.alignment}...")
        mapping = verifier.build_mapping_from_alignment(template_code)
        if not mapping:
            mapping = verifier.build_mapping_from_alignment(raw_code)

    if not mapping:
        sys.exit(f"❌ Error: No mapping generated. Check if '{template_code}' is in {args.alignment}.")

    verifier.run_rmsd_check(mapping, args.orig_chain, args.mod_chain)

if __name__ == "__main__":
    main()
