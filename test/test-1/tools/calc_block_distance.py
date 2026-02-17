#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM Tool: calc_block_distance
#
'''
Script to compute the distance between protein CA atoms and BLK gravity centers.
This matches the logic used in PRISM's add_hetatm_repulsion_shield.
'''

import sys
import argparse
import os
from typing import List, Tuple, Dict, Any

from modeller import *
from modeller.selection import Selection


def calculate_distances(pdb_file: str, protein_chain: str, blk_chain: str) -> List[Dict[str, Any]]:
    '''
    Load PDB, identify BLK residues, calculate gravity centers, and find min distances to CA atoms.
    '''
    env = Environ()
    env.io.atom_files_directory = ['.']
    env.io.hetatm = True
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    mdl = Model(env, file=pdb_file)
    
    # Identify BLK residues
    het_residues = [r for r in mdl.residues if r.chain.name == blk_chain and (r.hetatm or r.name == 'BLK') and r.name != 'HOH']
    if not het_residues:
        available_chains = sorted(list(set(r.chain.name for r in mdl.residues)))
        print(f"No HETATM/BLK residues found in chain {blk_chain}.")
        print(f"Available chains in PDB: {', '.join(available_chains)}")
        return []

    # Calculate gravity centers
    het_centers = []
    for res in het_residues:
        het_centers.append(res)

    # Get protein CA atoms
    try:
        prot_chain_obj = mdl.chains[protein_chain]
    except KeyError:
        available_chains = sorted(list(mdl.chains.keys()))
        print(f"Error: Chain {protein_chain} not found in {pdb_file}.")
        print(f"Available chains in PDB: {', '.join(available_chains)}")
        return []

    protein_ca = Selection(prot_chain_obj).only_atom_types('CA')
    if not protein_ca:
        print(f"No CA atoms found in chain {protein_chain}.")
        return []

    results = []
    
    def get_gravity_center(residue):
        atoms = residue.atoms
        x = sum(a.x for a in atoms) / len(atoms)
        y = sum(a.y for a in atoms) / len(atoms)
        z = sum(a.z for a in atoms) / len(atoms)
        return (x, y, z)

    for ca in protein_ca:
        ca_coord = (ca.x, ca.y, ca.z)
        min_dist = float('inf')
        closest_het = None
        
        for res in het_residues:
            het_coord = get_gravity_center(res)
            dist = ((ca_coord[0]-het_coord[0])**2 + (ca_coord[1]-het_coord[1])**2 + (ca_coord[2]-het_coord[2])**2)**0.5
            if dist < min_dist:
                min_dist = dist
                closest_het = res
        
        results.append({
            'residue_index': ca.residue.index,
            'residue_num': ca.residue.num,
            'residue_name': ca.residue.name,
            'min_distance': min_dist,
            'closest_het': f"{closest_het.name}:{closest_het.num}:{closest_het.chain.name}" if closest_het else "N/A"
        })

    return results

def main():
    parser = argparse.ArgumentParser(description="Calculate distances between protein CA atoms and BLK gravity centers.")
    parser.add_argument("pdb_file", help="Path to the PDB file.")
    parser.add_argument("--protein_chain", default="A", help="Chain ID for the protein (default: A).")
    parser.add_argument("--blk_chain", default="B", help="Chain ID for the BLK residues (default: B).")
    parser.add_argument("--threshold", type=float, default=10.0, help="Output only residues closer than this threshold (default: 10.0).")
    
    args = parser.parse_args()

    if not os.path.exists(args.pdb_file):
        print(f"Error: File {args.pdb_file} not found.")
        sys.exit(1)

    print(f"Analyzing {args.pdb_file}...")
    print(f"Protein Chain: {args.protein_chain}, BLK Chain: {args.blk_chain}")
    
    distances = calculate_distances(args.pdb_file, args.protein_chain, args.blk_chain)
    
    if not distances:
        print("No distances calculated.")
        return

    # Sort by distance
    distances.sort(key=lambda x: x['min_distance'])

    print("\n" + "="*80)
    print(f"{'Residue':<15} | {'Min Distance (A)':<18} | {'Closest BLK Group':<20}")
    print("-" * 80)
    
    min_all = distances[0]['min_distance']
    
    for d in distances:
        if d['min_distance'] <= args.threshold:
            res_str = f"{d['residue_name']} {d['residue_num']}"
            print(f"{res_str:<15} | {d['min_distance']:<18.3f} | {d['closest_het']}")

    print("="*80)
    print(f"\nGLOBAL MINIMUM DISTANCE: {min_all:.3f} A")
    print(f"SUGGESTED BLOCK_REPULSION_RADIUS: {min_all:.1f} A (Example suggestion)")
    if min_all < 2.0:
        print('''
        WARNING: The minimum distance is less than 2.0 A. 
        Make sure that the BLK residues are not too close to the protein.
        BLK is going to be fixed during the modeling process, but the protein could move its relative position to the BLK.
        ''')
    print("="*80)

if __name__ == "__main__":
    main()
