#!/usr/bin/env python3
"""
Renumbers residues in a PDB file sequentially, handles TER records,
and optionally converts ATOM records to HETATM for specified chains.
"""

import sys
from typing import Optional, Set, Union

def renumber_pdb_residues(pdb_file: str, hetatm_chains: Optional[Union[List[str], Set[str]]] = None):
    """
    Renumbers residues of a PDB file sequentially, handles TER records,
    and changes ATOM to HETATM in specified chains.

    Args:
        pdb_file (str): Path to the input PDB file.
        hetatm_chains (Set[str], optional): A set of chain identifiers (e.g., {'B', 'D'})
                                            that should be written as HETATM.
                                            All others will be ATOM.
                                            If None or empty, all are written as ATOM.
    """
    
    if hetatm_chains is None:
        hetatm_chains_set = set()
        hetatm_flag = ""
    elif isinstance(hetatm_chains, str):
        # Allow a single string to be passed
        hetatm_chains_set = {hetatm_chains}
        hetatm_flag = "_HETATM"
    else:
        hetatm_chains_set = set(hetatm_chains)
        hetatm_flag = "_HETATM" if hetatm_chains_set else ""
        
    # Construct the new filename: original_renum[_HETATM].pdb
    output_file = pdb_file.replace(".pdb", f"_renum{hetatm_flag}.pdb")

    with open(pdb_file, "r") as f_in, open(output_file, "w") as f_out:
        current_residue_id = None
        new_resnum = 0

        for line in f_in:
            record_name = line[0:6].strip()

            if record_name in ("ATOM", "HETATM"):
                # A unique residue ID is (resname, chain, old_resnum)
                res_id = (line[17:20].strip(), line[21].strip(), line[22:26].strip())  
                chain_id = res_id[1]

                # If this is a new residue, increment the counter
                if res_id != current_residue_id:
                    new_resnum += 1
                    current_residue_id = res_id

                # Change record type to HETATM if chain is in the specified set
                new_record_name = "HETATM" if chain_id in hetatm_chains_set else "ATOM  "
                
                # Replace record name (cols 1-6) and residue number (cols 23-26)
                new_line = new_record_name + line[6:22] + f"{new_resnum:4d}" + line[26:]
                f_out.write(new_line)

            elif record_name == "TER":
                # Include the TER record with the same index as the previous residue
                if new_resnum > 0:
                    new_line = line[:22] + f"{new_resnum:4d}" + line[26:]
                    f_out.write(new_line)
                else:
                    f_out.write(line)
            
            else:
                # Write other lines (HEADER, REMARK, etc.) as-is
                f_out.write(line)

    print(f"Renumbered file created: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python3 renumber_pdb.py <file.pdb> [HETATM_chains]")
        print("Example 1 (renumber only): python3 renumber_pdb.py 1abc.pdb")
        print("Example 2 (renumber & chain B to HETATM): python3 renumber_pdb.py 1abc.pdb B")
        print("Example 3 (renumber & chains B,D to HETATM): python3 renumber_pdb.py 1abc.pdb B,D")
        sys.exit(1)

    pdb_file_arg = sys.argv[1]
    
    if len(sys.argv) == 3:
        chains_input_arg = sys.argv[2]
        hetatm_chains_list = [chain.strip() for chain in chains_input_arg.split(',')]
        renumber_pdb_residues(pdb_file_arg, hetatm_chains_list)
    else:
        renumber_pdb_residues(pdb_file_arg)