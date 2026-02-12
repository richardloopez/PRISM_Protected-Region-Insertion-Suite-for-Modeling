#!/usr/bin/env python3
import argparse
import sys
import os
import json

# ====================================================================================
#                                   CLASSES
# ====================================================================================

class PDBAtom:
    def __init__(self, line):
        self.line = line
        self.record_type = line[0:6].strip()
        try:
            self.serial = int(line[6:11])
        except ValueError:
            self.serial = 0
        self.name = line[12:16] 
        self.alt_loc = line[16]
        self.res_name = line[17:20].strip()
        self.chain_id = line[21]
        try:
            self.res_seq = int(line[22:26])
        except ValueError:
            self.res_seq = 0
        self.i_code = line[26]
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occ = float(line[54:60]) if len(line) > 54 and line[54:60].strip() else 1.00
        self.temp = float(line[60:66]) if len(line) > 60 and line[60:66].strip() else 0.00
        self.element = line[76:78].strip() if len(line) > 76 else ""

    def to_pdb_line(self):
        '''
        Generates a PDB formatted line from current attributes.
        '''
        if len(self.name) == 4:
            name_str = f"{self.name}"
        else:
            if self.name[0].isdigit():
                name_str = f"{self.name:<4}"
            else:
                name_str = f" {self.name:<3}"

        return (f"{self.record_type:<6}{self.serial:>5} {name_str:4}{self.alt_loc}{self.res_name:>3} {self.chain_id}{self.res_seq:>4}{self.i_code}   "
                f"{self.x:>8.3f}{self.y:>8.3f}{self.z:>8.3f}{self.occ:>6.2f}{self.temp:>6.2f}          {self.element:>2}")

# ====================================================================================
#                                 HELPER FUNCTIONS
# ====================================================================================

def parse_args():
    parser = argparse.ArgumentParser(description="PRISM PDB Pre/Post-Processor.")
    subparsers = parser.add_subparsers(dest="mode", required=True, help="Mode of operation")

    # --- PREP MODE ---
    parser_prep = subparsers.add_parser("prep", help="Prepare PDB for PRISM (Generate Chain A/B and Log).")
    parser_prep.add_argument("input_pdb", help="Original input PDB file.")
    parser_prep.add_argument("protein_chains", help="Chains to merge into Chain A (e.g. 'A,C').")
    parser_prep.add_argument("ligand_chains", help="Chains to convert to BLK Chain B (e.g. 'B,D').")

    # --- RETRO MODE ---
    parser_retro = subparsers.add_parser("retro", help="Restore original ligand info to PRISM output.")
    parser_retro.add_argument("prism_output_pdb", help="The output PDB from PRISM (with BLK ligands).")
    parser_retro.add_argument("original_pdb", help="The original PDB file (used for reference/verification).")
    parser_retro.add_argument("log_file", help="The .json log file generated during the 'prep' stage.")

    return parser.parse_args()

def format_atom_name_blk(counter):
    '''
    Generates X1, X2... atom names.
    '''
    return f"X{counter}"

# ====================================================================================
#                                   PREP LOGIC
# ====================================================================================

def run_prep(input_path, prot_chains_str, lig_chains_str):
    '''
    Prepares a PDB file for PRISM by splitting chains and generating a log.
    '''
    if not os.path.exists(input_path):
        sys.exit(f"Error: File {input_path} not found.")

    prot_chains = [c.strip() for c in prot_chains_str.split(',')]
    lig_chains = [c.strip() for c in lig_chains_str.split(',')]

    print(f"[PREP] Processing {input_path}")
    print(f"[PREP] Protein Chains -> A: {prot_chains}")
    print(f"[PREP] Ligand Chains  -> B: {lig_chains}")

    protein_atoms = []
    ligand_atoms = []

    # 1. Read Atoms
    with open(input_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                chain = line[21]
                if chain in prot_chains:
                    protein_atoms.append(PDBAtom(line))
                elif chain in lig_chains:
                    ligand_atoms.append(PDBAtom(line))

    # 2. Prepare Data Structures
    new_chain_a_lines = []
    new_chain_b_lines = []
    
    # The Log Dictionary
    log_data = {
        "original_filename": os.path.basename(input_path),
        "protein_map": [], 
        "ligand_map": {}
    }

    current_serial = 1
    
    # --- PROCESS PROTEIN (CHAIN A) ---
    current_res_seq = 1
    prev_id = None
    
    for atom in protein_atoms:
        curr_id = (atom.chain_id, atom.res_seq, atom.i_code)
        if prev_id is not None and curr_id != prev_id:
            current_res_seq += 1
        
        new_atom = PDBAtom(atom.line)
        new_atom.serial = current_serial
        new_atom.chain_id = 'A'
        new_atom.res_seq = current_res_seq
        new_atom.record_type = "ATOM  "
        
        log_data["protein_map"].append({
            "new_chain": "A",
            "new_res_seq": current_res_seq,
            "new_atom_name": new_atom.name.strip(),
            "orig_chain": atom.chain_id,
            "orig_res_name": atom.res_name,
            "orig_res_seq": atom.res_seq,
            "orig_atom_name": atom.name.strip()
        })
        
        new_chain_a_lines.append(new_atom.to_pdb_line())
        current_serial += 1
        prev_id = curr_id

    # --- PROCESS LIGAND (CHAIN B) ---
    current_res_seq = 1
    prev_id = None
    atom_counter = 1
    
    for atom in ligand_atoms:
        curr_id = (atom.chain_id, atom.res_seq, atom.i_code)
        if prev_id is not None and curr_id != prev_id:
            current_res_seq += 1
            atom_counter = 1
        
        new_name = format_atom_name_blk(atom_counter)
        
        new_atom = PDBAtom(atom.line)
        new_atom.serial = current_serial
        new_atom.chain_id = 'B'
        new_atom.res_name = 'BLK'
        new_atom.res_seq = current_res_seq
        new_atom.name = new_name
        new_atom.element = "X"
        new_atom.temp = 99.99
        new_atom.record_type = "HETATM"
        
        key = f"{current_res_seq}_{new_name}"
        log_data["ligand_map"][key] = {
            "orig_chain": atom.chain_id,
            "orig_res_name": atom.res_name,
            "orig_res_seq": atom.res_seq,
            "orig_atom_name": atom.name.strip(),
            "orig_element": atom.element,
            "orig_record": atom.record_type.strip()
        }

        new_chain_b_lines.append(new_atom.to_pdb_line())
        current_serial += 1
        atom_counter += 1
        prev_id = curr_id

    # 3. Output
    base, _ = os.path.splitext(os.path.basename(input_path))
    out_pdb = f"{base}_prism_prep.pdb"
    out_log = f"{base}_prism_data.json"

    with open(out_pdb, 'w') as f:
        for l in new_chain_a_lines: f.write(l + "\n")
        if new_chain_a_lines: f.write("TER\n")
        for l in new_chain_b_lines: f.write(l + "\n")
        if new_chain_b_lines: f.write("TER\n")
        f.write("END\n")

    with open(out_log, 'w') as f:
        json.dump(log_data, f, indent=4)

    print(f"[PREP] Success!")
    print(f" > Generated PDB: {out_pdb}")
    print(f" > Generated Log: {out_log}")


# ====================================================================================
#                                   RETRO LOGIC
# ====================================================================================

def run_retro(model_path, original_path, log_path):
    '''
    Restores original ligand information to a PRISM output PDB.
    '''
    if not os.path.exists(model_path):
        sys.exit(f"Error: Model file {model_path} not found.")
    if not os.path.exists(log_path):
        sys.exit(f"Error: Log file {log_path} not found.")
    
    print(f"[RETRO] Restoring original ligand info...")
    print(f" > Model: {model_path}")
    print(f" > Log:   {log_path}")

    with open(log_path, 'r') as f:
        log_data = json.load(f)
    
    ligand_map = log_data.get("ligand_map", {})

    restored_lines = []
    
    with open(model_path, 'r') as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                restored_lines.append(line.strip())
                continue

            atom = PDBAtom(line)

            # --- RESTORE LIGAND (Chain B) ---
            if atom.chain_id == 'B':
                key = f"{atom.res_seq}_{atom.name.strip()}"
                
                if key in ligand_map:
                    info = ligand_map[key]
                    
                    # Restore Identity
                    atom.chain_id = info['orig_chain']
                    atom.res_name = info['orig_res_name']
                    atom.res_seq = info['orig_res_seq']
                    atom.name = info['orig_atom_name']
                    atom.element = info['orig_element']
                    atom.record_type = f"{info['orig_record']:<6}"
                    
                    atom.temp = 0.00 
                else:
                    print(f"[RETRO] Warning: No map found for Chain B atom {key}. Keeping as is.")

            # --- PROTEIN (Chain A) ---
            # Only for ligand_chains since protein chains could have changed
            # We keep Chain A as is (output from Modeller).
            
            restored_lines.append(atom.to_pdb_line())

    # Output
    base_model = os.path.splitext(os.path.basename(model_path))[0]
    out_name = f"{base_model}_restored.pdb"
    
    with open(out_name, 'w') as f:
        for line in restored_lines:
            f.write(line + "\n")
    
    print(f"[RETRO] Success! Restored file saved as: {out_name}")


# ====================================================================================
#                                      MAIN
# ====================================================================================

if __name__ == "__main__":
    args = parse_args()
    
    if args.mode == "prep":
        run_prep(args.input_pdb, args.protein_chains, args.ligand_chains)
    elif args.mode == "retro":
        run_retro(args.prism_output_pdb, args.original_pdb, args.log_file)