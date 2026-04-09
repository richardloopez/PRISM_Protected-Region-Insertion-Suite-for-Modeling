#!/usr/bin/env python3
import argparse
import sys
import os
import json
import math
import shutil

# ====================================================================================
#                                 MATH HELPERS
# ====================================================================================

def vec_sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def vec_dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
def vec_cross(a, b):
    return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])
def vec_norm(v):
    l = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    return (v[0]/l, v[1]/l, v[2]/l) if l > 1e-10 else (0,0,0)

def get_local_frame(n_pos, ca_pos, c_pos):
    '''Computes an orthonormal basis from N, CA, C coordinates.'''
    # x-axis along CA->C
    v_x = vec_norm(vec_sub(c_pos, ca_pos))
    v_nc = vec_norm(vec_sub(n_pos, ca_pos))
    # z-axis is normal to the plane
    v_z = vec_norm(vec_cross(v_x, v_nc))
    # y-axis is perpendicular to x and z
    v_y = vec_norm(vec_cross(v_z, v_x))
    return v_x, v_y, v_z

# ====================================================================================
#                                   CLASSES
# ====================================================================================

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from pdb_utils import PDBAtom

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
    parser_prep.add_argument("ptm_args", nargs="*", help="Optional PTM arguments (e.g. posttranslational=2 E-1,E-2:A-265 ...)")

    # --- RETRO MODE ---
    parser_retro = subparsers.add_parser("retro", help="Restore original ligand info to PRISM output.")
    parser_retro.add_argument("prism_output_pdb", help="The output PDB from PRISM (with BLK ligands).")
    parser_retro.add_argument("log_file", help="The .json log file generated during the 'prep' stage.")
    parser_retro.add_argument("ptm_args", nargs="*", help="Optional remapping arguments (e.g. posttranslational=2 A-209-new-A-232 ...)")

    return parser.parse_args()

def format_atom_name_blk(counter):
    '''
    Generates X1, X2... atom names.
    '''
    return f"X{counter}"

def parse_ptm_string(ptm_str):
    '''Parses strings like E-1,E-2:A-265 into ([(E, 1), (E, 2)], (A, 265))'''
    if ':' not in ptm_str: return None
    ptm_part, attach_part = ptm_str.split(':')
    def parse_id(s):
        if '-' not in s: return None
        c, seq = s.split('-', 1)
        return (c.strip(), int(seq.strip()))
    ptms = [parse_id(s) for s in ptm_part.split(',') if parse_id(s)]
    attach = parse_id(attach_part)
    if not ptms or not attach: return None
    return ptms, attach

def parse_remap_string(s):
    '''Parses strings like A-209-new-A-232 into ((A, 209), (A, 232))'''
    if '-new-' not in s: return None
    parts = s.split('-new-')
    if len(parts) != 2: return None
    orig_s, model_s = parts
    def parse_id(ss):
        if '-' not in ss: return None
        c, seq = ss.split('-', 1)
        return (c.strip(), int(seq.strip()))
    o_id, m_id = parse_id(orig_s), parse_id(model_s)
    if not o_id or not m_id: return None
    return o_id, m_id

def copy_to_input(filename):
    '''Copies a file to the ../input/ directory if it exists.'''
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_dir = os.path.join(os.path.dirname(script_dir), "input")
    if os.path.exists(input_dir) and os.path.isdir(input_dir):
        try:
            shutil.copy2(filename, os.path.join(input_dir, os.path.basename(filename)))
            print(f" [COPY] Copied {filename} to {input_dir}")
        except Exception as e:
            print(f" [COPY] Warning: Could not copy {filename} to {input_dir}: {e}")
    else:
        local_input = os.path.join(os.getcwd(), "input")
        if os.path.exists(local_input) and os.path.isdir(local_input):
            try:
                shutil.copy2(filename, os.path.join(local_input, os.path.basename(filename)))
                print(f" [COPY] Copied {filename} to {local_input}")
            except Exception as e:
                print(f" [COPY] Warning: Could not copy {filename} to {local_input}: {e}")

# ====================================================================================
#                                   PREP LOGIC
# ====================================================================================

def run_prep(input_path, prot_chains_str, lig_chains_str, ptm_args=[]):
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
                atom = PDBAtom(line)
                chain = atom.chain_id
                if chain in prot_chains:
                    protein_atoms.append(atom)
                elif chain in lig_chains:
                    ligand_atoms.append(atom)
                else: 
                    pass

    # 2. Handle PTMs
    ptm_specs = []
    for arg in ptm_args:
        if arg.startswith("posttranslational="): continue
        spec = parse_ptm_string(arg)
        if spec: ptm_specs.append(spec)

    ptm_map = []
    all_atoms = protein_atoms + ligand_atoms
    full_atom_list = []
    with open(input_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                full_atom_list.append(PDBAtom(line))

    res_dict = {}
    for atom in full_atom_list:
        key = (atom.chain_id, atom.res_seq)
        if key not in res_dict: res_dict[key] = []
        res_dict[key].append(atom)

    for ptms, attach in ptm_specs:
        if attach not in res_dict:
            print(f"[PREP] Warning: Attachment residue {attach} not found. Skipping PTM.")
            continue
        
        attach_atoms = res_dict[attach]
        n_at = next((a for a in attach_atoms if a.name.strip() == 'N'), None)
        ca_at = next((a for a in attach_atoms if a.name.strip() == 'CA'), None)
        c_at = next((a for a in attach_atoms if a.name.strip() == 'C'), None)

        if not (n_at and ca_at and c_at):
            print(f"[PREP] Warning: N, CA, or C missing for {attach}. Cannot compute local frame.")
            continue
        
        # Calc Frame
        n_p, ca_p, c_p = (n_at.x, n_at.y, n_at.z), (ca_at.x, ca_at.y, ca_at.z), (c_at.x, c_at.y, c_at.z)
        vx, vy, vz = get_local_frame(n_p, ca_p, c_p)

        ptm_residue_data = []
        for ptm_res_id in ptms:
            if ptm_res_id not in res_dict:
                print(f"[PREP] Warning: PTM residue {ptm_res_id} not found.")
                continue
            
            for p_at in res_dict[ptm_res_id]:
                # Calc Relative
                rel_p = vec_sub((p_at.x, p_at.y, p_at.z), ca_p)
                rx = vec_dot(rel_p, vx)
                ry = vec_dot(rel_p, vy)
                rz = vec_dot(rel_p, vz)

                ptm_residue_data.append({
                    "name": p_at.name.strip(),
                    "res_name": p_at.res_name,
                    "res_seq": p_at.res_seq,
                    "chain_id": p_at.chain_id,
                    "element": p_at.element,
                    "record_type": p_at.record_type,
                    "rel_coords": (rx, ry, rz)
                })
                
                # REMOVE from output lists if present
                protein_atoms = [a for a in protein_atoms if not (a.chain_id == p_at.chain_id and a.res_seq == p_at.res_seq)]
                ligand_atoms = [a for a in ligand_atoms if not (a.chain_id == p_at.chain_id and a.res_seq == p_at.res_seq)]

        ptm_map.append({
            "orig_attach_chain": attach[0],
            "orig_attach_res_seq": attach[1],
            "ptm_atoms": ptm_residue_data
        })
        print(f"[PREP] PTM: {ptms} attached to {attach} -> Processed {len(ptm_residue_data)} atoms.")

    # 3. Prepare Data Structures
    new_chain_a_lines = []
    new_chain_b_lines = []
    
    log_data = {
        "original_filename": os.path.basename(input_path),
        "protein_map": [], 
        "ligand_map": {},
        "ptm_map": ptm_map
    }

    current_serial = 1
    
    # --- PROCESS PROTEIN (CHAIN A) ---
    current_res_seq = 0
    prev_id = None
    
    for atom in protein_atoms:
        curr_id = (atom.chain_id, atom.res_seq, atom.i_code)
        if curr_id != prev_id:
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
    prev_id = None
    atom_counter = 1
    log_res_seq = 0
    
    for atom in ligand_atoms:
        curr_id = (atom.chain_id, atom.res_seq, atom.i_code)
        if curr_id != prev_id:
            current_res_seq += 1
            log_res_seq += 1
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
        
        key = f"{log_res_seq}_{new_name}"
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

    copy_to_input(out_pdb)


# ====================================================================================
#                                   RETRO LOGIC
# ====================================================================================

def run_retro(model_path, log_path, extra_args=[]):
    '''
    Restores original ligand and protein information to a PRISM output PDB.
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
    ptm_map_data = log_data.get("ptm_map", [])

    # 1. Parse Extra Args (Remappings)
    model_to_orig_explicit = {}
    for arg in extra_args:
        if arg.startswith("posttranslational="): continue
        remap = parse_remap_string(arg)
        if remap:
            o_id, m_id = remap
            model_to_orig_explicit[m_id] = o_id

    # 2. Build Lookups from Log
    prot_id_lookup = {}
    for p in log_data.get("protein_map", []):
        seq = p['new_res_seq']
        if seq not in prot_id_lookup:
            prot_id_lookup[seq] = (p['orig_chain'], p['orig_res_seq'], p['orig_res_name'])
            
    found_attachments = {}
    target_orig_ids = set((p['orig_attach_chain'], p['orig_attach_res_seq']) for p in ptm_map_data)

    restored_lines = []
    log_res_seq = 0
    prev_id = None
    current_serial = 1
    
    with open(model_path, 'r') as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                # Keep other header/footer lines but we will handle TER/END manually
                lstrip = line.strip()
                if lstrip and not lstrip.startswith(('TER', 'END', 'CONECT')):
                    restored_lines.append(lstrip)
                continue

            atom = PDBAtom(line)
            curr_id = (atom.chain_id, atom.res_seq, atom.i_code)
            
            # Identify the original ID of this residue for logic, but don't change modeled IDs
            orig_id = None
            m_id = (atom.chain_id, atom.res_seq)
            
            if m_id in model_to_orig_explicit:
                orig_id = model_to_orig_explicit[m_id]
            elif atom.chain_id == 'A':
                if atom.res_seq in prot_id_lookup:
                    info = prot_id_lookup[atom.res_seq]
                    orig_id = (info[0], info[1])
            elif atom.chain_id == 'B':
                # Restore ligands as they are usually treated as separate from the modeled protein backbone
                if curr_id != prev_id:
                    log_res_seq += 1
                key = f"{log_res_seq}_{atom.name.strip()}"
                if key in ligand_map:
                    info = ligand_map[key]
                    orig_id = (info['orig_chain'], info['orig_res_seq'])
                    atom.chain_id = info['orig_chain']
                    atom.res_name = info['orig_res_name']
                    atom.res_seq = info['orig_res_seq']
                    atom.name = info['orig_atom_name']
                    atom.element = info['orig_element']
                    atom.record_type = f"{info['orig_record']:<6}"
                    atom.temp = 0.00

            if False: pass
            if orig_id in target_orig_ids:
                if orig_id not in found_attachments: found_attachments[orig_id] = {}
                a_name = atom.name.strip()
                if a_name in ['N', 'CA', 'C']:
                    found_attachments[orig_id][a_name] = atom
            
            # Update serial
            atom.serial = current_serial
            restored_lines.append(atom.to_pdb_line())
            current_serial += 1
            prev_id = curr_id

    # --- RECONSTRUCT PDB ---
    final_output = []
    # Add headers
    for l in restored_lines:
        if not l.startswith(('ATOM', 'HETATM')):
            final_output.append(l)

    # Add Atoms by Chain
    all_atoms = [PDBAtom(l) for l in restored_lines if l.startswith(('ATOM', 'HETATM'))]
    
    # Process PTMs and add to all_atoms
    ptm_atoms = []
    if ptm_map_data:
        for ptm_spec in ptm_map_data:
            orig_attach_id = (ptm_spec['orig_attach_chain'], ptm_spec['orig_attach_res_seq'])
            if orig_attach_id not in found_attachments:
                print(f"[RETRO] Warning: Could not find attachment {orig_attach_id} in model.")
                continue
            
            bits = found_attachments[orig_attach_id]
            if not ('N' in bits and 'CA' in bits and 'C' in bits):
                continue
            
            n_at, ca_at, c_at = bits['N'], bits['CA'], bits['C']
            vx, vy, vz = get_local_frame((n_at.x, n_at.y, n_at.z), (ca_at.x, ca_at.y, ca_at.z), (c_at.x, c_at.y, c_at.z))
            
            for p_at_info in ptm_spec['ptm_atoms']:
                rx, ry, rz = p_at_info['rel_coords']
                new_x = ca_at.x + rx*vx[0] + ry*vy[0] + rz*vz[0]
                new_y = ca_at.y + rx*vx[1] + ry*vy[1] + rz*vz[1]
                new_z = ca_at.z + rx*vx[2] + ry*vy[2] + rz*vz[2]

                dummy_line = f"{'HETATM':<6}{0:>5} {'X':<4} {'PTM':>3} {'C'}{0:>4}    {0.0:>8.3f}{0.0:>8.3f}{0.0:>8.3f}{1.0:>6.2f}{0.0:>6.2f}"
                new_atom = PDBAtom(dummy_line)
                new_atom.record_type = f"{p_at_info['record_type']:<6}"
                new_atom.name = p_at_info['name']
                new_atom.res_name = p_at_info['res_name']
                new_atom.chain_id = p_at_info['chain_id']
                new_atom.res_seq = p_at_info['res_seq']
                new_atom.x, new_atom.y, new_atom.z = new_x, new_y, new_z
                new_atom.element = p_at_info['element']
                ptm_atoms.append(new_atom)

    full_atom_list = all_atoms + ptm_atoms
    
    chains_in_order = []
    for a in full_atom_list:
        if a.chain_id not in chains_in_order:
            chains_in_order.append(a.chain_id)
    
    serial = 1
    final_atom_lines = []
    for c in sorted(chains_in_order):
        for a in full_atom_list:
            if a.chain_id == c:
                a.serial = serial
                final_atom_lines.append(a.to_pdb_line())
                serial += 1
        final_atom_lines.append("TER")

    final_output.extend(final_atom_lines)
    final_output.append("END")

    # Output
    base_model = os.path.splitext(os.path.basename(model_path))[0]
    out_name = f"{base_model}_restored.pdb"
    
    with open(out_name, 'w') as f:
        for line in final_output:
            f.write(line + "\n")
    
    print(f"[RETRO] Success! Restored file saved as: {out_name}")


# ====================================================================================
#                                      MAIN
# ====================================================================================

if __name__ == "__main__":
    args = parse_args()
    
    if args.mode == "prep":
        run_prep(args.input_pdb, args.protein_chains, args.ligand_chains, args.ptm_args)
    elif args.mode == "retro":
        run_retro(args.prism_output_pdb, args.log_file, args.ptm_args)

