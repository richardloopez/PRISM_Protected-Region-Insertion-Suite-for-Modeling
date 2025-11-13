# Script to remove residues ("gaps", "-") from a PDB chain according to a reference sequence in one-letter code from a TXT file.
# Requires Biopython (`pip install biopython`).

from Bio.PDB import PDBParser, PDBIO
from Bio.SeqUtils import seq1

def get_chain_sequence(pdb_file, chain_id):
    """Returns the sequence (1-letter code) and a list of residue objects from a PDB chain."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb_file)
    chain = None
    for model in structure:
        chain = model[chain_id]
        break
    seq = ""
    residues = []
    for residue in chain.get_residues():
        if residue.get_id()[0] == " ":  # Standard residues only
            aa_one = seq1(residue.get_resname())
            seq += aa_one
            residues.append(residue)
    return seq, residues

def filter_pdb_by_txt(pdb_file, txt_file, chain_id="A", output_pdb="output.pdb"):
    # Read target sequence (TXT with "-" for gaps and one-letter codes for residues)
    with open(txt_file) as f:
        target_sequence = f.read().replace('\n', '').strip()

    pdb_sequence, pdb_residues = get_chain_sequence(pdb_file, chain_id)

    if len(pdb_sequence) != len(target_sequence):
        raise ValueError(f"PDB sequence length ({len(pdb_sequence)}) and TXT ({len(target_sequence)}) length differ.")

    # Compare residue at each position; ignore where TXT is "-"
    for i, (pdb_aa, txt_aa) in enumerate(zip(pdb_sequence, target_sequence)):
        if txt_aa != '-' and txt_aa != pdb_aa:
            raise ValueError(f"Mismatch at position {i+1}: PDB={pdb_aa}, TXT={txt_aa}")

    # Identify indices to keep (where txt_aa != "-")
    keep_indices = [i for i, aa in enumerate(target_sequence) if aa != "-"]

    # Build new structure with kept residues
    from Bio.PDB import Structure, Model, Chain, Residue, Atom

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb_file)
    new_structure = Structure.Structure("filtered")
    model = Model.Model(0)
    new_chain = Chain.Chain(chain_id)

    for i in keep_indices:
        residue = pdb_residues[i]
        # Clone residue to avoid references to the old structure
        new_residue = Residue.Residue(residue.get_id(), residue.get_resname(), residue.get_segid())
        for atom in residue.get_atoms():
            new_atom = Atom.Atom(atom.get_name(), atom.coord, atom.bfactor, atom.occupancy,
                                 atom.altloc, atom.fullname, atom.serial_number, atom.element)
            new_residue.add(new_atom)
        new_chain.add(new_residue)

    model.add(new_chain)
    new_structure.add(model)

    # Save new PDB file
    io = PDBIO()
    io.set_structure(new_structure)
    io.save(output_pdb)
    print(f"Filtered PDB saved as {output_pdb}")

# Example usage:
pdb_file = "original.pdb"      # Your original PDB file (change to your filename)
txt_file = "final_chain.txt"   # Your TXT file (with '-' and one-letter residue codes)
filter_pdb_by_txt(pdb_file, txt_file, chain_id="A", output_pdb="output.pdb")
