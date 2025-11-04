#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM (Protected-Region Insertion Suite for Modeling)
#
# Citation:
# If you use this software in your research, please cite:
# Lopez-Corbalan, R.

"""
Custom MODELLER model classes with fixed-region support.

This module defines specialized AutoModel and LoopModel classes that:
1. Keep experimental coordinates completely fixed during optimization
2. Apply HETATM repulsion shields to prevent loop clashes with bound molecules
3. Preserve coordinate integrity while allowing targeted loop refinement
"""

from modeller import *
from modeller.automodel import AutoModel, DOPEHRLoopModel
from modeller.selection import Selection
from typing import Set

from .config import MIN_DIST_FROM_NUCLEOTIDE_COM


class FixedRegionAutoModel(AutoModel):
    """
    Custom AutoModel that maintains fixed experimental coordinates from the main template.

    This class prevents optimization of experimental residues and applies a
    per-nucleotide repulsion shield to avoid clashes with bound molecules (DNA/RNA).

    Attributes:
        experimental_residues: Set of residue numbers that must remain fixed
        chain_id: Chain identifier for residue selection
    """

    def __init__(self, env, experimental_residues: Set[int], chain_id: str = 'A', **kwargs):
        """
        Initialize FixedRegionAutoModel with protected experimental residues.

        Args:
            env: MODELLER environment
            experimental_residues: Set of residue numbers (1-indexed) to keep fixed
            chain_id: Chain identifier (default: 'A')
            **kwargs: Additional arguments passed to AutoModel
        """
        super().__init__(env, **kwargs)
        self.experimental_residues = experimental_residues
        self.chain_id = chain_id
        print(f"\n[FixedRegionAutoModel] Initialized with {len(experimental_residues)} FIXED experimental residues")

    def select_atoms(self):
        """
        Define which atoms will be optimized during modeling.

        EXCLUDES experimental residues from the main template to prevent
        coordinate drift during optimization.

        Returns:
            Selection of atoms to optimize (all standard residues minus experimental)
        """
        if not self.experimental_residues:
            return Selection(self).only_std_residues()

        all_residues = Selection(self).only_std_residues()

        fixed_selection = Selection()
        for res_num in sorted(self.experimental_residues):
            try:
                res_range = self.residue_range(f'{res_num}:{self.chain_id}', f'{res_num}:{self.chain_id}')
                fixed_selection = fixed_selection | Selection(res_range)
            except Exception:
                continue

        optimizable_atoms = all_residues - fixed_selection
        print(f"[FixedRegionAutoModel] Optimizable atoms: {len(optimizable_atoms)} (excluding {len(self.experimental_residues)} experimental residues)")
        return optimizable_atoms

    def special_restraints(self, aln):
        """
        Add special restraints to the model.

        Note: Experimental core freezing is achieved via select_atoms() only.
        Default VdW restraints remain active.

        Args:
            aln: Alignment object
        """
        super().nonstd_restraints(aln)

        print("[FixedRegionAutoModel] Default VdW restraints active.")
        print("[FixedRegionAutoModel] Experimental core frozen exclusively via select_atoms().")

    def nonstd_restraints(self, aln):
        """
        Add HETATM repulsion shield restraints after standard restraints.

        Creates a per-nucleotide repulsion shield that prevents loop regions
        from clashing with bound molecules (DNA, RNA, ligands).

        Args:
            aln: Alignment object
        """
        super().nonstd_restraints(aln)
        rsr = self.restraints
        print("\n[FixedRegionAutoModel] Executing nonstd_restraints to add HETATM shield...")

        print("\n[FixedRegionAutoModel] Searching for HETATM (DNA/RNA) for PER-NUCLEOTIDE repulsion shield...")

        het_residue_list = [res for res in self.residues if res.hetatm and res.name != 'HOH']

        if not het_residue_list:
            print("[FixedRegionAutoModel] No HETATM residues (excl. HOH) found in nonstd_restraints. No repulsion restraints added.")
            return

        num_het_residues = len(het_residue_list)
        print(f"[FixedRegionAutoModel] {num_het_residues} HETATM residues detected (e.g., nucleotides) in nonstd_restraints.")

        all_std_residues = Selection(self).only_std_residues()
        fixed_selection = Selection()
        if self.experimental_residues:
            for res_num in self.experimental_residues:
                try:
                    res_range = self.residue_range(f'{res_num}:{self.chain_id}', f'{res_num}:{self.chain_id}')
                    fixed_selection = fixed_selection | Selection(res_range)
                except Exception:
                    continue

        optimizable_atoms = all_std_residues - fixed_selection
        optimizable_ca_atoms = optimizable_atoms.only_atom_types('CA')

        num_ca_atoms = len(optimizable_ca_atoms)
        if not optimizable_ca_atoms:
            print("[FixedRegionAutoModel] No optimizable CA atoms found in nonstd_restraints. No restraints added.")
            return

        print(f"[FixedRegionAutoModel] Preparing to add {num_ca_atoms} (C-alphas) x {num_het_residues} (HET res) = {num_ca_atoms * num_het_residues} restraints.")

        total_repulsion_count = 0
        stdev = 1.0

        het_pseudo_atoms = []
        for res in het_residue_list:
            het_atoms_in_res = Selection(res)
            pseudo = pseudo_atom.GravityCenter(het_atoms_in_res)
            rsr.pseudo_atoms.append(pseudo)
            het_pseudo_atoms.append(pseudo)
        print(f"[FixedRegionAutoModel] Created and added {len(het_pseudo_atoms)} GravityCenter pseudo atoms.")

        for ca_atom in optimizable_ca_atoms:
            for het_pseudo_atom in het_pseudo_atoms:
                rsr.add(forms.LowerBound(
                    group=physical.xy_distance,
                    feature=features.Distance(ca_atom, het_pseudo_atom),
                    mean=MIN_DIST_FROM_NUCLEOTIDE_COM,
                    stdev=stdev
                ))
                total_repulsion_count += 1

        print(f"[FixedRegionAutoModel] Added {total_repulsion_count} repulsion restraints (CA -> Nucleotide-PseudoAtom) in nonstd_restraints")
        print(f"                      (Minimum distance: {MIN_DIST_FROM_NUCLEOTIDE_COM} Å, σ={stdev})")


class FixedRegionLoopRefiner(DOPEHRLoopModel):
    """
    Custom DOPEHRLoopModel that NEVER refines experimental residues.

    This class ensures loop refinement only occurs in flexible regions
    while applying HETATM repulsion shields to prevent clashes.

    Attributes:
        loop_start: First residue of the loop to refine
        loop_end: Last residue of the loop to refine
        chain_id: Chain identifier
        experimental_residues: Set of residue numbers that must remain fixed
    """

    def __init__(self, env, inimodel, sequence, loop_start, loop_end, chain_id,
                 experimental_residues: Set[int], **kwargs):
        """
        Initialize FixedRegionLoopRefiner with validation.

        Args:
            env: MODELLER environment
            inimodel: Initial model PDB file
            sequence: Sequence code
            loop_start: First residue of loop (1-indexed)
            loop_end: Last residue of loop (1-indexed)
            chain_id: Chain identifier
            experimental_residues: Set of residue numbers to keep fixed
            **kwargs: Additional arguments passed to DOPEHRLoopModel

        Raises:
            ValueError: If loop contains any experimental residues
        """
        super().__init__(env,
                         inimodel=inimodel,
                         sequence=sequence,
                         **kwargs)
        self.loop_start = loop_start
        self.loop_end = loop_end
        self.chain_id = chain_id
        self.experimental_residues = experimental_residues

        loop_residues = set(range(loop_start, loop_end + 1))
        overlap = loop_residues.intersection(experimental_residues)

        if overlap:
            raise ValueError(
                f"[ERROR] Loop [{loop_start}-{loop_end}] contains experimental residues {sorted(overlap)}. "
                f"This loop must NOT be refined."
            )

        print(f"[FixedRegionLoopRefiner] Loop [{loop_start}-{loop_end}] validated: NO experimental residues")

    def select_loop_atoms(self):
        """
        Define the loop residues that will be refined.

        Returns:
            Selection of atoms in the loop region
        """
        range_start = f'{self.loop_start}:{self.chain_id}'
        range_end = f'{self.loop_end}:{self.chain_id}'
        return Selection(self.residue_range(range_start, range_end))

    def special_restraints(self, aln):
        """
        Special restraints placeholder.

        HETATM restraints are now handled in nonstd_restraints.

        Args:
            aln: Alignment object
        """
        pass

    def nonstd_restraints(self, aln):
        """
        Apply HETATM repulsion shield during loop refinement.

        Creates per-nucleotide repulsion restraints to prevent the loop
        from clashing with bound molecules.

        Args:
            aln: Alignment object
        """
        super().nonstd_restraints(aln)
        rsr = self.restraints

        print(f"\n[FixedRegionLoopRefiner] Searching for HETATM (DNA/RNA) for PER-NUCLEOTIDE shield (Loop {self.loop_start}-{self.loop_end})...")

        het_residue_list = [res for res in self.residues if res.hetatm and res.name != 'HOH']

        if not het_residue_list:
            print(f"[FixedRegionLoopRefiner] No HETATM residues (excl. HOH) found in nonstd_restraints.")
            return

        num_het_residues = len(het_residue_list)

        loop_selection = self.select_loop_atoms()
        loop_ca_atoms = loop_selection.only_atom_types('CA')

        num_ca_atoms = len(loop_ca_atoms)
        if not loop_ca_atoms:
            print(f"[FixedRegionLoopRefiner] No CA atoms found in loop.")
            return

        print(f"[FixedRegionLoopRefiner] Preparing to add {num_ca_atoms} (C-alphas) x {num_het_residues} (HET res) = {num_ca_atoms * num_het_residues} restraints.")

        het_pseudo_atoms = []
        for res in het_residue_list:
            het_atoms_in_res = Selection(res)
            pseudo = pseudo_atom.GravityCenter(het_atoms_in_res)
            rsr.pseudo_atoms.append(pseudo)
            het_pseudo_atoms.append(pseudo)
        print(f"[FixedRegionLoopRefiner] Created and added {len(het_pseudo_atoms)} GravityCenter pseudo atoms.")

        total_repulsion_count = 0
        stdev = 1.0

        for ca_atom in loop_ca_atoms:
            for het_pseudo_atom in het_pseudo_atoms:
                rsr.add(forms.LowerBound(
                    group=physical.xy_distance,
                    feature=features.Distance(ca_atom, het_pseudo_atom),
                    mean=MIN_DIST_FROM_NUCLEOTIDE_COM,
                    stdev=stdev
                ))
                total_repulsion_count += 1

        print(f"[FixedRegionLoopRefiner] Added {total_repulsion_count} repulsion restraints (CA -> Nucleotide-PseudoAtom) in nonstd_restraints")
        print(f"                       (Minimum distance: {MIN_DIST_FROM_NUCLEOTIDE_COM} Å, σ={stdev})")
