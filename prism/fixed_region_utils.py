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
Utilities for identifying and managing experimental (fixed) regions.

This module provides functions to identify residues that correspond to
experimental coordinates in the template structure, which must remain
fixed during modeling to prevent coordinate drift.
"""

from typing import List, Tuple, Set


def identify_experimental_residues(aligned_template_seq: str, aligned_target_seq: str) -> Set[int]:
    """
    Identify target residues that map to experimental regions in the main template.

    These residues correspond to positions where both the template and target
    have actual residues (not gaps), indicating experimental coordinates that
    must be preserved during modeling.

    Args:
        aligned_template_seq: Aligned sequence from the main template (MAIN_PDB)
        aligned_target_seq: Aligned sequence from the target protein

    Returns:
        Set of residue numbers (1-indexed) in the target sequence that map
        to experimental positions in the template (not gaps in either sequence)

    Raises:
        ValueError: If aligned sequences have different lengths
    """
    experimental_residues = set()
    target_res_num = 0

    if len(aligned_template_seq) != len(aligned_target_seq):
        print("[ERROR] Aligned sequence lengths do not match in identify_experimental_residues.")
        raise ValueError("Aligned template and target sequences must have equal lengths.")

    for template_char, target_char in zip(aligned_template_seq, aligned_target_seq):

        if template_char == '.' or template_char == '/':
            continue

        if target_char != '-':
            target_res_num += 1

        if template_char != '-' and target_char != '-':
            experimental_residues.add(target_res_num)

    print(f"\n[FIXED_REGION] Identified {len(experimental_residues)} experimental residues mapped from template")
    print(f"[FIXED_REGION] These residues will NOT be optimized or refined")

    return experimental_residues


def filter_loops_exclude_experimental(loop_ranges: List[Tuple[int, int]],
                                     experimental_residues: Set[int]) -> List[Tuple[int, int]]:
    """
    Filter loop ranges to exclude any loops containing experimental residues.

    This function ensures that no loop refinement is attempted on regions
    that contain experimental coordinates from the main template.

    Args:
        loop_ranges: List of (start, end) tuples defining loop regions
        experimental_residues: Set of residue numbers that must not be modified

    Returns:
        Filtered list of loop ranges that do not contain experimental residues
    """
    filtered_ranges = []

    for start, end in loop_ranges:
        loop_residues = set(range(start, end + 1))

        if loop_residues.intersection(experimental_residues):
            print(f"[FIXED_REGION] ⚠️  Loop [{start}-{end}] EXCLUDED: Contains experimental residues from MAIN_PDB")
            continue

        filtered_ranges.append((start, end))

    if filtered_ranges:
        range_strings = [f"[{s}-{e}]" if s != e else f"[{s}]" for s, e in filtered_ranges]
        print(f"\n[FIXED_REGION] ✓ Valid loops (without experimental residues): {range_strings}")
    else:
        print(f"\n[FIXED_REGION] ℹ️  No valid loops to refine (all contain experimental residues)")

    return filtered_ranges
