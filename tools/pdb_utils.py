#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM Tool: pdb_utils
#
'''
Shared tools module for PRISM. Includes PDBAtom class for parsing PDB records.
'''

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
