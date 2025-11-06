#!/usr/bin/env python3

"""
replicate_template.py

This script automates the "trick" of replicating the experimental template
(MAIN_PDB) N times to give it more weight in AutoModel's averaging.

It generates:

1. N copies of the PDB file in the 'input/' directory.
2. A 'config_update.txt' with the code block to update 'config.py'.
3. An updated 'manual_template_FullSeq.ali' in the current (tools/)
   directory, which you must then move to 'input/'.

Usage (from the 'tools/' folder):

python3 replicate_template.py ../input/9CB5_DS_renum_HETATM-B.pdb 1000

Arguments:

1. Path to the MAIN_PDB (e.g., ../input/my_pdb.pdb)
2. Number of replicas to generate (e.g., 1000)
"""

import sys
import os
import shutil
from typing import Tuple, List

def replicate_pdb_files(main_pdb_path: str, pdb_base_name: str, pdb_ext: str, num_replicas: int) -> Tuple[str, List[str]]:
    """
    Create N copies of the MAIN_PDB in its original directory.
    Returns the PDB directory and the list of replica names.
    """
    pdb_dir = os.path.dirname(main_pdb_path)
    replica_names = []
    print(f"\n--- 1. Replicating PDB files in {pdb_dir} ---")
    for i in range(1, num_replicas + 1):
        # Naming: base-1.pdb, base-2.pdb
        replica_name = f"{pdb_base_name}-{i}{pdb_ext}"
        replica_path = os.path.join(pdb_dir, replica_name)
        try:
            shutil.copy2(main_pdb_path, replica_path)
            replica_names.append(replica_name)
        except Exception as e:
            print(f" Error copying to {replica_path}: {e}")
    print(f" ✓ Created {len(replica_names)} PDB replicas.")
    return pdb_dir, replica_names

def update_config_file(replica_names: List[str]):
    """
    Reads ../prism/config.py and generates a 'config_update.txt' in the CWD
    with the updated PDB_TEMPLATE_FILES_NAMES list, ensuring proper commas.
    """
    config_path = os.path.join("..", "prism", "config.py")
    output_path = "config_update.txt"
    print(f"\n--- 2. Generating {output_path} ---")
    try:
        with open(config_path, 'r') as f:
            lines = f.readlines()

        new_config_block = []
        in_list = False
        list_lines = []
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("PDB_TEMPLATE_FILES_NAMES"):
                in_list = True
                new_config_block.append(line)
            elif in_list:
                if "]" in stripped:
                    # Process the collected list lines
                    # Remove trailing newlines and spaces
                    original_items = [l.rstrip('\n') for l in list_lines]

                    # Identify if the last item has a trailing comma or not
                    if original_items:
                        last_line = original_items[-1]
                        last_has_comma = last_line.endswith(",")
                    else:
                        last_has_comma = False

                    # Remove closing bracket line from list_lines
                    closing_bracket_line = line

                    # Clean and strip quotes from original items
                    cleaned_items = []
                    for item in original_items:
                        # Expect item format like "'filename'," or '"filename",'
                        cleaned = item.strip().rstrip(',')
                        # strip quotes
                        if (cleaned.startswith('"') and cleaned.endswith('"')) or (cleaned.startswith("'") and cleaned.endswith("'")):
                            cleaned = cleaned[1:-1]
                        cleaned_items.append(cleaned)

                    # Construct new lines: keep original items as is
                    # Append replicas after original items, with single commas except last item if originally with no comma, preserve that

                    # Prepare output lines list
                    output_lines = []
                    # Add original items back with their original commas except for last (handle carefully)
                    for i, orig_line in enumerate(original_items):
                        # Add original line as is to preserve commas
                        output_lines.append(orig_line + '\n')

                    # Append replicas
                    # If the original last item did not have comma, add comma to it before adding replicas
                    if not last_has_comma and output_lines:
                        output_lines[-1] = output_lines[-1].rstrip('\n') + ',\n'

                    for i, name in enumerate(replica_names):
                        # For last replica, comma optional: put comma anyway for safety except last overall
                        # But this code will put commas always, better to not put comma in last replica
                        is_last_replica = (i == len(replica_names) - 1)
                        comma = ',' if not is_last_replica else ''
                        output_lines.append(f'    "{name}"{comma}\n')

                    # Add the closing bracket line
                    output_lines.append(closing_bracket_line)

                    new_config_block.extend(output_lines)
                    in_list = False
                    list_lines.clear()
                else:
                    # Collect lines inside the list
                    list_lines.append(line.rstrip('\n'))
            else:
                new_config_block.append(line)

        if not new_config_block:
            print(f" Error: Could not find 'PDB_TEMPLATE_FILES_NAMES' in {config_path}")
            return

        with open(output_path, 'w') as f:
            f.write("# === Copy and paste this block into ../prism/config.py ===\n")
            f.writelines(new_config_block)

        print(f" ✓ '{output_path}' generated. Please review and copy the content to 'config.py'.")
    except FileNotFoundError:
        print(f" Error: Could not find {config_path}. Update file not generated.")
    except Exception as e:
        print(f" Error processing config.py: {e}")

def update_ali_files(pdb_dir: str, main_pdb_name: str, pdb_base_name: str, pdb_ext: str, num_replicas: int):
    """
    Copies and modifies 'manual_template_FullSeq.ali' to include replicas.
    The modified file is saved in the CWD (tools/).
    """
    ali_files_to_process = ["manual_template_FullSeq.ali"]
    cwd = os.getcwd()
    print(f"\n--- 3. Modifying .ali files in {cwd} ---")
    for ali_name in ali_files_to_process:
        src_ali_path = os.path.join(pdb_dir, ali_name)
        dst_ali_path = os.path.join(cwd, ali_name)  # Saved to ./tools/
        try:
            shutil.copy2(src_ali_path, dst_ali_path)
        except FileNotFoundError:
            print(f" Warning: Could not find {src_ali_path}. Skipping...")
            continue
        except Exception as e:
            print(f" Error copying {src_ali_path}: {e}")
            continue
        try:
            with open(dst_ali_path, 'r') as f:
                lines = f.readlines()

            new_ali_lines = []
            main_template_sequence_lines = []
            insert_index = -1
            in_main_template_block = False

            for i, line in enumerate(lines):
                new_ali_lines.append(line)
                if line.startswith(f">P1;{main_pdb_name}"):
                    in_main_template_block = True
                    main_template_sequence_lines = []  # Reset
                elif in_main_template_block:
                    if line.startswith("structureX:"):
                        continue
                    if line.strip() == "":
                        continue
                    main_template_sequence_lines.append(line)
                    if "*" in line:
                        in_main_template_block = False
                        insert_index = i + 1  # Mark to insert after this line

            if insert_index != -1 and main_template_sequence_lines:
                replica_blocks = []
                for j in range(1, num_replicas + 1):
                    replica_pdb_name = f"{pdb_base_name}-{j}{pdb_ext}"
                    replica_blocks.append(f">P1;{replica_pdb_name}\n")
                    header = f"structure:{replica_pdb_name}:FIRST:@:END:@:::-1.00:-1.00\n"
                    replica_blocks.append(header)
                    replica_blocks.extend(main_template_sequence_lines)

                new_ali_lines[insert_index:insert_index] = replica_blocks

                with open(dst_ali_path, 'w') as f:
                    f.writelines(new_ali_lines)
                print(f" ✓ {ali_name} modified and saved to {dst_ali_path}")
            else:
                print(f" Error: Could not find 'structureX' block for {main_pdb_name} in {ali_name}")
        except Exception as e:
            print(f" Error processing {dst_ali_path}: {e}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 replicate_template.py <path_to_main_pdb> <number_of_replicas>")
        print("Example (from tools/): python3 replicate_template.py ../input/9CB5_DS_renum_HETATM-B.pdb 50")
        sys.exit(1)

    main_pdb_path = sys.argv[1]

    try:
        num_replicas = int(sys.argv[2])
        if num_replicas < 1:
            raise ValueError
    except ValueError:
        print("Error: Number of replicas must be a positive integer (e.g., 50)")
        sys.exit(1)

    if not os.path.exists(main_pdb_path):
        print(f"Error: Main PDB file not found at: {main_pdb_path}")
        sys.exit(1)

    print(f"Starting replication for '{main_pdb_path}' ({num_replicas} replicas)...")

    main_pdb_name = os.path.basename(main_pdb_path)
    pdb_base_name, pdb_ext = os.path.splitext(main_pdb_name)

    pdb_dir, replica_names = replicate_pdb_files(main_pdb_path, pdb_base_name, pdb_ext, num_replicas)
    update_config_file(replica_names)
    update_ali_files(pdb_dir, main_pdb_name, pdb_base_name, pdb_ext, num_replicas)

    print("\n--- Process complete! ---")
    print("REMEMBER:")
    print("1. Copy the content of 'config_update.txt' into your '../prism/config.py' file.")
    print("2. Move the 'manual_template_FullSeq.ali' generated here (in tools/) to your 'input/' directory.")

if __name__ == "__main__":
    main()
