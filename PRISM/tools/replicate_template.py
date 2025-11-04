#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
replicate_template.py

Este script automatiza el "truco" de replicar la plantilla experimental (MAIN_PDB)
X veces para darle más peso en el promedio de AutoModel.

Genera 3 cosas:
1. Las N copias del archivo PDB en el directorio 'input'.
2. Un 'config_update.txt' con el bloque de código para actualizar 'config.py'.
3. Un 'manual_template_FullSeq.ali' actualizado en 'tools/' (que deberás mover a 'input').

Uso (desde la carpeta 'tools/'):
    python3 replicate_template.py ../input/9CB5_DS_renum_HETATM-B.pdb 50

Argumentos:
    1. Ruta al MAIN_PDB (ej: ../input/mi_pdb.pdb)
    2. Número de réplicas a generar (ej: 50)
"""

import sys
import os
import shutil
import re
from typing import Tuple, List 

def replicate_pdb_files(main_pdb_path: str, pdb_base_name: str, pdb_ext: str, num_replicas: int) -> Tuple[str, List[str]]:
    """
    Crea X copias del MAIN_PDB en el mismo directorio.
    Retorna el directorio de PDBs y la lista de nombres de réplicas.
    """
    pdb_dir = os.path.dirname(main_pdb_path)
    replica_names = []
    
    print(f"\n--- 1. Replicando archivos PDB en {pdb_dir} ---")
    
    for i in range(1, num_replicas + 1):
        # Nomenclatura: base-1.pdb, base-2.pdb
        replica_name = f"{pdb_base_name}-{i}{pdb_ext}"
        replica_path = os.path.join(pdb_dir, replica_name)
        
        try:
            shutil.copy2(main_pdb_path, replica_path)
            replica_names.append(replica_name)
        except Exception as e:
            print(f"  Error copiando a {replica_path}: {e}")
            
    print(f"  ✓ Creadas {len(replica_names)} réplicas de PDB.")
    return pdb_dir, replica_names

def update_config_file(replica_names: List[str]):
    """
    Lee ../src/config.py y genera un 'config_update.txt' en el CWD
    con la lista PDB_TEMPLATE_FILES_NAMES actualizada.
    """
    config_path = os.path.join("..", "src", "config.py")
    output_path = "config_update.txt"
    
    print(f"\n--- 2. Generando {output_path} ---")
    
    try:
        with open(config_path, 'r') as f:
            lines = f.readlines()

        new_config_block = []
        in_list = False

        for line in lines:
            if line.strip().startswith("PDB_TEMPLATE_FILES_NAMES"):
                in_list = True
                new_config_block.append(line)
            elif in_list:
                if "]" in line:
                    # Insertar las réplicas ANTES del ']'
                    for name in replica_names:
                        new_config_block.append(f'    "{name}",\n')
                    new_config_block.append(line)
                    in_list = False
                    break # Terminamos de procesar la lista
                else:
                    new_config_block.append(line)
        
        if not new_config_block:
            print("  Error: No se pudo encontrar 'PDB_TEMPLATE_FILES_NAMES' en config.py")
            return

        with open(output_path, 'w') as f:
            f.write("# === Copia y pega este bloque en ../src/config.py ===\n")
            f.writelines(new_config_block)
            
        print(f"  ✓ 'config_update.txt' generado. Por favor, revisa y copia el contenido a 'config.py'.")

    except FileNotFoundError:
        print(f"  Error: No se encontró {config_path}. No se pudo generar el update.")
    except Exception as e:
        print(f"  Error procesando config.py: {e}")

# --- MODIFICACIÓN ---
# Esta función ahora solo procesa el .ali limpio.
def update_ali_files(pdb_dir: str, main_pdb_name: str, pdb_base_name: str, pdb_ext: str, num_replicas: int):
    """
    Copia y modifica el archivo 'manual_template_FullSeq.ali' para incluir las réplicas.
    El archivo modificado se guarda en el CWD.
    """
    # --- MODIFICACIÓN ---
    # Ya no se procesa el archivo _cde.ali
    ali_files_to_process = ["manual_template_FullSeq.ali"]
    cwd = os.getcwd()

    print(f"\n--- 3. Modificando archivos .ali en {cwd} ---")

    for ali_name in ali_files_to_process:
        src_ali_path = os.path.join(pdb_dir, ali_name)
        dst_ali_path = os.path.join(cwd, ali_name) # Se guarda en ./tools/
        
        try:
            shutil.copy2(src_ali_path, dst_ali_path)
        except FileNotFoundError:
            print(f"  Advertencia: No se encontró {src_ali_path}. Saltando...")
            continue
        except Exception as e:
            print(f"  Error copiando {src_ali_path}: {e}")
            continue

        try:
            with open(dst_ali_path, 'r') as f:
                lines = f.readlines()

            new_ali_lines = []
            main_template_sequence_lines = []
            insert_index = -1
            in_main_template_block = False

            # Buscamos el bloque del MAIN_PDB (structureX)
            for i, line in enumerate(lines):
                new_ali_lines.append(line)
                
                if line.startswith(f">P1;{main_pdb_name}"):
                    in_main_template_block = True
                    main_template_sequence_lines = [] # Reset
                elif in_main_template_block:
                    if line.startswith("structureX:"):
                        continue
                    if line.strip() == "": # Ignorar líneas vacías
                        continue
                    
                    main_template_sequence_lines.append(line)
                    
                    if "*" in line:
                        in_main_template_block = False
                        insert_index = i + 1 # Marcar para insertar después de esta línea
            
            if insert_index != -1 and main_template_sequence_lines:
                # Construimos los bloques de réplica
                replica_blocks = []
                for j in range(1, num_replicas + 1):
                    replica_pdb_name = f"{pdb_base_name}-{j}{pdb_ext}"
                    
                    replica_blocks.append(f">P1;{replica_pdb_name}\n")
                    
                    # Usamos 'structure' y la cabecera genérica 'FIRST:@:END:@'
                    header = f"structure:{replica_pdb_name}:FIRST:@:END:@:::-1.00:-1.00\n"
                    replica_blocks.append(header)
                    replica_blocks.extend(main_template_sequence_lines)
                
                # Insertamos los nuevos bloques en la lista de líneas
                new_ali_lines[insert_index:insert_index] = replica_blocks
                
                # Sobrescribimos el archivo en el CWD
                with open(dst_ali_path, 'w') as f:
                    f.writelines(new_ali_lines)
                print(f"  ✓ {ali_name} modificado y guardado en {dst_ali_path}")
            else:
                print(f"  Error: No se pudo encontrar el bloque 'structureX' para {main_pdb_name} en {ali_name}")
                
        except Exception as e:
            print(f"  Error procesando {dst_ali_path}: {e}")

def main():
    if len(sys.argv) != 3:
        print("Uso: python3 replicate_template.py <ruta_al_main_pdb> <numero_de_replicas>")
        print("Ejemplo (desde tools/): python3 replicate_template.py ../input/9CB5_DS_renum_HETATM-B.pdb 50")
        sys.exit(1)

    main_pdb_path = sys.argv[1]
    
    try:
        num_replicas = int(sys.argv[2])
        if num_replicas < 1:
            raise ValueError
    except ValueError:
        print(f"Error: El número de réplicas debe ser un entero positivo (ej: 50)")
        sys.exit(1)

    if not os.path.exists(main_pdb_path):
        print(f"Error: No se encuentra el archivo PDB principal en: {main_pdb_path}")
        sys.exit(1)

    print(f"Iniciando replicación para '{main_pdb_path}' ({num_replicas} réplicas)...")
    
    # Extraer nombres
    main_pdb_name = os.path.basename(main_pdb_path)
    pdb_base_name, pdb_ext = os.path.splitext(main_pdb_name)

    # --- Tarea 1: Replicar PDBs ---
    pdb_dir, replica_names = replicate_pdb_files(main_pdb_path, pdb_base_name, pdb_ext, num_replicas)

    # --- Tarea 2: Adaptar Config ---
    update_config_file(replica_names)

    # --- Tarea 3: Adaptar .ali ---
    update_ali_files(pdb_dir, main_pdb_name, pdb_base_name, pdb_ext, num_replicas)

    print("\n--- ¡Proceso completado! ---")
    print("RECUERDA:")
    print("1. Copia el contenido de 'config_update.txt' en tu archivo '../src/config.py'.")
    # --- MODIFICACIÓN ---
    # Se elimina la mención al archivo _cde.ali
    print("2. Mueve el archivo 'manual_template_FullSeq.ali' generado aquí (en tools/) a tu carpeta 'input/' (o donde estén tus PDBs).")

if __name__ == "__main__":
    main()