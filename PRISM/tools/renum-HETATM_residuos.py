#!/usr/bin/env python3
import sys

def renumerar_residuos(pdb_file, hetatm_chains=None):
    """
    Renumera los residuos de un archivo PDB, maneja registros TER,
    y cambia ATOM a HETATM en las cadenas especificadas.

    :param pdb_file: Ruta del archivo PDB de entrada.
    :param hetatm_chains: Lista o conjunto de identificadores de cadena (ej: ['B', 'D'])
                          que deben ser escritos como HETATM. Las demás serán ATOM.
                          Si es None o vacío, todas se escriben como ATOM.
    """
    
    if hetatm_chains is None:
        hetatm_chains = set()
        hetatm_flag = ""
    elif isinstance(hetatm_chains, str):
        # Para facilitar la entrada de una única cadena como string
        hetatm_chains = {hetatm_chains}
        hetatm_flag = "_HETATM"
    else:
        hetatm_chains = set(hetatm_chains)
        hetatm_flag = "_HETATM" if hetatm_chains else ""
        
    # Construir el nuevo nombre de archivo: original_renum[_HETATM].pdb
    output_file = pdb_file.replace(".pdb", f"_renum{hetatm_flag}.pdb")

    with open(pdb_file, "r") as f_in, open(output_file, "w") as f_out:
        current_residue = None
        new_resnum = 0

        for line in f_in:
            record_name = line[0:6].strip()

            if record_name in ("ATOM", "HETATM"):
                # (resname, chain, old_resnum)
                res_id = (line[17:20].strip(), line[21].strip(), line[22:26].strip())  
                chain_id = res_id[1]

                if res_id != current_residue:
                    new_resnum += 1
                    current_residue = res_id

                # Cambiar a HETATM si la cadena está en la lista
                new_record_name = "HETATM" if chain_id in hetatm_chains else "ATOM  "
                
                # Reemplazar nombre del registro (col. 1-6) y número de residuo (col. 23-26)
                new_line = new_record_name + line[6:22] + f"{new_resnum:4d}" + line[26:]
                f_out.write(new_line)

            elif record_name == "TER":
                # Incluir el registro TER con el mismo índice del residuo anterior
                if new_resnum > 0:
                    new_line = line[:22] + f"{new_resnum:4d}" + line[26:]
                    f_out.write(new_line)
                else:
                    f_out.write(line)
            
            else:
                # Escribir otras líneas tal cual
                f_out.write(line)

    print(f"Archivo renumerado creado: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Uso: python3 renumerar_residuos.py archivo.pdb [cadenas_HETATM]")
        print("Ejemplo 1 (solo renumera): python3 renumerar_residuos.py 1abc.pdb")
        print("Ejemplo 2 (renumera y cadena B a HETATM): python3 renumerar_residuos.py 1abc.pdb B")
        print("Ejemplo 3 (renumera y cadenas B y D a HETATM): python3 renumerar_residuos.py 1abc.pdb B,D")
        sys.exit(1)

    pdb_file = sys.argv[1]
    
    if len(sys.argv) == 3:
        chains_input = sys.argv[2]
        hetatm_chains_list = [chain.strip() for chain in chains_input.split(',')]
        renumerar_residuos(pdb_file, hetatm_chains_list)
    else:
        renumerar_residuos(pdb_file)