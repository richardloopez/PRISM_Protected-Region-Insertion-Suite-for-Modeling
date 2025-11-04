#!/usr/bin/env python3
"""
Script mejorado para la detección de HETATM/BLK.
Detecta automáticamente todas las cadenas y las analiza individualmente.
No requiere Modeller, solo lee el archivo PDB directamente.
"""

import sys
import os
from typing import List, Dict, Any

# --- FUNCIONES DE EXTRACCIÓN (Modificadas para eliminar try/except) ---

def extract_hetatm_residues_simple(pdb_file: str, chain_id: str) -> List[Dict[str, Any]]:
    """
    Extrae información de residuos HETATM únicos para una CADENA específica.
    """
    hetatm_residues = []
    last_atom_resnum = 0
    seen_hetatm = set()
    
    with open(pdb_file, 'r') as f:
        for line in f:
            res_chain = line[21:22].strip()
            
            # Solo procesar líneas de la cadena solicitada
            if res_chain != chain_id:
                continue

            if line.startswith('ATOM'):
                try:
                    resnum = int(line[22:26].strip())
                    last_atom_resnum = max(last_atom_resnum, resnum)
                except ValueError:
                    continue
            
            elif line.startswith('HETATM'):
                resname = line[17:20].strip()
                try:
                    resnum = int(line[22:26].strip())
                except ValueError:
                    continue
                
                res_key = (resname, resnum, res_chain)
                if res_key not in seen_hetatm:
                    seen_hetatm.add(res_key)
                    hetatm_residues.append({
                        'resname': resname,
                        'resnum': resnum,
                        'chain': res_chain,
                        'position_after_atom_resnum': last_atom_resnum
                    })
    
    hetatm_residues.sort(key=lambda x: x['resnum'])
    return hetatm_residues

def get_all_chains(pdb_file: str) -> Dict[str, Dict[str, int]]:
    """Obtiene información sobre todas las cadenas en el PDB"""
    chains_info = {}
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain = line[21:22].strip()
                if not chain: # Omitir cadenas sin identificar
                    continue
                record_type = 'ATOM' if line.startswith('ATOM') else 'HETATM'
                
                if chain not in chains_info:
                    chains_info[chain] = {'ATOM': 0, 'HETATM': 0}
                chains_info[chain][record_type] += 1
    
    return chains_info

# --- FUNCIONES DE IMPRESIÓN AUXILIARES ---

def print_hetatm_summary_tables(hetatm_residues_list: List[Dict[str, Any]]):
    """Imprime las tablas de resumen de tipos y detalles de HETATM."""
    
    # Agrupar por tipo de residuo
    residue_types = {}
    for het in hetatm_residues_list:
        resname = het['resname']
        if resname not in residue_types:
            residue_types[resname] = 0
        residue_types[resname] += 1
    
    print(f"\n   Tipos de residuos HETATM:")
    for resname, count in sorted(residue_types.items()):
        print(f"     - {resname}: {count} residuo(s)")
    
    print(f"\n   Primeros 10 residuos HETATM detectados:")
    print(f"     {'Tipo':<8} {'ResNum':<8} {'Cadena':<8} {'Después de ATOM #':<20}")
    print(f"     {'-'*50}")
    for het in hetatm_residues_list[:10]:
        print(f"     {het['resname']:<8} {het['resnum']:<8} {het['chain']:<8} {het['position_after_atom_resnum']:<20}")
    
    if len(hetatm_residues_list) > 10:
        print(f"     ... y {len(hetatm_residues_list) - 10} más")

def print_modeller_explanation(hetatm_count: int, chain_id: str):
    """Imprime las instrucciones de modelado para una cadena específica."""
    print("\n   --- 💡 Instrucciones de Modelado (Cadena {chain_id}) ---")
    print(f"    a) env.io.hetatm = True está configurado ✓")
    print(f"    b) Modeller leerá {hetatm_count} residuos HETATM automáticamente para esta cadena.")
    print(f"    c) En el alineamiento PIR para la CADENA {chain_id} DEBERÁS AÑADIR:")
    print(f"       - Template: {hetatm_count} caracteres '.' (BLK) al final de la secuencia de la cadena {chain_id}")
    print(f"       - Target:   {hetatm_count} caracteres '.' (BLK) al final de la secuencia de la cadena {chain_id}")
    print(f"    d) Los residuos HETATM serán tratados como obstáculos durante el modelado.")
    print("   -----------------------------------------------------")


# --- FUNCIÓN PRINCIPAL ---

def main():
    """
    Ejecuta el análisis de detección de HETATM en todas las cadenas de un PDB.
    """
    
    # --- 1. Configuración y Carga de Archivo ---
    if len(sys.argv) > 1:
        PDB_FILE = sys.argv[1]
    else:
        print("Uso: python tu_script.py <archivo.pdb>")
        PDB_FILE = '8vx1_DS_renum_HETATM.pdb' # Archivo por defecto
        print(f"Usando archivo de prueba por defecto: '{PDB_FILE}'\n")

    if not os.path.exists(PDB_FILE):
        print(f"[ERROR] Archivo no encontrado: {PDB_FILE}")
        return 1

    # Asumimos que la cadena 'A' es siempre la peptídica principal
    PEPTIDIC_CHAIN_ID = 'A'
    modelled_hetatm_chains = {} # Para el resumen final

    print("="*70)
    print("PRUEBA DE DETECCIÓN DE RESIDUOS HETATM/BLK (v2.0 Automática)")
    print("="*70)
    print(f"\n1. Archivo PDB Template: {PDB_FILE}")

    # --- 2. Obtener Resumen de Cadenas ---
    try:
        chains_info = get_all_chains(PDB_FILE)
    except Exception as e:
        print(f"[ERROR] No se pudo leer el archivo PDB '{PDB_FILE}'. Error: {e}")
        return 1
        
    if not chains_info:
        print("[ERROR] Archivo PDB vacío o no se pudieron leer cadenas.")
        return 1

    print(f"\n2. Resumen de Cadenas en el PDB:")
    print(f"   {'Cadena':<10} {'ATOM (líneas)':<15} {'HETATM (líneas)':<15}")
    print(f"   {'-'*42}")
    for chain, counts in sorted(chains_info.items()):
        print(f"   {chain:<10} {counts['ATOM']:<15} {counts['HETATM']:<15}")

    # --- 3. Análisis Detallado por Cadena ---
    print("\n" + "="*70)
    print("3. ANÁLISIS DETALLADO POR CADENA")
    print("="*70)

    for chain_id in sorted(chains_info.keys()):
        print(f"\n### 🧬 Análisis Cadena: {chain_id} ###")
        
        try:
            hetatm_residues = extract_hetatm_residues_simple(PDB_FILE, chain_id)
        except Exception as e:
            print(f"  [ERROR] Ocurrió un error analizando la cadena {chain_id}: {e}")
            continue

        num_het_residues = len(hetatm_residues)
        atom_count = chains_info[chain_id]['ATOM']

        # Lógica para la cadena peptídica principal
        if chain_id == PEPTIDIC_CHAIN_ID:
            if num_het_residues > 0:
                print(f"  ❌ ¡ADVERTENCIA! Se detectaron {num_het_residues} residuos HETATM en la cadena peptídica '{chain_id}'.")
                print("     Esta cadena solo debería contener ATOM. Revise su archivo PDB.")
                print_hetatm_summary_tables(hetatm_residues)
            else:
                print(f"  ✓ CORRECTO: La cadena peptídica '{chain_id}' no contiene residuos HETATM.")
                print(f"     Contiene {atom_count} registros ATOM.")
        
        # Lógica para otras cadenas (ligandos, ADN, otras proteínas)
        else:
            if num_het_residues > 0:
                print(f"  ⓘ DETECTADOS: {num_het_residues} residuos HETATM (ligandos, ADN, etc.).")
                print(f"     Esta cadena también contiene {atom_count} registros ATOM.")
                modelled_hetatm_chains[chain_id] = num_het_residues
                
                # Mostrar tablas y explicación de modelado
                print_hetatm_summary_tables(hetatm_residues)
                print_modeller_explanation(num_het_residues, chain_id)
            else:
                print(f"  ⓘ No se detectaron residuos HETATM en esta cadena.")
                if atom_count > 0:
                    print(f"     Contiene {atom_count} registros ATOM (considerada cadena proteica).")
                else:
                    print(f"     Esta cadena no contiene residuos proteicos (0 ATOM).")

    # --- 4. Resumen Final ---
    print("\n" + "="*70)
    print("✓ RESUMEN DE PRUEBA COMPLETADA")
    print("="*70)
    
    if not modelled_hetatm_chains:
        print("\nEl sistema está listo. No se han detectado residuos HETATM")
        print("que requieran ser modelados como BLK (aparte de la cadena 'A').")
    else:
        # Calculamos el total sumando los valores (counts) del diccionario
        total_hetatm_count = sum(modelled_hetatm_chains.values())

        print("\nEl sistema está listo para modelar con residuos HETATM/BLK.")
        print("Resumen de cadenas que incluirán BLK en el alineamiento:")
            
        for chain_id, count in modelled_hetatm_chains.items():
            print(f"    - Cadena {chain_id}: {count} residuos HETATM (se deben añadir {count} caracteres '.' al final de los .ali por esta cadena)")
            
        print("---") # Separador para claridad
        print(f"TOTAL HETATM (BLK): {total_hetatm_count}")
        print(f"(Se deberán añadir un total de {total_hetatm_count} caracteres '.' al final de los .ali)")

    print("\n" + "="*70)
    return 0

if __name__ == '__main__':
    sys.exit(main())