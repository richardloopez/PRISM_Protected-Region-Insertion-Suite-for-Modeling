#!/usr/bin/env python3
"""
Script para verificar que las coordenadas experimentales NO se movieron.

Compara el MAIN_PDB original con un modelo generado.
PRIMERO alinea el modelo sobre el original (usando un mapeo de C-alfas)
y LUEGO calcula el RMSD de esos mismos residuos.
"""

import sys
import os

try:
    from modeller import *
    # No necesitamos complete_pdb
except ImportError:
    print("ERROR: MODELLER no está instalado o no está en el PYTHONPATH")
    sys.exit(1)

def verificar_coordenadas_fijas(pdb_original, pdb_modelo, mapeado_residuos, 
                               chain_id_orig='A', chain_id_mod='A'):
    """
    Verifica que los residuos experimentales no se hayan movido,
    usando un mapeo explícito (diccionario).
    """
    env = Environ()
    
    env.libs.topology.read(file='$(LIB)/top.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    
    env.io.atom_files_directory = ['.']
    
    # --- Cargar los modelos COMPLETOS ---
    print(f"\nCargando PDB original: {pdb_original} (Cadena: {chain_id_orig})")
    mdl_original = Model(env)
    mdl_original.read(file=pdb_original)
    
    print(f"Cargando modelo generado: {pdb_modelo} (Cadena: {chain_id_mod})")
    mdl_modelo = Model(env)
    mdl_modelo.read(file=pdb_modelo)

    # --- Crear alineamiento 1-a-1 solo con los FRAGMENTOS ---
    print("\nCreando alineamiento 1-a-1 para la superposición...")
    aln = Alignment(env)
    
    # Fragmento Original
    r_orig_min = min(mapeado_residuos.keys())
    r_orig_max = max(mapeado_residuos.keys())
    mdl_frag_orig = Model(env)
    mdl_frag_orig.read(file=pdb_original, 
                       model_segment=(f'{r_orig_min}:{chain_id_orig}', f'{r_orig_max}:{chain_id_orig}'))
    aln.append_model(mdl_frag_orig, align_codes='ORIG_FRAG', atom_files=pdb_original)

    # Fragmento Modelo
    r_mod_min = min(mapeado_residuos.values())
    r_mod_max = max(mapeado_residuos.values())
    mdl_frag_mod = Model(env)
    mdl_frag_mod.read(file=pdb_modelo, 
                      model_segment=(f'{r_mod_min}:{chain_id_mod}', f'{r_mod_max}:{chain_id_mod}'))
    aln.append_model(mdl_frag_mod, align_codes='MOD_FRAG', atom_files=pdb_modelo)
    
    # --- Selección de C-alfas del FRAGMENTO de REFERENCIA ---
    print(f"\nCreando selección de {len(mdl_frag_orig.residues)} C-alfas del fragmento original...")
    sel_frag_orig = Selection(mdl_frag_orig).only_atom_types('CA')
    
    if len(sel_frag_orig) == 0:
        print("❌ ERROR: La selección para el alineamiento está vacía.")
        return False

    print(f"Seleccionados {len(sel_frag_orig)} átomos C-alfa del PDB original.")
    
    # --- Superponer FRAGMENTO vs FRAGMENTO ---
    print("\nRealizando superposición estructural (fragmento vs fragmento)...")
    r = sel_frag_orig.superpose(mdl_frag_mod, aln)
    global_rmsd = r.rms
    
    print(f"✓ Superposición de fragmentos OK (RMSD={global_rmsd:.6f} Å).")
    
    # --- Aplicar esa misma transformación al MODELO COMPLETO ---
    print(f"Aplicando transformación al {pdb_modelo} completo...")
    sel_full_modelo = Selection(mdl_modelo)
    sel_full_modelo.transform(r.rotation)
    sel_full_modelo.translate(r.translation)
    
    # --- BUCLE DE CÁLCULO (SIN PRINTS) ---
    # Primero calculamos todo y lo guardamos en una lista
    
    max_rmsd = 0.0
    total_rmsd = 0.0
    lista_resultados = [] # Aquí guardaremos todos los resultados
    hubo_errores_calculo = False
    
    for res_num_orig, res_num_mod in sorted(mapeado_residuos.items()):
        try:
            res_orig = mdl_original.residues[f'{res_num_orig}:{chain_id_orig}']
            res_model = mdl_modelo.residues[f'{res_num_mod}:{chain_id_mod}']
            
            aa_orig = res_orig.pdb_name
            aa_mod = res_model.pdb_name
            
            ca_orig = res_orig.atoms['CA']
            ca_model = res_model.atoms['CA']
            
            rmsd = ((ca_orig.x - ca_model.x)**2 + 
                    (ca_orig.y - ca_model.y)**2 + 
                    (ca_orig.z - ca_model.z)**2)**0.5
            
            total_rmsd += rmsd
            max_rmsd = max(max_rmsd, rmsd)
            
            lista_resultados.append((res_num_orig, res_num_mod, aa_orig, aa_mod, rmsd))
                    
        except Exception as e:
            print(f"❌  Error al calcular {res_num_orig} (Orig) -> {res_num_mod} (Mod): {e}")
            lista_resultados.append((res_num_orig, res_num_mod, '???', '???', None))
            hubo_errores_calculo = True
    
    if not mapeado_residuos:
        print("\n❌ RESULTADO: El diccionario de mapeo está vacío.")
        return False

    # --- SECCIÓN 1: ORDENADOS POR POSICIÓN ---
    print(f"\n--------------------------------------------------------------------------------")
    print(f"RESIDUOS ORDENADOS POR POSICIÓN")
    print(f"--------------------------------------------------------------------------------")
    
    for res_o, res_m, aa_o, aa_m, rmsd in lista_resultados: # Ya está ordenado por posición
        if rmsd is not None:
            # Añadimos un indicador visual simple
            prefijo = "✓  " if rmsd <= 0.01 else "   "
            print(f"{prefijo} Res {aa_o:3s} {res_o:4d} (Orig) -> {aa_m:3s} {res_m:4d} (Mod): {rmsd:.6f} Å")
        else:
            print(f"❌   Res {aa_o:3s} {res_o:4d} (Orig) -> {aa_m:3s} {res_m:4d} (Mod): Error de comparación")

    # --- SECCIÓN 2: ORDENADOS POR RMSD ---
    print(f"\n--------------------------------------------------------------------------------")
    print(f"RESIDUOS ORDENADOS POR RMSD")
    print(f"--------------------------------------------------------------------------------")
    
    # Ordenar por RMSD (descendente), Nones al final
    resultados_sorted = sorted(lista_resultados, 
                             key=lambda x: x[4] if x[4] is not None else -float('inf'), 
                             reverse=True)
    
    for res_o, res_m, aa_o, aa_m, rmsd in resultados_sorted:
        if rmsd is not None:
            prefijo = "✓  " if rmsd <= 0.01 else "   "
            print(f"{prefijo} Res {aa_o:3s} {res_o:4d} (Orig) -> {aa_m:3s} {res_m:4d} (Mod): {rmsd:.6f} Å")
        else:
            print(f"❌   Res {aa_o:3s} {res_o:4d} (Orig) -> {aa_m:3s} {res_m:4d} (Mod): Error de comparación")

    # --- SECCIÓN 3: RESUMEN RMSD ---
    print(f"\n--------------------------------------------------------------------------------")
    print(f"RESUMEN RMSD")
    print(f"--------------------------------------------------------------------------------")
    print(f"  - RMSD promedio (individual): {total_rmsd / len(mapeado_residuos):.6f} Å")
    print(f"  - RMSD máximo (individual):   {max_rmsd:.6f} Å")
    print(f"  - RMSD global (de 'superpose'): {global_rmsd:.6f} Å")
    
    # Decidir el éxito o fracaso para el código de salida
    hubo_errores_rmsd = any(r[4] > 0.01 for r in lista_resultados if r[4] is not None)

    if not hubo_errores_rmsd and not hubo_errores_calculo:
        print(f"\n✓  ÉXITO: Todos los residuos experimentales mapeados se mantuvieron fijos (RMSD < 0.01 Å)")
        return True
    else:
        if hubo_errores_calculo:
            print("\nAVISO: Hubo errores al comparar algunos residuos.")
        if hubo_errores_rmsd:
             print("\nAVISO: Algunos residuos se movieron (RMSD > 0.01 Å).")
        return False

if __name__ == '__main__':
    print("="*80)
    print("VERIFICACIÓN DE COORDENADAS (CON ALINEAMIENTO)")
    print("="*80)
    
    PDB_ORIGINAL = 'pdb-del-pdb.pdb'
    PDB_MODELO = 'pdb-de-modeller.pdb'
    
    # --- CONFIGURA TU MAPEO AQUÍ ---
    # Mapea 1-166 (Orig) a 307-472 (Mod)
    
    orig_res = list(range(1, 167))      # 1 a 166
    mod_res = list(range(307, 473))   # 307 a 472
    
    if len(orig_res) != len(mod_res):
        print("ERROR en la configuración del script: los rangos del mapeo no coinciden.")
        print(f"Original: {len(orig_res)} residuos, Modelo: {len(mod_res)} residuos")
        sys.exit(1)
        
    MAPEADO_RESIDUOS = dict(zip(orig_res, mod_res))
    
    # -------------------------------------------------------------------
    
    if not os.path.exists(PDB_ORIGINAL):
        print(f"\n❌  ERROR: No se encuentra el archivo {PDB_ORIGINAL}")
        sys.exit(1)
    
    if not os.path.exists(PDB_MODELO):
        print(f"\n❌  ERROR: No se encuentra el archivo {PDB_MODELO}")
        sys.exit(1)
    
    print(f"\nPDB Original:  {PDB_ORIGINAL}")
    print(f"PDB Modelo:    {PDB_MODELO}")
    print(f"Residuos a mapear y alinear: {len(MAPEADO_RESIDUOS)}")
    
    exito = verificar_coordenadas_fijas(
        PDB_ORIGINAL, 
        PDB_MODELO, 
        MAPEADO_RESIDUOS,
        chain_id_orig='A', # Cadena a usar en el PDB Original
        chain_id_mod='A'   # Cadena a usar en el PDB Modelo
    )
    
    sys.exit(0 if exito else 1)