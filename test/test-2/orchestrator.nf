#!/usr/bin/env nextflow

/*
========================================================================================
    PRISM: Protected-Region Insertion Suite for Modeling
========================================================================================
    Author  : Richard Lopez Corbalan
    GitHub  : github.com/richardloopez
    License : MIT
----------------------------------------------------------------------------------------
    Orchestrator: Main DSL2 Workflow 
    Handles the parallelization and scheduling of the MODELLER pipeline.
========================================================================================
*/

import groovy.json.JsonSlurper

// --- 1. CONFIGURATION LOADING ---
def pythonPath = ['bash', '-c', 'which python3'].execute().text.trim()
def config_executor = [pythonPath, 'config.py'].execute(null, new File('PRISM'))
def stdout = new StringBuilder()
def stderr = new StringBuilder()
config_executor.waitForProcessOutput(stdout, stderr)

if (config_executor.exitValue() != 0 ) {
    error "Failed to load config.py: ${stderr.toString()}"          
}

def config_data = new JsonSlurper().parseText(stdout.toString())

nextflow.enable.dsl=2

params.input_dir = config_data.INPUT_DIR
params.results_dir = config_data.MODELING_RESULTS_DIR
params.num_auto    = config_data.TOTAL_HOMOLOGY_MODELS
params.num_split   = config_data.TOTAL_PARALLEL_JOBS
params.num_refine  = config_data.TOP_MODELS_FOR_REFINEMENT
params.execution_paradigm  = config_data.EXECUTION_PARADIGM
params.rsr_file = config_data.CUSTOM_RSRFILE_PATH
params.ini_file = config_data.CUSTOM_INIFILE_PATH

def flatPath = { filename -> file(filename).name }

// --- 2. PROCESS DEFINITIONS ---

process PREREQ_CDE {
    tag "S1-Prereq"
    publishDir "${params.input_dir}", mode: 'copy', saveAs: flatPath, pattern: "input/*.{ali,ss2}"
    publishDir "${params.results_dir}/logs", mode: 'copy', pattern: "*.log"
    
    output:
    path "input/*.ali", emit: alignments
    path "input/*.ss2", emit: ss2_file, optional: true 
    path "psipred_results", emit: psipred_results_dir, optional: true
    path "*.log", emit: log

    script:
    """
    cp -r ${projectDir}/PRISM .
    cp ${projectDir}/config.yaml .
    mkdir -p input modeling_results
    
    ln -sf ${projectDir}/input/* ./input/
    
    ROOT_DIR=\$(pwd)

    export PYTHONPATH=\$ROOT_DIR:\$PYTHONPATH
    python3 -m PRISM.controller --stage prereq-cde > S1_prereq_cde.log 2>&1
    """
}

process AUTOMODEL_PRECALC {
    tag "S2-AutoModel-Precalc:${job_id}"
    publishDir "${params.input_dir}", mode: 'copy', saveAs: flatPath, pattern: "input/*.{rsr,pdb,ali}"
    publishDir "${params.results_dir}/logs", mode: 'copy', pattern: "*.log"
    
    input:
    val job_id
    val input_mode
    path "modeling_results/*"

    output:
    path "input/*.rsr", emit: rsr_file
    path "input/*_ini.pdb", emit: ini_file 
    path "input/*precalculation.ali", emit: precalculation_ali
    path "*.log", emit: log

    script:
    """
    cp -r ${projectDir}/PRISM .
    cp ${projectDir}/config.yaml .
    mkdir -p input modeling_results
    
    ln -sf ${projectDir}/input/* ./input/
    
    ROOT_DIR=\$(pwd)

    cd modeling_results/
    export PYTHONPATH=\$ROOT_DIR:\$PYTHONPATH
    python3 -m PRISM.controller --stage automodel --job-id ${job_id} --input-mode ${input_mode} > ../S2_automodel_precalc_${job_id}.log 2>&1

    """
}

process AUTOMODEL {
    tag "S2-AutoModel:${job_id}"
    publishDir "${params.input_dir}", mode: 'copy', pattern: "input/*precomputed.ali", saveAs: flatPath
    publishDir "${params.results_dir}/logs", mode: 'copy', pattern: "*.log"

    input:
    val job_id
    val ali_file
    val input_mode
    val rsr_file
    val ini_file

    output:
    path "modeling_results/*.pdb", emit: models
    path "input/*precomputed.ali", emit: precomputed_ali, optional: true
    path "*.log", emit: log

    script:
    """
    cp -r ${projectDir}/PRISM .
    cp ${projectDir}/config.yaml .
    mkdir -p input modeling_results

    ln -sf ${projectDir}/modeling_results/* ./modeling_results/
    ln -sf ${projectDir}/input/* ./input/
    
    ROOT_DIR=\$(pwd)

    cd modeling_results/
    export PYTHONPATH=\$ROOT_DIR:\$PYTHONPATH
    python3 -m PRISM.controller --stage automodel --job-id ${job_id} --input-mode ${input_mode} > ../S2_automodel_${job_id}.log 2>&1

    """
}

process RANK_AUTOMODEL {
    tag "S3-Rank"
    publishDir "${params.results_dir}", mode: 'copy', saveAs: flatPath, pattern: "modeling_results/*.pdb"
    publishDir "${params.results_dir}/logs", mode: 'copy', pattern: "*.log"

    input:
    path "modeling_results/*"

    output:
    path "modeling_results/AUTO_*.pdb", emit: ranked_models
    path "*.log", emit: log

    script:
    """
    cp -r ${projectDir}/PRISM .
    cp ${projectDir}/config.yaml .
    mkdir -p input 
    
    ln -sf ${projectDir}/input/* ./input/
    ln -sf ${projectDir}/modeling_results/* ./modeling_results/
    
    ROOT_DIR=\$(pwd)

    cd modeling_results/
    export PYTHONPATH=\$ROOT_DIR:\$PYTHONPATH
    python3 -m PRISM.controller --stage rank-automodel > ../S3_rank_automodel.log 2>&1
    """
}

process LOOP_MODEL {
    tag "S4-Loop:${model_file.baseName}"
    publishDir "${params.results_dir}", mode: 'copy', saveAs: flatPath, pattern: "modeling_results/*.pdb"
    publishDir "${params.results_dir}/logs", mode: 'copy', pattern: "*.log"
    
    input:
    path model_file
    path "modeling_results/*"

    output:
    path "modeling_results/*LOOP*.pdb", emit: refined_pdbs, optional: true
    path "*.log", emit: log

    script:
    """
    cp -r ${projectDir}/PRISM .
    cp ${projectDir}/config.yaml .
    mkdir -p input
    
    ln -sf ${projectDir}/input/* ./input/
    
    ROOT_DIR=\$(pwd)

    mv ${model_file} modeling_results/
    cd modeling_results/
    export PYTHONPATH=\$ROOT_DIR:\$PYTHONPATH
    python3 -m PRISM.controller --stage loopmodel --input-model ${model_file.name} > ../S4_loop_model_${model_file.baseName}.log 2>&1
    """
}


process FINAL_RANKING {
    tag "S5-Final"
    publishDir "${params.results_dir}", mode: 'copy', saveAs: flatPath, pattern: "modeling_results/*.csv"
    publishDir "${params.results_dir}/logs", mode: 'copy', pattern: "*.log"

    input:
    path "S4_refined_models/*"
    path "S3_ranked_models/*"

    output:
    path "modeling_results/final_ranking.csv" 
    path "*.log", emit: log

    script:
    """
    cp -r ${projectDir}/PRISM .
    cp ${projectDir}/config.yaml .
    mkdir -p input modeling_results
    
    ln -sf ${projectDir}/input/* ./input/
    
    ROOT_DIR=\$(pwd)

    mv S3_ranked_models/*.pdb modeling_results/ 2>/dev/null || true
    mv S4_refined_models/*.pdb modeling_results/ 2>/dev/null || true
    cd modeling_results/
    export PYTHONPATH=\$ROOT_DIR:\$PYTHONPATH
    python3 -m PRISM.controller --stage final-rank > ../S5_final_ranking.log 2>&1
    """
}

// --- 3. WORKFLOW ---

workflow {
    PREREQ_CDE() // S1

    if (params.execution_paradigm == "precalculation") {  
        AUTOMODEL_PRECALC(1, params.execution_paradigm, PREREQ_CDE.out.alignments) // S2 (precalculation)
    } 
    else if (params.execution_paradigm == "prism-power") {  
        AUTOMODEL_PRECALC(1, "precalculation", PREREQ_CDE.out.alignments)
        job_ids = Channel.of(1..params.num_split)
        AUTOMODEL(job_ids, PREREQ_CDE.out.alignments.collect(), "precomputed", AUTOMODEL_PRECALC.out.rsr_file.collect(), AUTOMODEL_PRECALC.out.ini_file.collect()) // S2 (prism-power)
        
        RANK_AUTOMODEL(AUTOMODEL.out.models.collect()) // S3 (prism-power)
        
        top_models = RANK_AUTOMODEL.out.ranked_models
            .flatten()
            .filter { file -> 
                def match = (file.name =~ /AUTO_(\d+)\.pdb/)
                if(match) {
                    int rank = match[0][1].toInteger()
                    return rank <= params.num_refine
                }
                return false
            }
        
        LOOP_MODEL(top_models, PREREQ_CDE.out.alignments.collect()) // S4 (prism-power)
        
        FINAL_RANKING(
            LOOP_MODEL.out.refined_pdbs.collect().ifEmpty([]), 
            RANK_AUTOMODEL.out.ranked_models.collect()
        ) // S5
    } 
    else {
        job_ids = Channel.of(1..params.num_split)
        
        AUTOMODEL(job_ids, PREREQ_CDE.out.alignments.collect(), params.execution_paradigm, params.rsr_file, params.ini_file) // S2
        
        RANK_AUTOMODEL(AUTOMODEL.out.models.collect()) // S3
        
        top_models = RANK_AUTOMODEL.out.ranked_models
            .flatten()
            .filter { file -> 
                def match = (file.name =~ /AUTO_(\d+)\.pdb/)
                if(match) {
                    int rank = match[0][1].toInteger()
                    return rank <= params.num_refine
                }
                return false
            }
        
        LOOP_MODEL(top_models, PREREQ_CDE.out.alignments.collect()) // S4
        
        FINAL_RANKING(
            LOOP_MODEL.out.refined_pdbs.collect().ifEmpty([]), 
            RANK_AUTOMODEL.out.ranked_models.collect()
        ) // S5
    }
}