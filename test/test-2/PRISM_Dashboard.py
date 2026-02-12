#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM (PRISM_Dashboard)

'''
PRISM Dashboard: Main graphical user interface for the PRISM modeling suite.
Built with Streamlit for a user-friendly, high-fidelity modeling experience.
'''

import streamlit as st
import os
import yaml
import json
import sys
import pandas as pd
import plotly.express as px
from PRISM.config import settings, PrismConfig, load_settings
from PRISM.ui_utils import (
    run_nextflow, get_nextflow_progress, visualize_pdb, 
    save_uploaded_file, load_ranking_csv, list_files_in_dir,
    run_tool, get_score_distribution_data, delete_file
)
import streamlit.components.v1 as components

st.set_page_config(page_title="PRISM Dashboard", layout="wide")

# Dynamic Logo placement
logo_path = os.path.join("PRISM", "logo.png")
if os.path.exists(logo_path):
    st.logo(logo_path, icon_image=None) # Sidebar
    col_logo, col_title = st.columns([1, 10])
    with col_logo:
        st.image(logo_path, width=60) # Header logo
    with col_title:
        st.title("PRISM: Protected-Region Insertion Suite for Modeling")
else:
    st.title("🛡️ PRISM: Protected-Region Insertion Suite for Modeling")
st.markdown("---")

if 'config' not in st.session_state:
    try:
        st.session_state.config = load_settings()
    except:
        st.error("Could not load config.yaml. Please ensure it exists.")
        st.stop()

# Sidebar
# Theme
theme = st.sidebar.select_slider("Select GUI Theme", options=["Default", "Professional", "High Contrast", "Dark Modern"])

if theme == "Professional":
    st.markdown("""
        <style>
        .stApp { background-color: #f0f2f6 !important; }
        </style>
        """, unsafe_allow_html=True)
elif theme == "High Contrast":
    st.markdown("""
        <style>
        .stApp { background-color: #000000 !important; color: #FFFF00 !important;}
        .stMarkdown p, .stMarkdown h1, .stMarkdown h2, .stMarkdown h3, .stMarkdown h4, 
        label, .stWidgetLabel p { color: #FFFF00 !important; }
        [data-testid="stImage"], [data-testid="stArrowDataTable"], .stPlotlyChart {background-color: #111111 !important;
            border: 2px solid #FFFF00 !important;padding: 10px;border-radius: 5px;}
        button { background-color: #FFFF00 !important; 
            color: #000000 !important; 
            border: 2px solid #FFFFFF !important; 
            font-weight: bold !important;}
        button:hover {background-color: #FFFFFF !important;color: #000000 !important;}
        .stTextInput input, .stTextArea textarea, .stNumberInput input { 
            background-color: #222222 !important; 
            color: #FFFF00 !important; 
            border: 2px solid #FFFF00 !important; }
        section[data-testid="stSidebar"] {background-color: #000000 !important;
            border-right: 2px solid #FFFF00 !important;}
        section[data-testid="stSidebar"] * {color: #FFFF00 !important;}
        </style>
        """, unsafe_allow_html=True)
elif theme == "Dark Modern":
    st.markdown("""
        <style>
        .stApp { background-color: #0e1117 !important; color: white !important; }
        .stApp p, .stApp h1, .stApp h2, .stApp h3, .stApp h4, .stApp label, .stApp span { color: white !important; }
        div.stButton > button { background-color: #262730 !important; color: white !important; border: 1px solid #4B4B4B !important; }
        div.stButton > button:hover { border: 1px solid #FF4B4B !important; color: #FF4B4B !important; }
        .stTextInput input, .stTextArea textarea, .stNumberInput input { color: white !important; background-color: #262730 !important; }
        section[data-testid="stSidebar"] { background-color: #1a1c24 !important; }
        section[data-testid="stSidebar"] * { color: white !important; }
        </style>
        """, unsafe_allow_html=True)

st.sidebar.markdown("---")
st.sidebar.caption("PRISM v1.2.0 | Production Refactor")

# Main Navigation
tab_config, tab_files, tab_viz, tab_exec, tab_results, tab_tools = st.tabs([
    "⚙️ Config", "📂 Files", "🧬 Visualization", "🚀 Execution", "📊 Results", "🛠️ Tools"
])

with tab_config:
    st.header("Global Pipeline Configuration")
    
    # Group 1: Core Modeling Parameters
    with st.container(border=True):
        st.subheader("🧬 Modeling & Refinement")
        col1, col2 = st.columns(2)
        with col1:
            total_models = st.number_input("TOTAL_HOMOLOGY_MODELS", value=st.session_state.config.TOTAL_HOMOLOGY_MODELS, min_value=1)
            top_refine = st.number_input("TOP_MODELS_FOR_REFINEMENT", value=st.session_state.config.TOP_MODELS_FOR_REFINEMENT, min_value=0)
            loop_models = st.number_input("LOOP_MODELS_PER_TARGET", value=st.session_state.config.LOOP_MODELS_PER_TARGET, min_value=1)
            flank_size = st.number_input("MOBILE_FLANK_RESIDUES", value=st.session_state.config.MOBILE_FLANK_RESIDUES, min_value=0)
        with col2:
            repulsion = st.number_input("BLOCK_REPULSION_RADIUS (Å)", value=st.session_state.config.BLOCK_REPULSION_RADIUS, step=0.1)
            refine_flanks = st.checkbox("REFINE_FLANKS_DURING_AUTOMODEL", value=st.session_state.config.REFINE_FLANKS_DURING_AUTOMODEL)
            best_final = st.text_input("NUM_BEST_FINAL_MODELS (e.g., 5 or 'inf')", value=str(st.session_state.config.NUM_BEST_FINAL_MODELS))

        st.markdown("---")
        
    # Group 2: External Dependencies & API
    with st.container(border=True):
        st.subheader("🌐 External Tools & Predictions")
        col3, col4 = st.columns(2)
        with col3:
            psipred_pred = st.toggle("PERFORM_PSIPRED_PREDICTION", value=st.session_state.config.PERFORM_PSIPRED_PREDICTION)
        with col4:
            psipred_email = st.text_input("PSIPRED_EMAIL", value=st.session_state.config.PSIPRED_EMAIL)
            psipred_poll = st.number_input("PSIPRED_POLL_INTERVAL (sec)", value=st.session_state.config.PSIPRED_POLL_INTERVAL, min_value=10)

    # Group 3: Manual Overrides & Data Control
    with st.container(border=True):
        st.subheader("🛠️ Manual Overrides & Control")
        col_ov1, col_ov2 = st.columns(2)
        with col_ov1:
            use_manual_ali = st.toggle("USE_MANUAL_ALIGNMENT", value=st.session_state.config.USE_MANUAL_ALIGNMENT)
        with col_ov2:
            manual_opt = st.toggle("USE_MANUAL_OPTIMIZATION_SELECTION", value=st.session_state.config.USE_MANUAL_OPTIMIZATION_SELECTION)
            
        if manual_opt:
            opt_res_str = st.text_area("MANUAL_OPTIMIZATION_RESIDUES (Space separated)", 
                                     value=" ".join(map(str, st.session_state.config.MANUAL_OPTIMIZATION_RESIDUES)),
                                     help="Add residue numbers. This box appears instantly when the toggle is ON.")
        else:
            opt_res_str = ""

    # Group 4: Execution Engine
    with st.container(border=True):
        st.subheader("💻 Execution Engine")
        col5, col6 = st.columns(2)
        with col5:
            paradigm = st.selectbox("EXECUTION_PARADIGM", ["prism-power", "precalculation", "precomputed", "normal"], 
                                    index=["prism-power", "precalculation", "precomputed", "normal"].index(st.session_state.config.EXECUTION_PARADIGM))
            modeller_cores = st.number_input("MODELLER_CORES", value=st.session_state.config.MODELLER_CORES, min_value=1)
        with col6:
            parallel_jobs = st.number_input("TOTAL_PARALLEL_JOBS", value=st.session_state.config.TOTAL_PARALLEL_JOBS, min_value=1)

    # Group 5: Templates & Power Settings
    with st.container(border=True):
        st.subheader("🧬 Templates & PRISM Power")
        pdb_names_str = st.text_area("PDB_TEMPLATE_FILES_NAMES (One per line, first is MAIN)", 
                                     value="\n".join(st.session_state.config.PDB_TEMPLATE_FILES_NAMES))
        
        if paradigm == "prism-power":
            with st.expander("⚡ PRISM Power Settings (JSON format)"):
                # Handle cases where settings might be a dict or a Pydantic model
                power_settings = st.session_state.config.PRISM_POWER_SETTINGS
                if power_settings:
                    power_data = power_settings.model_dump() if hasattr(power_settings, 'model_dump') else power_settings
                else:
                    power_data = {}
                power_json = st.text_area("Config JSON", value=json.dumps(power_data, indent=2), height=200)
        else:
            power_json = "{}"

    # Advanced / Naming (Expander)
    with st.expander("📝 Naming & Path Constraints (Advanced)"):
        st.warning("Changing these may require re-organizing your input directory.")
        col_adv1, col_adv2 = st.columns(2)
        with col_adv1:
            seq_code = st.text_input("ALIGN_CODE_SEQUENCE", value=st.session_state.config.ALIGN_CODE_SEQUENCE)
            chain_id = st.text_input("CHAIN_ID", value=st.session_state.config.CHAIN_ID)
            blk_chain = st.text_input("BLK_CHAIN_ID", value=st.session_state.config.BLK_CHAIN_ID)
            fasta_base = st.text_input("FASTA_FILE_BASENAME", value=st.session_state.config.FASTA_FILE_BASENAME)
            ss2_base = st.text_input("SS2_FILE_BASENAME", value=st.session_state.config.SS2_FILE_BASENAME)
            manual_ali_base = st.text_input("MANUAL_ALIGNMENT_BASENAME", value=st.session_state.config.MANUAL_ALIGNMENT_BASENAME)
        with col_adv2:
            ini_base = st.text_input("CUSTOM_INIFILE_BASENAME", value=st.session_state.config.CUSTOM_INIFILE_BASENAME)
            rsr_base = st.text_input("CUSTOM_RSRFILE_BASENAME", value=st.session_state.config.CUSTOM_RSRFILE_BASENAME)
            input_dir_name = st.text_input("INPUT_DIR_NAME", value=st.session_state.config.INPUT_DIR_NAME)
            modeling_dir_name = st.text_input("MODELING_RESULTS_DIR_NAME", value=st.session_state.config.MODELING_RESULTS_DIR_NAME)
            psipred_dir_name = st.text_input("PSIPRED_RESULTS_DIR_NAME", value=st.session_state.config.PSIPRED_RESULTS_DIR_NAME)

    st.markdown("---")
    if st.button("💾 Save All Configuration", width='stretch'):
            try:
                # Update main params
                st.session_state.config.TOTAL_HOMOLOGY_MODELS = total_models
                st.session_state.config.TOP_MODELS_FOR_REFINEMENT = top_refine
                st.session_state.config.LOOP_MODELS_PER_TARGET = loop_models
                st.session_state.config.MOBILE_FLANK_RESIDUES = flank_size
                st.session_state.config.BLOCK_REPULSION_RADIUS = repulsion
                st.session_state.config.REFINE_FLANKS_DURING_AUTOMODEL = refine_flanks
                st.session_state.config.NUM_BEST_FINAL_MODELS = int(best_final) if best_final.isdigit() else best_final
                
                st.session_state.config.PERFORM_PSIPRED_PREDICTION = psipred_pred
                st.session_state.config.PSIPRED_EMAIL = psipred_email
                st.session_state.config.PSIPRED_POLL_INTERVAL = psipred_poll
                st.session_state.config.USE_MANUAL_ALIGNMENT = use_manual_ali
                
                st.session_state.config.EXECUTION_PARADIGM = paradigm
                st.session_state.config.MODELLER_CORES = modeller_cores
                st.session_state.config.TOTAL_PARALLEL_JOBS = parallel_jobs
                st.session_state.config.USE_MANUAL_OPTIMIZATION_SELECTION = manual_opt
                
                # Update complex params
                st.session_state.config.MANUAL_OPTIMIZATION_RESIDUES = [int(r) for r in opt_res_str.split()] if opt_res_str else []
                st.session_state.config.PDB_TEMPLATE_FILES_NAMES = [n.strip() for n in pdb_names_str.split("\n") if n.strip()]
                
                if paradigm == "prism-power":
                    power_data = json.loads(power_json)
                    st.session_state.config.PRISM_POWER_SETTINGS = power_data
                
                # Update advanced params
                st.session_state.config.ALIGN_CODE_SEQUENCE = seq_code
                st.session_state.config.CHAIN_ID = chain_id
                st.session_state.config.BLK_CHAIN_ID = blk_chain
                st.session_state.config.FASTA_FILE_BASENAME = fasta_base
                st.session_state.config.SS2_FILE_BASENAME = ss2_base
                st.session_state.config.MANUAL_ALIGNMENT_BASENAME = manual_ali_base
                
                # New advanced params
                st.session_state.config.CUSTOM_INIFILE_BASENAME = ini_base
                st.session_state.config.CUSTOM_RSRFILE_BASENAME = rsr_base
                st.session_state.config.INPUT_DIR_NAME = input_dir_name
                st.session_state.config.MODELING_RESULTS_DIR_NAME = modeling_dir_name
                st.session_state.config.PSIPRED_RESULTS_DIR_NAME = psipred_dir_name

                st.session_state.config.save_settings()
                st.success("✅ Configuration saved and synced to config.yaml")
                st.rerun()
            except Exception as e:
                st.error(f"❌ Error saving config: {e}")

with tab_files:
    st.header("Input Data Management")
    
    col_inv, col_up = st.columns([2, 1])
    
    with col_inv:
        st.subheader("📁 `input/` Directory Inventory")
        files = list_files_in_dir(st.session_state.config.INPUT_DIR)
        if files:
            for f in files:
                c1, c2, c3, c4 = st.columns([3, 1, 2, 1])
                with c1: st.text(f["name"])
                with c2: st.text(f["size"])
                with c3: st.text(f["modified"])
                with c4:
                    if st.button("🗑️", key=f"del_{f['name']}", help=f"Delete {f['name']}"):
                        if delete_file(st.session_state.config.INPUT_DIR, f["name"]):
                            st.success(f"Deleted {f['name']}")
                            st.rerun()
        else:
            st.info("No files found in input/ directory.")

    with col_up:
        st.subheader("📤 Upload Files")
        
        # Mandatory for most runs
        fasta = st.file_uploader("Upload FASTA Sequence", type=["fasta", "fa"])
        if fasta:
            save_uploaded_file(fasta, st.session_state.config.INPUT_DIR)
            st.success(f"Uploaded {fasta.name}")
            
        pdbs = st.file_uploader("Upload Template PDBs", type=["pdb"], accept_multiple_files=True)
        if pdbs:
            for pdb in pdbs:
                save_uploaded_file(pdb, st.session_state.config.INPUT_DIR)
            st.success(f"Uploaded {len(pdbs)} PDB(s)")

        st.divider()
        st.subheader("💡 Dynamic Requirements")
        # Conditional help/status based on config
        if not st.session_state.config.PERFORM_PSIPRED_PREDICTION:
            st.warning("⚠️ Prediction is OFF: A .ss2 file is REQUIRED in input/")
            ss2 = st.file_uploader("Upload SS2 File", type=["ss2"])
            if ss2:
                save_uploaded_file(ss2, st.session_state.config.INPUT_DIR)
                st.success(f"Uploaded {ss2.name}")
        else:
            st.info("ℹ️ Prediction is ON: .ss2 will be generated automatically.")

        if st.session_state.config.USE_MANUAL_ALIGNMENT:
            st.warning("⚠️ Manual Alignment is ON: A .ali file is REQUIRED in input/")
            ali = st.file_uploader("Upload Manual Alignment (.ali)", type=["ali", "pir"])
            if ali:
                save_uploaded_file(ali, st.session_state.config.INPUT_DIR)
                st.success(f"Uploaded {ali.name}")
        else:
            st.info("ℹ️ Manual Alignment is OFF: Alignment will be automated.")

        st.divider()
        other_files = st.file_uploader("Upload Additional Files (General)", accept_multiple_files=True, help="Use this for any other required files.")
        if other_files:
            for f in other_files:
                save_uploaded_file(f, st.session_state.config.INPUT_DIR)
            st.success(f"Uploaded {len(other_files)} additional file(s)")

with tab_viz:
    st.header("3D Structure Explorer")
    
    col_v1, col_v2 = st.columns([1, 3])
    
    with col_v1:
        st.subheader("Controls")
        folder = st.radio("Search Folder", ["input", "modeling_results"])
        target_dir = st.session_state.config.INPUT_DIR if folder == "input" else st.session_state.config.MODELING_RESULTS_DIR
        
        if os.path.exists(target_dir):
            file_list = [f for f in os.listdir(target_dir) if f.endswith(".pdb")]
            selected_pdb = st.selectbox("Select PDB File", file_list)
            
            st.divider()
            style = st.selectbox("Style", ["cartoon", "stick", "sphere", "line", "cross"])
            color = st.selectbox("Color Scheme", ["spectrum", "chain", "element", "ss"])
        else:
            st.error(f"Folder {folder} not found.")
            selected_pdb = None

    with col_v2:
        if selected_pdb:
            st.caption(f"Viewing: {selected_pdb} in {folder}/")
            pdb_path = os.path.join(target_dir, selected_pdb)
            html_data = visualize_pdb(pdb_path, style=style, color=color)
            components.html(html_data, height=600)

with tab_exec:
    st.header("Pipeline Execution")
    st.info("Ensure all files are uploaded and configuration is saved before running.")
    
    if st.button("🚀 Start PRISM Pipeline"):
        process = run_nextflow()
        log_area = st.empty()
        logs = ""
        for line in get_nextflow_progress(process):
            logs += line
            log_area.text_area("Nextflow Real-time Logs", logs, height=500)

with tab_results:
    st.header("Modeling Performance & Results")
    results_df = get_score_distribution_data()
    
    if results_df is not None:
        col_r1, col_r2 = st.columns([1, 1])
        
        with col_r1:
            st.subheader("Top Models Data")
            st.dataframe(results_df, width='stretch')
            
        with col_r2:
            st.subheader("DOPEHR Score Distribution")
            fig = px.bar(results_df, x="Model_Name", y="DOPEHR_score", 
                         title="Score by Model (Lower is Better)",
                         color="DOPEHR_score", color_continuous_scale="Viridis")
            st.plotly_chart(fig, width='stretch')
            
        st.markdown("---")
        st.subheader("Z-Score Analysis")
        fig_z = px.scatter(results_df, x="DOPEHR_score", y="DOPEHR_zscore", 
                           text="Model_Name", title="DOPEHR vs Z-Score")
        st.plotly_chart(fig_z, width='stretch')
    else:
        st.info("No ranking results found. Finish a pipeline run to see analytics.")

with tab_tools:
    st.header("PRISM Toolbox")
    st.markdown("""
        Run utility scripts from the `tools/` directory. 
        **Note:** All scripts are executed from the project root. 
        If a script requires a file argument, use paths relative to the project root (e.g., `input/my_file.pdb`).
    """)
    
    tools_dir = "tools"
    if os.path.exists(tools_dir):
        tools = [f for f in os.listdir(tools_dir) if f.endswith(".py")]
        selected_tool = st.selectbox("Select Tool to Run", tools, help="Scripts located in the /tools folder")
        
        tool_args = st.text_input("Arguments (space separated)", placeholder="e.g. --input input/template.pdb", 
                                 help="Provide arguments as if running from the command line")
        
        if st.button("🛠️ Execute Tool"):
            with st.spinner(f"Running {selected_tool}..."):
                args_list = tool_args.split() if tool_args else []
                output = run_tool(selected_tool, args_list)
                st.code(output, language="text")
    else:
        st.error("Tools directory not found.")
