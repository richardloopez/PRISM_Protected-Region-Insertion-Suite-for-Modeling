#!/usr/bin/env python3
"""PRISM Dashboard: Main graphical user interface for the PRISM modeling suite.

Built with Streamlit for a user-friendly, high-fidelity modeling
experience.  Provides configuration editing, file management, tool
execution, 3D structure visualization, pipeline orchestration, and
results analysis.
"""

from __future__ import annotations

import importlib.util
import json
import logging
import os
import sys
from pathlib import Path

import plotly.express as px
import streamlit as st


def _load_module(name: str, path: Path):
    """Load a module dynamically while ensuring package context for relative imports.

    Args:
        name: Full module name (e.g., 'PRISM.config').
        path: Path to the .py file.

    Returns:
        The loaded module object.
    """
    import types  # noqa: PLC0415

    if "." in name:
        pkg_name = name.rpartition(".")[0]
        if pkg_name not in sys.modules:
            pkg_mod = types.ModuleType(pkg_name)
            pkg_mod.__path__ = [str(path.parent)]
            sys.modules[pkg_name] = pkg_mod

    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load spec for {name} at {path}")

    mod = importlib.util.module_from_spec(spec)
    mod.__package__ = name.rpartition(".")[0] if "." in name else ""
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_PRISM_DIR = Path(__file__).resolve().parent
_config = _load_module("PRISM.config", _PRISM_DIR / "config.py")
_ui_utils = _load_module("PRISM.ui_utils", _PRISM_DIR / "ui_utils.py")
_icons = _load_module("PRISM.ui_icons", _PRISM_DIR / "ui_icons.py")

load_settings = _config.load_settings
(
    copy_path,
    create_directory,
    delete_path,
    get_all_project_files,
    get_nextflow_progress,
    get_score_distribution_data,
    list_files_in_dir,
    list_root_dirs,
    move_path,
    run_nextflow,
    run_tool,
    save_uploaded_file,
    visualize_pdb,
) = (
    _ui_utils.copy_path,
    _ui_utils.create_directory,
    _ui_utils.delete_path,
    _ui_utils.get_all_project_files,
    _ui_utils.get_nextflow_progress,
    _ui_utils.get_score_distribution_data,
    _ui_utils.list_files_in_dir,
    _ui_utils.list_root_dirs,
    _ui_utils.move_path,
    _ui_utils.run_nextflow,
    _ui_utils.run_tool,
    _ui_utils.save_uploaded_file,
    _ui_utils.visualize_pdb,
)

# Map icons
(
    ICON_CONFIG,
    ICON_FILE_MGMT,
    ICON_TOOLS,
    ICON_INPUT_FILES,
    ICON_VISUALIZATION,
    ICON_EXECUTION,
    ICON_RESULTS,
    ICON_PRISM,
    ICON_CORE,
    ICON_COMPUTING,
    ICON_MOLECULE,
    ICON_TOGGLES,
    ICON_POWER,
    ICON_WEB,
    ICON_SAVE,
    ICON_SUCCESS,
    ICON_ERROR,
    ICON_DELETE,
    ICON_MOVE,
    ICON_COPY,
    ICON_ADD,
    ICON_SEARCH,
    ICON_UPLOAD,
    ICON_FOLDER,
    ICON_FILE,
    ICON_IDEA,
    ICON_WARNING,
    ICON_INFO,
) = (
    _icons.ICON_CONFIG,
    _icons.ICON_FILE_MGMT,
    _icons.ICON_TOOLS,
    _icons.ICON_INPUT_FILES,
    _icons.ICON_VISUALIZATION,
    _icons.ICON_EXECUTION,
    _icons.ICON_RESULTS,
    _icons.ICON_PRISM,
    _icons.ICON_CORE,
    _icons.ICON_COMPUTING,
    _icons.ICON_MOLECULE,
    _icons.ICON_TOGGLES,
    _icons.ICON_POWER,
    _icons.ICON_WEB,
    _icons.ICON_SAVE,
    _icons.ICON_SUCCESS,
    _icons.ICON_ERROR,
    _icons.ICON_DELETE,
    _icons.ICON_MOVE,
    _icons.ICON_COPY,
    _icons.ICON_ADD,
    _icons.ICON_SEARCH,
    _icons.ICON_UPLOAD,
    _icons.ICON_FOLDER,
    _icons.ICON_FILE,
    _icons.ICON_IDEA,
    _icons.ICON_WARNING,
    _icons.ICON_INFO,
)

logger = logging.getLogger("PRISM_Dashboard")

# --- GUI TAB CONFIGURATION ---
GUI_TAB_ORDER = [
    {"id": "config", "label": f"{ICON_CONFIG} Config"},
    {"id": "file_mgmt", "label": f"{ICON_FILE_MGMT} File Management"},
    {"id": "tools", "label": f"{ICON_TOOLS} Tools"},
    {"id": "files", "label": f"{ICON_INPUT_FILES} Input Files"},
    {"id": "viz", "label": f"{ICON_VISUALIZATION} Visualization"},
    {"id": "exec", "label": f"{ICON_EXECUTION} Execution"},
    {"id": "results", "label": f"{ICON_RESULTS} Results"},
]

st.set_page_config(page_title="PRISM Dashboard", layout="wide")

# Dynamic Logo placement
LOGO_PATH = Path(__file__).parent / "logo.png"
if LOGO_PATH.exists():
    st.logo(str(LOGO_PATH), icon_image=None)
    col_logo, col_title = st.columns([1, 10])
    with col_logo:
        st.image(str(LOGO_PATH), width=60)
    with col_title:
        st.title("PRISM: Protected-Region Insertion Suite for Modeling")
else:
    st.title(f"{ICON_PRISM} PRISM: Protected-Region Insertion Suite for Modeling")

st.markdown("---")

if "config" not in st.session_state:
    try:
        st.session_state.config = load_settings()
    except (FileNotFoundError, ValueError, OSError) as exc:
        st.error(f"Could not load config.yaml. Please ensure it exists. Error: {exc}")
        st.stop()

# --- SIDEBAR (THEME SETTINGS) ---
theme = st.sidebar.select_slider(
    "Select GUI Theme",
    options=["Default", "Professional", "High Contrast", "Dark Modern"],
)

# Load external CSS themes
_css_path = Path(__file__).parent / "styles" / "themes.css"
_css_content = _css_path.read_text(encoding="utf-8") if _css_path.exists() else ""

THEME_CSS_CLASS = {
    "Professional": "theme-professional",
    "High Contrast": "theme-high-contrast",
    "Dark Modern": "theme-dark-modern",
}

if _css_content:
    if theme in THEME_CSS_CLASS:
        css_class = THEME_CSS_CLASS[theme]
        active_css = _css_content.replace(f".{css_class}.stApp", ".stApp")
        active_css = active_css.replace(f".{css_class} ", ".stApp ")
        active_css = active_css.replace(f".{css_class}:", ".stApp:")
        st.markdown(f"<style>{active_css}</style>", unsafe_allow_html=True)
    else:
        st.markdown(f"<style>{_css_content}</style>", unsafe_allow_html=True)

st.sidebar.markdown("---")
st.sidebar.caption("PRISM v1.0.0 | Professional Edition")

# --- MAIN NAVIGATION ---
tabs = st.tabs([tab["label"] for tab in GUI_TAB_ORDER])
tab_map = {tab["id"]: tabs[i] for i, tab in enumerate(GUI_TAB_ORDER)}

tab_config = tab_map["config"]
tab_files = tab_map["files"]
tab_viz = tab_map["viz"]
tab_exec = tab_map["exec"]
tab_results = tab_map["results"]
tab_tools = tab_map["tools"]
tab_file_mgmt = tab_map["file_mgmt"]

# ============================================================================
# TAB: CONFIGURATION  # noqa: ERA001
# ============================================================================
with tab_config:
    st.header("Global Pipeline Configuration")

    # 1. Core Settings
    with st.container(border=True):
        st.subheader(f"{ICON_CORE} Core Settings")
        col_c1, col_c2, col_c3 = st.columns(3)
        with col_c1:
            SEQ_CODE = st.text_input(
                "ALIGN_CODE_SEQUENCE", value=st.session_state.config.ALIGN_CODE_SEQUENCE
            )
        with col_c2:
            CHAIN_ID = st.text_input("CHAIN_ID", value=st.session_state.config.CHAIN_ID)
        with col_c3:
            BLK_CHAIN = st.text_input(
                "BLK_CHAIN_ID", value=st.session_state.config.BLK_CHAIN_ID
            )

    # 2. Execution & Parallelization
    with st.container(border=True):
        st.subheader(f"{ICON_COMPUTING} Execution & Parallelization")
        col_e1, col_e2 = st.columns(2)
        with col_e1:
            MODELLER_CORES = st.number_input(
                "MODELLER_CORES",
                value=st.session_state.config.MODELLER_CORES,
                min_value=1,
            )
        with col_e2:
            PARALLEL_JOBS = st.number_input(
                "TOTAL_PARALLEL_JOBS",
                value=st.session_state.config.TOTAL_PARALLEL_JOBS,
                min_value=1,
            )

    # 3. Modeling Parameters
    with st.container(border=True):
        st.subheader(f"{ICON_MOLECULE} Modeling Parameters")
        col_m1, col_m2 = st.columns(2)
        with col_m1:
            TOTAL_MODELS = st.number_input(
                "TOTAL_HOMOLOGY_MODELS",
                value=st.session_state.config.TOTAL_HOMOLOGY_MODELS,
                min_value=1,
            )
            LOOP_MODELS = st.number_input(
                "LOOP_MODELS_PER_TARGET",
                value=st.session_state.config.LOOP_MODELS_PER_TARGET,
                min_value=1,
            )
        with col_m2:
            TOP_REFINE = st.number_input(
                "TOP_MODELS_FOR_REFINEMENT",
                value=st.session_state.config.TOP_MODELS_FOR_REFINEMENT,
                min_value=0,
            )
            BEST_FINAL = st.text_input(
                "NUM_BEST_FINAL_MODELS",
                value=str(st.session_state.config.NUM_BEST_FINAL_MODELS),
            )
            REPULSION = st.text_input(
                "BLOCK_REPULSION_RADIUS",
                value=str(st.session_state.config.BLOCK_REPULSION_RADIUS),
            )

    # 4. Feature Toggles
    with st.container(border=True):
        st.subheader(f"{ICON_TOGGLES} Feature Toggles")
        col_t1, col_t2 = st.columns(2)
        with col_t1:
            USE_MANUAL_ALI = st.toggle(
                "USE_MANUAL_ALIGNMENT",
                value=st.session_state.config.USE_MANUAL_ALIGNMENT,
            )
            REFINE_FLANKS = st.toggle(
                "REFINE_FLANKS_DURING_AUTOMODEL",
                value=st.session_state.config.REFINE_FLANKS_DURING_AUTOMODEL,
            )
        with col_t2:
            MANUAL_OPT = st.toggle(
                "USE_MANUAL_OPTIMIZATION_SELECTION",
                value=st.session_state.config.USE_MANUAL_OPTIMIZATION_SELECTION,
            )
            MANUAL_FIX = st.toggle(
                "USE_MANUAL_FIXATION_SELECTION",
                value=st.session_state.config.USE_MANUAL_FIXATION_SELECTION,
            )

    # 5. Advanced Modifiers (Context-sensitive based on toggles)
    with st.container(border=True):
        st.subheader(f"{ICON_TOOLS} Advanced PRISM Modifiers")
        FLANK_SIZE = st.number_input(
            "MOBILE_FLANK_RESIDUES",
            value=st.session_state.config.MOBILE_FLANK_RESIDUES,
            min_value=0,
        )

        OPT_RES_STR = ""
        if MANUAL_OPT:
            OPT_RES_STR = st.text_area(
                "MANUAL_OPTIMIZATION_RESIDUES (Space separated)",
                value=" ".join(
                    map(str, st.session_state.config.MANUAL_OPTIMIZATION_RESIDUES)
                ),
            )

        FIX_RES_STR = ""
        if MANUAL_FIX:
            FIX_RES_STR = st.text_area(
                "MANUAL_FIXATION_RESIDUES (Space separated)",
                value=" ".join(
                    map(str, st.session_state.config.MANUAL_FIXATION_RESIDUES)
                ),
            )

    # 6. Paradigm Settings
    with st.container(border=True):
        st.subheader(f"{ICON_POWER} Paradigm Settings")
        PARADIGM = st.selectbox(
            "EXECUTION_PARADIGM",
            ["prism-power", "precalculation", "precomputed", "normal"],
            index=["prism-power", "precalculation", "precomputed", "normal"].index(
                st.session_state.config.EXECUTION_PARADIGM
            ),
        )
        PDB_NAMES_STR = st.text_area(
            "PDB_TEMPLATE_FILES_NAMES (One per line, first is MAIN)",
            value="\n".join(st.session_state.config.PDB_TEMPLATE_FILES_NAMES),
        )

        if PARADIGM == "prism-power":
            POWER_SETTINGS = st.session_state.config.PRISM_POWER_SETTINGS
            POWER_DATA = (
                POWER_SETTINGS.model_dump()
                if hasattr(POWER_SETTINGS, "model_dump")
                else POWER_SETTINGS
            )
            POWER_JSON = st.text_area(
                "PRISM_POWER_SETTINGS (JSON)",
                value=json.dumps(POWER_DATA, indent=2),
                height=200,
            )
        else:
            POWER_JSON = "{}"

    # 7. External Services
    with st.container(border=True):
        st.subheader(f"{ICON_WEB} External Services (PSIPRED)")
        PSIPRED_PRED = st.toggle(
            "PERFORM_PSIPRED_PREDICTION",
            value=st.session_state.config.PERFORM_PSIPRED_PREDICTION,
        )
        col_ps1, col_ps2 = st.columns(2)
        with col_ps1:
            PSIPRED_EMAIL = st.text_input(
                "PSIPRED_EMAIL", value=st.session_state.config.PSIPRED_EMAIL
            )
        with col_ps2:
            PSIPRED_POLL = st.number_input(
                "PSIPRED_POLL_INTERVAL",
                value=st.session_state.config.PSIPRED_POLL_INTERVAL,
                min_value=1,
            )

    # 8. I/O Directories & Files
    with st.expander(f"{ICON_FOLDER} I/O Directories & Files (Advanced)"):
        col_io1, col_io2 = st.columns(2)
        with col_io1:
            INPUT_DIR_NAME = st.text_input(
                "INPUT_DIR_NAME", value=st.session_state.config.INPUT_DIR_NAME
            )
            MODELING_DIR_NAME = st.text_input(
                "MODELING_RESULTS_DIR_NAME",
                value=st.session_state.config.MODELING_RESULTS_DIR_NAME,
            )
            PSIPRED_DIR_NAME = st.text_input(
                "PSIPRED_RESULTS_DIR_NAME",
                value=st.session_state.config.PSIPRED_RESULTS_DIR_NAME,
            )
        with col_io2:
            FASTA_BASE = st.text_input(
                "FASTA_FILE_BASENAME", value=st.session_state.config.FASTA_FILE_BASENAME
            )
            SS2_BASE = st.text_input(
                "SS2_FILE_BASENAME", value=st.session_state.config.SS2_FILE_BASENAME
            )
            MANUAL_ALI_BASE = st.text_input(
                "MANUAL_ALIGNMENT_BASENAME",
                value=st.session_state.config.MANUAL_ALIGNMENT_BASENAME,
            )
            INI_BASE = st.text_input(
                "CUSTOM_INIFILE_BASENAME",
                value=st.session_state.config.CUSTOM_INIFILE_BASENAME,
            )
            RSR_BASE = st.text_input(
                "CUSTOM_RSRFILE_BASENAME",
                value=st.session_state.config.CUSTOM_RSRFILE_BASENAME,
            )

    st.markdown("---")
    col_btn_reload, col_btn_save = st.columns(2)
    with col_btn_reload:
        reload_clicked = st.button(
            f"{ICON_CONFIG} Reload Config from Disk",
            width="stretch",
            help="Reload config from disk if modified externally.",
        )
    with col_btn_save:
        save_clicked = st.button(f"{ICON_SAVE} Save All Configuration", width="stretch")

    if reload_clicked:
        try:
            st.session_state.config = load_settings()
            st.rerun()
        except Exception as e:
            st.error(f"Error reloading config: {e}")

    if save_clicked:
        try:
            # Core
            st.session_state.config.ALIGN_CODE_SEQUENCE = SEQ_CODE
            st.session_state.config.CHAIN_ID = CHAIN_ID
            st.session_state.config.BLK_CHAIN_ID = BLK_CHAIN

            # Exec
            st.session_state.config.MODELLER_CORES = MODELLER_CORES
            st.session_state.config.TOTAL_PARALLEL_JOBS = PARALLEL_JOBS

            # Modeling
            st.session_state.config.TOTAL_HOMOLOGY_MODELS = TOTAL_MODELS
            st.session_state.config.TOP_MODELS_FOR_REFINEMENT = TOP_REFINE
            st.session_state.config.LOOP_MODELS_PER_TARGET = LOOP_MODELS

            # Parse BLOCK_REPULSION_RADIUS allowing strings ("remember-...") or floats
            try:
                rep_val = float(REPULSION)
            except ValueError:
                rep_val = REPULSION
            st.session_state.config.BLOCK_REPULSION_RADIUS = rep_val

            if BEST_FINAL.isdigit():
                st.session_state.config.NUM_BEST_FINAL_MODELS = int(BEST_FINAL)
            else:
                st.session_state.config.NUM_BEST_FINAL_MODELS = BEST_FINAL

            # Features
            st.session_state.config.USE_MANUAL_ALIGNMENT = USE_MANUAL_ALI
            st.session_state.config.REFINE_FLANKS_DURING_AUTOMODEL = REFINE_FLANKS
            st.session_state.config.USE_MANUAL_OPTIMIZATION_SELECTION = MANUAL_OPT
            st.session_state.config.USE_MANUAL_FIXATION_SELECTION = MANUAL_FIX

            # Advanced
            st.session_state.config.MOBILE_FLANK_RESIDUES = FLANK_SIZE
            st.session_state.config.MANUAL_OPTIMIZATION_RESIDUES = (
                [int(r) for r in OPT_RES_STR.split()] if OPT_RES_STR else []
            )
            st.session_state.config.MANUAL_FIXATION_RESIDUES = (
                [int(r) for r in FIX_RES_STR.split()] if FIX_RES_STR else []
            )

            # Paradigm
            st.session_state.config.EXECUTION_PARADIGM = PARADIGM
            st.session_state.config.PDB_TEMPLATE_FILES_NAMES = [
                n.strip() for n in PDB_NAMES_STR.split("\n") if n.strip()
            ]
            if PARADIGM == "prism-power":
                st.session_state.config.PRISM_POWER_SETTINGS = json.loads(POWER_JSON)

            # Ex Services
            st.session_state.config.PERFORM_PSIPRED_PREDICTION = PSIPRED_PRED
            st.session_state.config.PSIPRED_EMAIL = PSIPRED_EMAIL
            st.session_state.config.PSIPRED_POLL_INTERVAL = PSIPRED_POLL

            # I/O
            st.session_state.config.INPUT_DIR_NAME = INPUT_DIR_NAME
            st.session_state.config.MODELING_RESULTS_DIR_NAME = MODELING_DIR_NAME
            st.session_state.config.PSIPRED_RESULTS_DIR_NAME = PSIPRED_DIR_NAME
            st.session_state.config.FASTA_FILE_BASENAME = FASTA_BASE
            st.session_state.config.SS2_FILE_BASENAME = SS2_BASE
            st.session_state.config.MANUAL_ALIGNMENT_BASENAME = MANUAL_ALI_BASE
            st.session_state.config.CUSTOM_INIFILE_BASENAME = INI_BASE
            st.session_state.config.CUSTOM_RSRFILE_BASENAME = RSR_BASE

            st.session_state.config.save_settings()
            st.success(f"{ICON_SUCCESS} Configuration saved and synced to config.yaml")
            st.rerun()
        except Exception as exc:
            logger.exception("Error saving config")
            st.error(f"{ICON_ERROR} Error saving config: {exc}")


# ============================================================================
# TAB: FILE MANAGEMENT
# ============================================================================
with tab_file_mgmt:
    st.markdown(f"### {ICON_FILE_MGMT} File Management & Inputs")
    st.write(
        "Upload your structures (PDB) and sequence (FASTA) files here. "
        "These will be stored in the `input/` directory."
    )

    with st.expander(f"{ICON_ADD} Create New Folder", expanded=False):
        col_new1, col_new2 = st.columns([3, 1], vertical_alignment="bottom")
        with col_new1:
            new_folder_name = st.text_input(
                "Folder Name", placeholder="e.g. custom_results"
            )
        with col_new2:
            if st.button("Create", width="stretch"):
                if new_folder_name:
                    if create_directory(new_folder_name):
                        st.success(f"Folder '{new_folder_name}' created!")
                        st.rerun()
                    else:
                        st.error("Failed to create folder.")
                else:
                    st.warning("Please enter a name.")

    st.divider()

    root_dirs = list_root_dirs()

    if "selected_paths" not in st.session_state:
        st.session_state.selected_paths = []

    action_col1, action_col2, action_col3, action_col4 = st.columns(
        [1, 1, 1, 3], vertical_alignment="bottom"
    )

    with action_col4:
        destination = st.selectbox(
            "Target Folder for Move/Copy",
            [".", *root_dirs],
            help="Select destination for move or copy actions",
        )

    with action_col1:
        if st.button(
            f"{ICON_DELETE} Delete",
            type="primary",
            width="stretch",
            help="Delete selected items",
        ):
            if st.session_state.selected_paths:
                deleted_count = sum(
                    1 for p in st.session_state.selected_paths if delete_path(p)
                )
                st.session_state.selected_paths = []
                st.success(f"Deleted {deleted_count} items.")
                st.rerun()
            else:
                st.warning("No items selected.")

    with action_col2:
        if st.button(
            f"{ICON_MOVE} Move",
            width="stretch",
            help="Move selected items to destination",
        ):
            if st.session_state.selected_paths and destination:
                moved_count = 0
                for path in st.session_state.selected_paths:
                    if destination != "." and destination in path:
                        st.error(
                            f"Cannot move {path} into its own subfolder {destination}"
                        )
                        continue
                    target = str(Path(destination) / Path(path).name)
                    if move_path(path, target):
                        moved_count += 1
                st.session_state.selected_paths = []
                st.success(f"Moved {moved_count} items to {destination}.")
                st.rerun()
            else:
                st.warning("Select items and a destination.")

    with action_col3:
        if st.button(
            f"{ICON_COPY} Copy",
            width="stretch",
            help="Copy selected items to destination",
        ):
            if st.session_state.selected_paths and destination:
                copied_count = 0
                for path in st.session_state.selected_paths:
                    target_path = Path(destination) / Path(path).name
                    if target_path.exists():
                        target_path = Path(destination) / f"copy_{Path(path).name}"
                    if copy_path(path, target_path):
                        copied_count += 1
                st.session_state.selected_paths = []
                st.success(f"Copied {copied_count} items to {destination}.")
                st.rerun()
            else:
                st.warning("Select items and a destination.")

    if root_dirs:
        cols = st.columns(2)
        for i, directory in enumerate(root_dirs):
            with cols[i % 2], st.container(border=True):
                st.subheader(f"{ICON_FOLDER} {directory}")
                if st.checkbox(
                    f"Select folder: {directory}", key=f"sel_dir_{directory}"
                ):
                    if directory not in st.session_state.selected_paths:
                        st.session_state.selected_paths.append(directory)
                elif directory in st.session_state.selected_paths:
                    st.session_state.selected_paths.remove(directory)

                path_obj = Path(directory)
                dir_files = list(path_obj.iterdir()) if path_obj.exists() else []

                if dir_files:
                    for f in sorted(dir_files):
                        if f.name.startswith("."):
                            continue

                        is_dir = f.is_dir()
                        icon = ICON_FOLDER if is_dir else ICON_FILE

                        path_str = str(f)
                        if st.checkbox(f"{icon} {f.name}", key=f"sel_f_{path_str}"):
                            if path_str not in st.session_state.selected_paths:
                                st.session_state.selected_paths.append(path_str)
                        elif path_str in st.session_state.selected_paths:
                            st.session_state.selected_paths.remove(path_str)
                else:
                    st.caption("Empty folder")
    else:
        st.subheader("Project Root")
        st.info("No directories found in the project root.")


# ============================================================================
# TAB: TOOLS  # noqa: ERA001
# ============================================================================
with tab_tools:
    st.header("PRISM Toolbox")
    st.markdown(
        """
        Run utility scripts from the `tools/` directory.
        **Note:** All scripts are executed from the project root.
        If a script requires a file argument, use paths relative to the
        project root (e.g., `input/my_file.pdb`).
        """
    )

    tools_dir_path = Path("tools")
    if tools_dir_path.exists():
        tools = [f.name for f in tools_dir_path.glob("*.py")]
        selected_tool = st.selectbox(
            "Select Tool to Run", tools, help="Scripts located in the /tools folder"
        )

        if "tool_args_str" not in st.session_state:
            st.session_state.tool_args_str = ""
        if "do_autocomplete" not in st.session_state:
            st.session_state.do_autocomplete = False

        def autocomplete_args() -> None:
            """Autocomplete file paths in tool arguments."""
            typed = st.session_state.get("args_input_key", "")
            if not typed:
                return

            words = typed.split()
            if not words:
                return

            last_word = words[-1]
            proj_files = get_all_project_files()
            matches = [f for f in proj_files if f.startswith(last_word)]

            if len(matches) == 1:
                words[-1] = matches[0]
                new_val = " ".join(words) + " "
                st.session_state.args_input_key = new_val
                st.session_state.tool_args_str = new_val
                st.toast(
                    f"{ICON_SUCCESS} Path completed: {matches[0]}", icon=ICON_TOOLS
                )
            elif len(matches) > 1:
                common_prefix = os.path.commonprefix(matches)
                if common_prefix and len(common_prefix) > len(last_word):
                    words[-1] = common_prefix
                    new_val = " ".join(words)
                    st.session_state.args_input_key = new_val
                    st.session_state.tool_args_str = new_val
                st.info(
                    f"{ICON_IDEA} Multiple matches found: {', '.join(matches[:10])}..."
                )
            else:
                st.session_state.tool_args_str = typed

        if st.session_state.do_autocomplete:
            autocomplete_args()
            st.session_state.do_autocomplete = False

        col_args, col_comp = st.columns([4, 1], vertical_alignment="bottom")

        with col_args:
            st.text_input(
                "Tool Arguments",
                key="args_input_key",
                placeholder="e.g. prep input/my_file.pdb A B",
                help="Type your command. Click 'Complete' to autocomplete file paths.",
            )

        with col_comp:
            if st.button(
                f"{ICON_SEARCH} Complete",
                width="stretch",
                help="Autocomplete the last word",
            ):
                st.session_state.do_autocomplete = True
                st.rerun()

        current_args = st.session_state.get(
            "tool_args_str", ""
        ) or st.session_state.get("args_input_key", "")
        st.caption(f"{ICON_EXECUTION} Final Command: `{selected_tool} {current_args}`")

        if st.button(f"{ICON_TOOLS} Execute Tool", width="stretch", type="primary"):
            final_args = (
                st.session_state.args_input_key.split()
                if st.session_state.args_input_key
                else []
            )
            log_area = st.empty()
            full_output = ""
            for line in run_tool(selected_tool, final_args):
                full_output += line
                log_area.code(full_output, language="text")
    else:
        st.error("Tools directory not found.")


# ============================================================================
# TAB: INPUT FILES
# ============================================================================
with tab_files:
    st.header("Input Data Management")

    col_inv, col_up = st.columns([2, 1])

    with col_inv:
        st.subheader(f"{ICON_FOLDER} `input/` Directory Inventory")
        files = list_files_in_dir(st.session_state.config.INPUT_DIR)
        if files:
            for f in files:
                c1, c2, c3, c4 = st.columns([3, 1, 2, 1])
                with c1:
                    st.text(f["name"])
                with c2:
                    st.text(f["size"])
                with c3:
                    st.text(f["modified"])
                with c4:
                    if st.button(
                        ICON_DELETE, key=f"del_{f['name']}", help=f"Delete {f['name']}"
                    ):
                        target = Path(st.session_state.config.INPUT_DIR) / f["name"]
                        if delete_path(target):
                            st.success(f"Deleted {f['name']}")
                            st.rerun()
        else:
            st.info("No files found in input/ directory.")

    with col_up:
        st.subheader(f"{ICON_UPLOAD} Upload Files")

        fasta = st.file_uploader("Upload FASTA Sequence", type=["fasta", "fa"])
        if fasta:
            save_uploaded_file(fasta, st.session_state.config.INPUT_DIR)
            st.success(f"Uploaded {fasta.name}")

        pdbs = st.file_uploader(
            "Upload Template PDBs", type=["pdb"], accept_multiple_files=True
        )
        if pdbs:
            for pdb in pdbs:
                save_uploaded_file(pdb, st.session_state.config.INPUT_DIR)
            st.success(f"Uploaded {len(pdbs)} PDB(s)")

        st.divider()
        st.subheader(f"{ICON_IDEA} Dynamic Requirements")
        if not st.session_state.config.PERFORM_PSIPRED_PREDICTION:
            st.warning(
                f"{ICON_WARNING} Prediction is OFF: A .ss2 file is REQUIRED in input/"
            )
            ss2 = st.file_uploader("Upload SS2 File", type=["ss2"])
            if ss2:
                save_uploaded_file(ss2, st.session_state.config.INPUT_DIR)
                st.success(f"Uploaded {ss2.name}")
        else:
            st.info(
                f"{ICON_INFO} Prediction is ON: .ss2 will be generated automatically."
            )

        if st.session_state.config.USE_MANUAL_ALIGNMENT:
            st.warning(
                f"{ICON_WARNING} Manual Alignment is ON: A .ali file is "
                "REQUIRED in input/"
            )
            ali = st.file_uploader(
                "Upload Manual Alignment (.ali)", type=["ali", "pir"]
            )
            if ali:
                save_uploaded_file(ali, st.session_state.config.INPUT_DIR)
                st.success(f"Uploaded {ali.name}")
        else:
            st.info(
                f"{ICON_INFO} Manual Alignment is OFF: Alignment will be automated."
            )

        st.divider()
        other_files = st.file_uploader(
            "Upload Additional Files", accept_multiple_files=True
        )
        if other_files:
            for f in other_files:
                save_uploaded_file(f, st.session_state.config.INPUT_DIR)
            st.success(f"Uploaded {len(other_files)} additional file(s)")


# ============================================================================
# TAB: VISUALIZATION  # noqa: ERA001
# ============================================================================
with tab_viz:
    st.header("3D Structure Explorer")

    col_v1, col_v2 = st.columns([1, 3])

    with col_v1:
        st.subheader("Controls")
        folder = st.radio("Search Folder", ["input", "modeling_results"])

        target_dir = (
            st.session_state.config.INPUT_DIR
            if folder == "input"
            else st.session_state.config.MODELING_RESULTS_DIR
        )
        target_path = Path(target_dir)

        if target_path.exists():
            file_list = [f.name for f in target_path.glob("*.pdb")]
            selected_pdb = (
                st.selectbox("Select PDB File", file_list) if file_list else None
            )

            st.divider()
            style = st.selectbox(
                "Style", ["cartoon", "stick", "sphere", "line", "cross"]
            )
            color = st.selectbox("Color Scheme", ["spectrum", "chain", "element", "ss"])
        else:
            st.error(f"Folder {folder} not found.")
            selected_pdb = None

    with col_v2:
        if selected_pdb:
            st.caption(f"Viewing: {selected_pdb} in {folder}/")
            pdb_path = target_path / selected_pdb
            if pdb_path.exists():
                html_data = visualize_pdb(str(pdb_path), style=style, color=color)
                st.iframe(html_data, height=600)


# ============================================================================
# TAB: EXECUTION  # noqa: ERA001
# ============================================================================
with tab_exec:
    st.header("Pipeline Execution")
    st.info("Ensure all files are uploaded and configuration is saved before running.")

    if st.button(f"{ICON_EXECUTION} Start PRISM Pipeline"):
        process = run_nextflow()
        log_area = st.empty()
        logs = ""
        for line in get_nextflow_progress(process):
            logs += line
            log_area.text_area("Nextflow Real-time Logs", logs, height=500)


# ============================================================================
# TAB: RESULTS  # noqa: ERA001
# ============================================================================
with tab_results:
    st.header("Modeling Performance & Results")
    results_df = get_score_distribution_data()

    if results_df is not None and not results_df.empty:
        col_r1, col_r2 = st.columns([1, 1])

        with col_r1:
            st.subheader("Top Models Data")
            st.dataframe(results_df, width="stretch")

        with col_r2:
            st.subheader("DOPEHR Score Distribution")
            fig = px.bar(
                results_df,
                x="Model_Name",
                y="DOPEHR_score",
                title="Score by Model (Lower is Better)",
                color="DOPEHR_score",
                color_continuous_scale="Viridis",
            )
            st.plotly_chart(fig, width="stretch")

        st.markdown("---")
        st.subheader("Z-Score Analysis")
        fig_z = px.scatter(
            results_df,
            x="DOPEHR_score",
            y="DOPEHR_zscore",
            text="Model_Name",
            title="DOPEHR vs Z-Score",
        )
        st.plotly_chart(fig_z, width="stretch")
    else:
        st.info("No ranking results found. Finish a pipeline run to see analytics.")
