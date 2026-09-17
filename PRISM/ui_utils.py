#!/usr/bin/env python3
"""Utility functions for the PRISM Streamlit dashboard.

Provides helpers for Nextflow execution, PDB visualisation, file management,
and project-level file discovery.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import sys
from collections.abc import Generator
from pathlib import Path
from typing import Any

import pandas as pd
import py3Dmol

from PRISM.config import settings

logger = logging.getLogger("PRISM.ui_utils")


def run_nextflow() -> subprocess.Popen:
    """Launch the Nextflow pipeline as a subprocess.

    Returns:
        A ``Popen`` handle whose *stdout* streams combined stdout/stderr.
    """
    cmd = ["nextflow", "run", "pipeline/orchestrator.nf", "-params-file", "config.yaml"]
    return subprocess.Popen(  # noqa: S603
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )


def get_nextflow_progress(process: subprocess.Popen) -> Generator[str, None, None]:
    """Yield lines from a running Nextflow process.

    Args:
        process: The ``Popen`` handle returned by :func:`run_nextflow`.

    Yields:
        Each line of Nextflow output as it becomes available.
    """
    yield from process.stdout


def visualize_pdb(
    pdb_path: str | Path,
    style: str = "cartoon",
    color: str = "spectrum",
) -> str:
    """Render a PDB file as an interactive HTML viewer using py3Dmol.

    Args:
        pdb_path: Path to the PDB file.
        style: Rendering style (``cartoon``, ``sphere``, ``stick``,
            ``line``, or ``cross``).
        color: Colour scheme applied to the chosen style.

    Returns:
        An HTML string suitable for embedding via
        ``streamlit.iframe()``.
    """
    pdb_data = Path(pdb_path).read_text(encoding="utf-8")

    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_data, "pdb")

    style_dict: dict[str, Any] = {style: {"color": color}}
    view.setStyle(style_dict)

    view.zoomTo()
    return view._make_html()  # noqa: SLF001


def list_files_in_dir(directory: str | Path) -> list[dict[str, str]]:
    """List files in a directory with human-readable metadata.

    Args:
        directory: Path to the directory to scan.

    Returns:
        A list of dicts, each with keys ``name``, ``size``, and
        ``modified``.  Returns an empty list if the directory does not
        exist.
    """
    path = Path(directory)
    if not path.exists():
        return []

    return [
        {
            "name": f.name,
            "size": f"{f.stat().st_size / 1024:.1f} KB",
            "modified": pd.to_datetime(f.stat().st_mtime, unit="s").strftime(
                "%Y-%m-%d %H:%M:%S"
            ),
        }
        for f in path.iterdir()
        if f.is_file()
    ]


def run_tool(
    script_name: str, args: list[str] | None = None
) -> Generator[str, None, None]:
    """Execute a utility script from the ``tools/`` directory and stream output.

    Output is streamed line-by-line. Standard error is merged into
    standard output.

    Args:
        script_name: Filename of the script inside ``tools/``.
        args: Optional list of CLI arguments to pass to the script.

    Yields:
        Each line of the tool's output as it becomes available.
    """
    project_root = Path(__file__).resolve().parent.parent

    tool_path = project_root / "tools" / script_name
    if not tool_path.exists():
        yield f"Error: Tool {script_name} not found."
        return

    processed_args: list[str] = []
    if args:
        for arg in args:
            is_path_like = arg.startswith("./") or "/" in arg or Path(arg).exists()
            is_flag = arg.startswith(("-", "--"))
            if is_path_like and not is_flag:
                processed_args.append(str(Path(arg).resolve()))
            else:
                processed_args.append(arg)

    cmd = [sys.executable, str(tool_path), *processed_args]

    try:
        process = subprocess.Popen(  # noqa: S603
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=str(project_root),
            bufsize=1,
        )

        if process.stdout:
            yield from process.stdout

        process.wait()

        if process.returncode != 0:
            yield (
                f"\n[ERROR] Tool {script_name} exited with "
                f"return code {process.returncode}\n"
            )
    except Exception as exc:
        logger.exception("Unexpected error running tool %s", script_name)
        yield f"\n[CRITICAL ERROR] Failed to execute tool: {exc}\n"


def get_score_distribution_data() -> pd.DataFrame | None:
    """Read the final ranking CSV and return it as a DataFrame.

    Returns:
        A ``DataFrame`` with model ranking data, or *None* if the CSV
        does not exist.
    """
    csv_path = Path(settings.FINAL_RANKING_CSV)
    if csv_path.exists():
        return pd.read_csv(csv_path)
    return None


def _validate_project_path(path: str | Path) -> Path:
    """Resolve and validate that a path is within the project workspace.

    Args:
        path: The path to validate.

    Returns:
        The resolved Path.

    Raises:
        ValueError: If the path resides outside the project root directory.
    """
    project_root = Path(__file__).resolve().parent.parent
    resolved = Path(path).resolve()
    if project_root not in resolved.parents and resolved != project_root:
        raise ValueError(f"Path traversal detected: {path} resolves outside workspace.")
    return resolved


def delete_path(path: str | Path) -> bool:
    """Delete a file or directory.

    Args:
        path: File or directory to remove.

    Returns:
        *True* if deletion succeeded, *False* otherwise.
    """
    try:
        p = _validate_project_path(path)
    except ValueError:
        logger.exception("Path validation failed")
        return False

    if not p.exists():
        return False

    try:
        if p.is_file():
            p.unlink()
        elif p.is_dir():
            shutil.rmtree(p)
    except Exception:
        logger.exception("Failed to delete %s", path)
        return False
    else:
        return True


def create_directory(path: str | Path) -> bool:
    """Create a directory (including parents) if it does not exist.

    Args:
        path: Directory path to create.

    Returns:
        *True* if creation succeeded, *False* otherwise.
    """
    try:
        Path(path).mkdir(parents=True, exist_ok=True)
    except Exception:
        logger.exception("Failed to create directory %s", path)
        return False
    else:
        return True


def move_path(src: str | Path, dst: str | Path) -> bool:
    """Move a file or directory to a new location.

    Args:
        src: Source path.
        dst: Destination path.

    Returns:
        *True* if the move succeeded, *False* otherwise.
    """
    try:
        s = _validate_project_path(src)
        d = _validate_project_path(dst)
    except ValueError:
        logger.exception("Path validation failed")
        return False

    try:
        shutil.move(str(s), str(d))
    except Exception:
        logger.exception("Failed to move %s -> %s", src, dst)
        return False
    else:
        return True


def copy_path(src: str | Path, dst: str | Path) -> bool:
    """Copy a file or directory to a new location.

    Args:
        src: Source path.
        dst: Destination path.

    Returns:
        *True* if the copy succeeded, *False* otherwise.
    """
    try:
        s = _validate_project_path(src)
        d = _validate_project_path(dst)
    except ValueError:
        logger.exception("Path validation failed")
        return False

    try:
        if s.is_dir():
            shutil.copytree(s, d)
        else:
            shutil.copy2(s, d)
    except Exception:
        logger.exception("Failed to copy %s -> %s", src, dst)
        return False
    else:
        return True


def list_root_dirs() -> list[str]:
    """List non-hidden directories in the project root.

    Returns:
        Sorted list of directory names.
    """
    root = Path()
    dirs = [d.name for d in root.iterdir() if d.is_dir() and not d.name.startswith(".")]
    return sorted(dirs)


def get_all_project_files() -> list[str]:
    """Recursively list all project files for autocompletion.

    Excludes hidden directories and common internal directories
    (``.pixi``, ``__pycache__``, ``work``, etc.).

    Returns:
        Sorted list of relative file paths.
    """
    exclude_dirs = {".pixi", ".git", "__pycache__", ".snakemake", ".nextflow", "work"}
    all_files: list[str] = []
    root = Path()

    for item in root.rglob("*"):
        if item.is_dir():
            continue
        parts = item.parts
        if any(p in exclude_dirs or p.startswith(".") for p in parts):
            continue
        if item.name.startswith("."):
            continue
        all_files.append(str(item))

    return sorted(all_files)


def save_uploaded_file(uploaded_file: Any, target_dir: str | Path) -> str:
    """Save a Streamlit uploaded file to the target directory.

    Args:
        uploaded_file: A Streamlit ``UploadedFile`` object.
        target_dir: Directory in which to save the file.

    Returns:
        The full path to the saved file.
    """
    target = _validate_project_path(target_dir)
    target.mkdir(parents=True, exist_ok=True)
    file_path = _validate_project_path(target / uploaded_file.name)
    file_path.write_bytes(uploaded_file.getbuffer())
    return str(file_path)
