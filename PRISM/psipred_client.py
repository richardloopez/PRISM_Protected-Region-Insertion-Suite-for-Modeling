#!/usr/bin/env python3
"""PSIPRED REST API client.

Handles submission of FASTA sequences to the UCL PSIPRED web server,
polls for job completion, and retrieves secondary structure (``.ss2``) files.
"""

from __future__ import annotations

import logging
import shutil
import time
from pathlib import Path
from typing import Any

import requests

from . import config

logger = logging.getLogger("PRISM.psipred_client")

URL_SUBMIT = "https://bioinf.cs.ucl.ac.uk/psipred/api/submission/"
URL_DOWNLOAD_BASE = "https://bioinf.cs.ucl.ac.uk/psipred/api/submissions/"
HEADERS = {"Accept": "application/json"}

MAX_POLL_SECONDS: int = 5 * 3600
"""Maximum time (in seconds) to wait for a PSIPRED job before timing out."""


def submit_job(fasta_path: Path, email: str) -> str:
    """Submit a FASTA sequence to the UCL PSIPRED web server.

    Args:
        fasta_path: Path to the FASTA file to submit.
        email: Email address for job notifications.

    Returns:
        UUID assigned to the submitted job.

    Raises:
        FileNotFoundError: If the FASTA file does not exist.
        RuntimeError: If the server response lacks a UUID.
    """
    logger.info("Submitting %s to server...", fasta_path.name)

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    fasta_content = fasta_path.read_text(encoding="utf-8")

    payload = {
        "job": "psipred",
        "submission_name": fasta_path.stem,
        "email": email,
    }

    files = {"input_data": (fasta_path.name, fasta_content)}

    response = requests.post(
        URL_SUBMIT, data=payload, files=files, headers=HEADERS, timeout=120
    )
    response.raise_for_status()
    result = response.json()

    uuid = result.get("UUID")
    if not uuid:
        raise RuntimeError(f"Server accepted request but returned no UUID: {result}")

    logger.info("Job submitted successfully. UUID: %s", uuid)
    return uuid


def poll_job(uuid: str, interval: int) -> dict[str, Any]:
    """Poll the PSIPRED server until the job completes or times out.

    Args:
        uuid: UUID of the submitted job.
        interval: Polling interval in seconds.

    Returns:
        Dictionary containing the final job status and results.

    Raises:
        RuntimeError: If the job fails on the server.
        TimeoutError: If polling exceeds ``MAX_POLL_SECONDS``.
    """
    elapsed_time = 0
    while elapsed_time < MAX_POLL_SECONDS:
        time.sleep(interval)
        elapsed_time += interval
        response = requests.get(f"{URL_SUBMIT}{uuid}", headers=HEADERS, timeout=60)
        response.raise_for_status()
        status_data = response.json()

        state = status_data.get("state", "Unknown")
        msg = status_data.get("last_message", "")

        if state.lower() == "complete":
            logger.info("Job completed successfully")
            return status_data
        if state.lower() == "error":
            raise RuntimeError(f"Job failed with message: {msg}")

    raise TimeoutError("Job polling timed out after 5 hours.")


def download_results(status_data: dict[str, Any], output_dir: Path) -> Path | None:
    """Download result files from the PSIPRED server.

    Args:
        status_data: Dictionary returned by ``poll_job``.
        output_dir: Local directory to save downloaded files.

    Returns:
        Path to the downloaded ``.ss2`` file, or *None* if no SS2 file
        was found.

    Raises:
        RuntimeError: If no submissions or results are present in the
            status data.
    """
    submissions = status_data.get("submissions", [])
    if not submissions:
        raise RuntimeError("No submissions found in status data")

    results_list = submissions[0].get("results", [])
    if not results_list:
        raise RuntimeError("No results found in status data")

    data_paths = [r["data_path"] for r in results_list if "data_path" in r]
    logger.info("[PSIPRED] Downloading %d files to %s", len(data_paths), output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded_ss2: Path | None = None

    for path_fragment in data_paths:
        filename = Path(path_fragment).name
        url = URL_DOWNLOAD_BASE + filename
        local_path = output_dir / filename

        r = requests.get(url, timeout=120)
        r.raise_for_status()
        local_path.write_bytes(r.content)
        logger.info("[PSIPRED] Downloaded %s to %s", filename, local_path)

        if filename.endswith(".ss2"):
            downloaded_ss2 = local_path

    return downloaded_ss2


def run_psipred_request() -> None:
    """Orchestrate the full PSIPRED workflow.

    Submits the FASTA sequence, polls until complete, downloads results,
    and copies the SS2 file to the pipeline input directory.

    Raises:
        RuntimeError: If the job finishes without producing an SS2 file.
    """
    fasta_path = Path(config.FASTA_FILE_PATH)
    results_dir = Path(config.PSIPRED_RESULTS_DIR)
    target_ss2_path = Path(config.INPUT_DIR) / config.SS2_FILE_BASENAME

    uuid = submit_job(fasta_path, config.PSIPRED_EMAIL)
    final_status = poll_job(uuid, config.PSIPRED_POLL_INTERVAL)
    ss2_path = download_results(final_status, results_dir)

    if ss2_path:
        shutil.copy2(ss2_path, target_ss2_path)
        logger.info("Copied %s to %s", ss2_path.name, target_ss2_path)
    else:
        raise RuntimeError("Job finished, but no SS2 file was found in results")
