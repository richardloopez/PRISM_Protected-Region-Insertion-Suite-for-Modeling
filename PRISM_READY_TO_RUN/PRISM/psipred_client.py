#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM (psipred client)

'''
PSIPRED API Client.

Handles submission of FASTA sequences to the UCL PSIPRED web server,
polls for job completion, and retrieves secondary structure (.ss2) files.
'''

import time 
import shutil
import requests
from pathlib import Path 
from typing import Dict, Any, Optional

from . import config

URL_SUBMIT = "https://bioinf.cs.ucl.ac.uk/psipred/api/submission/"
URL_CHECK_BASE = "https://bioinf.cs.ucl.ac.uk/psipred/api/submission/"
URL_DOWNLOAD_BASE = "https://bioinf.cs.ucl.ac.uk/psipred/api/submissions/"
HEADERS = {"Accept": "application/json"}

def submit_job(fasta_path: Path, email: str) -> str:
    '''
    Submit a job to the UCL PSIPRED web server.

    Args:
        fasta_path: Path to the FASTA file to submit.
        email: Email address to receive job notifications.

    Returns:
        UUID of the submitted job.
    '''
    print(f"\n[PSIPRED] Submitting {fasta_path.name} to server...")

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    try:
        fasta_content = fasta_path.read_text()
    except Exception as e:
        raise IOError(f"Failed to read FASTA file: {e}")
    
    payload = {
        'job': 'psipred',
        'submission_name': fasta_path.stem,
        'email': email
    }

    files = {
        'input_data': (fasta_path.name, fasta_content)
    }

    try:
        response = requests.post(URL_SUBMIT, data=payload, files=files, headers=HEADERS)
        response.raise_for_status()
        result = response.json()
    except requests.RequestException as e:
        raise RuntimeError(f"Failed to contact PSIPRED server: {e}")
    except ValueError:
        raise RuntimeError(f"Invalid JSON response from server: {response.text}")
    
    uuid = result.get('UUID')
    if not uuid:
        raise RuntimeError(f"Server accepted request but returned no UUID: {result}")
    
    print(f"\n[PSIPRED] Job submitted successfully. UUID: {uuid}")
    return uuid


def poll_job(uuid: str, interval: int) -> Dict[str, Any]:
    '''
    Poll the UCL PSIPRED web server for job completion.

    Args:
        uuid: UUID of the submitted job.
        interval: Polling interval in seconds.

    Returns:
        Dictionary containing job status and results.
    '''
    last_state = None
    while True:
        try:
            time.sleep(interval)
            response = requests.get(f'{URL_CHECK_BASE}{uuid}', headers=HEADERS)
            response.raise_for_status()
            status_data = response.json()
            
            state = status_data.get('state', 'Unknown')
            msg = status_data.get('last_message', '')

            if state != last_state:
                log_msg = f"[PSIPRED] Job {state}: {msg}"
                print(log_msg)
                last_state = state
            if state.lower() == 'complete':
                print("\n[PSIPRED] Job completed successfully")
                return status_data
            if state.lower() == 'error':
                raise RuntimeError(f"Job failed with message: {msg}")
            
        except requests.RequestException as e: 
            print(f"[PSIPRED] Network error checking status (retrying): {e}")
            time.sleep(interval)


def download_results(status_data: Dict[str, Any], output_dir: Path) -> Optional[Path]:
    '''
    Download results from the UCL PSIPRED web server.

    Args:
        status_data: Dictionary containing job status and results.
        output_dir: Directory to save downloaded files.

    Returns:
        Optional path to downloaded .ss2 file.
    '''
    submissions = status_data.get('submissions', [])
    if not submissions:
        raise RuntimeError("No submissions found in status data")
    
    results_list = submissions[0].get('results', [])
    if not results_list:
        raise RuntimeError("No results found in status data")
    
    data_paths = [r['data_path'] for r in results_list if 'data_path' in r]
    print(f"\n[PSIPRED] Downloading {len(data_paths)} files to {output_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded_ss2 = None

    for path_fragment in data_paths:
        filename = Path(path_fragment).name
        url = URL_DOWNLOAD_BASE + filename
        local_path = output_dir / filename

        try:
            r = requests.get(url)
            r.raise_for_status()
            local_path.write_bytes(r.content)
            print(f"[PSIPRED] Downloaded {filename} to {local_path}")
            
            if filename.endswith('.ss2'):
                downloaded_ss2 = local_path

        except Exception as e:
            print(f"[PSIPRED] Failed to download {filename}: {e}")
            
    return downloaded_ss2


def run_psipred_request():
    """
    Orchestrates the PSIPRED pipeline: Submit job > Poll for completion > Download results > Copy SS2.
    """
    fasta_path = Path(config.FASTA_FILE_PATH)
    results_dir = Path(config.PSIPRED_RESULTS_DIR)
    target_ss2_path = Path(config.INPUT_DIR) / config.SS2_FILE_BASENAME

    uuid = submit_job(fasta_path, config.PSIPRED_EMAIL)
    final_status = poll_job(uuid, config.PSIPRED_POLL_INTERVAL)
    ss2_path = download_results(final_status, results_dir)

    if ss2_path:
        try:
            target_dir = Path(config.INPUT_DIR)
            target_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ss2_path, target_ss2_path)
            print(f"[PSIPRED] Copied {ss2_path.name} to {target_ss2_path}")
        except Exception as e:
            raise IOError(f"[PSIPRED] Failed to copy {ss2_path.name} to {target_ss2_path}: {e}")
    else:
        raise RuntimeError("Job finished, but no SS2 file was found in results")







