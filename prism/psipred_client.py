#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM (Protected-Region Insertion Suite for Modeling)
#
# Citation:
# If you use this software in your research, please cite:
# Lopez-Corbalan, R.

"""
PSIPRED API Client.

This module handles the submission of a FASTA sequence to the
UCL PSIPRED web server, polls for results, and downloads
the final secondary structure (.ss2) file.
"""

import requests
import time
import os
import sys
import shutil

# Use relative imports, as this is now part of the 'prism' package
from .config import (
    FASTA_FILE_PATH,
    PSIPRED_EMAIL,
    PSIPRED_RESULTS_DIR,
    SS2_FILE_BASENAME,
    INPUT_DIR,
    PSIPRED_POLL_INTERVAL
)

def run_psipred_request():
    """
    Main function to submit, poll, and retrieve PSIPRED results.
    
    Reads all configuration from config.py.
    - Submits FASTA_FILE_PATH.
    - Downloads all results to PSIPRED_RESULTS_DIR.
    - Copies the final .ss2 file to INPUT_DIR with the name SS2_FILE_BASENAME.

    Raises:
        FileNotFoundError: If the FASTA file is missing.
        RuntimeError: If the server job fails or no .ss2 file is returned.
        IOError: If the final .ss2 file cannot be copied.
    """

    print(f"[PSIPRED] Submitting {FASTA_FILE_PATH} to PSIPRED server...")

    # 1️⃣ Read the FASTA file
    try:
        with open(FASTA_FILE_PATH, "r") as f:
            fasta_content = f.read()
    except FileNotFoundError:
        print(f"❌ Error: FASTA file not found at: {FASTA_FILE_PATH}")
        raise

    # 2️⃣ Define endpoints
    url_submit = "https://bioinf.cs.ucl.ac.uk/psipred/api/submission/"
    url_check_base = "https://bioinf.cs.ucl.ac.uk/psipred/api/submission/"

    # 3️⃣ Submission data
    submission_name = os.path.splitext(os.path.basename(FASTA_FILE_PATH))[0]
    data = {
        "job": "psipred",
        "submission_name": submission_name,
        "email": PSIPRED_EMAIL
    }
    files = {
        "input_data": (os.path.basename(FASTA_FILE_PATH), fasta_content)
    }
    headers = {"Accept": "application/json"}

    # 4️⃣ Submit request to server
    print("📤 Submitting job to PSIPRED server...")
    try:
        response = requests.post(url_submit, data=data, files=files, headers=headers)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"❌ Error contacting PSIPRED server: {e}")
        raise

    # 5️⃣ Parse response
    try:
        result = response.json()
    except Exception:
        print("❌ Server response is not valid JSON. Text returned:")
        print(response.text)
        raise

    uuid = result.get("UUID")
    if not uuid:
        raise RuntimeError(f"❌ Did not receive a valid UUID: {result}")

    print(f"✅ Job accepted (UUID = {uuid})")

    # 6️⃣ Poll for progress
    last_state = None
    status_data = {}
    while True:
        try:
            time.sleep(PSIPRED_POLL_INTERVAL)
            r = requests.get(f"{url_check_base}{uuid}", headers=headers)
            r.raise_for_status()
            status_data = r.json()

            state = status_data.get("state", "Unknown")
            msg = status_data.get("last_message", "")
            progress_msg = f"🔁 Status: {state}"
            if msg:
                progress_msg += f" | Message: {msg}"
            if state != last_state:
                print(progress_msg)
                last_state = state

            if state.lower() == "complete":
                print("✅ Job completed successfully.")
                break
            elif state.lower() == "error":
                raise RuntimeError(f"❌ Job failed: {status_data}")
        
        except requests.exceptions.RequestException as e:
            print(f"⚠️ Error during status check (will retry): {e}")
            time.sleep(PSIPRED_POLL_INTERVAL) # Extra wait before retry
        except KeyboardInterrupt:
            print("\nInterrupted by user.")
            raise

    # 7️⃣ Download ALL results to psipred-results/
    results_list = status_data.get("submissions", [])[0].get("results", [])
    data_paths = [result["data_path"] for result in results_list if "data_path" in result]
    if not data_paths:
        raise RuntimeError("⚠️ No result paths found in submission data.")

    print(f"📥 Downloading {len(data_paths)} result file(s) to: {PSIPRED_RESULTS_DIR}")
    os.makedirs(PSIPRED_RESULTS_DIR, exist_ok=True)
    
    downloaded_ss2_path = None

    for path in data_paths:
        filename_on_server = os.path.basename(path)
        url = "https://bioinf.cs.ucl.ac.uk/psipred/api/submissions/" + filename_on_server
        local_out_path = os.path.join(PSIPRED_RESULTS_DIR, filename_on_server)

        try:
            r_file = requests.get(url)
            r_file.raise_for_status()
            
            with open(local_out_path, "wb") as f_out:
                f_out.write(r_file.content)

            print(f"  > Saved: {local_out_path}")
            
            # Save the path to the .ss2 file if we find it
            if local_out_path.endswith(".ss2"):
                downloaded_ss2_path = local_out_path
                
        except Exception as e:
            print(f"  > ⚠️ Failed to download {filename_on_server}: {e}")

    print(f"🎉 All PSIPRED results downloaded to {PSIPRED_RESULTS_DIR}.")

    # 8️⃣ Copy the .ss2 file to the /input folder
    if downloaded_ss2_path:
        final_destination_path = os.path.join(INPUT_DIR, SS2_FILE_BASENAME)
        print(f"📥 Copying {downloaded_ss2_path} to {final_destination_path}...")
        
        try:
            # Ensure the /input directory exists
            os.makedirs(INPUT_DIR, exist_ok=True)
            shutil.copy2(downloaded_ss2_path, final_destination_path)
            print(f"✅ .ss2 file successfully copied to: {final_destination_path}")
        except Exception as e:
            raise IOError(f"❌ Error copying .ss2 file to {INPUT_DIR}: {e}")
    else:
        raise RuntimeError("❌ No .ss2 file was found in the PSIPRED download results.")

