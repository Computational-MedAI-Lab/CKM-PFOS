#!/usr/bin/env python3
"""
multiple-docking.py

Parallel AutoDock Vina docking with resume support.

Usage examples:
    python multiple-docking.py --workers 8 --timeout 300 --retries 2
    python multiple-docking.py --workers 16 --resume-file resume.json

Notes:
 - Make sure 'vina' is in PATH (or give full path via --vina-bin).
 - This script only checks for existing output files to decide resume.
"""

import os
import json
import csv
import time
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count

# -----------------------------
# Worker function (runs in processes)
# -----------------------------
def dock_single(ligand_path,
                receptor_path,
                out_path,
                log_path,
                center,
                grid_size,
                exhaustiveness,
                vina_bin,
                timeout,
                retries):
    """
    Runs vina for a single ligand, with retries and timeout.
    Returns a dict with status and metadata.
    """
    ligand_name = os.path.basename(ligand_path)
    attempt = 0
    last_returncode = None
    start_time = time.time()

    # Quick pre-check (file size)
    try:
        if not os.path.exists(ligand_path) or os.path.getsize(ligand_path) < 200:
            return {
                "ligand": ligand_name,
                "status": "bad_file",
                "msg": "ligand missing or too small",
                "time": 0.0
            }
    except Exception as e:
        return {
            "ligand": ligand_name,
            "status": "error",
            "msg": f"stat error: {e}",
            "time": 0.0
        }

    cmd_base = [
        vina_bin,
        "--receptor", receptor_path,
        "--ligand", ligand_path,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(grid_size),
        "--size_y", str(grid_size),
        "--size_z", str(grid_size),
        "--out", out_path,
        "--exhaustiveness", str(exhaustiveness)
    ]

    while attempt <= retries:
        attempt += 1
        try:
            # Open log file per attempt (overwrite previous attempt log for clarity)
            with open(log_path, "w") as lf:
                proc = subprocess.run(
                    cmd_base,
                    stdout=lf,
                    stderr=lf,
                    timeout=timeout
                )
                last_returncode = proc.returncode

            # If Vina returned success, confirm output exists and is reasonable
            if last_returncode == 0 and os.path.exists(out_path) and os.path.getsize(out_path) > 200:
                elapsed = time.time() - start_time
                return {
                    "ligand": ligand_name,
                    "status": "success",
                    "attempts": attempt,
                    "time": round(elapsed, 2),
                    "msg": "",
                    "returncode": 0
                }
            else:
                # Non-zero return or broken output - will retry if attempts remain
                msg = f"vina returned {last_returncode}"
                if attempt > retries:
                    elapsed = time.time() - start_time
                    return {
                        "ligand": ligand_name,
                        "status": "failed",
                        "attempts": attempt,
                        "time": round(elapsed, 2),
                        "msg": msg,
                        "returncode": last_returncode
                    }
                # small backoff before retry
                time.sleep(1.0)

        except subprocess.TimeoutExpired:
            msg = f"timeout after {timeout}s"
            if attempt > retries:
                elapsed = time.time() - start_time
                return {
                    "ligand": ligand_name,
                    "status": "timeout",
                    "attempts": attempt,
                    "time": round(elapsed, 2),
                    "msg": msg,
                    "returncode": None
                }
            # else retry
            time.sleep(1.0)
        except Exception as e:
            # Unexpected exception in worker
            elapsed = time.time() - start_time
            return {
                "ligand": ligand_name,
                "status": "error",
                "attempts": attempt,
                "time": round(elapsed, 2),
                "msg": f"exception: {e}",
                "returncode": None
            }

    # Fallback (shouldn't reach)
    return {
        "ligand": ligand_name,
        "status": "error",
        "attempts": attempt,
        "time": 0.0,
        "msg": "unknown error",
        "returncode": last_returncode
    }


# -----------------------------
# Main
# -----------------------------
def main():
    p = argparse.ArgumentParser(description="Parallel Vina docking with resume support")
    p.add_argument("--receptor", default="receptor/3EJH.pdbqt", help="receptor pdbqt path")
    p.add_argument("--ligand-folder", default="pdbqt_ligand_harbdatabase", help="folder with ligand .pdbqt files")
    p.add_argument("--output-folder", default="COL1A1_3EJH_output", help="folder to write outputs")
    p.add_argument("--center", nargs=3, type=float, default=[-25.22,24.11,-3.85], help="center_x center_y center_z")
    p.add_argument("--grid-size", type=float, default=20.0, help="grid box size (x=y=z)")
    p.add_argument("--exhaustiveness", type=int, default=16, help="vina exhaustiveness")
    p.add_argument("--workers", type=int, default=4, help="number of parallel workers (max 64)")
    p.add_argument("--timeout", type=int, default=300, help="per-ligand timeout (seconds)")
    p.add_argument("--retries", type=int, default=1, help="number of retries on failure/timeout")
    p.add_argument("--vina-bin", default="vina", help="path to vina binary")
    p.add_argument("--resume-file", default="resume_3EJH.json", help="resume progress JSON file")
    p.add_argument("--summary-csv", default="summary.csv", help="summary CSV output")
    p.add_argument("--max-workers-allowed", type=int, default=64, help="hard cap on workers")
    args = p.parse_args()

    # sanitize workers
    max_allowed = min(args.max_workers_allowed, 64)
    workers = max(1, min(args.workers, max_allowed, cpu_count() * 2))
    # make folders
    os.makedirs(args.output_folder, exist_ok=True)

    # load ligands
    ligands = sorted([f for f in os.listdir(args.ligand_folder) if f.endswith(".pdbqt")])
    if not ligands:
        print("No ligands found in", args.ligand_folder)
        return

    # load resume state if exists
    resume_state = {"done": {}, "failed": {}}
    if os.path.exists(args.resume_file):
        try:
            with open(args.resume_file, "r") as rf:
                resume_state = json.load(rf)
        except Exception:
            print("Warning: could not read resume file, starting fresh.")

    # build task list skipping already-done (basic check: output exists and non-trivial size)
    tasks = []
    for lig in ligands:
        out_name = f"{os.path.splitext(lig)[0]}_out.pdbqt"
        out_path = os.path.join(args.output_folder, f"{os.path.splitext(lig)[0]}_out.pdbqt")
        # consider done if recorded in resume_state done OR output file exists and >=200 bytes
        already_done = False
        if lig in resume_state.get("done", {}):
            already_done = True
        elif os.path.exists(out_path) and os.path.getsize(out_path) > 200:
            already_done = True

        if already_done:
            # ensure entry in resume_state
            resume_state["done"].setdefault(lig, {"status": "skipped_existing"})
            continue

        tasks.append(lig)

    total_to_run = len(tasks)
    print(f"Workers: {workers} | Ligands total: {len(ligands)} | To run: {total_to_run}")

    # prepare summary CSV header if new
    summary_csv_path = os.path.join(args.output_folder, args.summary_csv)
    write_header = not os.path.exists(summary_csv_path)
    if write_header:
        with open(summary_csv_path, "w", newline="") as csvf:
            writer = csv.DictWriter(csvf, fieldnames=["ligand", "status", "attempts", "time", "msg", "returncode"])
            writer.writeheader()

    # main executor
    with ProcessPoolExecutor(max_workers=workers) as exe:
        future_to_lig = {}
        for lig in tasks:
            ligand_path = os.path.join(args.ligand_folder, lig)
            out_basename = f"{os.path.splitext(lig)[0]}_out.pdbqt"
            out_path = os.path.join(args.output_folder, out_basename)
            log_path = os.path.join(args.output_folder, f"{os.path.splitext(lig)[0]}.log")

            fut = exe.submit(
                dock_single,
                ligand_path,
                args.receptor,
                out_path,
                log_path,
                args.center,
                args.grid_size,
                args.exhaustiveness,
                args.vina_bin,
                args.timeout,
                args.retries
            )
            future_to_lig[fut] = lig

        # collect results as completed
        for fut in as_completed(future_to_lig):
            lig = future_to_lig[fut]
            try:
                res = fut.result()
            except Exception as e:
                res = {
                    "ligand": lig,
                    "status": "error",
                    "attempts": 0,
                    "time": 0.0,
                    "msg": f"worker exception: {e}",
                    "returncode": None
                }

            # record into resume_state and summary CSV
            ligand_name = res.get("ligand", lig)
            status = res.get("status", "error")
            # write resume_state
            if status == "success":
                resume_state.setdefault("done", {})[ligand_name] = {"status": "success", "attempts": res.get("attempts", 1)}
                # ensure failed cleanup if existed
                resume_state.get("failed", {}).pop(ligand_name, None)
            else:
                resume_state.setdefault("failed", {})[ligand_name] = {"status": status, "msg": res.get("msg", "")}

            # append to summary CSV
            with open(summary_csv_path, "a", newline="") as csvf:
                writer = csv.DictWriter(csvf, fieldnames=["ligand", "status", "attempts", "time", "msg", "returncode"])
                writer.writerow({
                    "ligand": ligand_name,
                    "status": status,
                    "attempts": res.get("attempts", 0),
                    "time": res.get("time", 0.0),
                    "msg": res.get("msg", ""),
                    "returncode": res.get("returncode", "")
                })

            # persist resume file after every completed job (safe small writes)
            try:
                with open(args.resume_file, "w") as rf:
                    json.dump(resume_state, rf, indent=2)
            except Exception as e:
                print("Warning: cannot write resume file:", e)

            print(f"[{status.upper()}] {ligand_name} attempts={res.get('attempts', 0)} time={res.get('time', 0.0)} msg={res.get('msg','')}")

    # finished
    print("\nAll submitted jobs finished.")
    print("Summary written to:", summary_csv_path)
    print("Resume state written to:", args.resume_file)
    done_count = len(resume_state.get("done", {}))
    failed_count = len(resume_state.get("failed", {}))
    print(f"Done: {done_count}  Failed: {failed_count}")

if __name__ == "__main__":
    main()
