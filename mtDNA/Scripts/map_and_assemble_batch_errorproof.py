#!/usr/bin/env python3

import argparse
import subprocess
import os
import csv
from pathlib import Path
import shutil
import sys


def check_tool(tool):
    """Check if a tool exists on PATH."""
    result = shutil.which(tool)
    if result is None:
        print(f"[ERROR] Required tool '{tool}' not found in PATH.")
        return False
    return True


def run_cmd(cmd):
    """Run a shell command and return success/failure + output."""
    process = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    return process.returncode, process.stdout, process.stderr


def process_sample(reads1, reads2, reference, sample, threads, error_log):
    sample_dir = Path(sample)
    mapped1 = f"{sample}_mapped_1.fastq.gz"
    mapped2 = f"{sample}_mapped_2.fastq.gz"
    contigs_path = sample_dir / "contigs.fasta"
    final_contigs = f"{sample}_contigs.fasta"

    # ---- Skip if output already exists ----
    if contigs_path.exists() or Path(final_contigs).exists():
        msg = f"[SKIP] Sample '{sample}' already has output. Skipping.\n"
        print(msg)
        error_log.write(msg)
        return

    print(f"\n=== Processing sample: {sample} ===")

    # ---- Step 1: BBMap ----
    print("Running BBMap...")
    cmd = [
        "bbmap.sh",
        f"in1={reads1}",
        f"in2={reads2}",
        f"ref={reference}",
        f"outm1={mapped1}",
        f"outm2={mapped2}",
        f"threads={threads}"
    ]

    code, out, err = run_cmd(cmd)
    if code != 0:
        msg = f"[ERROR] BBMap failed for sample '{sample}': {err}\n"
        print(msg)
        error_log.write(msg)
        return

    # ---- Step 2: SPAdes ----
    print("Running SPAdes...")
    cmd = [
        "spades.py",
        "--pe-1", "1", mapped1,
        "--pe-2", "1", mapped2,
        "-t", threads,
        "-o", sample
    ]

    code, out, err = run_cmd(cmd)
    if code != 0:
        msg = f"[ERROR] SPAdes failed for sample '{sample}': {err}\n"
        print(msg)
        error_log.write(msg)
        return

    # ---- Step 3: Rename contigs ----
    if not contigs_path.exists():
        msg = f"[ERROR] SPAdes produced no contigs.fasta for sample '{sample}'\n"
        print(msg)
        error_log.write(msg)
        return

    shutil.copy(contigs_path, final_contigs)
    print(f"[OK] Finished sample '{sample}'. Output: {final_contigs}\n")


def main():
    parser = argparse.ArgumentParser(description="Batch BBMap + SPAdes pipeline.")
    parser.add_argument("--csv", required=True,
                        help="CSV file with columns: sample,reads1,reads2")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--threads", default="8")

    args = parser.parse_args()

    # ---- Check required tools ----
    if not check_tool("bbmap.sh") or not check_tool("spades.py"):
        sys.exit(1)

    # ---- Open error log (append mode) ----
    with open("error.log", "a") as error_log:
        error_log.write("\n=== NEW RUN ===\n")

        # ---- Read CSV ----
        with open(args.csv) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                sample = row["sample"]
                reads1 = row["reads1"]
                reads2 = row["reads2"]

                process_sample(
                    reads1, reads2, args.reference, sample, args.threads, error_log
                )

    print("\n=== ALL SAMPLES COMPLETE ===")
    print("Check error.log for any failed samples.")


if __name__ == "__main__":
    main()
