#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys

def run_cmd(cmd):
    print(f"\n[RUNNING] {' '.join(cmd)}\n")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print("[ERROR] Command failed:")
        print(result.stderr)
        sys.exit(result.returncode)

    print(result.stdout)


def main():
    parser = argparse.ArgumentParser(description="Map reads with BBMap and assemble with SPAdes.")

    parser.add_argument("--reads1", required=True, help="Input FASTQ.gz file (read 1)")
    parser.add_argument("--reads2", required=True, help="Input FASTQ.gz file (read 2)")
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--sample", required=True, help="Sample name (used for output directory)")
    parser.add_argument("--threads", default="8", help="Number of threads (default: 8)")

    args = parser.parse_args()

    # Output mapped reads
    mapped1 = f"{args.sample}_mapped_1.fastq.gz"
    mapped2 = f"{args.sample}_mapped_2.fastq.gz"

    # ---- Step 1: BBMap ----
    bbmap_cmd = [
        "bbmap.sh",
        f"in1={args.reads1}",
        f"in2={args.reads2}",
        f"ref={args.reference}",
        f"outm1={mapped1}",
        f"outm2={mapped2}",
        f"threads={args.threads}"
    ]

    print("=== Running BBMap mapping step ===")
    run_cmd(bbmap_cmd)

    # ---- Step 2: SPAdes ----
    spades_cmd = [
        "spades.py",
        "--pe-1", "1", mapped1,
        "--pe-2", "1", mapped2,
        "-t", args.threads,
        "-o", args.sample
    ]

    print("=== Running SPAdes assembly step ===")
    run_cmd(spades_cmd)

    print("\n=== DONE ===")
    print(f"Mapped reads: {mapped1}, {mapped2}")
    print(f"Assembly output directory: {args.sample}")


if __name__ == "__main__":
    main()
