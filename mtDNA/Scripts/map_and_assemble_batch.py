#!/usr/bin/env python3

import argparse
import subprocess
import shutil
import sys
from pathlib import Path
import csv
import shutil as sh


def check_dependencies():
    """Verify that required external programs are available on PATH."""
    required_tools = ["bbmap.sh", "spades.py"]
    missing = [tool for tool in required_tools if sh.which(tool) is None]

    if missing:
        print("\n[ERROR] Missing required tools in PATH:")
        for tool in missing:
            print(f"  - {tool}")
        print("\nPlease install or add them to your PATH.\n")
        sys.exit(1)

    print("=== Dependency check passed: bbmap.sh and spades.py found ===\n")


def run_cmd(cmd):
    print(f"\n[RUNNING] {' '.join(cmd)}\n")
    result = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if result.returncode != 0:
        print("[ERROR] Command failed:")
        print(result.stderr)
        sys.exit(result.returncode)
    print(result.stdout)


def process_sample(sample, reads1, reads2, reference, threads):
    """Run BBMap + SPAdes + rename contigs for a single sample."""

    print(f"\n==============================")
    print(f"=== Processing sample: {sample} ===")
    print(f"==============================\n")

    sample_dir = Path(sample)
    contigs_source = sample_dir / "contigs.fasta"
    contigs_target = f"{sample}_contigs.fasta"

    mapped1 = f"{sample}_mapped_1.fastq.gz"
    mapped2 = f"{sample}_mapped_2.fastq.gz"

    # ---- Step 1: BBMap ----
    run_cmd([
        "bbmap.sh",
        f"in1={reads1}",
        f"in2={reads2}",
        f"ref={reference}",
        f"outm1={mapped1}",
        f"outm2={mapped2}",
        f"threads={threads}"
    ])

    # ---- Step 2: SPAdes ----
    run_cmd([
        "spades.py",
        "--pe-1", "1", mapped1,
        "--pe-2", "1", mapped2,
        "-t", threads,
        "-o", sample
    ])

    # ---- Step 3: Copy contigs to sample-named output ----
    if not contigs_source.exists():
        print(f"[ERROR] Missing contigs.fasta for sample: {sample}")
        sys.exit(1)

    shutil.copy(contigs_source, contigs_target)

    print(f"\n✔ Sample {sample} complete!")
    print(f"→ Contigs saved as: {contigs_target}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Batch mapping + assembly pipeline using BBMap and SPAdes."
    )

    parser.add_argument("--csv", required=True,
                        help="CSV file with columns: sample,reads1,reads2,reference")
    parser.add_argument("--threads", default="8",
                        help="Threads to use per sample (default: 8)")

    args = parser.parse_args()

    # ---- Check dependencies ----
    check_dependencies()

    # ---- Read CSV ----
    with open(args.csv) as f:
        reader = csv.DictReader(f)
        required_cols = {"sample", "reads1", "reads2", "reference"}

        if not required_cols.issubset(reader.fieldnames):
            print(f"\n[ERROR] CSV must contain columns: {', '.join(required_cols)}\n")
            sys.exit(1)

        samples = list(reader)

    print(f"=== Loaded {len(samples)} samples from CSV ===\n")

    # ---- Process samples ----
    for row in samples:
        process_sample(
            sample=row["sample"],
            reads1=row["reads1"],
            reads2=row["reads2"],
            reference=row["reference"],
            threads=args.threads
        )

    print("\n==============================")
    print("=== ALL SAMPLES COMPLETED! ===")
    print("==============================\n")


if __name__ == "__main__":
    main()
