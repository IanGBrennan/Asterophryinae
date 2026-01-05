#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path
import shutil
import sys
import glob
from datetime import datetime


def check_dependency(tool):
    if shutil.which(tool) is None:
        print(f"[ERROR] Required tool '{tool}' not found in PATH.")
        sys.exit(1)


def run_cmd(cmd):
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    return proc.returncode, proc.stdout, proc.stderr


def assemble_sample(sample, r1, r2, threads, error_log):
    sample_dir = Path(sample)
    contigs_source = sample_dir / "contigs.fasta"
    contigs_target = f"{sample}_contigs.fasta"

    # ---- Skip if output already exists ----
    if contigs_source.exists() or Path(contigs_target).exists():
        msg = f"[SKIP] {sample}: assembly already exists\n"
        print(msg.strip())
        error_log.write(msg)
        return False

    print(f"\n=== Assembling sample: {sample} ===")

    cmd = [
        "spades.py",
        "--pe-1", "1", r1,
        "--pe-2", "1", r2,
        "-t", str(threads),
        "-o", sample
    ]

    code, out, err = run_cmd(cmd)
    if code != 0:
        msg = f"[ERROR] {sample}: SPAdes failed\n{err}\n"
        print(msg.strip())
        error_log.write(msg)
        return False

    if not contigs_source.exists():
        msg = f"[ERROR] {sample}: contigs.fasta not found after SPAdes\n"
        print(msg.strip())
        error_log.write(msg)
        return False

    shutil.copy(contigs_source, contigs_target)
    print(f"[OK] {sample}: contigs written to {contigs_target}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Batch SPAdes assembly from mapped read pairs"
    )
    parser.add_argument(
        "--input_dir", default=".",
        help="Directory containing *_mapped_1.fastq.gz files (default: current directory)"
    )
    parser.add_argument(
        "--threads", default=8, type=int,
        help="Threads per SPAdes run (default: 8)"
    )

    args = parser.parse_args()

    # ---- Dependency check ----
    check_dependency("spades.py")

    input_dir = Path(args.input_dir)

    r1_files = sorted(input_dir.glob("*_mapped_1.fastq.gz"))

    if not r1_files:
        print("[ERROR] No *_mapped_1.fastq.gz files found.")
        sys.exit(1)

    print(f"=== Found {len(r1_files)} samples to process ===")

    successes = []
    failures = []

    with open("spades_errors.log", "a") as error_log:
        error_log.write(f"\n=== NEW RUN {datetime.now()} ===\n")

        for r1 in r1_files:
            sample = r1.name.replace("_mapped_1.fastq.gz", "")
            r2 = input_dir / f"{sample}_mapped_2.fastq.gz"

            if not r2.exists():
                msg = f"[ERROR] {sample}: missing {r2.name}\n"
                print(msg.strip())
                error_log.write(msg)
                failures.append(sample)
                continue

            ok = assemble_sample(
                sample=sample,
                r1=str(r1),
                r2=str(r2),
                threads=args.threads,
                error_log=error_log
            )

            if ok:
                successes.append(sample)
            else:
                failures.append(sample)

    # ---- Summary ----
    print("\n==============================")
    print("=== SPADES BATCH COMPLETE ===")
    print("==============================\n")

    print(f"Successful assemblies ({len(successes)}):")
    for s in successes:
        print(f"  ✔ {s}")

    print(f"\nFailed assemblies ({len(failures)}):")
    for s in failures:
        print(f"  ✖ {s}")

    if failures:
        print("\nSee spades_errors.log for details.\n")


if __name__ == "__main__":
    main()
