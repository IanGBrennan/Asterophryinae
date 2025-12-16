#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys
from pathlib import Path
from Bio import SeqIO

def run_cmd(cmd):
    print(f"\n[RUNNING] {' '.join(cmd)}\n")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print("[ERROR] Command failed:")
        print(result.stderr)
        sys.exit(result.returncode)
    return result.stdout


def blast_and_rename_contigs(contig_file, reference_fasta, sample_name, threads):
    renamed_output = f"{sample_name}_contigs_renamed.fasta"

    # Make BLAST database
    print("=== Creating BLAST database ===")
    run_cmd(["makeblastdb", "-in", reference_fasta, "-dbtype", "nucl"])

    # Run BLASTN
    print("=== Running BLASTN on SPAdes contigs ===")
    blast_output = f"{sample_name}_blast_results.tsv"
    run_cmd([
        "blastn",
        "-query", contig_file,
        "-db", reference_fasta,
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
        "-num_threads", threads,
        "-out", blast_output
    ])

    # Load BLAST hits
    best_hits = {}
    with open(blast_output) as f:
        for line in f:
            qseqid, sseqid, pident, length, evalue, bitscore = line.strip().split("\t")
            bitscore = float(bitscore)

            # Save best hit per contig
            if qseqid not in best_hits or bitscore > best_hits[qseqid][1]:
                best_hits[qseqid] = (sseqid, bitscore)

    print(f"=== Renaming contigs based on best BLAST hits ===")
    renamed_records = []
    for record in SeqIO.parse(contig_file, "fasta"):
        if record.id in best_hits:
            hit, score = best_hits[record.id]
            record.id = f"{hit}__BLASTscore_{int(score)}"
        else:
            record.id = f"{record.id}__NO_BLAST_HIT"
        record.description = ""
        renamed_records.append(record)

    # Write renamed FASTA
    SeqIO.write(renamed_records, renamed_output, "fasta")
    print(f"Renamed contigs written to: {renamed_output}")

    return renamed_output


def main():
    parser = argparse.ArgumentParser(description="Map reads, assemble with SPAdes, and rename contigs via BLAST.")

    parser.add_argument("--reads1", required=True, help="Input FASTQ.gz file (read 1)")
    parser.add_argument("--reads2", required=True, help="Input FASTQ.gz file (read 2)")
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--sample", required=True, help="Sample name (used for output directory)")
    parser.add_argument("--threads", default="8", help="Number of threads (default: 8)")

    args = parser.parse_args()

    sample = args.sample
    contig_file = os.path.join(sample, "contigs.fasta")

    # Output mapped reads
    mapped1 = f"{sample}_mapped_1.fastq.gz"
    mapped2 = f"{sample}_mapped_2.fastq.gz"

    # ---- Step 1: BBMap ----
    print("=== Running BBMap mapping step ===")
    run_cmd([
        "bbmap.sh",
        f"in1={args.reads1}",
        f"in2={args.reads2}",
        f"ref={args.reference}",
        f"outm1={mapped1}",
        f"outm2={mapped2}",
        f"threads={args.threads}"
    ])

    # ---- Step 2: SPAdes ----
    print("=== Running SPAdes assembly step ===")
    run_cmd([
        "spades.py",
        "--pe-1", "1", mapped1,
        "--pe-2", "1", mapped2,
        "-t", args.threads,
        "-o", sample
    ])

    # Check contigs output
    if not Path(contig_file).exists():
        print(f"[ERROR] SPAdes did not produce contigs.fasta at: {contig_file}")
        sys.exit(1)

    # ---- Step 3: BLAST + rename contigs ----
    renamed = blast_and_rename_contigs(contig_file, args.reference, sample, args.threads)

    print("\n=== DONE ===")
    print(f"Mapped reads: {mapped1}, {mapped2}")
    print(f"Assembly directory: {sample}")
    print(f"Renamed contigs: {renamed}")


if __name__ == "__main__":
    main()
