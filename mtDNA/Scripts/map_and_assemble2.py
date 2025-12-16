#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys
import shutil
from pathlib import Path
import shutil as sh


def check_dependencies():
    """Verify that required external programs are available on PATH."""
    required_tools = ["bbmap.sh", "spades.py"]

    missing = []

    for tool in required_tools:
        if sh.which(tool) is None:
            missing.append(tool)

    if missing:
        print("\n[ERROR] Missing required tools in your PATH:")
        for tool in missing:
            print(f"  - {tool}")
        print("\nPlease install them or add them to your PATH before running this script.\n")
        sys.exit(1)

    print("=== Dependency check passed: bbmap.sh and spades.py found ===")


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


def main():
    parser = argparse.ArgumentParser(
        description="Map reads with BBMap and assemble with SPAdes."
    )

    parser.add_argument("--reads1", required=True_
