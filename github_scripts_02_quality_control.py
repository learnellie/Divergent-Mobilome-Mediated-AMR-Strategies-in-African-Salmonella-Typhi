#!/usr/bin/env python3
"""
Step 2: Quality Control with QUAST
====================================
Run QUAST v5.2.0 on all assemblies, identify contaminated assemblies,
and apply QC filters.

Contaminated assemblies identified (>5.0 Mb total length):
  - 10060_6#72 (Tanzania, 5.20 Mb) — Peribacillus frigoritolerans kit contaminant
  - 16399_4#37 (Tanzania, 5.14 Mb) — S. Bareilly plasmid contaminant
  - 16404_6#46 (Ghana, 5.27 Mb) — S. Typhimurium plasmid contaminant

Usage:
  python 02_quality_control.py

Requirements:
  conda install -c bioconda quast
"""

import os
import sys
import subprocess
import pandas as pd
from pathlib import Path

# ---- Configuration ----
DATA_DIR = Path("data")
ASSEMBLIES_DIR = DATA_DIR / "assemblies"
METADATA_DIR = DATA_DIR / "metadata"
RESULTS_DIR = Path("results")
QUAST_DIR = RESULTS_DIR / "quast"

# Contaminated assemblies to exclude
CONTAMINATED = ["10060_6#72", "16399_4#37", "16404_6#46"]

# QC thresholds
MIN_N50 = 50000  # Below this = poor assembly (except ERR4451552 which is kept)
KEEP_DESPITE_LOW_N50 = ["ERR4451552"]  # Good contig count despite low N50


def run_quast():
    """Run QUAST on all assemblies."""
    print("=== Running QUAST v5.2.0 ===")

    QUAST_DIR.mkdir(parents=True, exist_ok=True)

    # Decompress assemblies for QUAST (it can't read .gz directly in some versions)
    decompressed_dir = RESULTS_DIR / "decompressed_tmp"
    decompressed_dir.mkdir(exist_ok=True)

    fasta_files = list(ASSEMBLIES_DIR.glob("*.fasta.gz")) + list(ASSEMBLIES_DIR.glob("*.fna.gz"))
    print(f"  Found {len(fasta_files)} compressed assemblies")

    # Decompress
    for f in fasta_files:
        out_name = f.name.replace(".gz", "")
        out_path = decompressed_dir / out_name
        if not out_path.exists():
            subprocess.run(f"gunzip -c {f} > {out_path}", shell=True, check=True)

    # Run QUAST
    input_files = list(decompressed_dir.glob("*.fasta")) + list(decompressed_dir.glob("*.fna"))
    input_str = " ".join(str(f) for f in input_files)

    cmd = f"""quast.py {input_str} \\
      -o {QUAST_DIR} \\
      --min-contig 0 \\
      --threads 4 \\
      --fast"""

    print(f"  Running QUAST on {len(input_files)} assemblies...")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  QUAST error: {result.stderr[-500:]}")
    else:
        print("  QUAST completed successfully")

    # Parse results
    report_path = QUAST_DIR / "report.tsv"
    if report_path.exists():
        quast_df = pd.read_csv(report_path, sep="\t", index_col=0)
        print(f"\n  QUAST report: {quast_df.shape}")
        print(f"  Median N50: {quast_df.loc['N50'].median():,.0f}")
        print(f"  Median contigs: {quast_df.loc['# contigs'].median():.0f}")

    return quast_df if report_path.exists() else None


def identify_contamination(quast_df):
    """Flag oversized assemblies (>5.0 Mb) as potentially contaminated."""
    print("\n=== Contamination Check ===")

    if quast_df is None:
        print("  No QUAST results to check")
        return

    total_lengths = quast_df.loc["Total length"]
    oversized = total_lengths[total_lengths.astype(float) > 5_000_000]

    print(f"  Assemblies >5.0 Mb (potential contamination):")
    for asm, length in oversized.items():
        status = "KNOWN CONTAMINANT" if asm in CONTAMINATED else "REVIEW NEEDED"
        print(f"    {asm}: {float(length)/1e6:.2f} Mb — {status}")

    return oversized.index.tolist()


def apply_qc_filters(quast_df):
    """Apply QC filters and update metadata."""
    print("\n=== Applying QC Filters ===")

    # Load metadata
    meta_path = METADATA_DIR / "study_subset_with_filenames.csv"
    if not meta_path.exists():
        print(f"  ERROR: Metadata not found at {meta_path}")
        return

    df = pd.read_csv(meta_path)
    print(f"  Starting with {len(df)} isolates")

    # Remove contaminated assemblies
    df["assembly_name"] = df["fasta_filename"].str.replace(r"\.fasta\.gz$|\.fna\.gz$", "", regex=True)
    contaminated_mask = df["assembly_name"].isin(CONTAMINATED)
    n_contam = contaminated_mask.sum()
    df = df[~contaminated_mask].copy()
    print(f"  Removed {n_contam} contaminated assemblies")

    # Remove poor N50 assemblies (except ERR4451552)
    if quast_df is not None:
        n50 = quast_df.loc["N50"]
        poor_n50 = [asm for asm in n50.index if float(n50[asm]) < MIN_N50 and asm not in KEEP_DESPITE_LOW_N50]
        poor_mask = df["assembly_name"].isin(poor_n50)
        n_poor = poor_mask.sum()
        df = df[~poor_mask].copy()
        print(f"  Removed {n_poor} poor-N50 assemblies (<{MIN_N50:,})")

    # Save filtered metadata
    df.to_csv(METADATA_DIR / "study_subset_qc.csv", index=False)
    print(f"\n  Final dataset: {len(df)} isolates")
    print(f"    East 4.3.1: {len(df[df['lineage']=='East_4.3.1'])}")
    print(f"    West 3.1.1: {len(df[df['lineage']=='West_3.1.1'])}")

    return df


def main():
    quast_df = run_quast()
    identify_contamination(quast_df)
    filtered_df = apply_qc_filters(quast_df)
    print("\n=== QUALITY CONTROL COMPLETE ===")


if __name__ == "__main__":
    main()
