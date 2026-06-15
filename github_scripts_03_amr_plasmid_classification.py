#!/usr/bin/env python3
"""
Step 3: AMR Gene Detection and Plasmid Classification
=======================================================
Run AMRFinderPlus and PlasmidFinder on all QC-passing assemblies,
then cross-reference to classify each AMR gene as plasmid-borne or chromosomal.

AMRFinderPlus parameters:
  --organism Salmonella --ident_min 0.9 --coverage_min 0.5

PlasmidFinder parameters:
  -d enterobacteriaceae -t 0.95 -l 0.60 -x

Cross-referencing logic:
  1. Build lookup: (assembly, contig) → plasmid replicon from PlasmidFinder
  2. For each AMR hit, check if its contig has a PlasmidFinder hit → plasmid-borne
  3. Refine: short contigs (<15 kb) with ≥2 MDR genes in IncHI1+ isolates → "likely_plasmid"

Usage:
  python 03_amr_plasmid_classification.py

Requirements:
  conda install -c bioconda amrfinderplus plasmidfinder
  amrfinder_update  (to download latest database)
"""

import os
import sys
import subprocess
import pandas as pd
import gzip
import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

# ---- Configuration ----
DATA_DIR = Path("data")
ASSEMBLIES_DIR = DATA_DIR / "assemblies"
METADATA_DIR = DATA_DIR / "metadata"
RESULTS_DIR = Path("results")
AMR_DIR = RESULTS_DIR / "amrfinderplus"
PLASMID_DIR = RESULTS_DIR / "plasmidfinder"

# AMRFinderPlus parameters
AMR_ORGANISM = "Salmonella"
AMR_IDENT_MIN = 0.9
AMR_COVERAGE_MIN = 0.5
AMR_DB = None  # Use default; set path if custom DB needed

# PlasmidFinder parameters
PLASMID_DB = "enterobacteriaceae"
PLASMID_IDENTITY = 0.95
PLASMID_COVERAGE = 0.60

# MDR cassette contig heuristic
MDR_GENES = {"blaTEM-1", "catA1", "sul1", "dfrA15", "dfrA7", "tet(B)", "qacEdelta1"}
MAX_CASSETTE_CONTIG_LEN = 15000  # bp
MIN_MDR_GENES_ON_CONTIG = 2

# Region classification
EAST_COUNTRIES = ["Kenya", "Malawi", "Tanzania", "Uganda", "Zimbabwe",
                  "Ethiopia", "Rwanda", "Mozambique", "South Africa"]


def get_assembly_list():
    """Load QC-filtered assembly list."""
    meta_path = METADATA_DIR / "study_subset_qc.csv"
    df = pd.read_csv(meta_path)
    df["assembly_name"] = df["fasta_filename"].str.replace(
        r"\.fasta\.gz$|\.fna\.gz$", "", regex=True
    )
    print(f"  Loaded {len(df)} QC-passing isolates")
    return df


def run_amrfinderplus(assembly_name, fasta_path, output_dir):
    """Run AMRFinderPlus on a single assembly."""
    out_file = output_dir / f"{assembly_name}.tsv"

    if out_file.exists():
        return out_file

    # Decompress if needed
    decompressed = None
    if str(fasta_path).endswith(".gz"):
        decompressed = Path(str(fasta_path).replace(".gz", ""))
        if not decompressed.exists():
            with gzip.open(fasta_path, "rt") as f_in:
                with open(decompressed, "w") as f_out:
                    shutil.copyfileobj(f_in, f_out)
        fasta_path = decompressed

    cmd = f"""amrfinder \\
      --nucleotide {fasta_path} \\
      --organism {AMR_ORGANISM} \\
      --ident_min {AMR_IDENT_MIN} \\
      --coverage_min {AMR_COVERAGE_MIN} \\
      --name {assembly_name} \\
      --output {out_file}"""

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  AMRFinderPlus ERROR for {assembly_name}: {result.stderr[:200]}")
        return None

    # Cleanup decompressed file
    if decompressed and decompressed.exists():
        decompressed.unlink()

    return out_file


def run_plasmidfinder(assembly_name, fasta_path, output_dir):
    """Run PlasmidFinder on a single assembly."""
    out_dir = output_dir / assembly_name

    if (out_dir / "results_tab.tsv").exists():
        return out_dir

    # Decompress if needed (PlasmidFinder cannot read gzipped FASTA)
    decompressed = None
    if str(fasta_path).endswith(".gz"):
        decompressed = Path(str(fasta_path).replace(".gz", ""))
        if not decompressed.exists():
            with gzip.open(fasta_path, "rt") as f_in:
                with open(decompressed, "w") as f_out:
                    shutil.copyfileobj(f_in, f_out)
        fasta_path = decompressed

    # PlasmidFinder database path (conda installation)
    db_path = None
    for candidate in [
        "/opt/conda/share/plasmidfinder-2.1.6/database",
        os.path.expanduser("~/.conda/envs/default/share/plasmidfinder-2.1.6/database"),
    ]:
        if os.path.exists(candidate):
            db_path = candidate
            break

    cmd_parts = [
        "plasmidfinder.py",
        f"-i {fasta_path}",
        f"-o {out_dir}",
        f"-d {PLASMID_DB}",
        f"-t {PLASMID_IDENTITY}",
        f"-l {PLASMID_COVERAGE}",
        "-x",  # Extended output (do NOT use -q, it suppresses results_tab.tsv)
    ]
    if db_path:
        cmd_parts.append(f"-p {db_path}")

    cmd = " ".join(cmd_parts)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  PlasmidFinder ERROR for {assembly_name}: {result.stderr[:200]}")
        return None

    # Cleanup decompressed file
    if decompressed and decompressed.exists():
        decompressed.unlink()

    return out_dir


def run_all_tools(metadata_df):
    """Run AMRFinderPlus and PlasmidFinder on all assemblies."""
    AMR_DIR.mkdir(parents=True, exist_ok=True)
    PLASMID_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\n=== Running AMRFinderPlus + PlasmidFinder on {len(metadata_df)} assemblies ===")

    amr_results = {}
    plasmid_results = {}

    for idx, row in metadata_df.iterrows():
        asm = row["assembly_name"]
        fn = row["fasta_filename"]
        fasta_path = ASSEMBLIES_DIR / fn

        if not fasta_path.exists():
            # Try without .gz
            fasta_path = ASSEMBLIES_DIR / fn.replace(".gz", "")

        if not fasta_path.exists():
            print(f"  SKIP {asm}: FASTA not found")
            continue

        if idx % 20 == 0:
            print(f"  Processing {idx+1}/{len(metadata_df)}: {asm}")

        # AMRFinderPlus
        amr_out = run_amrfinderplus(asm, fasta_path, AMR_DIR)
        if amr_out:
            amr_results[asm] = amr_out

        # PlasmidFinder
        plasmid_out = run_plasmidfinder(asm, fasta_path, PLASMID_DIR)
        if plasmid_out:
            plasmid_results[asm] = plasmid_out

    print(f"  AMRFinderPlus: {len(amr_results)} results")
    print(f"  PlasmidFinder: {len(plasmid_results)} results")
    return amr_results, plasmid_results


def parse_amrfinderplus_results(amr_results):
    """Parse all AMRFinderPlus output files into one DataFrame."""
    all_hits = []
    for asm, path in amr_results.items():
        try:
            df = pd.read_csv(path, sep="\t")
            all_hits.append(df)
        except Exception as e:
            print(f"  Error parsing {path}: {e}")

    if not all_hits:
        return pd.DataFrame()

    combined = pd.concat(all_hits, ignore_index=True)

    # Filter to acquired AMR only (Type == "AMR", not "STRESS")
    amr_only = combined[combined["Type"] == "AMR"].copy()
    print(f"  Total AMRFinderPlus hits: {len(combined)}")
    print(f"  Acquired AMR hits (Type=AMR): {len(amr_only)}")

    return amr_only


def parse_plasmidfinder_results(plasmid_results):
    """Parse all PlasmidFinder output files into one DataFrame."""
    all_hits = []
    for asm, out_dir in plasmid_results.items():
        results_tab = out_dir / "results_tab.tsv"
        if results_tab.exists():
            try:
                df = pd.read_csv(results_tab, sep="\t")
                if len(df) > 0 and not df.iloc[0].astype(str).eq("-").all():
                    df["Assembly"] = asm
                    all_hits.append(df)
            except Exception as e:
                print(f"  Error parsing {results_tab}: {e}")

    if not all_hits:
        return pd.DataFrame()

    combined = pd.concat(all_hits, ignore_index=True)
    print(f"  Total PlasmidFinder hits: {len(combined)}")
    return combined


def cross_reference(amr_df, plasmid_df, metadata_df):
    """Cross-reference AMR genes with plasmid replicons to classify location."""
    print("\n=== Cross-Referencing AMR ↔ PlasmidFinder ===")

    # Build plasmid contig lookup: (assembly, contig) → replicon info
    plasmid_contigs = {}
    if len(plasmid_df) > 0:
        for _, row in plasmid_df.iterrows():
            asm = row.get("Assembly", "")
            contig = row.get("Contig", "")
            repl_type = row.get("Plasmid", "")
            if asm and contig:
                plasmid_contigs[(asm, contig)] = repl_type

    print(f"  Plasmid contig entries: {len(plasmid_contigs)}")

    # Classify each AMR hit
    classifications = []
    for _, hit in amr_df.iterrows():
        asm = hit.get("Name", hit.get("Assembly", ""))
        contig = hit.get("Contig id", "")
        gene = hit.get("Element symbol", "")
        gene_name = hit.get("Element name", "")
        hit_type = hit.get("Type", "")
        amr_class = hit.get("Class", "")
        subclass = hit.get("Subclass", "")
        start = hit.get("Start", "")
        stop = hit.get("Stop", "")
        strand = hit.get("Strand", "")

        # Primary classification: is the contig in PlasmidFinder results?
        if (asm, contig) in plasmid_contigs:
            location = "plasmid-borne"
            replicon = plasmid_contigs[(asm, contig)]
        else:
            location = "chromosomal"
            replicon = ""

        classifications.append({
            "Assembly": asm,
            "Contig": contig,
            "Gene": gene,
            "Gene_name": gene_name,
            "Type": hit_type,
            "AMR_Class": amr_class,
            "Subclass": subclass,
            "Start": start,
            "Stop": stop,
            "Strand": strand,
            "Location": location,
            "Plasmid_replicon": replicon,
        })

    class_df = pd.DataFrame(classifications)

    # Add metadata (Country, Region, Genotype)
    meta_lookup = metadata_df.set_index("assembly_name")
    class_df["Country"] = class_df["Assembly"].map(
        lambda a: meta_lookup.loc[a, "Country"] if a in meta_lookup.index else "Unknown"
    )
    class_df["Region"] = class_df["Country"].apply(
        lambda c: "East_4.3.1" if c in EAST_COUNTRIES else "West_3.1.1"
    )
    class_df["Genotype"] = class_df["Assembly"].map(
        lambda a: meta_lookup.loc[a, "Final_genotype"] if a in meta_lookup.index else "Unknown"
    )

    # ---- Refinement: "likely_plasmid" for MDR cassette contigs ----
    # Short contigs (<15 kb) with ≥2 MDR genes in IncHI1-positive isolates
    # are likely fragments of the IncHI1 MDR cassette that split across contigs
    # during short-read assembly

    # Identify IncHI1-positive isolates
    inchi1_isolates = set()
    if len(plasmid_df) > 0:
        for _, row in plasmid_df.iterrows():
            repl = str(row.get("Plasmid", ""))
            if "IncHI1" in repl:
                inchi1_isolates.add(row.get("Assembly", ""))

    # Find MDR cassette contigs
    contig_mdr_counts = {}
    for _, row in class_df.iterrows():
        asm = row["Assembly"]
        contig = row["Contig"]
        gene = row["Gene"]
        if gene in MDR_GENES and asm in inchi1_isolates:
            key = (asm, contig)
            if key not in contig_mdr_counts:
                contig_mdr_counts[key] = {"genes": set(), "contig": contig}
            contig_mdr_counts[key]["genes"].add(gene)

    # Get contig lengths from assembly files (approximate from AMR hit positions)
    # For simplicity, use the Stop position as a proxy for contig length
    # (more accurate: read from FASTA headers)
    contig_lengths = {}
    for _, row in class_df.iterrows():
        key = (row["Assembly"], row["Contig"])
        stop = row["Stop"]
        if key not in contig_lengths or (pd.notna(stop) and int(stop) > contig_lengths[key]):
            contig_lengths[key] = int(stop) if pd.notna(stop) else 0

    # Mark likely_plasmid
    mdr_cassette_contigs = set()
    for (asm, contig), info in contig_mdr_counts.items():
        if len(info["genes"]) >= MIN_MDR_GENES_ON_CONTIG:
            contig_len = contig_lengths.get((asm, contig), 0)
            if contig_len < MAX_CASSETTE_CONTIG_LEN:
                mdr_cassette_contigs.add((asm, contig))

    print(f"  MDR cassette contigs identified: {len(mdr_cassette_contigs)}")

    # Apply refinement
    class_df["Location_refined"] = class_df["Location"]
    for idx, row in class_df.iterrows():
        key = (row["Assembly"], row["Contig"])
        if key in mdr_cassette_contigs and row["Location"] == "chromosomal":
            class_df.at[idx, "Location_refined"] = "likely_plasmid"

    # Summary
    loc_counts = class_df["Location_refined"].value_counts()
    print(f"\n  Refined location classification:")
    for loc, count in loc_counts.items():
        print(f"    {loc}: {count}")

    # Save
    class_df.to_csv(RESULTS_DIR / "amr_plasmid_classification.tsv", sep="\t", index=False)
    print(f"\n  Saved: {RESULTS_DIR / 'amr_plasmid_classification.tsv'}")

    return class_df


def generate_summary_tables(class_df, metadata_df):
    """Generate regional AMR gene prevalence and plasmid replicon tables."""
    print("\n=== Generating Summary Tables ===")

    key_genes = ["blaTEM-1", "catA1", "sul1", "dfrA7", "dfrA15", "tet(B)", "sul2"]

    # AMR gene prevalence by region
    prevalence_data = []
    for gene in key_genes:
        for region in ["East_4.3.1", "West_3.1.1"]:
            region_df = class_df[class_df["Region"] == region]
            total_isolates = metadata_df[
                metadata_df["Country"].apply(
                    lambda c: ("East_4.3.1" if c in EAST_COUNTRIES else "West_3.1.1")
                ) == region
            ].shape[0]

            gene_hits = region_df[region_df["Gene"] == gene]
            n_isolates_with_gene = gene_hits["Assembly"].nunique()
            pct = n_isolates_with_gene / total_isolates * 100 if total_isolates > 0 else 0

            # Location breakdown
            loc_breakdown = gene_hits["Location_refined"].value_counts().to_dict()
            n_plasmid = loc_breakdown.get("plasmid-borne", 0)
            n_likely = loc_breakdown.get("likely_plasmid", 0)
            n_chrom = loc_breakdown.get("chromosomal", 0)

            prevalence_data.append({
                "Gene": gene,
                "Region": region,
                "N_isolates_with_gene": n_isolates_with_gene,
                "Total_isolates": total_isolates,
                "Prevalence_pct": round(pct, 1),
                "Plasmid_borne": n_plasmid,
                "Likely_plasmid": n_likely,
                "Chromosomal": n_chrom,
            })

    prev_df = pd.DataFrame(prevalence_data)
    prev_df.to_csv(RESULTS_DIR / "amr_gene_prevalence_by_region.tsv", sep="\t", index=False)
    print(f"  Saved: amr_gene_prevalence_by_region.tsv")

    # Plasmid replicon prevalence by region
    plasmid_df_parsed = parse_plasmidfinder_results(
        {asm: PLASMID_DIR / asm for asm in metadata_df["assembly_name"]
         if (PLASMID_DIR / asm / "results_tab.tsv").exists()}
    )
    if len(plasmid_df_parsed) > 0:
        repl_data = []
        for repl in ["IncHI1", "IncQ1"]:
            for region in ["East_4.3.1", "West_3.1.1"]:
                total = metadata_df[
                    metadata_df["Country"].apply(
                        lambda c: ("East_4.3.1" if c in EAST_COUNTRIES else "West_3.1.1")
                    ) == region
                ].shape[0]
                has_repl = plasmid_df_parsed[
                    plasmid_df_parsed["Plasmid"].str.contains(repl, na=False)
                ]["Assembly"].nunique()
                repl_data.append({
                    "Replicon": repl,
                    "Region": region,
                    "N_isolates": has_repl,
                    "Total_isolates": total,
                    "Prevalence_pct": round(has_repl / total * 100, 1) if total > 0 else 0,
                })

        repl_df = pd.DataFrame(repl_data)
        repl_df.to_csv(RESULTS_DIR / "plasmid_replicon_prevalence_by_region.tsv", sep="\t", index=False)
        print(f"  Saved: plasmid_replicon_prevalence_by_region.tsv")


def main():
    metadata_df = get_assembly_list()
    amr_results, plasmid_results = run_all_tools(metadata_df)
    amr_df = parse_amrfinderplus_results(amr_results)
    plasmid_df = parse_plasmidfinder_results(plasmid_results)
    class_df = cross_reference(amr_df, plasmid_df, metadata_df)
    generate_summary_tables(class_df, metadata_df)
    print("\n=== AMR + PLASMID CLASSIFICATION COMPLETE ===")


if __name__ == "__main__":
    main()
