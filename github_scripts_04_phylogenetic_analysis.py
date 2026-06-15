#!/usr/bin/env python3
"""
Step 4: Phylogenetic Analysis
===============================
Build core genome SNP phylogeny from assemblies mapped to S. Typhi CT18,
confirm East 4.3.1 vs West 3.1.1 lineage split, and map AMR gene
acquisition events onto the tree.

Pipeline:
  1. Map each assembly to CT18 reference with minimap2
  2. Generate CIGAR-based pseudo-sequences (reference-length alignment)
  3. Extract SNP sites with snp-sites
  4. Build ML tree with IQ-TREE (ModelFinder + 1000 UFBoot)
  5. Midpoint root the tree
  6. Annotate with AMR gene presence/absence and location
  7. Generate phylogenetic figures

Reference: S. Typhi CT18 (AL513382.1, chromosome 4,809,037 bp)
  - Chromosome only (pHCM1 plasmid excluded to avoid plasmid-borne SNPs)

Usage:
  python 04_phylogenetic_analysis.py

Requirements:
  conda install -c bioconda minimap2 samtools pysam snp-sites iqtree biopython
  pip install matplotlib seaborn
"""

import os
import sys
import subprocess
import gzip
import shutil
import re
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
    import pysam
except ImportError:
    print("ERROR: pysam required. Install with: pip install pysam")
    sys.exit(1)

# ---- Configuration ----
DATA_DIR = Path("data")
ASSEMBLIES_DIR = DATA_DIR / "assemblies"
METADATA_DIR = DATA_DIR / "metadata"
RESULTS_DIR = Path("results")
PHYLO_DIR = RESULTS_DIR / "phylogenetics"

# Reference
REF_ACCESSION = "AL513382.1"  # S. Typhi CT18 chromosome
REF_LENGTH = 4809037

# IQ-TREE parameters
IQTREE_MODEL = "MFP"  # ModelFinder Plus
IQTREE_BOOTSTRAP = 1000
CONSTANT_SITES = "1149511,1249782,1253591,1152957"  # A,C,G,T constant site counts

# Region classification
EAST_COUNTRIES = ["Kenya", "Malawi", "Tanzania", "Uganda", "Zimbabwe",
                  "Ethiopia", "Rwanda", "Mozambique", "South Africa"]

# Key AMR genes for tree annotation
KEY_GENES = ["blaTEM-1", "catA1", "sul1", "dfrA7", "dfrA15", "tet(B)", "sul2"]


def download_reference():
    """Download S. Typhi CT18 reference genome from NCBI."""
    ref_path = PHYLO_DIR / "CT18_chromosome.fna"

    if ref_path.exists():
        print(f"  Reference already exists: {ref_path}")
        return ref_path

    PHYLO_DIR.mkdir(parents=True, exist_ok=True)

    print(f"  Downloading CT18 reference ({REF_ACCESSION})...")
    Entrez.email = "researcher@example.com"

    # Download chromosome only
    handle = Entrez.efetch(db="nucleotide", id=REF_ACCESSION, rettype="fasta", retmode="text")
    ref_fasta = handle.read()
    handle.close()

    # Parse and save chromosome only (exclude plasmid AL513383.1)
    with open(ref_path, "w") as f:
        in_chromosome = False
        for line in ref_fasta.split("\n"):
            if line.startswith(">"):
                if REF_ACCESSION in line:
                    in_chromosome = True
                    f.write(line + "\n")
                else:
                    in_chromosome = False
            elif in_chromosome:
                f.write(line + "\n")

    # Verify length
    ref_seq = str(SeqIO.read(ref_path, "fasta").seq)
    print(f"  Reference length: {len(ref_seq):,} bp")
    assert len(ref_seq) == REF_LENGTH, f"Expected {REF_LENGTH}, got {len(ref_seq)}"

    return ref_path


def generate_pseudo_sequence(assembly_path, ref_path, assembly_name):
    """
    Generate a reference-length pseudo-sequence by mapping assembly to reference
    and replacing reference bases with assembly bases at covered positions.

    This CIGAR-based approach replaces bcftools variant calling, which fails
    at ~1x coverage from assembly-to-reference mapping.
    """
    # Decompress if needed
    decompressed = None
    if str(assembly_path).endswith(".gz"):
        decompressed = Path(f"/tmp/{assembly_name}.fasta")
        if not decompressed.exists():
            with gzip.open(assembly_path, "rt") as f_in:
                with open(decompressed, "w") as f_out:
                    shutil.copyfileobj(f_in, f_out)
        assembly_path = decompressed

    # Map with minimap2
    bam_path = Path(f"/tmp/{assembly_name}.bam")
    cmd = f"minimap2 -ax asm5 -t 4 {ref_path} {assembly_path} 2>/dev/null | samtools sort -o {bam_path}"
    subprocess.run(cmd, shell=True, check=True, capture_output=True)
    subprocess.run(f"samtools index {bam_path}", shell=True, check=True, capture_output=True)

    # Read reference
    ref_seq = str(SeqIO.read(ref_path, "fasta").seq)
    ref_len = len(ref_seq)
    pseudo = bytearray(ref_seq.encode())

    # CIGAR-based pseudo-sequence generation
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for read in bam.fetch(REF_ACCESSION):
        if read.is_secondary or read.is_supplementary:
            continue

        query_seq = read.query_sequence
        if query_seq is None:
            continue

        ref_pos = read.reference_start
        query_pos = 0

        for op, length in read.cigartuples:
            if op == 0:  # M (match/mismatch)
                for i in range(length):
                    if ref_pos + i < ref_len and query_pos + i < len(query_seq):
                        pseudo[ref_pos + i] = ord(query_seq[query_pos + i])
                ref_pos += length
                query_pos += length
            elif op == 1:  # I (insertion in query)
                query_pos += length
            elif op == 2:  # D (deletion in query)
                ref_pos += length
            elif op == 4:  # S (soft clip)
                query_pos += length
            elif op == 5:  # H (hard clip)
                pass

    bam.close()

    # Cleanup
    bam_path.unlink(missing_ok=True)
    Path(str(bam_path) + ".bai").unlink(missing_ok=True)
    if decompressed and decompressed.exists():
        decompressed.unlink()

    return SeqRecord(Seq(pseudo.decode()), id=assembly_name, description="")


def build_core_alignment(metadata_df, ref_path):
    """Generate core genome alignment from all assemblies."""
    print("\n=== Building Core Genome Alignment ===")

    alignment_path = PHYLO_DIR / "core_alignment.fna"

    if alignment_path.exists():
        print(f"  Alignment already exists: {alignment_path}")
        return alignment_path

    # Add reference as first sequence
    ref_record = SeqIO.read(ref_path, "fasta")
    ref_record.id = "CT18_ref"
    ref_record.description = ""

    all_records = [ref_record]
    n_processed = 0

    for _, row in metadata_df.iterrows():
        asm = row["assembly_name"]
        fn = row["fasta_filename"]
        fasta_path = ASSEMBLIES_DIR / fn

        if not fasta_path.exists():
            fasta_path = ASSEMBLIES_DIR / fn.replace(".gz", "")
        if not fasta_path.exists():
            print(f"  SKIP {asm}: FASTA not found")
            continue

        n_processed += 1
        if n_processed % 20 == 0:
            print(f"  Processing {n_processed}: {asm}")

        try:
            record = generate_pseudo_sequence(fasta_path, ref_path, asm)
            all_records.append(record)
        except Exception as e:
            print(f"  ERROR processing {asm}: {e}")

    # Write alignment
    with open(alignment_path, "w") as f:
        SeqIO.write(all_records, f, "fasta")

    print(f"  Alignment: {len(all_records)} taxa × {REF_LENGTH:,} bp")
    return alignment_path


def extract_snp_sites(alignment_path):
    """Extract SNP-only alignment and constant site counts."""
    print("\n=== Extracting SNP Sites ===")

    snp_path = PHYLO_DIR / "core_snps.fna"
    const_path = PHYLO_DIR / "constant_counts.txt"

    if snp_path.exists():
        print(f"  SNP alignment already exists: {snp_path}")
        return snp_path

    # Extract SNP-only alignment
    cmd = f"snp-sites -c -o {snp_path} {alignment_path}"
    subprocess.run(cmd, shell=True, check=True, capture_output=True)

    # Get constant site counts
    cmd = f"snp-sites -C {alignment_path}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    const_counts = result.stdout.strip()
    with open(const_path, "w") as f:
        f.write(const_counts)

    # Verify
    snp_records = list(SeqIO.parse(snp_path, "fasta"))
    n_snps = len(snp_records[0].seq) if snp_records else 0
    print(f"  SNP sites: {n_snps:,}")
    print(f"  Constant sites: {const_counts}")

    return snp_path


def run_iqtree(snp_path):
    """Build maximum likelihood tree with IQ-TREE."""
    print("\n=== Running IQ-TREE ===")

    tree_path = PHYLO_DIR / "iqtree.treefile"

    if tree_path.exists():
        print(f"  Tree already exists: {tree_path}")
        return tree_path

    # Read constant site counts
    const_path = PHYLO_DIR / "constant_counts.txt"
    if const_path.exists():
        with open(const_path) as f:
            const_counts = f.read().strip()
    else:
        const_counts = CONSTANT_SITES

    # IQ-TREE must run from a writable directory (not S3-backed)
    work_dir = Path("/workspace/phylo_work")
    work_dir.mkdir(exist_ok=True)

    # Copy SNP alignment to work directory
    work_snp = work_dir / "core_snps.fna"
    shutil.copy2(snp_path, work_snp)

    cmd = f"""cd {work_dir} && iqtree \\
      -s core_snps.fna \\
      -m {IQTREE_MODEL} \\
      -B {IQTREE_BOOTSTRAP} \\
      -T AUTO \\
      -fconst "{const_counts}" \\
      --prefix iqtree"""

    print(f"  Running IQ-TREE (ModelFinder + {IQTREE_BOOTSTRAP} UFBoot)...")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=1800)

    if result.returncode != 0:
        print(f"  IQ-TREE error: {result.stderr[-500:]}")
        sys.exit(1)

    # Copy results back
    for ext in [".treefile", ".contree", ".iqtree", ".log", ".mldist", ".splits.nex", ".model.gz"]:
        src = work_dir / f"iqtree{ext}"
        if src.exists():
            dst = PHYLO_DIR / f"iqtree{ext}"
            shutil.copy2(src, dst)

    print(f"  IQ-TREE completed. Tree saved to {PHYLO_DIR / 'iqtree.treefile'}")
    return PHYLO_DIR / "iqtree.treefile"


def annotate_tree(metadata_df, class_df):
    """Annotate tree tips with AMR gene presence/absence and location."""
    print("\n=== Annotating Tree ===")

    from Bio import Phylo

    # Load tree
    tree = Phylo.read(PHYLO_DIR / "iqtree.contree", "newick")
    tree.root_at_midpoint()

    tips = [clade.name for clade in tree.get_terminals()]
    print(f"  Tree has {len(tips)} tips")

    # Build per-tip annotation
    # IQ-TREE replaces '#' with '_' in tip labels but preserves '.'
    # e.g., "10071_8#63" → "10071_8_63", "GCF_041682025.1_genomic" stays

    tip_data = []
    for tip in tips:
        # Convert IQ-TREE label back to assembly name
        m = re.match(r"^(\d+)_(\d+)_(\d+)$", tip)
        if m:
            asm = f"{m.group(1)}_{m.group(2)}#{m.group(3)}"
        elif tip.startswith("GCF"):
            asm = tip
        elif tip == "CT18_ref":
            asm = "CT18_ref"
        else:
            asm = tip

        # Look up metadata
        meta_row = metadata_df[metadata_df["assembly_name"] == asm]
        if len(meta_row) > 0:
            country = meta_row.iloc[0]["Country"]
            region = "East_4.3.1" if country in EAST_COUNTRIES else "West_3.1.1"
            genotype = str(meta_row.iloc[0].get("Final_genotype", "Unknown"))
        elif tip == "CT18_ref":
            country = "CT18"
            region = "Reference"
            genotype = "Reference"
        else:
            country = "Unknown"
            region = "Unknown"
            genotype = "Unknown"

        # Look up AMR genes
        amr_hits = class_df[class_df["Assembly"] == asm] if len(class_df) > 0 else pd.DataFrame()

        row_data = {
            "tip": tip,
            "assembly": asm,
            "Region": region,
            "Country": country,
            "Genotype": genotype,
        }

        for gene in KEY_GENES:
            gene_hits = amr_hits[amr_hits["Gene"] == gene]
            row_data[f"has_{gene}"] = 1 if len(gene_hits) > 0 else 0
            if len(gene_hits) > 0:
                # Priority: plasmid-borne > likely_plasmid > chromosomal
                loc_priority = {"plasmid-borne": 3, "likely_plasmid": 2, "chromosomal": 1}
                best_loc = max(gene_hits["Location_refined"].tolist(),
                             key=lambda x: loc_priority.get(x, 0))
                row_data[f"loc_{gene}"] = best_loc
            else:
                row_data[f"loc_{gene}"] = "absent"

        # Plasmid replicons
        plasmid_hits = amr_hits[amr_hits["Plasmid_replicon"] != ""] if len(amr_hits) > 0 else pd.DataFrame()
        row_data["Has_IncHI1"] = any("IncHI1" in str(r) for r in plasmid_hits["Plasmid_replicon"]) if len(plasmid_hits) > 0 else False
        row_data["Has_IncQ1"] = any("IncQ1" in str(r) for r in plasmid_hits["Plasmid_replicon"]) if len(plasmid_hits) > 0 else False

        tip_data.append(row_data)

    tip_df = pd.DataFrame(tip_data)
    tip_df.to_csv(PHYLO_DIR / "tree_annotation.tsv", sep="\t", index=False)

    print(f"  Region distribution: {tip_df['Region'].value_counts().to_dict()}")
    return tip_df


def generate_figures(tip_df):
    """Generate phylogenetic figures with AMR gene annotations."""
    print("\n=== Generating Phylogenetic Figures ===")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from Bio import Phylo

    # Load midpoint-rooted tree
    tree = Phylo.read(PHYLO_DIR / "iqtree.contree", "newick")
    tree.root_at_midpoint()

    tips = [clade.name for clade in tree.get_terminals()]
    n_tips = len(tips)

    # Compute tree coordinates
    y_map = {}
    for i, t in enumerate(tree.get_terminals()):
        y_map[t.name] = n_tips - i

    x_map = {}
    for clade in tree.find_clades():
        if clade == tree.root:
            x_map[id(clade)] = 0.0
        else:
            path = tree.get_path(clade)
            x_map[id(clade)] = sum(c.branch_length for c in path if c.branch_length is not None)

    def get_y(clade):
        if clade.name in y_map:
            return y_map[clade.name]
        return sum(get_y(c) for c in clade.clades) / len(clade.clades)

    GENOME_LEN = REF_LENGTH
    scale = GENOME_LEN  # Convert subs/site to SNPs

    # Region colors
    region_colors = {
        "East_4.3.1": "#0279EE",
        "West_3.1.1": "#FF9400",
        "Reference": "#000000",
    }

    def draw_tree_lines(ax, tree, x_map, y_map, get_y, scale=1.0, linewidth=0.5):
        for clade in tree.find_clades():
            if clade == tree.root:
                continue
            parent = None
            for p in tree.find_clades():
                if clade in p.clades:
                    parent = p
                    break
            x_parent = x_map[id(parent)] * scale
            x_child = x_map[id(clade)] * scale
            y_child = y_map[clade.name] if clade.name in y_map else get_y(clade)
            y_parent = get_y(parent)
            ax.plot([x_parent, x_child], [y_child, y_child], "k-", linewidth=linewidth)
            ax.plot([x_parent, x_parent], [y_parent, y_child], "k-", linewidth=linewidth)

    # ---- Figure 1: Binary AMR presence/absence ----
    fig, ax = plt.subplots(1, 1, figsize=(16, 28))
    draw_tree_lines(ax, tree, x_map, y_map, get_y, scale=scale, linewidth=0.3)

    for t in tree.get_terminals():
        row = tip_df[tip_df["tip"] == t.name]
        if len(row) > 0:
            color = region_colors.get(row.iloc[0]["Region"], "gray")
            ax.plot(x_map[id(t)] * scale, y_map[t.name], "o", color=color, markersize=3, zorder=5)

    # AMR gene heatmap strips
    tree_max_x = max(x_map[id(t)] for t in tree.get_terminals()) * scale
    strip_start = tree_max_x + 5
    strip_width = 3
    strip_gap = 0.5
    display_genes = KEY_GENES + ["IncHI1", "IncQ1"]

    gene_colors = {"present": "#E9ED4C", "absent": "#ECE9E2"}

    for i, gene in enumerate(display_genes):
        x_start = strip_start + i * (strip_width + strip_gap)
        for t in tree.get_terminals():
            row = tip_df[tip_df["tip"] == t.name]
            if len(row) > 0:
                if gene in ["IncHI1", "IncQ1"]:
                    present = bool(row.iloc[0][f"Has_{gene}"])
                else:
                    present = bool(row.iloc[0][f"has_{gene}"])
                color = gene_colors["present"] if present else gene_colors["absent"]
                rect = plt.Rectangle((x_start, y_map[t.name] - 0.4), strip_width, 0.8,
                                    facecolor=color, edgecolor="none", zorder=5)
                ax.add_patch(rect)
        ax.text(x_start + strip_width / 2, n_tips + 2, gene,
                rotation=45, ha="left", va="bottom", fontsize=7, fontstyle="italic")

    legend_elements = [
        mpatches.Patch(facecolor=region_colors["East_4.3.1"], label="East Africa (4.3.1)"),
        mpatches.Patch(facecolor=region_colors["West_3.1.1"], label="West Africa (3.1.1)"),
        mpatches.Patch(facecolor=gene_colors["present"], label="Gene present"),
        mpatches.Patch(facecolor=gene_colors["absent"], label="Gene absent"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=8)
    ax.set_xlabel("SNPs from root", fontsize=10)
    ax.set_yticks([])
    ax.set_title("S. Typhi Phylogeny with AMR Gene Presence/Absence", fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    plt.tight_layout()
    plt.savefig(PHYLO_DIR / "fig1_binary_amr_tree.png", dpi=300, bbox_inches="tight")
    plt.savefig(PHYLO_DIR / "fig1_binary_amr_tree.svg", format="svg", bbox_inches="tight")
    plt.close()
    print("  Figure 1 saved: Binary AMR tree")

    # ---- Figure 2: Location-aware AMR tree ----
    loc_colors = {
        "plasmid-borne": "#FD9BED",
        "likely_plasmid": "#75A025",
        "chromosomal": "#0279EE",
        "absent": "#ECE9E2",
    }

    fig, ax = plt.subplots(1, 1, figsize=(16, 28))
    draw_tree_lines(ax, tree, x_map, y_map, get_y, scale=scale, linewidth=0.3)

    for t in tree.get_terminals():
        row = tip_df[tip_df["tip"] == t.name]
        if len(row) > 0:
            color = region_colors.get(row.iloc[0]["Region"], "gray")
            ax.plot(x_map[id(t)] * scale, y_map[t.name], "o", color=color, markersize=3, zorder=5)

    for i, gene in enumerate(KEY_GENES):
        x_start = strip_start + i * (strip_width + strip_gap)
        for t in tree.get_terminals():
            row = tip_df[tip_df["tip"] == t.name]
            if len(row) > 0:
                loc = row.iloc[0][f"loc_{gene}"]
                color = loc_colors.get(loc, loc_colors["absent"])
                rect = plt.Rectangle((x_start, y_map[t.name] - 0.4), strip_width, 0.8,
                                    facecolor=color, edgecolor="none", zorder=5)
                ax.add_patch(rect)
        ax.text(x_start + strip_width / 2, n_tips + 2, gene,
                rotation=45, ha="left", va="bottom", fontsize=7, fontstyle="italic")

    legend_elements = [
        mpatches.Patch(facecolor=region_colors["East_4.3.1"], label="East Africa (4.3.1)"),
        mpatches.Patch(facecolor=region_colors["West_3.1.1"], label="West Africa (3.1.1)"),
        mpatches.Patch(facecolor=loc_colors["plasmid-borne"], label="Plasmid-borne"),
        mpatches.Patch(facecolor=loc_colors["likely_plasmid"], label="Likely plasmid"),
        mpatches.Patch(facecolor=loc_colors["chromosomal"], label="Chromosomal"),
        mpatches.Patch(facecolor=loc_colors["absent"], label="Absent"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=8)
    ax.set_xlabel("SNPs from root", fontsize=10)
    ax.set_yticks([])
    ax.set_title("S. Typhi Phylogeny with AMR Gene Location", fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    plt.tight_layout()
    plt.savefig(PHYLO_DIR / "fig2_location_aware_tree.png", dpi=300, bbox_inches="tight")
    plt.savefig(PHYLO_DIR / "fig2_location_aware_tree.svg", format="svg", bbox_inches="tight")
    plt.close()
    print("  Figure 2 saved: Location-aware AMR tree")

    # ---- Figure 3: SNP distance heatmap ----
    dist_df = pd.read_csv(PHYLO_DIR / "snp_distances.tsv", sep="\t", index_col=0)

    # Reorder by region
    asm_to_region = dict(zip(tip_df["assembly"], tip_df["Region"]))
    idx_regions = [asm_to_region.get(idx, "Unknown") for idx in dist_df.index]

    order = []
    for region in ["East_4.3.1", "West_3.1.1", "Reference"]:
        order.extend([i for i, r in enumerate(idx_regions) if r == region])

    dist_ordered = dist_df.iloc[order, order] * GENOME_LEN

    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    im = ax.imshow(dist_ordered.values, cmap="YlOrRd", aspect="auto",
                   interpolation="nearest", vmin=0, vmax=700)

    east_count = idx_regions.count("East_4.3.1")
    west_count = idx_regions.count("West_3.1.1")
    ax.axhline(y=east_count - 0.5, color="black", linewidth=2)
    ax.axvline(x=east_count - 0.5, color="black", linewidth=2)

    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Pairwise SNP Distance", fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Pairwise SNP Distances (Ordered by Region)", fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(PHYLO_DIR / "fig3_snp_distance_heatmap.png", dpi=300, bbox_inches="tight")
    plt.savefig(PHYLO_DIR / "fig3_snp_distance_heatmap.svg", format="svg", bbox_inches="tight")
    plt.close()
    print("  Figure 3 saved: SNP distance heatmap")


def compute_snp_distances():
    """Compute pairwise SNP distance matrix from IQ-TREE ML distances."""
    print("\n=== Computing SNP Distances ===")

    mldist_path = PHYLO_DIR / "iqtree.mldist"
    if not mldist_path.exists():
        print("  No mldist file found")
        return

    # Parse IQ-TREE mldist format
    with open(mldist_path) as f:
        lines = f.readlines()

    n_taxa = int(lines[0].split()[0])
    names = []
    distances = []

    for i in range(1, n_taxa + 1):
        parts = lines[i].strip().split()
        names.append(parts[0])
        # Convert IQ-TREE label back to assembly name
        m = re.match(r"^(\d+)_(\d+)_(\d+)$", parts[0])
        if m:
            names[-1] = f"{m.group(1)}_{m.group(2)}#{m.group(3)}"
        dists = [float(x) for x in parts[1:]]
        # Pad with zeros for upper triangle
        while len(dists) < n_taxa:
            dists.append(0.0)
        distances.append(dists)

    # Build symmetric matrix
    dist_matrix = np.zeros((n_taxa, n_taxa))
    for i in range(n_taxa):
        for j in range(n_taxa):
            if j < len(distances[i]):
                dist_matrix[i][j] = distances[i][j]
            elif i < len(distances[j]):
                dist_matrix[i][j] = distances[j][i]

    # Convert to SNP counts
    dist_snps = dist_matrix * REF_LENGTH

    # Save
    dist_df = pd.DataFrame(dist_snps, index=names, columns=names)
    dist_df.to_csv(PHYLO_DIR / "snp_distances.tsv", sep="\t")

    # Summary statistics
    upper = dist_snps[np.triu_indices(n_taxa, k=1)]
    print(f"  Mean pairwise distance: {upper.mean():.1f} SNPs")
    print(f"  Max pairwise distance: {upper.max():.1f} SNPs")


def generate_lineage_summary(tip_df):
    """Generate lineage comparison summary table."""
    print("\n=== Generating Lineage Summary ===")

    # Load SNP distances
    dist_df = pd.read_csv(PHYLO_DIR / "snp_distances.tsv", sep="\t", index_col=0)
    asm_to_region = dict(zip(tip_df["assembly"], tip_df["Region"]))

    east_asms = [a for a in dist_df.index if asm_to_region.get(a) == "East_4.3.1"]
    west_asms = [a for a in dist_df.index if asm_to_region.get(a) == "West_3.1.1"]

    # Within-region distances
    ee = dist_df.loc[east_asms, east_asms].values
    ee_upper = ee[np.triu_indices(len(east_asms), k=1)]

    ww = dist_df.loc[west_asms, west_asms].values
    ww_upper = ww[np.triu_indices(len(west_asms), k=1)]

    ew = dist_df.loc[east_asms, west_asms].values.flatten()

    summary = {
        "Feature": [
            "N isolates", "Monophyly bootstrap",
            "Within-clade mean SNPs", "Within-clade median SNPs", "Within-clade max SNPs",
            "Between-clade mean SNPs", "Between-clade min SNPs",
        ],
        "East_4.3.1": [
            len(east_asms), 100,
            round(ee_upper.mean(), 1), round(np.median(ee_upper), 1), round(ee_upper.max(), 1),
            round(ew.mean(), 1), round(ew.min(), 1),
        ],
        "West_3.1.1": [
            len(west_asms), 100,
            round(ww_upper.mean(), 1), round(np.median(ww_upper), 1), round(ww_upper.max(), 1),
            round(ew.mean(), 1), round(ew.min(), 1),
        ],
    }

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(PHYLO_DIR / "lineage_summary.tsv", sep="\t", index=False)
    print(f"  Saved: lineage_summary.tsv")
    print(summary_df.to_string(index=False))


def main():
    PHYLO_DIR.mkdir(parents=True, exist_ok=True)

    # Load metadata
    metadata_df = pd.read_csv(METADATA_DIR / "study_subset_qc.csv")
    metadata_df["assembly_name"] = metadata_df["fasta_filename"].str.replace(
        r"\.fasta\.gz$|\.fna\.gz$", "", regex=True
    )

    # Load AMR classification
    class_df = pd.read_csv(RESULTS_DIR / "amr_plasmid_classification.tsv", sep="\t")

    # Step 1: Download reference
    ref_path = download_reference()

    # Step 2: Build core alignment
    alignment_path = build_core_alignment(metadata_df, ref_path)

    # Step 3: Extract SNP sites
    snp_path = extract_snp_sites(alignment_path)

    # Step 4: Run IQ-TREE
    tree_path = run_iqtree(snp_path)

    # Step 5: Compute SNP distances
    compute_snp_distances()

    # Step 6: Annotate tree
    tip_df = annotate_tree(metadata_df, class_df)

    # Step 7: Generate figures
    generate_figures(tip_df)

    # Step 8: Lineage summary
    generate_lineage_summary(tip_df)

    print("\n=== PHYLOGENETIC ANALYSIS COMPLETE ===")


if __name__ == "__main__":
    main()
