#!/usr/bin/env python3
"""
Step 5: Generate Summary Figures and Tables
=============================================
Create publication-quality figures and summary tables from the
AMR classification and phylogenetic analysis results.

Figures generated:
  1. AMR gene location by region (stacked bar chart)
  2. AMR gene prevalence heatmap
  3. Plasmid replicon distribution by region (bar chart)
  4. AMR architecture comparison schematic

Usage:
  python 05_summary_figures.py

Requirements:
  pip install pandas matplotlib seaborn
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path

# ---- Configuration ----
RESULTS_DIR = Path("results")
KEY_GENES = ["blaTEM-1", "catA1", "sul1", "dfrA7", "dfrA15", "tet(B)", "sul2"]
EAST_COUNTRIES = ["Kenya", "Malawi", "Tanzania", "Uganda", "Zimbabwe",
                  "Ethiopia", "Rwanda", "Mozambique", "South Africa"]

# Phylo color palette
COLORS = {
    "plasmid-borne": "#FD9BED",
    "likely_plasmid": "#75A025",
    "chromosomal": "#0279EE",
    "absent": "#ECE9E2",
    "East_4.3.1": "#0279EE",
    "West_3.1.1": "#FF9400",
    "present": "#E9ED4C",
}


def load_data():
    """Load classification and prevalence data."""
    class_df = pd.read_csv(RESULTS_DIR / "amr_plasmid_classification.tsv", sep="\t")
    prev_df = pd.read_csv(RESULTS_DIR / "amr_gene_prevalence_by_region.tsv", sep="\t")
    repl_df = pd.read_csv(RESULTS_DIR / "plasmid_replicon_prevalence_by_region.tsv", sep="\t")
    return class_df, prev_df, repl_df


def fig_amr_location_by_region(class_df, metadata_df):
    """Figure: Stacked bar chart of AMR gene location by region."""
    print("  Generating: AMR gene location by region...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    for ax_idx, region in enumerate(["East_4.3.1", "West_3.1.1"]):
        ax = axes[ax_idx]
        region_df = class_df[class_df["Region"] == region]

        # Count isolates with each gene by location
        data = []
        for gene in KEY_GENES:
            gene_hits = region_df[region_df["Gene"] == gene]
            n_plasmid = gene_hits[gene_hits["Location_refined"] == "plasmid-borne"]["Assembly"].nunique()
            n_likely = gene_hits[gene_hits["Location_refined"] == "likely_plasmid"]["Assembly"].nunique()
            n_chrom = gene_hits[gene_hits["Location_refined"] == "chromosomal"]["Assembly"].nunique()
            total_isolates = metadata_df[
                metadata_df["Country"].apply(
                    lambda c: ("East_4.3.1" if c in EAST_COUNTRIES else "West_3.1.1")
                ) == region
            ].shape[0]
            n_absent = total_isolates - gene_hits["Assembly"].nunique()

            data.append({
                "Gene": gene,
                "Plasmid-borne": n_plasmid,
                "Likely plasmid": n_likely,
                "Chromosomal": n_chrom,
                "Absent": n_absent,
            })

        plot_df = pd.DataFrame(data).set_index("Gene")
        plot_df.plot(kind="bar", stacked=True, ax=ax,
                    color=[COLORS["plasmid-borne"], COLORS["likely_plasmid"],
                           COLORS["chromosomal"], COLORS["absent"]])
        ax.set_title(f"{region}", fontsize=12, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("Number of isolates" if ax_idx == 0 else "")
        ax.legend(fontsize=8)
        ax.tick_params(axis="x", rotation=45)

    plt.suptitle("AMR Gene Location by Region", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "amr_gene_location_by_region.png", dpi=300, bbox_inches="tight")
    plt.savefig(RESULTS_DIR / "amr_gene_location_by_region.svg", format="svg", bbox_inches="tight")
    plt.close()


def fig_amr_prevalence_heatmap(prev_df):
    """Figure: Heatmap of AMR gene prevalence by region."""
    print("  Generating: AMR gene prevalence heatmap...")

    # Pivot for heatmap
    pivot = prev_df.pivot(index="Gene", columns="Region", values="Prevalence_pct")
    pivot = pivot.reindex(KEY_GENES)

    fig, ax = plt.subplots(figsize=(6, 8))
    sns.heatmap(pivot, annot=True, fmt=".1f", cmap="YlOrRd",
                vmin=0, vmax=100, linewidths=0.5, ax=ax,
                cbar_kws={"label": "Prevalence (%)"})
    ax.set_title("AMR Gene Prevalence by Region (%)", fontsize=12, fontweight="bold")
    ax.set_ylabel("")
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "amr_gene_prevalence_heatmap.png", dpi=300, bbox_inches="tight")
    plt.savefig(RESULTS_DIR / "amr_gene_prevalence_heatmap.svg", format="svg", bbox_inches="tight")
    plt.close()


def fig_plasmid_replicon_by_region(repl_df):
    """Figure: Bar chart of plasmid replicon distribution by region."""
    print("  Generating: Plasmid replicon distribution...")

    fig, ax = plt.subplots(figsize=(8, 5))
    repl_df.pivot(index="Replicon", columns="Region", values="Prevalence_pct").plot(
        kind="bar", ax=ax, color=[COLORS["East_4.3.1"], COLORS["West_3.1.1"]]
    )
    ax.set_ylabel("Prevalence (%)")
    ax.set_xlabel("")
    ax.set_title("Plasmid Replicon Prevalence by Region", fontsize=12, fontweight="bold")
    ax.legend(fontsize=10)
    ax.tick_params(axis="x", rotation=0)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "plasmid_replicon_by_region.png", dpi=300, bbox_inches="tight")
    plt.savefig(RESULTS_DIR / "plasmid_replicon_by_region.svg", format="svg", bbox_inches="tight")
    plt.close()


def fig_amr_architecture_comparison():
    """Figure: Schematic comparing East vs West AMR architecture."""
    print("  Generating: AMR architecture comparison...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # East 4.3.1: Chromosomal MDR island
    ax = axes[0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title("East Africa — Genotype 4.3.1\nChromosomal MDR Island", fontsize=11, fontweight="bold")

    # Chromosome
    ax.barh(5, 8, height=1.5, left=1, color="#0279EE", alpha=0.3, edgecolor="black")
    ax.text(5, 5, "Chromosome", ha="center", va="center", fontsize=10, fontweight="bold")

    # MDR island (integrated)
    ax.barh(5, 2, height=1.2, left=3, color="#75A025", alpha=0.7, edgecolor="black")
    ax.text(4, 5, "MDR\nisland", ha="center", va="center", fontsize=8, color="white", fontweight="bold")

    # IncQ1 plasmid
    ax.barh(2, 1.5, height=0.8, left=4, color="#FD9BED", alpha=0.7, edgecolor="black")
    ax.text(4.75, 2, "IncQ1\n(sul2)", ha="center", va="center", fontsize=7)

    # Labels
    ax.text(1, 8.5, "blaTEM-1: chromosomal (100%)", fontsize=8)
    ax.text(1, 7.8, "dfrA7: chromosomal/likely_plasmid (80%/20%)", fontsize=8)
    ax.text(1, 7.1, "IncHI1: degraded (17%)", fontsize=8)
    ax.text(1, 6.4, "IncQ1: intact (65%)", fontsize=8)
    ax.axis("off")

    # West 3.1.1: Plasmid-borne MDR
    ax = axes[1]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title("West Africa — Genotype 3.1.1\nPlasmid-Borne MDR", fontsize=11, fontweight="bold")

    # Chromosome
    ax.barh(7, 8, height=1.5, left=1, color="#0279EE", alpha=0.3, edgecolor="black")
    ax.text(5, 7, "Chromosome", ha="center", va="center", fontsize=10, fontweight="bold")

    # IncHI1 plasmid
    ax.barh(3.5, 6, height=1.5, left=2, color="#FD9BED", alpha=0.7, edgecolor="black")
    ax.text(5, 3.5, "IncHI1 Plasmid", ha="center", va="center", fontsize=10, fontweight="bold")

    # MDR cassette on plasmid
    ax.barh(3.5, 2.5, height=1.2, left=3.5, color="#75A025", alpha=0.7, edgecolor="black")
    ax.text(4.75, 3.5, "MDR\ncassette", ha="center", va="center", fontsize=8, color="white", fontweight="bold")

    # Labels
    ax.text(1, 1.5, "blaTEM-1: plasmid-borne (61%)", fontsize=8)
    ax.text(1, 0.8, "dfrA15: likely_plasmid (89%)", fontsize=8)
    ax.text(1, 0.1, "tet(B): plasmid-borne (47%)", fontsize=8)
    ax.text(1, -0.6, "IncHI1: intact (67%)", fontsize=8)
    ax.axis("off")

    plt.suptitle("Divergent AMR Architectures in African S. Typhi", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "amr_architecture_comparison.png", dpi=300, bbox_inches="tight")
    plt.savefig(RESULTS_DIR / "amr_architecture_comparison.svg", format="svg", bbox_inches="tight")
    plt.close()


def main():
    class_df, prev_df, repl_df = load_data()
    metadata_df = pd.read_csv(Path("data") / "metadata" / "study_subset_qc.csv")

    print("=== Generating Summary Figures ===")
    fig_amr_location_by_region(class_df, metadata_df)
    fig_amr_prevalence_heatmap(prev_df)
    fig_plasmid_replicon_by_region(repl_df)
    fig_amr_architecture_comparison()

    print("\n=== SUMMARY FIGURES COMPLETE ===")


if __name__ == "__main__":
    main()
