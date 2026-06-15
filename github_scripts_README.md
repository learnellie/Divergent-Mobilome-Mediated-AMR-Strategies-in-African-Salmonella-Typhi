# Divergent Mobilome-Mediated AMR Strategies in African Salmonella Typhi

## Chromosomal Integration in East African Genotype 4.3.1 vs Plasmid-Borne MDR in West African Genotype 3.1.1

---

### Overview

This repository contains the complete analysis pipeline for comparing antimicrobial resistance (AMR) gene mobilization strategies between two dominant African *Salmonella* Typhi lineages:

- **East African genotype 4.3.1 (H58)**: MDR genes are chromosomally integrated via a degraded IncHI1 island, with IncQ1 small plasmids carrying sul2
- **West African genotype 3.1.1**: MDR genes are carried on intact, transmissible IncHI1 plasmids with a characteristic MDR cassette

The analysis confirms these lineages form distinct monophyletic clades (bootstrap 100) separated by ~467 mean SNPs, with lineage-specific AMR gene repertoires (dfrA7/dfrA15, IncQ1/IncHI1) and divergent plasmid architectures.

### Key Findings

| Feature | East 4.3.1 | West 3.1.1 |
|---------|-----------|-----------|
| Isolates | 135 | 48 |
| Monophyly bootstrap | 100 | 100 |
| Within-clade mean SNPs | 104.6 | 56.2 |
| Between-clade mean SNPs | 466.5 | 466.5 |
| blaTEM-1 location | 100% chromosomal | 61% plasmid-borne |
| dfrA7/dfrA15 | dfrA7 (70%) | dfrA15 (73%) |
| IncHI1 prevalence | 17% (degraded) | 67% (intact) |
| IncQ1 prevalence | 65% | 0% |
| MDR architecture | Chromosomal island | Plasmid-borne cassette |

### Data Sources

1. **GTGC (Global Typhoid Genomics Consortium)** — Carey et al. 2023, eLife
   - FigShare DOI: [10.26180/21431883](https://doi.org/10.26180/21431883)
   - 183 draft assemblies from 7 East + 4 West African countries

2. **Zuza et al. 2024** — Gates Open Research
   - NCBI BioProject: PRJNA1127853
   - 4 complete genome assemblies (Malawi, long-read)

3. **Reference**: *S*. Typhi CT18 (AL513382.1, 4,809,037 bp chromosome)

### Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_dataset_curation.py` | Download metadata, filter to 184-isolate study subset, download assemblies |
| 2 | `02_quality_control.py` | QUAST quality assessment, contamination detection, QC filtering |
| 3 | `03_amr_plasmid_classification.py` | AMRFinderPlus + PlasmidFinder, cross-reference for plasmid vs chromosomal classification |
| 4 | `04_phylogenetic_analysis.py` | Core genome SNP phylogeny, IQ-TREE, tree annotation with AMR genes |
| 5 | `05_summary_figures.py` | Publication-quality figures and summary tables |

Run all steps:
```bash
bash run_pipeline.sh
```

Run a single step:
```bash
bash run_pipeline.sh 4    # Phylogenetic analysis only
```

### Tool Parameters

| Tool | Version | Parameters |
|------|---------|------------|
| AMRFinderPlus | 4.2.7 (DB 2026-05-15.1) | `--organism Salmonella --ident_min 0.9 --coverage_min 0.5` |
| PlasmidFinder | 2.1.6 | `-d enterobacteriaceae -t 0.95 -l 0.60 -x` |
| QUAST | 5.2.0 | `--min-contig 0 --fast` |
| minimap2 | 2.31 | `-ax asm5 -t 4` |
| IQ-TREE | 3.1.2 | `-m MFP -B 1000 -T AUTO -fconst "1149511,1249782,1253591,1152957"` |
| snp-sites | 2.5.1 | `-c` (ACGT-only SNP extraction) |

### Classification Methodology

**Primary classification**: AMR gene contig cross-referenced with PlasmidFinder hits
- If AMR gene is on a contig with a PlasmidFinder replicon hit → **plasmid-borne**
- Otherwise → **chromosomal**

**Refinement — "likely_plasmid"**: Short-read assemblies fragment large plasmids across multiple contigs. MDR cassette genes (blaTEM-1, catA1, sul1, dfrA15, dfrA7, tet(B), qacEdelta1) on short contigs (<15 kb) with ≥2 MDR genes in IncHI1-positive isolates are classified as **likely_plasmid**, representing fragments of the IncHI1 MDR cassette.

**Final classification counts**: 245 plasmid-borne / 196 likely_plasmid / 722 chromosomal (1,163 total AMR hits)

### Output Files

```
results/
├── amr_plasmid_classification.tsv          # Full 1,163-hit classification with metadata
├── amr_gene_prevalence_by_region.tsv       # Gene prevalence by region with location breakdown
├── plasmid_replicon_prevalence_by_region.tsv  # Replicon prevalence by region
├── isolate_amr_summary.tsv                 # Per-isolate AMR summary
├── amr_summary_by_country.tsv              # Per-country AMR summary
├── amr_gene_location_by_region.png/svg     # Stacked bar chart
├── amr_gene_prevalence_heatmap.png/svg     # Heatmap
├── plasmid_replicon_by_region.png/svg      # Bar chart
├── amr_architecture_comparison.png/svg     # Schematic
└── phylogenetics/
    ├── core_alignment.fna                  # Full core genome alignment (184 taxa × 4.8 Mb)
    ├── core_snps.fna                       # SNP-only alignment (3,196 sites)
    ├── iqtree.treefile                     # ML tree with bootstrap support
    ├── iqtree.contree                      # Bootstrap consensus tree
    ├── iqtree.iqtree                       # Full IQ-TREE report
    ├── snp_distances.tsv                   # Pairwise SNP distance matrix
    ├── tree_annotation.tsv                 # Tree annotation with AMR + location
    ├── lineage_summary.tsv                 # Lineage comparison table
    ├── fig1_binary_amr_tree.png/svg        # Binary AMR presence/absence tree
    ├── fig2_location_aware_tree.png/svg    # Location-aware AMR tree
    └── fig3_snp_distance_heatmap.png/svg   # SNP distance heatmap
```

### Requirements

```bash
# Conda (bioconda channel)
conda install -c bioconda quast amrfinderplus plasmidfinder minimap2 samtools snp-sites iqtree

# Python packages
pip install pandas matplotlib seaborn biopython pysam requests

# AMRFinderPlus database
amrfinder_update
```

### Citation

If you use this pipeline, please cite the data sources:

- Carey et al. 2023. "Global phylogenomics of *Salmonella* Typhi reveals a new H58 clade replacing the classical H58 clade." eLife.
- Zuza et al. 2024. Gates Open Research.
- Parkhill et al. 2001. "Complete genome sequence of a multiple drug resistant *Salmonella* enterica serovar Typhi CT18." Nature. (CT18 reference)

### License

MIT
