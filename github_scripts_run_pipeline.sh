#!/usr/bin/env bash
#
# run_pipeline.sh — Master pipeline for African S. Typhi AMR analysis
#
# Runs all steps in sequence:
#   01. Dataset curation and download
#   02. Quality control (QUAST)
#   03. AMR gene detection + plasmid classification
#   04. Phylogenetic analysis
#   05. Summary figures and tables
#
# Usage:
#   bash run_pipeline.sh [STEP]
#
#   If STEP is provided, runs only that step (1-5).
#   Otherwise, runs all steps.
#
# Requirements (install before running):
#   conda install -c bioconda quast amrfinderplus plasmidfinder minimap2 samtools snp-sites iqtree
#   pip install pandas matplotlib seaborn biopython pysam requests
#   amrfinder_update  (download AMRFinderPlus database)
#
# Project: Divergent Mobilome-Mediated AMR Strategies in African Salmonella Typhi
# Reference: S. Typhi CT18 (AL513382.1)

set -euo pipefail

STEP="${1:-all}"

echo "============================================"
echo "African S. Typhi AMR Analysis Pipeline"
echo "============================================"
echo "Date: $(date)"
echo "Step requested: ${STEP}"
echo ""

# Create output directories
mkdir -p data/metadata data/assemblies data/raw_tarballs data/logs
mkdir -p results/quast results/amrfinderplus results/plasmidfinder results/phylogenetics

run_step() {
    local step_num=$1
    local step_name=$2
    local step_script=$3

    echo ""
    echo "============================================"
    echo "Step ${step_num}: ${step_name}"
    echo "============================================"

    if [ ! -f "${step_script}" ]; then
        echo "ERROR: Script not found: ${step_script}"
        exit 1
    fi

    time python3 "${step_script}" 2>&1 | tee "data/logs/step${step_num}_$(date +%Y%m%d_%H%M%S).log"

    if [ $? -ne 0 ]; then
        echo "ERROR: Step ${step_num} failed. Check log above."
        exit 1
    fi

    echo "Step ${step_num} completed successfully."
}

case "${STEP}" in
    1)
        run_step 1 "Dataset Curation and Download" "01_dataset_curation.py"
        ;;
    2)
        run_step 2 "Quality Control (QUAST)" "02_quality_control.py"
        ;;
    3)
        run_step 3 "AMR Gene Detection + Plasmid Classification" "03_amr_plasmid_classification.py"
        ;;
    4)
        run_step 4 "Phylogenetic Analysis" "04_phylogenetic_analysis.py"
        ;;
    5)
        run_step 5 "Summary Figures and Tables" "05_summary_figures.py"
        ;;
    all)
        run_step 1 "Dataset Curation and Download" "01_dataset_curation.py"
        run_step 2 "Quality Control (QUAST)" "02_quality_control.py"
        run_step 3 "AMR Gene Detection + Plasmid Classification" "03_amr_plasmid_classification.py"
        run_step 4 "Phylogenetic Analysis" "04_phylogenetic_analysis.py"
        run_step 5 "Summary Figures and Tables" "05_summary_figures.py"
        ;;
    *)
        echo "Unknown step: ${STEP}. Use 1-5 or 'all'."
        exit 1
        ;;
esac

echo ""
echo "============================================"
echo "Pipeline complete!"
echo "============================================"
echo "Results saved in: results/"
echo "  - AMR classification: results/amr_plasmid_classification.tsv"
echo "  - Phylogenetics:      results/phylogenetics/"
echo "  - Summary figures:    results/*.png + results/*.svg"
echo ""
