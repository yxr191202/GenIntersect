# CUT&Tag Data Analysis Workflow

This folder contains a complete CUT&Tag processing and downstream analysis workflow customized for three paired-end NovaSeq PE150 samples (EC1, EC2, EU1) and closely following the official Henikoff lab protocol.

## Contents

- `cut_tag_pipeline.sh` – Shell pipeline that processes raw FASTQ files through QC, trimming, alignment, filtering, peak calling, consensus peak generation, and count matrix creation.
- `cut_tag_downstream_analysis.R` – R script that ingests the consensus peak counts, performs DESeq2-based differential analysis, annotates peaks, and generates publication-ready plots.

## Usage Overview

1. **Edit user-defined parameters** in both scripts to match your directory structure, reference genome, and organism-specific annotation packages.
2. Run the shell pipeline: `bash cut_tag_pipeline.sh`.
3. After the count matrix is produced, execute the R downstream script: `Rscript cut_tag_downstream_analysis.R`.

Refer to the inline comments for detailed explanations of each step and relevant tool parameters.
