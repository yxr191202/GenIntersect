#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# CUT&Tag data processing pipeline inspired by the official Henikoff Lab workflow
# https://www.protocols.io/view/bench-top-cut-tag-tutorial-bcuhiw6w
#
# This script assumes paired-end 150 bp FASTQ files generated on an Illumina
# NovaSeq platform and three experimental samples: EC1, EC2, and EU1.
# Update the variables in the "USER-DEFINED PARAMETERS" section as needed.
# -----------------------------------------------------------------------------
set -euo pipefail

# ============================ USER-DEFINED PARAMETERS =========================
# Directory containing raw FASTQ files (*.fastq.gz). Each sample should have
# two files, suffixed with _R1.fastq.gz and _R2.fastq.gz.
RAW_DIR="/path/to/raw_fastq"

# Output base directory for all results.
OUT_DIR="/path/to/cut_tag_results"

# Path to Bowtie2 index basename for the reference genome (e.g. hg38, mm10).
BOWTIE2_INDEX="/path/to/genome/index/bowtie2_index"

# Path to a blacklisted regions BED file (ENCODE-recommended) for filtering.
BLACKLIST_BED="/path/to/blacklist.bed"

# Effective genome size for MACS2 (e.g. 2.7e9 for human hg38).
EFFECTIVE_GENOME_SIZE="2.7e9"

# Samples to process (space-separated list matching FASTQ prefixes).
SAMPLES=(EC1 EC2 EU1)

# Number of threads to allocate to multithreaded steps.
THREADS=16

# =============================== PREPARATIONS ================================
mkdir -p "${OUT_DIR}"/{00_logs,01_fastqc_raw,02_trimmed,03_fastqc_trimmed,04_align,05_filtered,06_bigwig,07_peaks,08_consensus,09_counts}
LOG_DIR="${OUT_DIR}/00_logs"

# Helper function for logging
log(){
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_DIR}/pipeline.log"
}

log "CUT&Tag pipeline started."
log "Raw data directory: ${RAW_DIR}"
log "Output directory: ${OUT_DIR}"

# ============================ 1. QUALITY CONTROL =============================
log "Running FastQC on raw FASTQ files."
fastqc \
  --threads "${THREADS}" \
  --outdir "${OUT_DIR}/01_fastqc_raw" \
  "${RAW_DIR}"/*.fastq.gz

log "Summarizing QC reports with MultiQC."
multiqc \
  --title "CUT&Tag Raw FASTQ QC" \
  --outdir "${OUT_DIR}/01_fastqc_raw" \
  "${OUT_DIR}/01_fastqc_raw"

# ======================== 2. ADAPTOR & QUALITY TRIMMING ======================
log "Running Cutadapt for adaptor/quality trimming (official Henikoff settings)."
for SAMPLE in "${SAMPLES[@]}"; do
  R1="${RAW_DIR}/${SAMPLE}_R1.fastq.gz"
  R2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"
  TRIM_R1="${OUT_DIR}/02_trimmed/${SAMPLE}_trimmed_R1.fastq.gz"
  TRIM_R2="${OUT_DIR}/02_trimmed/${SAMPLE}_trimmed_R2.fastq.gz"

  cutadapt \
    -j "${THREADS}" \
    -q 10,10 \
    -m 25 \
    -a CTGTCTCTTATACACATCT \
    -A CTGTCTCTTATACACATCT \
    -o "${TRIM_R1}" \
    -p "${TRIM_R2}" \
    "${R1}" "${R2}" \
    2> "${LOG_DIR}/${SAMPLE}_cutadapt.log"

done

log "Running FastQC on trimmed reads."
fastqc \
  --threads "${THREADS}" \
  --outdir "${OUT_DIR}/03_fastqc_trimmed" \
  "${OUT_DIR}/02_trimmed"/*.fastq.gz

log "Summarizing trimmed QC with MultiQC."
multiqc \
  --title "CUT&Tag Trimmed FASTQ QC" \
  --outdir "${OUT_DIR}/03_fastqc_trimmed" \
  "${OUT_DIR}/03_fastqc_trimmed"

# ====================== 3. ALIGNMENT TO REFERENCE GENOME =====================
log "Aligning reads with Bowtie2 using recommended CUT&Tag parameters."
for SAMPLE in "${SAMPLES[@]}"; do
  TRIM_R1="${OUT_DIR}/02_trimmed/${SAMPLE}_trimmed_R1.fastq.gz"
  TRIM_R2="${OUT_DIR}/02_trimmed/${SAMPLE}_trimmed_R2.fastq.gz"
  SAM_OUT="${OUT_DIR}/04_align/${SAMPLE}.sam"
  BAM_OUT="${OUT_DIR}/04_align/${SAMPLE}.bam"

  bowtie2 \
    --local \
    --very-sensitive-local \
    --no-unal \
    --no-mixed \
    --no-discordant \
    --phred33 \
    -I 10 -X 700 \
    -x "${BOWTIE2_INDEX}" \
    -1 "${TRIM_R1}" \
    -2 "${TRIM_R2}" \
    -p "${THREADS}" \
    -S "${SAM_OUT}" \
    2> "${LOG_DIR}/${SAMPLE}_bowtie2.log"

  samtools view -@ "${THREADS}" -bS "${SAM_OUT}" | samtools sort -@ "${THREADS}" -o "${BAM_OUT}"
  samtools index "${BAM_OUT}"
  rm "${SAM_OUT}"
done

# ========================= 4. DEDUPLICATION & FILTERS ========================
log "Removing duplicates and mitochondrial reads."
for SAMPLE in "${SAMPLES[@]}"; do
  BAM_IN="${OUT_DIR}/04_align/${SAMPLE}.bam"
  DEDUP_BAM="${OUT_DIR}/05_filtered/${SAMPLE}.dedup.bam"
  FILTERED_BAM="${OUT_DIR}/05_filtered/${SAMPLE}.filtered.bam"

  picard MarkDuplicates \
    I="${BAM_IN}" \
    O="${DEDUP_BAM}" \
    M="${LOG_DIR}/${SAMPLE}_dedup_metrics.txt" \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=SILENT

  samtools view -@ "${THREADS}" -h "${DEDUP_BAM}" |
    awk '$3!="chrM" && $3!="MT"' |
    samtools view -@ "${THREADS}" -b -o "${FILTERED_BAM}" -

  samtools index "${FILTERED_BAM}"
  rm "${DEDUP_BAM}" "${DEDUP_BAM}.bai"
done

log "Filtering low-quality alignments and blacklisted regions."
for SAMPLE in "${SAMPLES[@]}"; do
  FILTERED_BAM="${OUT_DIR}/05_filtered/${SAMPLE}.filtered.bam"
  CLEAN_BAM="${OUT_DIR}/05_filtered/${SAMPLE}.clean.bam"

  samtools view \
    -@ "${THREADS}" \
    -b \
    -q 30 \
    "${FILTERED_BAM}" |
    bedtools intersect -v -abam stdin -b "${BLACKLIST_BED}" \
    > "${CLEAN_BAM}"

  samtools index "${CLEAN_BAM}"
done

# ====================== 5. FRAGMENT SIZE & INSERT METRICS ====================
log "Calculating fragment size distributions."
for SAMPLE in "${SAMPLES[@]}"; do
  bamPEFragmentSize \
    --bamfiles "${OUT_DIR}/05_filtered/${SAMPLE}.clean.bam" \
    --histogram "${OUT_DIR}/05_filtered/${SAMPLE}_fragment_size.pdf" \
    --maxFragmentLength 1000 \
    --numberOfProcessors "${THREADS}"
done

# =========================== 6. BIGWIG GENERATION ============================
log "Generating normalized coverage tracks (RPKM)."
for SAMPLE in "${SAMPLES[@]}"; do
  bamCoverage \
    --bam "${OUT_DIR}/05_filtered/${SAMPLE}.clean.bam" \
    --outFileName "${OUT_DIR}/06_bigwig/${SAMPLE}.rpkm.bw" \
    --outFileFormat bigwig \
    --normalizeUsing RPKM \
    --binSize 10 \
    --extendReads 200 \
    --ignoreDuplicates \
    --numberOfProcessors "${THREADS}"
done

# =========================== 7. PEAK CALLING (MACS2) =========================
log "Calling peaks with MACS2 (narrow peaks)."
for SAMPLE in "${SAMPLES[@]}"; do
  BAM_INPUT="${OUT_DIR}/05_filtered/${SAMPLE}.clean.bam"
  macs2 callpeak \
    -t "${BAM_INPUT}" \
    -f BAMPE \
    -g "${EFFECTIVE_GENOME_SIZE}" \
    --keep-dup all \
    --nomodel \
    --shift 0 \
    --extsize 200 \
    -B \
    --call-summits \
    -n "${SAMPLE}" \
    --outdir "${OUT_DIR}/07_peaks" \
    2> "${LOG_DIR}/${SAMPLE}_macs2.log"
done

# ============================ 8. CONSENSUS PEAK SET ==========================
log "Building consensus peak set across samples."
ls "${OUT_DIR}/07_peaks"/*_peaks.narrowPeak > "${OUT_DIR}/08_consensus/peak_list.txt"
mergePeaks \
  -d 200 \
  -peak "${OUT_DIR}/08_consensus/peak_list.txt" \
  > "${OUT_DIR}/08_consensus/consensus_peaks.bed"

# Convert consensus BED to SAF format required by featureCounts.
awk 'BEGIN{OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"} {print $4==""?"peak_"NR:$4,$1,$2,$3,"."}' \
  "${OUT_DIR}/08_consensus/consensus_peaks.bed" \
  > "${OUT_DIR}/08_consensus/consensus_peaks.saf"

# ======================== 9. READ COUNT MATRIX (FEATURECOUNTS) ===============
log "Quantifying reads per consensus peak."
featureCounts \
  -a "${OUT_DIR}/08_consensus/consensus_peaks.saf" \
  -o "${OUT_DIR}/09_counts/consensus_peak_counts.txt" \
  -F SAF \
  -T "${THREADS}" \
  -p \
  -B \
  -C \
  "${OUT_DIR}/05_filtered"/*.clean.bam

log "CUT&Tag pipeline completed successfully."
