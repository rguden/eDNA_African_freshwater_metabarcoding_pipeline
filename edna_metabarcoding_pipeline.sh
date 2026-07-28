#!/bin/bash -l

# ===============================================================================
# eDNA Illumina workflow: metabarcoding pipeline from African freshwater systems
# Author: Rodgee Mae Guden
# ===============================================================================
#
# Markers supported: COI, 12S, 16S, 18S
#
# Main steps:
#   concat          Concatenate lane-level reads per sample
#   qc_raw          FastQC on concatenated raw reads
#   fastp           Quality filtering with fastp
#   qc_filtered     FastQC/MultiQC on filtered reads
#   demux           Marker demultiplexing with cutadapt
#   asv             ASV inference with QIIME2/DADA2
#   otu             OTU clustering with VSEARCH
#   taxonomy_asv    BLAST + TaxonKit taxonomy for ASVs
#   taxonomy_otu    BLAST + TaxonKit taxonomy for OTUs
#   multiqc         MultiQC reports
#
# Example:
#   ./edna_metabarcoding_pipeline.sh -s concat,qc_raw,fastp,qc_filtered,demux,asv,taxonomy_asv
#
# ===============================================================================

set -euo pipefail

# ----------------- CONFIG -----------------

READS_DIR="${READS_DIR:-/staging/leuven/stg_00100/Rodgee/metabarcoding_data/raw_reads}"
METADATA="${METADATA:-/staging/leuven/stg_00100/Rodgee/metabarcoding_data/metadata/samples.tsv}"
PRIMER_DIR="${PRIMER_DIR:-/staging/leuven/stg_00100/Rodgee/metabarcoding_data/primers}"
OUTPUT_DIR="${OUTPUT_DIR:-/staging/leuven/stg_00100/Rodgee/metabarcoding_data/output}"

THREADS="${THREADS:-32}"
FORCE="${FORCE:-false}"

# Analysis mode used only when running the "all" shortcut
# Options: asv, otu, both
MODE="${MODE:-both}"

# ----------------- FASTP CONFIG -----------------

AVERAGE_QUAL="${AVERAGE_QUAL:-35}"
LENGTH_REQUIRED="${LENGTH_REQUIRED:-50}"
CUT_FRONT_WINDOW_SIZE="${CUT_FRONT_WINDOW_SIZE:-20}"
CUT_FRONT_MEAN_QUAL="${CUT_FRONT_MEAN_QUAL:-30}"

# ----------------- CUTADAPT CONFIG -----------------

CUTADAPT_ERROR_RATE="${CUTADAPT_ERROR_RATE:-0.10}"
CUTADAPT_OVERLAP="${CUTADAPT_OVERLAP:-20}"

# If cutadapt >= 4 is available and reads may occur in both orientations, keep --revcomp.
# If your cutadapt version does not support --revcomp, set CUTADAPT_EXTRA="" when running.
CUTADAPT_EXTRA="${CUTADAPT_EXTRA:---revcomp}"

# ----------------- OTU CONFIG -----------------

OTU_ID="${OTU_ID:-0.97}"

# ----------------- TAXONOMY CONFIG -----------------

# BLAST_DB should be the BLAST database prefix, not the FASTA file.
BLAST_DB="${BLAST_DB:-/staging/leuven/stg_00100/Rodgee/databases/metabarcoding/custom_silva_ncbi_metabarcoding}"

# Optional: if provided, the pipeline will build a BLAST database from this FASTA.
CUSTOM_DB_FASTA="${CUSTOM_DB_FASTA:-}"

BLAST_MIN_ID="${BLAST_MIN_ID:-80}"
BLAST_EVALUE="${BLAST_EVALUE:-0.001}"
BLAST_MAX_TARGETS="${BLAST_MAX_TARGETS:-10}"

# TaxonKit database directory containing nodes.dmp and names.dmp
TAXONKIT_DB="${TAXONKIT_DB:-/staging/leuven/stg_00100/Rodgee/databases/taxdump}"

# ----------------- MARKERS -----------------

MARKERS=("COI" "12S" "16S" "18S")

# Optional ASV length thresholds for downstream filtering, if needed.
# These are not automatically applied in this version, but are documented here
# so that the filtering step can be added consistently if required.
declare -A MIN_ASV_LEN
MIN_ASV_LEN["COI"]=160
MIN_ASV_LEN["12S"]=135
MIN_ASV_LEN["16S"]=170
MIN_ASV_LEN["18S"]=167

# ----------------- LOGGING -----------------

LOG_PATH="$OUTPUT_DIR/logs"
LOG_FILE="$LOG_PATH/edna_metabarcoding_pipeline.log"

mkdir -p "$OUTPUT_DIR" "$LOG_PATH"
touch "$LOG_FILE"

log() {
    command echo "[$(date '+%F %T')] $*" | tee -a "$LOG_FILE"
}

trap 'log "ERROR: Pipeline failed at line $LINENO"; exit 1' ERR

# ----------------- HELP -----------------

usage() {
cat <<EOF
Usage:
  $0 -s <step1,step2,...>

Available steps:
  check
  concat
  qc_raw
  fastp
  qc_filtered
  demux
  asv
  otu
  taxonomy_asv
  taxonomy_otu
  multiqc
  all

Example:
  $0 -s concat,qc_raw,fastp,qc_filtered,demux,asv,taxonomy_asv

Important environment variables:
  READS_DIR       Directory containing raw FASTQ files
  METADATA        Tab-delimited sample sheet
  PRIMER_DIR      Directory containing marker primer FASTA files
  OUTPUT_DIR      Output directory
  THREADS         Number of threads
  MODE            asv, otu, or both, used by all
  BLAST_DB        Local BLAST database prefix
  TAXONKIT_DB     Directory containing NCBI taxdump files
  FORCE           true/false; overwrite outputs when possible

EOF
}

# ------------------------------------------------------------
# Parse steps
# ------------------------------------------------------------

echo "RAW ARGS: $@"

STEPS=()
while getopts "s:h" opt; do
    case ${opt} in
        s) IFS=',' read -r -a STEPS <<< "$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [ ${#STEPS[@]} -eq 0 ]; then
    STEPS=(check)
fi

log "Steps to run: ${STEPS[*]}"

# ============================================================
# HELPER FUNCTIONS
# ============================================================

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || {
        log "ERROR: Required command not found: $1"
        exit 1
    }
}

make_dirs() {
    mkdir -p \
        "$OUTPUT_DIR/concat" \
        "$OUTPUT_DIR/qc/raw" \
        "$OUTPUT_DIR/qc/filtered" \
        "$OUTPUT_DIR/fastp" \
        "$OUTPUT_DIR/demux" \
        "$OUTPUT_DIR/qiime/manifests" \
        "$OUTPUT_DIR/qiime/imported" \
        "$OUTPUT_DIR/qiime/dada2" \
        "$OUTPUT_DIR/qiime/exported" \
        "$OUTPUT_DIR/otu" \
        "$OUTPUT_DIR/taxonomy" \
        "$OUTPUT_DIR/multiqc" \
        "$OUTPUT_DIR/logs/fastqc" \
        "$OUTPUT_DIR/logs/fastp" \
        "$OUTPUT_DIR/logs/cutadapt" \
        "$OUTPUT_DIR/logs/qiime" \
        "$OUTPUT_DIR/logs/vsearch" \
        "$OUTPUT_DIR/logs/blast"
}

sample_ids() {
    awk -F'\t' '
        NR == 1 && $1 == "sample_id" {next}
        NF > 0 && $1 !~ /^#/ {print $1}
    ' "$METADATA"
}

get_sample_row() {
    local SAMPLE="$1"
    awk -v s="$SAMPLE" -F'\t' '
        NR == 1 && $1 == "sample_id" {next}
        $1 == s {print; exit}
    ' "$METADATA"
}

get_read_files() {
    local SAMPLE="$1"
    local READ="$2"

    local ROW
    ROW="$(get_sample_row "$SAMPLE")"

    local R1_GLOB=""
    local R2_GLOB=""

    if [[ -n "$ROW" ]]; then
        R1_GLOB="$(echo "$ROW" | awk -F'\t' '{print $2}')"
        R2_GLOB="$(echo "$ROW" | awk -F'\t' '{print $3}')"
    fi

    if [[ "$READ" == "R1" && -n "${R1_GLOB:-}" ]]; then
        compgen -G "$READS_DIR/$R1_GLOB" | sort
    elif [[ "$READ" == "R2" && -n "${R2_GLOB:-}" ]]; then
        compgen -G "$READS_DIR/$R2_GLOB" | sort
    else
        find "$READS_DIR" -type f \( \
            -name "*${SAMPLE}*_${READ}_*.fastq.gz" -o \
            -name "*${SAMPLE}*_${READ}.fastq.gz" -o \
            -name "*${SAMPLE}*${READ}*.fq.gz" -o \
            -name "*${SAMPLE}*${READ}*.fastq.gz" \
        \) | sort
    fi
}

select_samples_for_array_or_sequential() {
    mapfile -t ALL_SAMPLES < <(sample_ids)
    N_SAMPLES=${#ALL_SAMPLES[@]}

    if [[ "$N_SAMPLES" -eq 0 ]]; then
        log "ERROR: No samples found in metadata file: $METADATA"
        exit 1
    fi

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        IDX=$((SLURM_ARRAY_TASK_ID - 1))

        if [[ "$IDX" -ge "$N_SAMPLES" ]]; then
            log "ERROR: SLURM_ARRAY_TASK_ID exceeds number of samples"
            exit 1
        fi

        SAMPLES_TO_RUN=("${ALL_SAMPLES[$IDX]}")
        log "Array mode detected: processing sample ${SAMPLES_TO_RUN[0]}"
    else
        SAMPLES_TO_RUN=("${ALL_SAMPLES[@]}")
        log "No SLURM array detected → running sequentially over $N_SAMPLES samples"
    fi
}

check_file_nonempty() {
    local FILE="$1"
    [[ -s "$FILE" ]]
}

# ============================================================
# STEP LOOP
# ============================================================

for STEP in "${STEPS[@]}"; do

log "Starting step: $STEP"

case "$STEP" in

# ============================================================
### CHECK
# ============================================================
check)

set -euo pipefail

log "Checking input files and dependencies"

[[ -d "$READS_DIR" ]] || { log "ERROR: READS_DIR not found: $READS_DIR"; exit 1; }
[[ -f "$METADATA" ]] || { log "ERROR: METADATA file not found: $METADATA"; exit 1; }
[[ -d "$PRIMER_DIR" ]] || { log "ERROR: PRIMER_DIR not found: $PRIMER_DIR"; exit 1; }

for MARKER in "${MARKERS[@]}"; do
    [[ -f "$PRIMER_DIR/${MARKER}_F.fasta" ]] || { log "ERROR: Missing primer file: $PRIMER_DIR/${MARKER}_F.fasta"; exit 1; }
    [[ -f "$PRIMER_DIR/${MARKER}_R.fasta" ]] || { log "ERROR: Missing primer file: $PRIMER_DIR/${MARKER}_R.fasta"; exit 1; }
done

require_cmd fastqc
require_cmd multiqc
require_cmd fastp
require_cmd cutadapt
require_cmd vsearch
require_cmd blastn
require_cmd makeblastdb

if [[ "$MODE" == "asv" || "$MODE" == "both" ]]; then
    require_cmd qiime
fi

if command -v taxonkit >/dev/null 2>&1; then
    log "TaxonKit found"
else
    log "WARNING: TaxonKit not found. Lineage formatting will be skipped."
fi

make_dirs

log "Input and dependency check completed"
;;

# ============================================================
### CONCATENATE LANES
# ============================================================
concat)

set -euo pipefail
shopt -s nullglob

make_dirs
select_samples_for_array_or_sequential

log "Concatenating lane-level reads per sample"

for SAMPLE in "${SAMPLES_TO_RUN[@]}"; do

    OUT_R1="$OUTPUT_DIR/concat/${SAMPLE}_R1.fastq.gz"
    OUT_R2="$OUTPUT_DIR/concat/${SAMPLE}_R2.fastq.gz"

    if [[ -s "$OUT_R1" && -s "$OUT_R2" && "$FORCE" != "true" ]]; then
        log "Skipping $SAMPLE concat: outputs already exist"
        continue
    fi

    mapfile -t R1_FILES < <(get_read_files "$SAMPLE" "R1")
    mapfile -t R2_FILES < <(get_read_files "$SAMPLE" "R2")

    if [[ "${#R1_FILES[@]}" -eq 0 || "${#R2_FILES[@]}" -eq 0 ]]; then
        log "ERROR: Missing R1 or R2 files for sample $SAMPLE"
        exit 1
    fi

    if [[ "${#R1_FILES[@]}" -ne "${#R2_FILES[@]}" ]]; then
        log "ERROR: Unequal R1/R2 file counts for sample $SAMPLE"
        log "R1 count: ${#R1_FILES[@]}"
        log "R2 count: ${#R2_FILES[@]}"
        exit 1
    fi

    log "$SAMPLE: concatenating ${#R1_FILES[@]} R1 files and ${#R2_FILES[@]} R2 files"

    cat "${R1_FILES[@]}" > "$OUT_R1"
    cat "${R2_FILES[@]}" > "$OUT_R2"

done

log "Concatenation completed"
;;

# ============================================================
### RAW QC
# ============================================================
qc_raw)

set -euo pipefail
shopt -s nullglob

THREADS="${THREADS:-8}"

make_dirs
select_samples_for_array_or_sequential

QC_DIR="$OUTPUT_DIR/qc/raw"
LOG_DIR="$OUTPUT_DIR/logs/fastqc"

log "Starting FastQC on concatenated raw reads"

for SAMPLE in "${SAMPLES_TO_RUN[@]}"; do

    READS1="$OUTPUT_DIR/concat/${SAMPLE}_R1.fastq.gz"
    READS2="$OUTPUT_DIR/concat/${SAMPLE}_R2.fastq.gz"

    [[ -s "$READS1" ]] || { log "ERROR: Missing $READS1. Run concat first."; exit 1; }
    [[ -s "$READS2" ]] || { log "ERROR: Missing $READS2. Run concat first."; exit 1; }

    log "FastQC raw reads: $SAMPLE"

    fastqc \
        -t "$THREADS" \
        -o "$QC_DIR" \
        "$READS1" "$READS2" \
        &>> "$LOG_DIR/${SAMPLE}.raw.fastqc.log"

done

log "FastQC raw reads completed"
;;

# ============================================================
### FASTP
# ============================================================
fastp)

set -euo pipefail
shopt -s nullglob

make_dirs
select_samples_for_array_or_sequential

FASTP_DIR="$OUTPUT_DIR/fastp"
LOG_DIR="$OUTPUT_DIR/logs/fastp"

log "Starting fastp quality filtering"
log "AVERAGE_QUAL=$AVERAGE_QUAL; LENGTH_REQUIRED=$LENGTH_REQUIRED"

for SAMPLE in "${SAMPLES_TO_RUN[@]}"; do

    IN_R1="$OUTPUT_DIR/concat/${SAMPLE}_R1.fastq.gz"
    IN_R2="$OUTPUT_DIR/concat/${SAMPLE}_R2.fastq.gz"

    OUT_R1="$FASTP_DIR/${SAMPLE}_filtered_R1.fastq.gz"
    OUT_R2="$FASTP_DIR/${SAMPLE}_filtered_R2.fastq.gz"

    HTML="$FASTP_DIR/${SAMPLE}_fastp.html"
    JSON="$FASTP_DIR/${SAMPLE}_fastp.json"

    [[ -s "$IN_R1" ]] || { log "ERROR: Missing $IN_R1. Run concat first."; exit 1; }
    [[ -s "$IN_R2" ]] || { log "ERROR: Missing $IN_R2. Run concat first."; exit 1; }

    if [[ -s "$OUT_R1" && -s "$OUT_R2" && "$FORCE" != "true" ]]; then
        log "Skipping $SAMPLE fastp: outputs already exist"
        continue
    fi

    log "Running fastp: $SAMPLE"

    fastp \
        --in1 "$IN_R1" \
        --in2 "$IN_R2" \
        --out1 "$OUT_R1" \
        --out2 "$OUT_R2" \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --cut_front \
        --cut_front_window_size "$CUT_FRONT_WINDOW_SIZE" \
        --cut_front_mean_quality "$CUT_FRONT_MEAN_QUAL" \
        --average_qual "$AVERAGE_QUAL" \
        --length_required "$LENGTH_REQUIRED" \
        --low_complexity_filter \
        --overrepresentation_analysis \
        --thread "$THREADS" \
        --html "$HTML" \
        --json "$JSON" \
        &>> "$LOG_DIR/${SAMPLE}.fastp.log"

done

log "fastp completed"
;;

# ============================================================
### FILTERED QC
# ============================================================
qc_filtered)

set -euo pipefail
shopt -s nullglob

make_dirs
select_samples_for_array_or_sequential

QC_DIR="$OUTPUT_DIR/qc/filtered"
LOG_DIR="$OUTPUT_DIR/logs/fastqc"

log "Starting FastQC on filtered reads"

for SAMPLE in "${SAMPLES_TO_RUN[@]}"; do

    READS1="$OUTPUT_DIR/fastp/${SAMPLE}_filtered_R1.fastq.gz"
    READS2="$OUTPUT_DIR/fastp/${SAMPLE}_filtered_R2.fastq.gz"

    [[ -s "$READS1" ]] || { log "ERROR: Missing $READS1. Run fastp first."; exit 1; }
    [[ -s "$READS2" ]] || { log "ERROR: Missing $READS2. Run fastp first."; exit 1; }

    log "FastQC filtered reads: $SAMPLE"

    fastqc \
        -t "$THREADS" \
        -o "$QC_DIR" \
        "$READS1" "$READS2" \
        &>> "$LOG_DIR/${SAMPLE}.filtered.fastqc.log"

done

log "FastQC filtered reads completed"
;;

# ============================================================
### MARKER DEMULTIPLEXING
# ============================================================
demux)

set -euo pipefail
shopt -s nullglob

make_dirs
select_samples_for_array_or_sequential

LOG_DIR="$OUTPUT_DIR/logs/cutadapt"

log "Starting marker demultiplexing with cutadapt"
log "CUTADAPT_ERROR_RATE=$CUTADAPT_ERROR_RATE; CUTADAPT_OVERLAP=$CUTADAPT_OVERLAP"

for MARKER in "${MARKERS[@]}"; do

    mkdir -p "$OUTPUT_DIR/demux/$MARKER"

    FWD="$PRIMER_DIR/${MARKER}_F.fasta"
    REV="$PRIMER_DIR/${MARKER}_R.fasta"

    for SAMPLE in "${SAMPLES_TO_RUN[@]}"; do

        IN_R1="$OUTPUT_DIR/fastp/${SAMPLE}_filtered_R1.fastq.gz"
        IN_R2="$OUTPUT_DIR/fastp/${SAMPLE}_filtered_R2.fastq.gz"

        OUT_R1="$OUTPUT_DIR/demux/$MARKER/${SAMPLE}_${MARKER}_R1.fastq.gz"
        OUT_R2="$OUTPUT_DIR/demux/$MARKER/${SAMPLE}_${MARKER}_R2.fastq.gz"

        [[ -s "$IN_R1" ]] || { log "ERROR: Missing $IN_R1. Run fastp first."; exit 1; }
        [[ -s "$IN_R2" ]] || { log "ERROR: Missing $IN_R2. Run fastp first."; exit 1; }

        if [[ -s "$OUT_R1" && -s "$OUT_R2" && "$FORCE" != "true" ]]; then
            log "Skipping $SAMPLE $MARKER demux: outputs already exist"
            continue
        fi

        log "Demultiplexing sample $SAMPLE marker $MARKER"

        cutadapt \
            $CUTADAPT_EXTRA \
            -e "$CUTADAPT_ERROR_RATE" \
            -O "$CUTADAPT_OVERLAP" \
            -g "^file:$FWD" \
            -G "^file:$REV" \
            --discard-untrimmed \
            --cores "$THREADS" \
            -o "$OUT_R1" \
            -p "$OUT_R2" \
            "$IN_R1" \
            "$IN_R2" \
            &> "$LOG_DIR/${SAMPLE}.${MARKER}.cutadapt.log"

    done
done

log "Marker demultiplexing completed"
;;

# ============================================================
### ASV INFERENCE WITH QIIME2 / DADA2
# ============================================================
asv)

set -euo pipefail
shopt -s nullglob

make_dirs

log "Starting ASV inference with QIIME2/DADA2"

for MARKER in "${MARKERS[@]}"; do

    MANIFEST="$OUTPUT_DIR/qiime/manifests/${MARKER}_manifest.csv"
    IMPORTED="$OUTPUT_DIR/qiime/imported/${MARKER}_demux.qza"
    TABLE="$OUTPUT_DIR/qiime/dada2/${MARKER}_table.qza"
    REPSEQS="$OUTPUT_DIR/qiime/dada2/${MARKER}_rep-seqs.qza"
    STATS="$OUTPUT_DIR/qiime/dada2/${MARKER}_stats.qza"
    STATS_QZV="$OUTPUT_DIR/qiime/dada2/${MARKER}_stats.qzv"

    log "Preparing QIIME2 manifest for $MARKER"

    echo "sample-id,forward-absolute-filepath,reverse-absolute-filepath" > "$MANIFEST"

    while read -r SAMPLE; do

        R1="$OUTPUT_DIR/demux/$MARKER/${SAMPLE}_${MARKER}_R1.fastq.gz"
        R2="$OUTPUT_DIR/demux/$MARKER/${SAMPLE}_${MARKER}_R2.fastq.gz"

        if [[ -s "$R1" && -s "$R2" ]]; then
            echo "$SAMPLE,$(realpath "$R1"),$(realpath "$R2")" >> "$MANIFEST"
        else
            log "WARNING: Missing demux reads for $SAMPLE $MARKER; skipping in manifest"
        fi

    done < <(sample_ids)

    if [[ -s "$TABLE" && -s "$REPSEQS" && "$FORCE" != "true" ]]; then
        log "Skipping DADA2 for $MARKER: outputs already exist"
        continue
    fi

    log "Importing $MARKER reads into QIIME2"

    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$MANIFEST" \
        --output-path "$IMPORTED" \
        --input-format PairedEndFastqManifestPhred33V2 \
        &>> "$OUTPUT_DIR/logs/qiime/${MARKER}.import.log"

    log "Running DADA2 for $MARKER"

    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "$IMPORTED" \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --p-n-threads "$THREADS" \
        --o-table "$TABLE" \
        --o-representative-sequences "$REPSEQS" \
        --o-denoising-stats "$STATS" \
        &>> "$OUTPUT_DIR/logs/qiime/${MARKER}.dada2.log"

    qiime metadata tabulate \
        --m-input-file "$STATS" \
        --o-visualization "$STATS_QZV" \
        &>> "$OUTPUT_DIR/logs/qiime/${MARKER}.stats.log"

    mkdir -p "$OUTPUT_DIR/qiime/exported/$MARKER"

    qiime tools export \
        --input-path "$REPSEQS" \
        --output-path "$OUTPUT_DIR/qiime/exported/$MARKER" \
        &>> "$OUTPUT_DIR/logs/qiime/${MARKER}.export_repseqs.log"

done

log "ASV inference completed"
;;

# ============================================================
### OTU CLUSTERING WITH VSEARCH
# ============================================================
otu)

set -euo pipefail
shopt -s nullglob

make_dirs

log "Starting OTU workflow with VSEARCH"

for MARKER in "${MARKERS[@]}"; do

    log "Processing marker: $MARKER"

    mkdir -p "$OUTPUT_DIR/otu/$MARKER/merged"
    mkdir -p "$OUTPUT_DIR/otu/$MARKER/fasta"

    COMBINED_FASTA="$OUTPUT_DIR/otu/$MARKER/${MARKER}_all_samples.fasta"
    UNIQUES="$OUTPUT_DIR/otu/$MARKER/${MARKER}_uniques.fasta"
    OTUS="$OUTPUT_DIR/otu/$MARKER/${MARKER}_otus97.fasta"
    UC="$OUTPUT_DIR/otu/$MARKER/${MARKER}_clusters.uc"
    OTU_TABLE="$OUTPUT_DIR/otu/$MARKER/${MARKER}_otu_table.tsv"

    : > "$COMBINED_FASTA"

    while read -r SAMPLE; do

        R1="$OUTPUT_DIR/demux/$MARKER/${SAMPLE}_${MARKER}_R1.fastq.gz"
        R2="$OUTPUT_DIR/demux/$MARKER/${SAMPLE}_${MARKER}_R2.fastq.gz"

        if [[ ! -s "$R1" || ! -s "$R2" ]]; then
            log "WARNING: Missing demux reads for $SAMPLE $MARKER; skipping"
            continue
        fi

        MERGED="$OUTPUT_DIR/otu/$MARKER/merged/${SAMPLE}_${MARKER}_merged.fastq"
        FASTA="$OUTPUT_DIR/otu/$MARKER/fasta/${SAMPLE}_${MARKER}.fasta"

        log "$SAMPLE $MARKER: merging paired reads"

        vsearch \
            --fastq_mergepairs "$R1" \
            --reverse "$R2" \
            --fastqout "$MERGED" \
            --threads "$THREADS" \
            &> "$OUTPUT_DIR/logs/vsearch/${SAMPLE}.${MARKER}.merge.log"

        log "$SAMPLE $MARKER: converting merged reads to FASTA"

        vsearch \
            --fastq_filter "$MERGED" \
            --fastaout "$FASTA" \
            --fastq_maxee 1.0 \
            --threads "$THREADS" \
            &> "$OUTPUT_DIR/logs/vsearch/${SAMPLE}.${MARKER}.filter.log"

        awk -v s="$SAMPLE" '
            /^>/ {
                sub(/^>/, ">" s "_")
                print $0 ";sample=" s ";"
                next
            }
            {print}
        ' "$FASTA" >> "$COMBINED_FASTA"

    done < <(sample_ids)

    log "$MARKER: dereplicating full-length sequences"

    vsearch \
        --derep_fulllength "$COMBINED_FASTA" \
        --sizeout \
        --output "$UNIQUES" \
        --threads "$THREADS" \
        &> "$OUTPUT_DIR/logs/vsearch/${MARKER}.derep.log"

    log "$MARKER: clustering OTUs at ${OTU_ID} identity"

    vsearch \
        --cluster_fast "$UNIQUES" \
        --id "$OTU_ID" \
        --centroids "$OTUS" \
        --uc "$UC" \
        --relabel "${MARKER}_OTU_" \
        --threads "$THREADS" \
        &> "$OUTPUT_DIR/logs/vsearch/${MARKER}.cluster.log"

    log "$MARKER: generating OTU table"

    vsearch \
        --usearch_global "$COMBINED_FASTA" \
        --db "$OTUS" \
        --id "$OTU_ID" \
        --otutabout "$OTU_TABLE" \
        --threads "$THREADS" \
        &> "$OUTPUT_DIR/logs/vsearch/${MARKER}.otutable.log"

done

log "OTU workflow completed"
;;

# ============================================================
### TAXONOMY FOR ASVs
# ============================================================
taxonomy_asv)

set -euo pipefail
shopt -s nullglob

make_dirs

log "Starting taxonomic assignment for ASVs"

if [[ -n "$CUSTOM_DB_FASTA" ]]; then
    log "Building BLAST database from $CUSTOM_DB_FASTA"
    mkdir -p "$(dirname "$BLAST_DB")"
    makeblastdb \
        -in "$CUSTOM_DB_FASTA" \
        -dbtype nucl \
        -parse_seqids \
        -out "$BLAST_DB" \
        &>> "$OUTPUT_DIR/logs/blast/makeblastdb.log"
fi

for MARKER in "${MARKERS[@]}"; do

    QUERY="$OUTPUT_DIR/qiime/exported/$MARKER/dna-sequences.fasta"

    if [[ ! -s "$QUERY" ]]; then
        log "WARNING: Missing ASV representative sequences for $MARKER: $QUERY"
        continue
    fi

    PREFIX="$OUTPUT_DIR/taxonomy/${MARKER}_ASV"
    BLAST_OUT="${PREFIX}_blast.tsv"
    BEST_OUT="${PREFIX}_best_hits.tsv"
    TAXIDS="${PREFIX}_taxids.txt"
    LINEAGES="${PREFIX}_lineages.tsv"
    FINAL="${PREFIX}_taxonomy.tsv"

    log "$MARKER ASV: running BLAST"

    blastn \
        -query "$QUERY" \
        -db "$BLAST_DB" \
        -perc_identity "$BLAST_MIN_ID" \
        -evalue "$BLAST_EVALUE" \
        -max_target_seqs "$BLAST_MAX_TARGETS" \
        -num_threads "$THREADS" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms stitle" \
        -out "$BLAST_OUT" \
        &> "$OUTPUT_DIR/logs/blast/${MARKER}.ASV.blast.log"

    log "$MARKER ASV: selecting best hits"

    sort -k1,1 -k12,12gr "$BLAST_OUT" \
        | awk -F'\t' '!seen[$1]++' \
        > "$BEST_OUT"

    if command -v taxonkit >/dev/null 2>&1; then

        log "$MARKER ASV: formatting lineages with TaxonKit"

        awk -F'\t' '{
            split($13, a, ";")
            print a[1]
        }' "$BEST_OUT" > "$TAXIDS"

        taxonkit lineage \
            --data-dir "$TAXONKIT_DB" \
            "$TAXIDS" \
            | taxonkit reformat \
                --data-dir "$TAXONKIT_DB" \
                -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" \
            > "$LINEAGES"

        paste "$BEST_OUT" "$LINEAGES" > "$FINAL"

    else
        cp "$BEST_OUT" "$FINAL"
    fi

done

log "ASV taxonomy completed"
;;

# ============================================================
### TAXONOMY FOR OTUs
# ============================================================
taxonomy_otu)

set -euo pipefail
shopt -s nullglob

make_dirs

log "Starting taxonomic assignment for OTUs"

if [[ -n "$CUSTOM_DB_FASTA" ]]; then
    log "Building BLAST database from $CUSTOM_DB_FASTA"
    mkdir -p "$(dirname "$BLAST_DB")"
    makeblastdb \
        -in "$CUSTOM_DB_FASTA" \
        -dbtype nucl \
        -parse_seqids \
        -out "$BLAST_DB" \
        &>> "$OUTPUT_DIR/logs/blast/makeblastdb.log"
fi

for MARKER in "${MARKERS[@]}"; do

    QUERY="$OUTPUT_DIR/otu/$MARKER/${MARKER}_otus97.fasta"

    if [[ ! -s "$QUERY" ]]; then
        log "WARNING: Missing OTU representative sequences for $MARKER: $QUERY"
        continue
    fi

    PREFIX="$OUTPUT_DIR/taxonomy/${MARKER}_OTU97"
    BLAST_OUT="${PREFIX}_blast.tsv"
    BEST_OUT="${PREFIX}_best_hits.tsv"
    TAXIDS="${PREFIX}_taxids.txt"
    LINEAGES="${PREFIX}_lineages.tsv"
    FINAL="${PREFIX}_taxonomy.tsv"

    log "$MARKER OTU: running BLAST"

    blastn \
        -query "$QUERY" \
        -db "$BLAST_DB" \
        -perc_identity "$BLAST_MIN_ID" \
        -evalue "$BLAST_EVALUE" \
        -max_target_seqs "$BLAST_MAX_TARGETS" \
        -num_threads "$THREADS" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms stitle" \
        -out "$BLAST_OUT" \
        &> "$OUTPUT_DIR/logs/blast/${MARKER}.OTU97.blast.log"

    log "$MARKER OTU: selecting best hits"

    sort -k1,1 -k12,12gr "$BLAST_OUT" \
        | awk -F'\t' '!seen[$1]++' \
        > "$BEST_OUT"

    if command -v taxonkit >/dev/null 2>&1; then

        log "$MARKER OTU: formatting lineages with TaxonKit"

        awk -F'\t' '{
            split($13, a, ";")
            print a[1]
        }' "$BEST_OUT" > "$TAXIDS"

        taxonkit lineage \
            --data-dir "$TAXONKIT_DB" \
            "$TAXIDS" \
            | taxonkit reformat \
                --data-dir "$TAXONKIT_DB" \
                -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" \
            > "$LINEAGES"

        paste "$BEST_OUT" "$LINEAGES" > "$FINAL"

    else
        cp "$BEST_OUT" "$FINAL"
    fi

done

log "OTU taxonomy completed"
;;

# ============================================================
### MULTIQC
# ============================================================
multiqc)

set -euo pipefail

make_dirs

log "Running MultiQC"

multiqc \
    "$OUTPUT_DIR/qc" \
    "$OUTPUT_DIR/fastp" \
    "$OUTPUT_DIR/logs" \
    -o "$OUTPUT_DIR/multiqc" \
    &>> "$OUTPUT_DIR/logs/multiqc.log" || true

log "MultiQC completed"
;;

# ============================================================
### ALL
# ============================================================
all)

log "Running full pipeline"

"$0" -s check
"$0" -s concat
"$0" -s qc_raw
"$0" -s fastp
"$0" -s qc_filtered
"$0" -s demux

if [[ "$MODE" == "asv" || "$MODE" == "both" ]]; then
    "$0" -s asv
    "$0" -s taxonomy_asv
fi

if [[ "$MODE" == "otu" || "$MODE" == "both" ]]; then
    "$0" -s otu
    "$0" -s taxonomy_otu
fi

"$0" -s multiqc

log "Full pipeline completed"
;;

# ============================================================
### UNKNOWN STEP
# ============================================================
*)

log "ERROR: Unknown step: $STEP"
usage
exit 1
;;

esac

log "Finished step: $STEP"

done

log "Requested steps completed successfully"
