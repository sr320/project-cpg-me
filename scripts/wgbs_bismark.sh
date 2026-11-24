#!/bin/bash
# WGBS validation for one assembly: SRA runs -> fastp -> Bismark -> dedup -> CpG coverage
# -> `cpgme validate`. Written for Raven; tool paths are overridable.
#
# Usage: scripts/wgbs_bismark.sh <assembly accession> <label> "<description>" <SRR> [SRR ...]
#   Runs from one library are aligned separately and deduplicated together.
#   Needs `cpgme run <accession> --genome` first (uses data/<accession>/ genome and GFF).
#
# Environment (defaults for Raven):
#   WGBS_ROOT  work directory (large; FASTQ and BAMs)   PARALLEL  Bismark instances
#   NON_DIRECTIONAL=1  align all four strands (PBAT-like libraries, e.g. Zymo Pico Methyl-Seq).
#     Check first: a directional run that maps at about half the expected rate is the sign.
set -euo pipefail

ACC=$1; LABEL=$2; DESC=$3; shift 3; RUNS=("$@")
REPO=$(cd "$(dirname "$0")/.." && pwd)
WGBS_ROOT=${WGBS_ROOT:-/home/shared/16TB_HDD_01/sr320/cpgme-wgbs}
PARALLEL=${PARALLEL:-6}
THREADS=${THREADS:-16}
NON_DIRECTIONAL=${NON_DIRECTIONAL:-0}
DIR_OPT=(); MODE=directional
[ "$NON_DIRECTIONAL" = 1 ] && { DIR_OPT=(--non_directional); MODE=non-directional; }
export PATH=/home/shared/Bismark-0.24.0:/home/shared/bowtie2-2.4.4-linux-x86_64:/home/shared/sratoolkit.3.0.2-ubuntu64/bin:/home/shared/fastp-v0.24.0:/home/shared/pigz-2.8/bin:$HOME/.local/bin:$PATH

GDIR=$WGBS_ROOT/$ACC/genome
WORK=$WGBS_ROOT/$ACC/$LABEL
mkdir -p "$GDIR" "$WORK"/{sra,fastq,trimmed,aln,meth}
log() { echo "[$(date '+%F %T')] $*"; }

for t in prefetch fasterq-dump fastp bismark bowtie2 samtools deduplicate_bismark bismark_methylation_extractor; do
  command -v $t >/dev/null || { echo "missing tool: $t" >&2; exit 1; }
done
log "tools: $(bismark --version | grep -o 'v[0-9.]*' | head -1), bowtie2 $(bowtie2 --version | head -1 | awk '{print $NF}'), fastp $(fastp --version 2>&1 | awk '{print $2}')"

# 1. Bismark genome (once per assembly)
if [ ! -d "$GDIR/Bisulfite_Genome" ]; then
  FNA=$(ls "$REPO"/data/"$ACC"/*_genomic.fna | grep -v cds_from | head -1)
  [ -s "$FNA" ] || { echo "no genome FASTA in $REPO/data/$ACC; run cpgme run $ACC --genome" >&2; exit 1; }
  ln -sf "$FNA" "$GDIR/genome.fa"
  log "preparing Bismark genome from $FNA"
  bismark_genome_preparation --parallel 8 "$GDIR" > "$WORK/genome_preparation.log" 2>&1
fi

# 2. Per run: download, trim, align
BAMS=()
for R in "${RUNS[@]}"; do
  BAM=$WORK/aln/${R}_1_bismark_bt2_pe.bam
  if [ ! -s "$BAM" ]; then
    log "$R: download"
    prefetch --max-size 100G -O "$WORK/sra" "$R" > /dev/null
    fasterq-dump --split-files -e 8 -t "$WORK/fastq" -O "$WORK/fastq" "$WORK/sra/$R" > /dev/null
    log "$R: fastp"
    fastp -i "$WORK/fastq/${R}_1.fastq" -I "$WORK/fastq/${R}_2.fastq" \
          -o "$WORK/trimmed/${R}_1.fq.gz" -O "$WORK/trimmed/${R}_2.fq.gz" \
          --detect_adapter_for_pe --trim_poly_g -w 16 \
          -j "$WORK/trimmed/${R}.fastp.json" -h "$WORK/trimmed/${R}.fastp.html" 2> "$WORK/trimmed/${R}.fastp.log"
    rm -f "$WORK/fastq/${R}"_*.fastq
    log "$R: bismark $MODE (parallel $PARALLEL)"
    bismark --genome "$GDIR" "${DIR_OPT[@]}" -1 "$WORK/trimmed/${R}_1.fq.gz" -2 "$WORK/trimmed/${R}_2.fq.gz" \
            --parallel "$PARALLEL" -p 2 --temp_dir "$WORK/aln" -o "$WORK/aln" \
            > "$WORK/aln/${R}.bismark.log" 2>&1
    rm -f "$WORK/trimmed/${R}"_*.fq.gz "$WORK/sra/$R"/*.sra
  fi
  grep -h "Mapping efficiency" "$WORK/aln/${R}_1_bismark_bt2_PE_report.txt" | sed "s/^/  $R: /" || true
  BAMS+=("$BAM")
done

# 3. Deduplicate (all runs of the library together) and extract CpG methylation
DEDUP=$WORK/aln/$LABEL.deduplicated.bam
if [ ! -s "$DEDUP" ]; then
  log "deduplicate"
  MULTI=(); [ ${#BAMS[@]} -gt 1 ] && MULTI=(--multiple)
  deduplicate_bismark --paired --bam "${MULTI[@]}" --outfile "$LABEL" --output_dir "$WORK/aln" \
      "${BAMS[@]}" > "$WORK/aln/dedup.log" 2>&1
  # output naming differs between single and --multiple; take the newest deduplicated BAM
  NEW=$(ls -t "$WORK"/aln/*deduplicated.bam | head -1)
  [ "$NEW" = "$DEDUP" ] || mv "$NEW" "$DEDUP"
fi
COV=$WORK/meth/$LABEL.deduplicated.bismark.cov.gz
if [ ! -s "$COV" ]; then
  log "methylation extractor"
  bismark_methylation_extractor --paired-end --no_overlap --comprehensive --bedGraph --gzip \
      --ignore_r2 2 --parallel 8 --buffer_size 20G -o "$WORK/meth" "$DEDUP" > "$WORK/meth/extractor.log" 2>&1
  rm -f "$WORK"/meth/C??_context_*.txt.gz   # per-read calls; the coverage file is what we keep
fi

# 4. Compare with CpG O/E
log "cpgme validate"
cd "$REPO"
uv run cpgme validate "$ACC" --cov "$COV" --label "$LABEL" --source "$DESC; runs ${RUNS[*]}; Bismark $MODE"
mkdir -p "results/$ACC/validation/$LABEL/reports"
cp "$WORK"/aln/*_PE_report.txt "$WORK"/aln/*deduplication_report.txt "$WORK"/meth/*splitting_report.txt \
   "$WORK"/meth/*M-bias.txt "$WORK"/trimmed/*.fastp.json "results/$ACC/validation/$LABEL/reports/" 2>/dev/null || true
log "WGBS_DONE $ACC $LABEL"
