#!/usr/bin/env bash
# eggNOG-mapper + GO Slim (`cpgme annotate`) for every manifest row that has results but no
# goslim_test.json yet, one species at a time. Safe to rerun: finished species are skipped.
#
# Usage:  scripts/annotate_all.sh
# Env:    EMAPPER_DATA (eggNOG DB dir), CPU (default 40)
# Logs:   logs/annotate/<accession>.log; progress lines in the caller's stdout.
set -uo pipefail
export PATH=$HOME/.local/bin:$PATH
export PYTHONUNBUFFERED=1
cd "$(dirname "$0")/.."

EMAPPER_DATA=${EMAPPER_DATA:-/home/shared/16TB_HDD_01/sr320/eggnog-data/emapperdb-5.0.2}
CPU=${CPU:-40}
mkdir -p logs/annotate

mapfile -t accs < <(grep -v '^#' manifest/species.tsv |
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} {print $h["accession"]}')
todo=()
for a in "${accs[@]}"; do
  [ -f "results/$a/summary.json" ] && [ ! -f "results/$a/goslim_test.json" ] && todo+=("$a")
done
echo "[$(date '+%F %T')] ${#todo[@]} species to annotate"

failed=0
for a in "${todo[@]}"; do
  echo "[$(date '+%F %T')] $a start"
  if uv run cpgme annotate "$a" --emapper-data "$EMAPPER_DATA" --cpu "$CPU" \
       > "logs/annotate/$a.log" 2>&1; then
    echo "[$(date '+%F %T')] $a done"
  else
    echo "[$(date '+%F %T')] $a FAILED"; failed=$((failed + 1))
  fi
done
echo "ANNOTATE_DONE ($failed failed)"
