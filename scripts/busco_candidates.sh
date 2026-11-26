#!/usr/bin/env bash
# BUSCO (protein mode, metazoa_odb10, as NCBI reports for RefSeq annotations) for candidate
# genomes whose NCBI annotation has no BUSCO score. One protein per gene (longest CDS).
# TSA project accessions (e.g. GCWQ01000000) are scored on the proteins of the longest ORF
# per contig, as `cpgme run` calls them (contigs are not grouped into genes, so duplicated
# BUSCOs run higher than for genomes).
#
# Usage:  scripts/busco_candidates.sh [ACCESSION ...]
#         with no arguments: every include=check row in manifest/candidates.tsv
# Env:    BUSCO_SIF (container image), CONTAINER (singularity|apptainer), CPU (per BUSCO run), PARALLEL (concurrent runs),
#         WORK (scratch, default data/busco)
# Writes  manifest/busco.tsv rows (accession, lineage, version, complete, single, duplicated,
#         fragmented, missing, n), which `cpgme candidates` uses to fill the BUSCO columns.
set -euo pipefail
export PATH=$HOME/.local/bin:$PATH
cd "$(dirname "$0")/.."

BUSCO_SIF=${BUSCO_SIF:-/home/shared/16TB_HDD_01/sr320/containers/busco_v5.8.2.sif}
CONTAINER=${CONTAINER:-singularity}  # apptainer fails on Raven (user namespaces)
CPU=${CPU:-8}
PARALLEL=${PARALLEL:-5}
WORK=${WORK:-data/busco}
OUT=manifest/busco.tsv
mkdir -p "$WORK"

if [ $# -gt 0 ]; then
  accs=("$@")
else
  mapfile -t accs < <(grep -v '^#' manifest/candidates.tsv |
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["include"]=="check"{print $h["accession"]}')
fi
[ -f "$OUT" ] || printf 'accession\tlineage\tversion\tcomplete\tsingle\tduplicated\tfragmented\tmissing\tn\n' > "$OUT"

one() {
  acc=$1
  if grep -q "^$acc	" "$OUT"; then echo "$acc: already in $OUT"; return 0; fi
  faa="$WORK/$acc.faa"
  uv run python - "$acc" "$faa" <<'EOF'
import gzip, shutil, sys
from pathlib import Path
from cpgme.fetch import fetch
from cpgme.annotate import write_gene_proteins
from cpgme.tsa import fetch_tsa, is_tsa, tsa_cds, write_proteins
acc, faa = sys.argv[1], Path(sys.argv[2])
gz = faa.with_suffix(".faa.gz")
if is_tsa(acc):
    f = fetch_tsa(acc, Path("data") / acc)
    n = len(write_proteins(tsa_cds(f["TSA_FASTA"]), gz))
else:
    f = fetch(acc, Path("data") / acc)
    n = len(write_gene_proteins(f["CDS_FASTA"], f["PROT_FASTA"], gz, f["GENOME_GFF"]))
with gzip.open(gz, "rb") as s, open(faa, "wb") as d:
    shutil.copyfileobj(s, d)
gz.unlink()
print(f"{acc}: {n} proteins")
EOF
  "$CONTAINER" exec --bind "$PWD" --pwd "$PWD/$WORK" "$BUSCO_SIF" \
    busco -i "$acc.faa" -m proteins -l metazoa_odb10 -c "$CPU" -o "$acc" -f \
          --download_path "$PWD/$WORK/busco_downloads" > "$WORK/$acc.log" 2>&1
  js=$(ls "$WORK/$acc"/short_summary*.json)
  python3 - "$acc" "$js" >> "$OUT" <<'EOF'
import json, sys
acc, js = sys.argv[1], json.load(open(sys.argv[2]))
r, n = js["results"], js["results"]["n_markers"]
print("\t".join(map(str, [acc, js["lineage_dataset"]["name"], js["versions"]["busco"],
      round(r["Complete percentage"] / 100, 4), round(r["Single copy percentage"] / 100, 4),
      round(r["Multi copy percentage"] / 100, 4), round(r["Fragmented percentage"] / 100, 4),
      round(r["Missing percentage"] / 100, 4), n])))
EOF
  rm -f "$faa"
  echo "$acc: $(tail -1 "$OUT")"
}
export -f one
export OUT WORK CPU BUSCO_SIF CONTAINER

# Download the lineage once before running in parallel.
"$CONTAINER" exec --bind "$PWD" "$BUSCO_SIF" busco --download_path "$PWD/$WORK/busco_downloads" \
  --download metazoa_odb10 > /dev/null 2>&1 || true
printf '%s\n' "${accs[@]}" | xargs -P "$PARALLEL" -I{} bash -c 'one {}'
echo BUSCO_DONE
