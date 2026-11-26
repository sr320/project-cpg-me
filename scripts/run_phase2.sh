#!/usr/bin/env bash
# Phase 2 on one machine: run N shards of the manifest in parallel, one log per shard.
# Usage: scripts/run_phase2.sh [N]   (default 8; more shards hit NCBI rate limits, HTTP 429)
set -uo pipefail
export PATH=$HOME/.local/bin:$PATH
export PYTHONUNBUFFERED=1   # logs show progress as it happens
# Cap numpy/BLAS threads: N shards x default threads oversubscribes a shared machine.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2} OPENBLAS_NUM_THREADS=${OMP_NUM_THREADS:-2} MKL_NUM_THREADS=${OMP_NUM_THREADS:-2}
cd "$(dirname "$0")/.."
N=${1:-8}
mkdir -p logs/phase2
for k in $(seq 1 "$N"); do
  uv run cpgme run-shard --shard "$k/$N" > "logs/phase2/shard_${k}of${N}.log" 2>&1 &
done
wait
echo PHASE2_DONE
uv run cpgme status | tail -1
