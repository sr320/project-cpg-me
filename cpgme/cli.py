"""Command line: cpgme run | run-shard | status | annotate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

from .run import is_done, run

DEFAULT_MANIFEST = Path("manifest/species.tsv")


def read_manifest(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader(
            (line for line in fh if not line.startswith("#")), delimiter="\t")]
    return [r for r in rows if r.get("accession")]


def in_shard(accession: str, k: int, n: int) -> bool:
    """Stable assignment of an accession to shard k of n (1-based)."""
    h = int(hashlib.md5(accession.encode()).hexdigest(), 16)
    return h % n == k - 1


def _genome_wanted(row: dict) -> bool:
    return row.get("genome_features", "").strip().lower() in {"1", "yes", "true"}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="cpgme", description=__doc__)
    ap.add_argument("--results", type=Path, default=Path("results"))
    ap.add_argument("--data", type=Path, default=Path("data"))
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("run", help="process one accession")
    p.add_argument("accession")
    p.add_argument("--genome", action="store_true",
                   help="also download the genome and compute gene-body/upstream features")
    p.add_argument("--force", action="store_true")

    p = sub.add_parser("run-shard", help="process this machine's share of the manifest")
    p.add_argument("--shard", required=True, help="k/N, e.g. 2/5")
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--force", action="store_true")

    p = sub.add_parser("status", help="compare manifest with results")
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)

    p = sub.add_parser("annotate", help="eggNOG/GO Slim annotation and the 2015 GO Slim test "
                                        "(after `run`)")
    p.add_argument("accession")
    p.add_argument("--emapper-annotations", type=Path,
                   help="existing *.emapper.annotations file (skip running eggNOG-mapper)")
    p.add_argument("--emapper-data", type=Path, help="eggNOG-mapper --data_dir")
    p.add_argument("--cpu", type=int, default=4)
    p.add_argument("--proteins-only", action="store_true",
                   help="only write results/<acc>/proteins.faa.gz (to annotate elsewhere)")

    args = ap.parse_args(argv)

    if args.cmd == "annotate":
        from .annotate import annotate
        annotate(args.accession, args.results, args.data,
                 emapper_annotations=args.emapper_annotations,
                 emapper_data=args.emapper_data, cpu=args.cpu,
                 proteins_only=args.proteins_only)
        return 0

    if args.cmd == "run":
        run(args.accession, args.results, args.data, genome=args.genome, force=args.force)
        return 0

    rows = read_manifest(args.manifest)

    if args.cmd == "run-shard":
        k, n = (int(x) for x in args.shard.split("/"))
        if not 1 <= k <= n:
            ap.error("--shard must be k/N with 1 <= k <= N")
        mine = [r for r in rows if in_shard(r["accession"], k, n)]
        print(f"shard {k}/{n}: {len(mine)} of {len(rows)} accessions")
        failed = []
        for r in mine:
            try:
                run(r["accession"], args.results, args.data,
                    genome=_genome_wanted(r), force=args.force)
            except Exception as exc:  # keep going; report at the end
                print(f"{r['accession']}: FAILED: {exc}", file=sys.stderr)
                failed.append(r["accession"])
        if failed:
            print(f"{len(failed)} failed: {' '.join(failed)}", file=sys.stderr)
            return 1
        return 0

    if args.cmd == "status":
        done = 0
        for r in rows:
            ok = is_done(args.results / r["accession"], _genome_wanted(r))
            done += ok
            annotated = (args.results / r["accession"] / "goslim_test.json").exists()
            print(f"{'done   ' if ok else 'pending'}  {'annotated' if annotated else '         '}  "
                  f"{r['accession']:<18} {r.get('group', ''):<10} {r.get('species', '')}")
        print(f"{done}/{len(rows)} done")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
