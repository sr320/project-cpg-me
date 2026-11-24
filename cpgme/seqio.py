"""Minimal FASTA and GFF3 readers (gzip-aware, streaming)."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterator
from urllib.parse import unquote


def _open(path: Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (header_without_>, sequence) one record at a time."""
    header = None
    chunks: list[str] = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].rstrip("\n")
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)


def parse_gff_attributes(field: str) -> dict[str, str]:
    attrs = {}
    for part in field.strip().split(";"):
        if not part:
            continue
        key, _, value = part.partition("=")
        attrs[key] = unquote(value)
    return attrs


def read_gff(path: Path, feature_types: set[str] | None = None) -> Iterator[dict]:
    """Yield GFF3 records as dicts. Coordinates stay 1-based, closed."""
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            if feature_types is not None and cols[2] not in feature_types:
                continue
            yield {
                "seqid": cols[0],
                "type": cols[2],
                "start": int(cols[3]),
                "end": int(cols[4]),
                "strand": cols[6],
                "attrs": parse_gff_attributes(cols[8]),
            }
