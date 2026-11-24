"""Download annotated genome packages from the NCBI Datasets v2 REST API.

No NCBI CLI is required: the package is a zip with a fixed layout, fetched with urllib.
Only the expected members are extracted, and each is checked against the package's
md5sum.txt.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import threading
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
USER_AGENT = "cpgme/0.1 (comparative CpG O/E pipeline)"

# Annotation types we may request -> filename inside ncbi_dataset/data/<acc>/
FILES = {
    "CDS_FASTA": "cds_from_genomic.fna",
    "GENOME_GFF": "genomic.gff",
    "PROT_FASTA": "protein.faa",
}


def _get(url: str, retries: int = 4) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=600) as resp:
                return resp.read()
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(5 * 2**attempt)
    raise RuntimeError("unreachable")


STALL_SECONDS = 120  # abort and retry a download that makes no progress for this long


def _copy_with_watchdog(resp, out) -> None:
    """Copy resp -> out; close the connection if no bytes arrive for STALL_SECONDS.

    A socket timeout alone did not catch a stalled NCBI transfer in the pilot (a download
    hung for 30+ minutes), so a separate thread watches progress and closes the socket.
    """
    last = [time.monotonic()]
    finished = threading.Event()
    stalled = threading.Event()

    def watch():
        while not finished.wait(5):
            if time.monotonic() - last[0] > STALL_SECONDS:
                stalled.set()
                try:
                    resp.close()
                finally:
                    return

    t = threading.Thread(target=watch, daemon=True)
    t.start()
    try:
        while True:
            block = resp.read(1 << 20)
            if not block:
                break
            out.write(block)
            last[0] = time.monotonic()
    except Exception:
        if stalled.is_set():
            raise TimeoutError(f"download stalled for {STALL_SECONDS}s") from None
        raise
    finally:
        finished.set()
    if stalled.is_set():
        raise TimeoutError(f"download stalled for {STALL_SECONDS}s")


def _download_to(url: str, dest: Path, retries: int = 4) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp, open(tmp, "wb") as out:
                _copy_with_watchdog(resp, out)
            tmp.replace(dest)
            return
        except Exception as exc:
            if attempt == retries - 1:
                raise
            print(f"  download attempt {attempt + 1} failed ({exc}); retrying", flush=True)
            time.sleep(5 * 2**attempt)


def assembly_report(accession: str) -> dict:
    url = f"{API}/genome/accession/{urllib.parse.quote(accession)}/dataset_report"
    data = json.loads(_get(url))
    reports = data.get("reports", [])
    if not reports:
        raise ValueError(f"No NCBI assembly report for {accession}")
    return reports[0]


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch(accession: str, outdir: Path, genome: bool = False) -> dict[str, Path]:
    """Download CDS, GFF and protein FASTA (plus genome FASTA if requested).

    Returns {annotation_type: local path}. Re-uses files already present and verified.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    types = list(FILES) + (["GENOME_FASTA"] if genome else [])

    done_marker = outdir / f"fetched_{'genome' if genome else 'annot'}.json"
    if done_marker.exists():
        paths = {k: outdir / v for k, v in json.loads(done_marker.read_text()).items()}
        if all(p.exists() for p in paths.values()):
            return paths

    query = urllib.parse.urlencode([("include_annotation_type", t) for t in types])
    url = f"{API}/genome/accession/{urllib.parse.quote(accession)}/download?{query}"
    zpath = outdir / "package.zip"
    _download_to(url, zpath)

    prefix = f"ncbi_dataset/data/{accession}/"
    paths: dict[str, Path] = {}
    with zipfile.ZipFile(zpath) as zf:
        names = set(zf.namelist())
        checksums = {}
        if "md5sum.txt" in names:
            for line in zf.read("md5sum.txt").decode().splitlines():
                digest, _, member = line.partition("  ")
                checksums[member.strip()] = digest.strip()

        wanted: dict[str, str] = {}
        for t, fname in FILES.items():
            if prefix + fname in names:
                wanted[t] = prefix + fname
        if genome:
            fasta = [n for n in names if n.startswith(prefix) and n.endswith("_genomic.fna")
                     and "cds_from" not in n and "rna" not in n]
            if len(fasta) != 1:
                raise ValueError(f"Expected one genome FASTA in package, found {fasta}")
            wanted["GENOME_FASTA"] = fasta[0]
        missing = [t for t in types if t not in wanted]
        if "CDS_FASTA" in missing or "GENOME_GFF" in missing:
            raise ValueError(f"{accession}: package lacks {missing}; is it annotated?")

        for t, member in wanted.items():
            # Write to a fixed local name; never trust paths inside the archive.
            dest = outdir / Path(member).name
            with zf.open(member) as src, open(dest, "wb") as out:
                shutil.copyfileobj(src, out, length=1 << 20)
            expected = checksums.get(member)
            if expected and _md5(dest) != expected:
                raise IOError(f"md5 mismatch for {member}")
            paths[t] = dest

    zpath.unlink()
    done_marker.write_text(json.dumps({k: p.name for k, p in paths.items()}))
    return paths
