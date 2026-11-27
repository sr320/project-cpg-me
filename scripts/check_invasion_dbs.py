"""Check invasion status of manifest species against GISD and WRiMS, and relabel manifest/traits.tsv.

Usage:  python scripts/check_invasion_dbs.py                  # query both databases, write evidence + labels
        python scripts/check_invasion_dbs.py --from-json X    # reuse a saved query result (JSON)

Sources
  GISD   IUCN Global Invasive Species Database, https://www.iucngisd.org/gisd/speciesname/<name>
         "account" = a species account with a summary; "stub" = listed with no account text.
  WRiMS  World Register of Introduced Marine Species, via the WoRMS REST API
         (AphiaDistributionsByAphiaID): per-locality establishmentMeans and invasiveness.
         Covers marine/brackish taxa; its few insect records are from sub-Antarctic islands.

Rule (applied in that order)
  invasive    GISD account describing the species as an alien, or >= 1 WRiMS locality with
              invasiveness "Invasive" or "Of concern"
  introduced  WRiMS alien records but none flagged invasive, or a GISD stub only
  native      marine/brackish taxon (WRiMS scope) with no alien records and not in GISD
  (not covered by either database, including terrestrial and freshwater taxa WoRMS lists: the earlier general-knowledge label is kept and marked so)
Overrides (OVERRIDE below) record cases where a database hit does not mean "alien".

Writes manifest/invasion_db.tsv (evidence per species) and rewrites the invasion/source/basis
columns of manifest/traits.tsv.
"""
from __future__ import annotations

import argparse, csv, json, re, time, urllib.error, urllib.parse, urllib.request
from pathlib import Path

SYN = {  # names a database may file the species under
    "Magallana gigas": ["Crassostrea gigas"], "Magallana angulata": ["Crassostrea angulata"],
    "Dreissena rostriformis bugensis": ["Dreissena bugensis", "Dreissena rostriformis"],
    "Mizuhopecten yessoensis": ["Patinopecten yessoensis"],
    "Sinanodonta woodiana": ["Anodonta woodiana"], "Sinohyriopsis cumingii": ["Hyriopsis cumingii"],
    "Ylistrum balloti": ["Amusium balloti"], "Desmophyllum pertusum": ["Lophelia pertusa"],
    "Pinctada imbricata": ["Pinctada radiata", "Pinctada imbricata radiata"],
    "Ruditapes philippinarum": ["Venerupis philippinarum", "Tapes philippinarum"],
}
OVERRIDE = {
    "Acanthaster planci": ("native", "GISD account is for population outbreaks within the native "
                           "Indo-Pacific range, not an alien introduction; WRiMS: 2 alien records, none invasive"),
    "Pinctada imbricata": ("introduced", "WRiMS: no alien records as P. imbricata; the Lessepsian invader "
                           "P. radiata (WRiMS invasive in Egypt, Tunisia) may or may not be the sequenced lineage"),
}
UA = {"User-Agent": "cpgme/0.1 (research script)"}


def get(url: str) -> tuple[int | None, str]:
    for i in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return r.status, r.read().decode("utf-8", "ignore")
        except urllib.error.HTTPError as e:
            if e.code in (204, 404):
                return e.code, ""
        except Exception:
            pass
        time.sleep(2 * (i + 1))
    return None, ""


def gisd(name: str) -> tuple[str, str]:
    st, html = get("https://www.iucngisd.org/gisd/speciesname/" + urllib.parse.quote_plus(name))
    if st is None:
        return "error", ""
    if "not present yet in our archive" in html:
        return "absent", ""
    txt = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ",
                 re.sub(r"<script.*?</script>|<style.*?</style>", "", html, flags=re.S)))
    m = re.search(r"Summary (.{0,400})", txt)
    return ("account", m.group(1).strip()) if m and m.group(1).strip() else ("stub", "")


def wrims(name: str) -> dict:
    st, body = get("https://www.marinespecies.org/rest/AphiaRecordsByName/"
                   + urllib.parse.quote(name) + "?like=false&marine_only=false")
    if st is None:
        return {"wrims": "error"}
    if not body.strip():
        return {"wrims": "not_in_worms"}
    recs = json.loads(body)
    r = next((x for x in recs if x.get("status") == "accepted"), recs[0])
    aid = r.get("valid_AphiaID") or r["AphiaID"]
    _, body = get(f"https://www.marinespecies.org/rest/AphiaDistributionsByAphiaID/{aid}")
    d = json.loads(body) if body.strip() else []
    alien = [x for x in d if x.get("establishmentMeans") == "Alien"]
    inv = [x for x in alien if x.get("invasiveness") in ("Invasive", "Of concern")]
    return {"wrims": "in_worms", "aphia": aid, "n_alien": len(alien), "n_invasive": len(inv),
            "inv_places": sorted({x.get("locality") for x in inv}), "marine": habitat(aid)}


def habitat(aid: int) -> bool | None:
    """True if WoRMS flags the taxon marine or brackish, i.e. within WRiMS scope."""
    _, body = get(f"https://www.marinespecies.org/rest/AphiaRecordByAphiaID/{aid}")
    if not body.strip():
        return None
    r = json.loads(body)
    return bool(r.get("isMarine") or r.get("isBrackish"))


def query(species: list[str]) -> list[dict]:
    out = []
    for sp in species:
        names = [sp] + SYN.get(sp, [])
        g = ("absent", "", "")
        for n in names:
            st, summ = gisd(n)
            time.sleep(0.4)
            if st in ("account", "stub") or st == "error":
                g = (st, summ, n)
                if st != "error":
                    break
        w = {"wrims": "not_in_worms"}
        for n in names:  # first name WoRMS knows
            w = wrims(n)
            time.sleep(0.4)
            if w["wrims"] in ("in_worms", "error"):
                break
        out.append({"species": sp, "gisd": g[0], "gisd_summary": g[1], "gisd_name": g[2], **w})
        print(sp, g[0], w.get("wrims"), w.get("n_alien"), w.get("n_invasive"), flush=True)
    return out


def label(e: dict, prior: str) -> tuple[str, str]:
    """(invasion, source) from the database evidence."""
    if e["species"] in OVERRIDE:
        return OVERRIDE[e["species"]][0], "override"
    g, w = e["gisd"], e.get("wrims")
    n_alien, n_inv = e.get("n_alien") or 0, e.get("n_invasive") or 0
    covered = g in ("account", "stub") or w == "in_worms"
    src = "+".join(x for x, ok in (("GISD", g in ("account", "stub")), ("WRiMS", n_alien > 0)) if ok)
    if g == "account" or n_inv > 0:
        return "invasive", src
    if n_alien > 0 or g == "stub":
        return "introduced", src
    if w == "in_worms" and e.get("marine"):
        return "native", "WRiMS (no alien records)"
    return prior, "general knowledge (not in GISD; outside WRiMS scope)"


def evidence_text(e: dict) -> str:
    parts = []
    if e["gisd"] == "account":
        parts.append("GISD account" + (f" (as {e['gisd_name']})" if e["gisd_name"] != e["species"] else ""))
    elif e["gisd"] == "stub":
        parts.append("GISD listed, no account")
    if e.get("wrims") == "in_worms" and e.get("n_alien"):
        s = f"WRiMS {e['n_alien']} alien localities, {e['n_invasive']} invasive"
        if e["n_invasive"]:
            s += " (" + "; ".join(e["inv_places"][:3]) + ("; ..." if len(e["inv_places"]) > 3 else "") + ")"
        parts.append(s)
    return "; ".join(parts)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--traits", type=Path, default=Path("manifest/traits.tsv"))
    ap.add_argument("--evidence", type=Path, default=Path("manifest/invasion_db.tsv"))
    ap.add_argument("--from-json", type=Path, help="reuse saved query results instead of querying")
    a = ap.parse_args()
    lines = a.traits.read_text().splitlines()
    header = [l for l in lines if l.startswith("#")]
    rows = list(csv.DictReader([l for l in lines if not l.startswith("#")], delimiter="\t"))
    ev = json.loads(a.from_json.read_text()) if a.from_json else query([r["species"] for r in rows])
    ev = {e["species"]: e for e in ev}

    with open(a.evidence, "w", newline="") as fh:
        fh.write(f"# GISD and WRiMS evidence per species, written by scripts/check_invasion_dbs.py "
                 f"on {time.strftime('%Y-%m-%d')}\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "gisd", "gisd_name", "worms", "aphia_id", "marine_or_brackish", "wrims_alien_localities",
                    "wrims_invasive_localities", "wrims_invasive_places"])
        for r in rows:
            e = ev[r["species"]]
            w.writerow([r["species"], e["gisd"], e.get("gisd_name") or "", e.get("wrims"), e.get("aphia", ""),
                        {True: "yes", False: "no"}.get(e.get("marine"), ""), e.get("n_alien", ""), e.get("n_invasive", ""), "; ".join(e.get("inv_places", []))])

    changes = []
    for r in rows:
        e = ev[r["species"]]
        prior = r.get("label_general_knowledge") or r["invasion"]
        new, src = label(e, prior)
        note = OVERRIDE[r["species"]][1] if r["species"] in OVERRIDE else evidence_text(e)
        old_basis = r.get("notes") if "notes" in r else r.get("basis", "")
        r.update({"invasion": new, "source": src, "evidence": note,
                  "label_general_knowledge": prior, "notes": old_basis or ""})
        if new != prior:
            changes.append((r["species"], prior, new, src))
    header = [h for h in header if not h.startswith("# PROVISIONAL")] + [
        "# source: which database(s) set the label (rule in scripts/check_invasion_dbs.py); "
        "'general knowledge' = species outside GISD and outside WRiMS scope (terrestrial/freshwater).",
        "# evidence: database hits; label_general_knowledge: the earlier label drafted 2026-11-27; notes: common name and context.",
        f"# Checked against GISD and WRiMS on {time.strftime('%Y-%m-%d')}; full evidence in manifest/invasion_db.tsv.",
    ]
    cols = ["species", "invasion", "source", "evidence", "label_general_knowledge", "notes"]
    with open(a.traits, "w", newline="") as fh:
        fh.write("\n".join(header) + "\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols])
    print(f"\n{len(changes)} labels changed:")
    for sp, o, n, s in changes:
        print(f"  {sp}: {o} -> {n}  [{s}]")


if __name__ == "__main__":
    main()
