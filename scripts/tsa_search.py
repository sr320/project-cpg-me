"""List NCBI TSA (transcriptome shotgun assembly) projects for species, for manifest gaps.

Usage:  python scripts/tsa_search.py "Genus species" ["Genus species" ...]
Prints  species, TSA master accession, organism, date, title. A project like GFYR00000000 is
        written GFYR01000000 in the manifest (version 01; check the FTP area
        https://ftp.ncbi.nlm.nih.gov/genbank/tsa/<first letter>/ for later versions).
"""
import json, sys, time, urllib.parse, urllib.request
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u):
    for i in range(4):
        try: return urllib.request.urlopen(urllib.request.Request(u,headers={"User-Agent":"cpgme/0.1"}),timeout=60).read()
        except Exception: time.sleep(3*(i+1))
    raise
for sp in sys.argv[1:]:
    term=f'"{sp}"[Organism] AND tsa-master[prop]'
    ids=json.loads(get(E+"esearch.fcgi?"+urllib.parse.urlencode({"db":"nuccore","term":term,"retmode":"json","retmax":50})))["esearchresult"]["idlist"]
    time.sleep(0.4)
    if not ids: print(f"{sp}\t-\t0 TSA"); continue
    s=json.loads(get(E+"esummary.fcgi?"+urllib.parse.urlencode({"db":"nuccore","id":",".join(ids),"retmode":"json"})))["result"]
    for i in ids:
        r=s[i]
        print(f"{sp}\t{r['caption']}\t{r.get('organism','')}\t{r.get('createdate','')}\t{r['title'][:110]}")
    time.sleep(0.4)
