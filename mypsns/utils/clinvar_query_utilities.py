"""
clinvar_query_utilities.py — query ClinVar by gene (+ optional genomic window)
and fetch RefSeq protein sequences, through NCBI E-utilities only.

Why E-utilities and not the web page: NCBI asks programmatic users to use the
API (the web UI rate-limits/blocks scrapers and its HTML changes without
notice), and the API returns structured records you don't have to scrape.

  fetch_clinvar_variants()   esearch → esummary, one row per VCV record
                             (used by services/compare_variants.py)
  fetch_reference_protein()  efetch a protein FASTA → sequence
  suggest_refseq_protein()   most-cited NM_ transcript in a variants table
  clinvar_via_eutilities()   older standalone variant of fetch_clinvar_variants
                             with a richer, assembly-aware parser
                             (parse_clinvar_summary); not used by the CLI.

The ClinVar `[chrpos]` search field is GRCh38 (GRCh37 would be `[chrpos37]`),
so the optional chrom/start/stop window is always in GRCh38 coordinates.

NCBI rate limits: 3 requests/s without an API key, 10/s with one. Pass
api_key=... (free, from your NCBI account) and email=... to be polite.

Deps: requests (optional `variants` extra: pip install -e .[variants])
"""
from __future__ import annotations
import time
import requests
import pandas as pd
import re

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# ── ClinVar via E-utilities ─────────────────────────────────────────────
def _eutils(endpoint: str, params: dict, email=None, api_key=None,
            timeout=60) -> requests.Response:
    p = {"tool": "mypsns_compare_variants", **params}
    if email:
        p["email"] = email
    if api_key:
        p["api_key"] = api_key
    r = requests.get(f"{EUTILS}/{endpoint}", params=p, timeout=timeout)
    r.raise_for_status()
    return r

def fetch_clinvar_variants(gene, chrom=None, start=None, stop=None,
                           assembly="GRCh38", email=None, api_key=None,
                           retmax=2000, pause=0.34) -> pd.DataFrame:
    """esearch → esummary (JSON), one row per ClinVar VCV record.

    `gene[gene]` alone already returns every record touching the gene; the
    optional chrom/start/stop window (all three or none) is a GRCh38 filter on
    top. `assembly` is not used by the search (`[chrpos]` is always GRCh38) and
    is kept only for signature compatibility.
    """
    term = f"{gene}[gene]"
    if chrom and start and stop:
        c = str(chrom).lower().replace("chr", "")
        term += f" AND {c}[chr] AND {int(start)}:{int(stop)}[chrpos]"
    js = _eutils("esearch.fcgi",
                 {"db": "clinvar", "term": term, "retmax": retmax,
                  "retmode": "json"}, email, api_key, timeout=30).json()
    ids = js["esearchresult"].get("idlist", [])
    rows = []
    for i in range(0, len(ids), 200):
        time.sleep(pause)
        res = _eutils("esummary.fcgi",
                      {"db": "clinvar", "id": ",".join(ids[i:i + 200]),
                       "retmode": "json"}, email, api_key).json().get("result", {})
        for uid in res.get("uids", []):
            rows.append(_parse_summary(res[uid]))
    return pd.DataFrame(rows)

def fetch_reference_protein(accession, email=None, api_key=None) -> str:
    """efetch a protein FASTA (e.g. 'NP_002894.1') → sequence."""
    txt = _eutils("efetch.fcgi",
                  {"db": "protein", "id": accession, "rettype": "fasta",
                   "retmode": "text"}, email, api_key, timeout=30).text
    return "".join(l.strip() for l in txt.splitlines()
                   if l and not l.startswith(">"))

def suggest_refseq_protein(variants_df) -> str | None:
    """Best-effort: the RefSeq transcript (NM_) most cited in the ClinVar
    'name' column. The matching protein is its NP_; pass that NP_ (or a UniProt
    sequence) to the aligner. Returns the most common NM_ accession or None."""
    accs = [m.group(1) for nm in variants_df.get("name", pd.Series([], dtype=str)).fillna("")
            for m in [re.search(r"\b(NM_\d+\.\d+)", nm)] if m]
    return pd.Series(accs).mode().iloc[0] if accs else None

def _parse_summary(rec: dict) -> dict:
    clin = rec.get("germline_classification") or rec.get("clinical_significance") or {}
    traits = clin.get("trait_set") or rec.get("trait_set") or []
    v0 = (rec.get("variation_set") or [{}])[0]
    return {
        "accession": rec.get("accession"),
        "name": rec.get("title"),
        "protein_change": rec.get("protein_change"),
        "molecular_consequence": "; ".join(rec.get("molecular_consequence_list") or []),
        "germline_classification": clin.get("description"),
        "review_status": clin.get("review_status"),
        "conditions": "; ".join(t.get("trait_name", "") for t in traits if t.get("trait_name")),
        "variant_type": v0.get("variant_type") or rec.get("obj_type"),
    }

def _norm_chr(chrom) -> str:
    """'chr17' / 17 / 'Chr17' -> '17'."""
    return str(chrom).lower().replace("chr", "").strip()

def parse_clinvar_summary(rec: dict, assembly="GRCh38") -> dict:
    """Extract the useful fields from one esummary record (schema-tolerant).

    NCBI renamed the classification block from 'clinical_significance' to
    'germline_classification' in the 2024 schema; read whichever is present.
    """
    clin = (rec.get("germline_classification")
            or rec.get("clinical_significance") or {})
    traits = clin.get("trait_set") or rec.get("trait_set") or []
    conditions = "; ".join(t.get("trait_name", "") for t in traits if t.get("trait_name"))

    genes = ", ".join(g.get("symbol", "") for g in rec.get("genes", []) if g.get("symbol"))

    vset = rec.get("variation_set") or [{}]
    v0 = vset[0] if vset else {}
    loc = {}
    for l in v0.get("variation_loc", []) or []:
        if l.get("assembly_name") == assembly:
            loc = l
            break

    mc = rec.get("molecular_consequence_list") or []

    return {
        "uid":                    rec.get("uid"),
        "accession":              rec.get("accession"),
        "name":                   rec.get("title"),
        "gene(s)":                genes,
        "variant_type":           v0.get("variant_type") or rec.get("obj_type"),
        "protein_change":         rec.get("protein_change"),
        "molecular_consequence":  "; ".join(mc),
        "germline_classification": clin.get("description"),
        "review_status":          clin.get("review_status"),
        "last_evaluated":         clin.get("last_evaluated"),
        "conditions":             conditions,
        "assembly":               loc.get("assembly_name"),
        "chr":                    loc.get("chr"),
        "start":                  loc.get("start") or loc.get("display_start"),
        "stop":                   loc.get("stop")  or loc.get("display_stop"),
        "cdna_change":            v0.get("cdna_change"),
        "canonical_spdi":         v0.get("canonical_spdi"),
    }

# ─────────────────────────────────────────────────────────────────────────────
# E-utilities clinvar API
# ─────────────────────────────────────────────────────────────────────────────
def clinvar_via_eutilities(gene, chrom=None, start=None, stop=None,
                           assembly="GRCh38", api_key=None, email=None,
                           retmax=2000, pause=0.34) -> pd.DataFrame:
    """esearch -> esummary (JSON). Returns one row per ClinVar VCV record.

    For a single gene, gene[gene] already returns every variant in that gene;
    the chr/start/stop window is an extra filter (useful for intergenic /
    multi-gene regions).
    """
    term = f'{gene}[gene]'
    if chrom and start and stop:
        term += f' AND {_norm_chr(chrom)}[chr] AND {int(start)}:{int(stop)}[chrpos]'

    common = {"tool": "clinvar_query"}
    if api_key:
        common["api_key"] = api_key
    if email:
        common["email"] = email

    # esearch -> UID list
    r = requests.get(f"{EUTILS}/esearch.fcgi",
                     params={"db": "clinvar", "term": term, "retmax": retmax,
                             "retmode": "json", **common}, timeout=30)
    r.raise_for_status()
    idlist = r.json()["esearchresult"].get("idlist", [])
    if not idlist:
        return pd.DataFrame()

    # esummary in batches of 200 (JSON)
    rows = []
    for i in range(0, len(idlist), 200):
        batch = idlist[i:i + 200]
        time.sleep(pause)
        s = requests.get(f"{EUTILS}/esummary.fcgi",
                         params={"db": "clinvar", "id": ",".join(batch),
                                 "retmode": "json", **common}, timeout=60)
        s.raise_for_status()
        res = s.json().get("result", {})
        for uid in res.get("uids", []):
            rows.append(parse_clinvar_summary(res[uid], assembly))
    return pd.DataFrame(rows)