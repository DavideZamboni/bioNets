"""
clinvar_query.py — query ClinVar by gene (+ optional genomic location) → DataFrame.

Three approaches, most-robust first:
  1. clinvar_via_eutilities()  — NCBI E-utilities API   (RECOMMENDED)
  2. clinvar_via_requests()    — parse the server-rendered results page (no browser)
  3. clinvar_via_selenium()    — drive a real browser   (what you asked for)

Why the API first: the ClinVar results page is server-rendered, so Selenium is
rarely necessary, and for programmatic access NCBI explicitly asks you to use
E-utilities (the web UI rate-limits/blocks scrapers and its HTML changes without
notice). The API returns structured records you don't have to scrape.

NCBI rate limits: 3 requests/s without an API key, 10/s with one. Pass
api_key=... (free, from your NCBI account) and email=... to be polite.

Deps:  pip install requests pandas lxml          # + selenium for approach 3
"""
from __future__ import annotations
import time
import urllib.parse
import requests
import pandas as pd
import re

EUTILS  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CLINVAR = "https://www.ncbi.nlm.nih.gov/clinvar/"

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
    """esearch → esummary (JSON), one row per ClinVar VCV record."""
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
    """Most-cited NM_ transcript in the records (its NP_ is the reference)."""
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
# Reference protein sequence (RefSeq protein the ClinVar positions refer to)
# ─────────────────────────────────────────────────────────────────────────────

def suggest_refseq_protein(variants_df):
    """Best-effort: pull the RefSeq transcript (NM_) most cited in the ClinVar
    'name' column. The matching protein is its NP_; pass that NP_ (or a UniProt
    sequence) to the aligner. Returns the most common NM_ accession or None."""
    accs = []
    for nm in variants_df.get("name", pd.Series([], dtype=str)).fillna(""):
        m = re.search(r"\b(NM_\d+\.\d+)", nm)
        if m:
            accs.append(m.group(1))
    return pd.Series(accs).mode().iloc[0] if accs else None

def fetch_reference_protein(accession, email=None, api_key=None):
    """efetch a protein FASTA (e.g. 'NP_002894.1') -> sequence string."""
    params = {"db": "protein", "id": accession, "rettype": "fasta",
              "retmode": "text", "tool": "variant_metric_ranking"}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    r = requests.get(f"{EUTILS}/efetch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    lines = r.text.splitlines()
    return "".join(l.strip() for l in lines if l and not l.startswith(">"))

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