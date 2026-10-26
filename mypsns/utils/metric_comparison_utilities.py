from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import networkx as nx
import re

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# ── residue → 1-letter coding (inlined; Biopython is used only for alignment) ─
# A 26-entry lookup doesn't justify importing Biopython's Data submodule (and
# coupling to its internals). The standard 20 + IUPAC ambiguity codes (B/X/Z/J/
# U/O) are copied from Bio.Data.IUPACData; forcefield protonation/modification
# variants that no library carries are merged in. Unknown → 'X'.
_THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
    "ASX": "B", "XAA": "X", "GLX": "Z", "XLE": "J", "SEC": "U", "PYL": "O",
    # forcefield His tautomers/protomers + modified residues
    "HSD": "H", "HSE": "H", "HSP": "H", "HID": "H", "HIE": "H", "HIP": "H",
    "GLM": "G", "CYX": "C", "CYM": "C", "LYN": "K", "ASH": "D", "GLH": "E",
    "ARN": "R", "LYK": "K",
    # nucleic acids — only needed if a nucleic chain is ever aligned (PSN work
    # rarely is; the Entity column + FF data already identify residues). Extend
    # here if a use case appears:
    "DA": "A", "DG": "G", "DC": "C", "DT": "T", "DU": "U",
    "A": "A", "G": "G", "C": "C", "U": "U", "T": "T",
}

def res3to1(resname: str) -> str:
    return _THREE2ONE.get((resname or "").strip().upper(), "X")

# ── structure sequences from the interactor table (Entity-driven) ────────
def chain_sequences(int_top: pd.DataFrame, resindex_to_pdb: Dict[int, Tuple[str, int]]) -> Dict[str, List[int]]:
    """{chain: [(resindex, resid, resname, aa1)]} for polymeric residues
    only (Entity == 'P'), ordered along the chain. No ion/water blocklist."""
    res = int_top.drop_duplicates("Resindex")
    res = res[res["Entity"] == "P"]
    chains = {}
    for ch, grp in res.groupby("Chain"):
        grp = grp.sort_values("Resindex")
        out = []
        for _, r in grp.iterrows():
            rix = int(r["Resindex"])
            resid = resindex_to_pdb.get(rix, (ch, r.get("resnum", rix)))[1]
            out.append((rix, int(resid), r["resname"], res3to1(r["resname"])))
        chains[ch] = out
    return chains

# ── alignment: reference position → resindex (auto-pick chain) ───────────
def map_reference_positions(int_top, resindex_to_pdb, ref_seq, min_identity=70.0,
                            open_gap=-10.0, extend_gap=-0.5):
    from Bio.Align import PairwiseAligner, substitution_matrices
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap
    aligner.mode = "global"
    best = (None, {}, -1.0)
    chain_seqs = chain_sequences(int_top, resindex_to_pdb)
    for ch, residues in chain_seqs:
        s_seq = "".join(x[3] for x in residues)
        if not s_seq:
            continue
        aln = aligner.align(ref_seq, s_seq)[0]
        a_ref, a_str = aln[0], aln[1]
        mapping, ri, si, ident = {}, 0, 0, 0
        for rc, sc in zip(a_ref, a_str):
            if rc != "-" and sc != "-":
                mapping[ri + 1] = residues[si]
                ident += (rc == sc)
            ri += rc != "-"
            si += sc != "-"
        pct = ident / max(len(s_seq), 1) * 100
        if pct > best[2]:
            best = (ch, mapping, pct)
    if best[2] < min_identity:
        raise ValueError(f"No chain aligns to the reference above "
                         f"{min_identity}% (best {best[0]} = {best[2]:.1f}%).")
    return best

# ── per-position truth table ─────────────────────────────────────────────────
_PROT_RE = re.compile(r"(?:p\.)?([A-Za-z]{1,3})(\d+)([A-Za-z]{1,3})")

def parse_protein_pos(text):
    if not isinstance(text, str):
        return None
    m = _PROT_RE.search(text)
    return int(m.group(2)) if m else None

def classify_significance(s):
    s = (s or "").lower()
    if "conflict" in s:
        return "vus"
    if "pathogenic" in s:
        return "pathogenic"
    if "benign" in s:
        return "benign"
    return "vus"

# ── scoring helpers ──────────────────────────────────────────────────────────
def _auroc(x, y):
    from scipy.stats import rankdata
    x = np.asarray(x, float); y = np.asarray(y, int)
    n1, n0 = int(y.sum()), int((y == 0).sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(x)
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

def _auprc(x:float, y: float) -> float:
    x = np.asarray(x, float); y = np.asarray(y, int)[np.argsort(-x)]
    tp, fp = np.cumsum(y), np.cumsum(1 - y)
    prec = tp / np.maximum(tp + fp, 1)
    rec = tp / max(int(y.sum()), 1)
    return float(np.sum((rec - np.concatenate([[0.0], rec[:-1]])) * prec))

# ── in-memory node/edge metrics from inherited network data ──────────────
def node_metric_table(G: nx.Graph, cdf: pd.DataFrame) -> pd.DataFrame:
    """Combine inherited centrality_df + graph degree + edge-derived features
    for one (inttype, source) graph → DataFrame keyed by resindex."""
    cdf["degree"] = cdf["resindex"].map(dict(G.degree())).fillna(0).astype(int)
    # edge-derived node features (edges evaluated on equal footing)
    rows = []
    for u, v, d in G.edges(data=True):
        w = float(d.get("weight", 0.0))
        rows += [(u, w), (v, w)]
    if rows:
        ef = pd.DataFrame(rows, columns=["resindex", "w"]).groupby("resindex")["w"]
        feats = pd.DataFrame({"edge_strength": ef.sum(), "edge_count": ef.size(),
                              "edge_weight_mean": ef.mean(),
                              "edge_weight_max": ef.max()}).reset_index()
        cdf = cdf.merge(feats, on="resindex", how="left")
    return cdf

# missense truth table builder
def missense_truth_table(variants_df: pd.DataFrame) -> pd.DataFrame:
    df = variants_df.copy()
    if df.empty:
        return pd.DataFrame()
    mc = df.get("molecular_consequence", pd.Series("", index=df.index)).fillna("")
    df = df[mc.str.contains("missense", case=False)].copy()
    pc = df.get("protein_change", pd.Series("", index=df.index)).fillna("")
    nm = df.get("name", pd.Series("", index=df.index)).fillna("")
    df["ref_pos"] = [parse_protein_pos(a) or parse_protein_pos(b) for a, b in zip(pc, nm)]
    df = df.dropna(subset=["ref_pos"])
    df["ref_pos"] = df["ref_pos"].astype(int)
    df["sig_class"] = df["germline_classification"].map(classify_significance)
    rows = []
    for pos, g in df.groupby("ref_pos"):
        vc = g["sig_class"].value_counts()
        npath, nben = int(vc.get("pathogenic", 0)), int(vc.get("benign", 0))
        rows.append({"ref_pos": pos, "n_pathogenic": npath, "n_benign": nben,
                     "n_vus": int(vc.get("vus", 0)),
                     "label": "non_benign" if npath else "benign" if nben else "vus"})
    return pd.DataFrame(rows).set_index("ref_pos").sort_index()

def score(merged, negative, NODE_METRICS: List[str], EDGE_FEATURES: List[str]) -> pd.DataFrame:
    from scipy.stats import mannwhitneyu
    pos = merged["label"].eq("non_benign")
    neg = merged["label"].eq("benign") if negative == "benign" else ~pos
    cols = [c for c in NODE_METRICS + EDGE_FEATURES if c in merged.columns]
    out = []
    for m in cols:
        sub = merged.loc[pos | neg, [m]].assign(y=pos[pos | neg].astype(int)).dropna()
        y, x = sub["y"].to_numpy(), sub[m].to_numpy()
        if int(y.sum()) == 0 or int((y == 0).sum()) == 0:
            continue
        auroc = _auroc(x, y)
        xs = x if auroc >= 0.5 else -x
        try:
            p = mannwhitneyu(x[y == 1], x[y == 0], alternative="two-sided").pvalue
        except ValueError:
            p = np.nan
        out.append({"metric": m, "kind": "edge" if m in EDGE_FEATURES else "node",
                    "auroc": round(auroc, 3), "auroc_oriented": round(max(auroc, 1 - auroc), 3),
                    "direction": "higher" if auroc >= 0.5 else "lower",
                    "auprc": round(_auprc(xs, y), 3), "mwu_p": p,
                    "n_pos": int(y.sum()), "n_neg": int((y == 0).sum())})
    return (pd.DataFrame(out).sort_values("auroc_oriented", ascending=False)
            .reset_index(drop=True))