"""
Configuration system for the myPSNs CLI.

All configuration flows through Pydantic models. Each stage owns its own
config; PipelineConfig composes them into the full run.

Loading precedence (lowest → highest):
    1. Pydantic defaults
    2. YAML file (--config)
    3. CLI argument overrides
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Literal, Optional

import yaml
from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator, model_validator


# Default NA strings shared across CSV reads in the pipeline.
# Note: "NA" is deliberately omitted to preserve sodium atom names.
DEFAULT_NA_VALUES: List[str] = [
    "", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN",
    "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A",
    "NULL", "NaN", "n/a", "nan", "null",
]

# ─────────────────────────────────────────────────────────────────────────────
# Per-stage configurations
# ─────────────────────────────────────────────────────────────────────────────

class FFReaderConfig(BaseModel):
    """Configuration for the forcefield-reader stage."""

    model_config = ConfigDict(extra="forbid")

    nonbonded_name: str = Field(
        default="ffnonbonded.itp",
        description="Name of the .itp file holding nonbonded parameters.",
    )
    additional_residues: bool = Field(
        default=True,
        description="Merge in extra residue definitions from resources/.",
    )
    na_values: List[str] = Field(default_factory=lambda: DEFAULT_NA_VALUES.copy())


class TopologyConfig(BaseModel):
    """Configuration for the topology-recognition stage."""

    model_config = ConfigDict(extra="forbid")

    max_to_check: int = Field(
        default=30,
        description=(
            "Max number of forcefield entries to inspect per residue topology. "
            "-1 = all (sorted by name similarity)."
        ),
    )
    na_values: List[str] = Field(default_factory=lambda: DEFAULT_NA_VALUES.copy())


# Default named-interactor (PyInteraph-style) role definitions, in the
# YAML-native form recorded in run_config.yaml. Each role maps to a list of
# RULES; a rule is a dict {attribute: values} whose conditions are AND-combined,
# and the rules in a role are OR-combined. A key prefixed with "not " negates.
# This mirrors NAMED_INTERACTOR_DEFAULTS in interactors_utilities (kept in sync
# by a test) so config.py need not import the heavy utilities module.
NAMED_INTERACTOR_DEFAULTS: Dict[str, list] = {
    # Hydrophobic side chains (PyInteraph hc_reslist); side-chain atoms only.
    # PyInteraph's EXACT exclusion set: the CHARMM amide HN is intentionally NOT
    # excluded (PyInteraph includes it in the side-chain COM), which is what
    # matches HC to ~0.1%.
    "HC": [
        {
            "resname": ["ALA", "VAL", "LEU", "ILE", "PHE", "PRO", "TRP", "MET"],
            "not atom_name": ["N", "CA", "C", "O", "H", "H1", "H2", "H3",
                              "O1", "O2", "OXT"],
        },
    ],
    # hydrogen_bonds.ini DONORS
    "HDon": [
        {"atom_name": ["OH", "OG", "NE2", "OG1", "NE", "N", "ND1", "NZ",
                       "NH1", "NH2", "OW", "ND2", "SG", "OH2", "NE1"]},
    ],
    # hydrogen_bonds.ini ACCEPTORS
    "HAcc": [
        {"atom_name": ["OH", "OG", "OD1", "OD2", "OG1", "O", "ND1", "NE2",
                       "OE2", "OW", "SG", "OE1", "OH2", "SD"]},
    ],
    # charged_groups.ini POSITIVE groups (*p). His is included but only acts as
    # basic when protonated — enforced by the net-charge gate (min_residue_charge),
    # so neutral His does not form salt bridges.
    "Basic": [
        {"resname": ["NH3+", "GLY-NH3+", "NHT3+", "PRO-NH2+"],
         "atom_name": ["N", "H1", "H2", "H3", "HT1", "HT2", "HT3"]},
        {"resname": ["ARG"],
         "atom_name": ["CD", "NE", "CZ", "NH1", "HH11", "HH12",
                       "NH2", "HH21", "HH22"]},
        {"resname": ["HIS", "HSD", "HSE", "HSP", "HSA", "HSB",
                     "HIS1", "HISH", "HID", "HIE", "HIP"],
         "atom_name": ["CG", "ND1", "HD1", "CD2", "CE1", "NE2", "HE2"]},
        {"resname": ["LYS", "LYSH"],
         "atom_name": ["CE", "NZ", "HZ1", "HZ2", "HZ3"]},
    ],
    # charged_groups.ini NEGATIVE groups (*n) + C-terminus (all oxygen variants)
    "Acidic": [
        {"resname": ["ASP"], "atom_name": ["CG", "OD1", "OD2"]},
        {"resname": ["GLU"], "atom_name": ["CD", "OE1", "OE2"]},
        {"resname": ["COO-"],
         "atom_name": ["C", "O1", "O2", "O", "OXT", "OT1", "OT2"]},
    ],
    "hb_masses": [1],
}


class InteractorsConfig(BaseModel):
    """Configuration for the interactors stage."""

    model_config = ConfigDict(extra="forbid")

    modes: List[Literal["q_interactors", "named_interactors"]] = Field(
        default_factory=lambda: ["named_interactors"],
        description="Interactor-typing schemes to run.",
    )

    # HBond
    max_DA_charge: float = Field(
        default=-0.3,
        description="Maximum (least negative) charge for an HBond donor/acceptor.",
    )
    min_H_charge: float = Field(
        default=0.3,
        description="Minimum H-atom charge for HBond participation.",
    )
    DA_masses: List[int] = Field(
        default_factory=lambda: [14, 16, 19],
        description="Rounded masses recognised as D/A elements (N/O/F).",
    )

    # Salt bridge
    min_abs_charge: float = Field(
        default=0.75,
        description="Minimum |local charge| for Acidic/Basic classification.",
    )

    # Hydrophobic cluster
    max_abs_charge: float = Field(
        default=0.27,
        description="Maximum |charge| for an atom to be HC-eligible.",
    )
    max_num_heavy_ats: int = Field(
        default=3, ge=1,
        description="Max heavy atoms per HC sub-cluster.",
    )
    min_bonds_to_charged: int = Field(
        default=3, ge=0,
        description="Min bond-graph distance from a charged atom for HC eligibility.",
    )
    min_bonds_to_polar: int = Field(
        default=1, ge=0,
        description="Min bond-graph distance from a polar atom for HC eligibility.",
    )

    # Net-charge gate for named_interactors salt-bridge roles (Basic/Acidic):
    # a residue is treated as basic/acidic only if |sum of its partial charges|
    # exceeds this value. Makes salt-bridge participation protonation-state aware
    # (neutral His ≈ 0 excluded; protonated His ≈ +1 included) regardless of
    # resname. Set to 0 to disable the gate (pure resname/atom-name selection).
    min_residue_charge: float = Field(
        default=0.8, ge=0.0,
        description="Min |net residue charge| for Basic/Acidic (salt-bridge) roles.",
    )

    # ── named_interactors role definitions ──────────────────────────────────
    # Full PyInteraph-style role rule set, recorded in run_config.yaml and used
    # verbatim (it REPLACES the built-in defaults when set). Each role → list of
    # {attribute: values} rules; conditions AND-combine within a rule, rules
    # OR-combine across a role, and a key prefixed "not " negates. To customise:
    # copy this block from a run_config.yaml into your config, edit it, re-run
    # the interactors stage.
    named_interactors_rules: Dict[str, Any] = Field(
        default_factory=lambda: {k: [dict(r) if isinstance(r, dict) else r
                                      for r in v] if isinstance(v, list) else v
                                 for k, v in NAMED_INTERACTOR_DEFAULTS.items()},
        description="named_interactors role definitions (PyInteraph-style); "
                    "replaces the defaults when set.",
    )

    na_values: List[str] = Field(default_factory=lambda: DEFAULT_NA_VALUES.copy())

class ContactsConfig(BaseModel):
    """Configuration for the contacts stage."""

    model_config = ConfigDict(extra="forbid")

    # Parallelization
    nc: int = Field(default=2, ge=1, description="Number of parallel workers.")
    ram_thresh: int = Field(
        default=80, ge=10, le=99,
        description="RAM-usage %% that freezes HBond batch growth.",
    )
    nframes: int = Field(
        default=500, ge=1,
        description="Initial frame batch size for adaptive HBond batching.",
    )

    # Trajectory window
    beg_frm: int = Field(default=0, ge=0)
    end_frm: int = Field(default=-1)

    # Interaction selectors
    interactors_type: List[Literal["q_interactors", "named_interactors"]] = Field(
        default_factory=lambda: ["q_interactors", "named_interactors"],
    )
    targets: List[Literal["HBonds", "Salt_bridges", "HCs"]] = Field(
        default_factory=lambda: ["HBonds", "Salt_bridges", "HCs"],
    )

    # Chain selectors (None → auto-detect)
    subject_chains: Optional[List[str]] = None
    surrounding_chains: Optional[List[str]] = None

    # Geometric cutoffs
    d_a_dist: float = Field(default=4.5, gt=0, description="D–A cutoff [Å]. PyInteraph-closest (≈ its 3.5 Å H–A).")
    angle: float = Field(default=120.0, ge=0, le=180, description="D–H–A angle cutoff [°]. PyInteraph 1 value.")
    sb_cutoff: float = Field(default=4.5, gt=0, description="Salt-bridge cutoff [Å]. PyInteraph 1 --sb-co.")
    hc_cutoff: float = Field(default=5.0, gt=0, description="Hydrophobic CoM cutoff [Å]. PyInteraph 1 --hc-co.")

    # Preliminary network analysis
    matrix_type: List[Literal["persistency", "energy"]] = Field(
        default_factory=lambda: ["persistency", "energy"],
        description="Which post-contact matrices to write.",
    )
    energy_use_abs: bool = Field(
        default=False,
        description="Whether to retain repulsive forces present in energy network analysis.",
    )
    thr_step: float = Field(
        default=0.1, gt=0.0, le=50.0,
        description="Step of % persistency or energy used when computing network statistics.",
    )


    
    @field_validator("matrix_type")
    @classmethod
    def _matrix_type_non_empty(cls, v):
        if not v:
            raise ValueError("matrix_type must select at least one of: persistency, energy")
        return list(dict.fromkeys(v))   # dedupe, keep order


    na_values: List[str] = Field(default_factory=lambda: DEFAULT_NA_VALUES.copy())

class NetworkConfig(BaseModel):
    """Configuration for the network stage (standalone — not in `run`)."""

    model_config = ConfigDict(extra="forbid")
    # ── Frame count (informational — logged only, not required) ──
    n_frames: Optional[int] = Field(
        default=None, ge=1,
        description="Trajectory frames analysed during Contacts. Optional; logged only.",
    )
    # ── Which contact tables to consume ──
    interactors_type: List[Literal["q_interactors", "named_interactors"]] = Field(
        default_factory=lambda: ["q_interactors", "named_interactors"],
        description="Adjacency/ sub-folders to scan.",
    )
    targets: List[Literal["HBonds", "Salt_bridges", "HCs"]] = Field(
        default_factory=lambda: ["HBonds", "Salt_bridges", "HCs"],
        description="Contact types contributing to the network.",
    )
    # ── Edge-weight source(s): persistency-based, energy-based, or both ──
    edge_weight_sources: List[Literal["persistency", "energy"]] = Field(
        default_factory=lambda: ["persistency", "energy"],
        description=("Build one static network per entry. 'persistency' uses the "
                     "*_adjacency.dat matrices; 'energy' uses *_total_energy.dat."),
    )
    energy_use_abs: bool = Field(
        default=False,
        description=("Energy → weight transform. False: weight = max(0, -E) "
                     "(attractive only). True: weight = |E| (also repulsive)."),
    )
    # ── Time-resolved energy network (optional) ──
    window_size: Optional[int] = Field(
        default=None, ge=1,
        description=("If set, also build one network per window of this many frames "
                     "from the per-frame energy CSVs and write community_evolution_*.csv."),
    )
    window_stride: Optional[int] = Field(
        default=None, ge=1,
        description="Step between windows. Defaults to window_size (non-overlapping).",
    )
    # ── Frame range (applies to BOTH static and time-resolved) ──
    frame_start: Optional[int] = Field(
        default=None, ge=0,
        description=("First frame (inclusive). When set, static networks are "
                     "re-aggregated from the per-frame energy CSVs over the range."),
    )
    frame_stop: Optional[int] = Field(
        default=None, ge=1,
        description="Last frame (exclusive). None → end of trajectory.",
    )
    # ── Output metrics (one .pml per node × community combination) ──
    node_metrics: List[Literal[
        "degree", "weighted_degree",
        "degree_centrality", "betweenness_centrality", "closeness_centrality",
    ]] = Field(
        default_factory=lambda: ["betweenness_centrality"],
        description="Node centrality metrics; each drives one CA-sphere scaling.",
    )
    community_methods: List[Literal[
        "louvain", "greedy_modularity", "label_propagation", "connected_components",
    ]] = Field(
        default_factory=lambda: ["louvain"],
        description="Community algorithms; each produces one PyMOL view.",
    )
    community_seed: int = Field(
        default=42,
        description="Seed for stochastic community algorithms.",
    )
    community_min_members: int = Field(
        default=2, ge=1,
        description="Minimum community size to keep; smaller communities are "
                    "dissolved (nodes marked solitary) and omitted from community "
                    "renders. 1 keeps every community including singletons.",
    )
    # ── Filtering ──
    # Edge thresholds are per-source and both expressed as a percentage on the
    # 0..100 scale the elbow / degree-distribution plots are scanned on:
    #   persistency % = 100 · contact-frequency
    #   energy %      = 100 · strength / max(strength)  (strength = max(0,-E) or |E|)
    # An edge is kept when its weight (in %) is strictly greater than the threshold.
    min_persistency: float = Field(
        default=0.0, ge=0.0, le=100.0,
        description="Drop persistency edges with weight ≤ this (percent, 0..100).",
    )
    min_energy: float = Field(
        default=0.0, ge=0.0, le=100.0,
        description="Drop energy edges with weight ≤ this (percent of strongest edge, 0..100).",
    )
    min_degree: int = Field(
        default=0, ge=0,
        description="Drop nodes (and edges) with degree below this.",
    )
    # ── PyMOL rendering ──
    reference_pdb: Optional[Path] = Field(
        default=None,
        description="PDB loaded in PyMOL; defaults to pipeline topology.",
    )
    sphere_scale: float = Field(
        default=2.0, gt=0.0,
        description="Max CA-sphere vdw radius [Å] (largest node metric).",
    )
    stick_scale: float = Field(
        default=0.5, gt=0.0,
        description="Max edge stick (dash) radius [Å] (largest edge weight).",
    )
    # ── Optional auxiliary outputs ──
    write_graphml: bool = Field(
        default=False,
        description="Also write a GraphML for Gephi/Cytoscape.",
    )
    shortest_paths: bool = Field(
        default=False,
        description="Compute and dump all-pairs shortest paths.",
    )
    na_values: List[str] = Field(default_factory=lambda: DEFAULT_NA_VALUES.copy())
 
    # ── Validators: non-empty lists, drop duplicates ──
    @field_validator(
        "interactors_type", "targets", "edge_weight_sources",
        "node_metrics", "community_methods",
    )
    @classmethod
    def _non_empty_unique(cls, v):
        if not v:
            raise ValueError("must select at least one entry")
        return list(dict.fromkeys(v))   # dedupe, preserve order
    @field_validator("reference_pdb", mode="before")
    @classmethod
    def _coerce_pdb_path(cls, v):
        from pathlib import Path
        if v is None or isinstance(v, Path):
            return v
        return Path(str(v)).expanduser()

    @model_validator(mode="after")
    def _check_frame_range(self):
        if (self.frame_start is not None and self.frame_stop is not None
                and self.frame_stop <= self.frame_start):
            raise ValueError(
                f"frame_stop ({self.frame_stop}) must be > frame_start "
                f"({self.frame_start})"
            )
        return self

# ─────────────────────────────────────────────────────────────────────────────
# Pipeline (top-level) configuration
# ─────────────────────────────────────────────────────────────────────────────

class PipelineConfig(BaseModel):
    """Full pipeline configuration."""

    model_config = ConfigDict(extra="forbid")

    # General
    topology: Path = Field(description="Topology file (PDB, GRO, …).")
    trajectory: Optional[Path] = Field(
        default=None, description="MD trajectory (only required for `contacts`).",
    )
    forcefield: Path = Field(description="Forcefield folder (GROMACS-style).")
    outpath: Path = Field(description="Root output folder for the run.")

    # Per-stage
    ff_reader: FFReaderConfig = Field(default_factory=FFReaderConfig)
    topology_recognition: TopologyConfig = Field(default_factory=TopologyConfig)
    interactors: InteractorsConfig = Field(default_factory=InteractorsConfig)
    contacts: ContactsConfig = Field(default_factory=ContactsConfig)
    network: NetworkConfig = Field(default_factory=NetworkConfig)

    @field_validator(
        "topology", "trajectory", "forcefield", "outpath",
        mode="before",
    )
    @classmethod
    def _coerce_path(cls, v: Any) -> Any:
        if v is None or isinstance(v, Path):
            return v
        return Path(str(v)).expanduser()

# ─────────────────────────────────────────────────────────────────────────────
# Loading / merging / saving
# ─────────────────────────────────────────────────────────────────────────────

class ConfigError(Exception):
    """Raised when configuration loading or validation fails."""

def load_yaml(path: Path) -> dict:
    """Load a YAML file into a dict (empty dict if file is empty)."""
    path = Path(path).expanduser()
    if not path.exists():
        raise ConfigError(f"Config file does not exist: {path}")
    try:
        with path.open("r") as f:
            data = yaml.safe_load(f) or {}
    except yaml.YAMLError as e:
        raise ConfigError(f"Failed to parse YAML {path}: {e}") from e
    if not isinstance(data, dict):
        raise ConfigError(
            f"Top-level YAML in {path} must be a mapping; got {type(data).__name__}."
        )
    return data

def _deep_clean(d: dict) -> dict:
    """Strip None leaves recursively. None ≡ unset CLI option."""
    out: dict = {}
    for k, v in d.items():
        if v is None:
            continue
        if isinstance(v, dict):
            cleaned = _deep_clean(v)
            if cleaned:
                out[k] = cleaned
        else:
            out[k] = v
    return out

def _deep_merge(base: dict, override: dict) -> dict:
    """Recursive dict merge. override values win; nested dicts merge."""
    out = dict(base)
    for k, v in override.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out

def build_pipeline_config(
    *,
    config_path: Optional[Path] = None,
    overrides: Optional[dict] = None,
) -> PipelineConfig:
    """
    Construct a PipelineConfig from optional YAML + overrides dict.

    `overrides` must be structured like the YAML (nested per-stage keys).
    """
    base: dict = {}
    if config_path is not None:
        base = load_yaml(config_path)

    if overrides:
        base = _deep_merge(base, _deep_clean(overrides))

    try:
        return PipelineConfig(**base)
    except ValidationError as e:
        msg = "Configuration validation failed:\n"
        for err in e.errors():
            loc = ".".join(str(x) for x in err["loc"])
            msg += f"  - {loc}: {err['msg']}\n"
        raise ConfigError(msg.rstrip()) from e

class _FlowList(list):
    """A list that YAML renders inline (flow style), e.g. [OH, OG, NE2]."""


def _flow_named_rules(rules):
    """Return a copy of a named_interactors_rules dict whose leaf value-lists are
    wrapped in _FlowList so they dump inline, keeping each rule compact and
    easy to copy/edit from run_config.yaml."""
    out = {}
    for role, rule_list in (rules or {}).items():
        if role == 'hb_masses' and isinstance(rule_list, list):
            out[role] = _FlowList(rule_list)
            continue
        new_rules = []
        if isinstance(rule_list, list):
            for rule in rule_list:
                if isinstance(rule, dict):
                    new_rules.append({k: (_FlowList(v) if isinstance(v, list) else v)
                                      for k, v in rule.items()})
                else:
                    new_rules.append(rule)
            out[role] = new_rules
        else:
            out[role] = rule_list
    return out


def save_pipeline_config(cfg: PipelineConfig, path: Path) -> None:
    """Serialise a PipelineConfig to YAML.

    The whole document is block-style, but the leaf value-lists inside
    ``interactors.named_interactors_rules`` are rendered inline so each role's
    rules stay compact and copy/paste-editable.
    """
    data = cfg.model_dump(mode="json")
    if isinstance(data.get("interactors"), dict) and \
            "named_interactors_rules" in data["interactors"]:
        data["interactors"]["named_interactors_rules"] = _flow_named_rules(
            data["interactors"]["named_interactors_rules"]
        )

    dumper = yaml.SafeDumper
    dumper.add_representer(
        _FlowList,
        lambda d, seq: d.represent_sequence(
            "tag:yaml.org,2002:seq", seq, flow_style=True),
    )

    path = Path(path).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        yaml.dump(data, f, Dumper=dumper, default_flow_style=False, sort_keys=False)