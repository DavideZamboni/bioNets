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
from typing import Any, List, Literal, Optional

import yaml
from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator


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
    d_a_dist: float = Field(default=4.5, gt=0, description="D–A cutoff [Å].")
    angle: float = Field(default=150.0, ge=0, le=180, description="D–H–A angle cutoff [°].")
    sb_cutoff: float = Field(default=3.5, gt=0, description="Salt-bridge cutoff [Å].")
    hc_cutoff: float = Field(default=5.5, gt=0, description="Hydrophobic CoM cutoff [Å].")

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
    # ── Frame count (needed to normalise frequencies) ──
    n_frames: Optional[int] = Field(
        default=None, ge=1,
        description="Trajectory frames analysed during Contacts. Required.",
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
    # ── Filtering ──
    min_frequency: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Drop edges with weight (frequency) below this.",
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
        "interactors_type", "targets", "node_metrics", "community_methods",
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

def save_pipeline_config(cfg: PipelineConfig, path: Path) -> None:
    """Serialise a PipelineConfig to YAML."""
    data = cfg.model_dump(mode="json")
    path = Path(path).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        yaml.safe_dump(data, f, default_flow_style=False, sort_keys=False)