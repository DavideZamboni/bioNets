"""
contacts.py – Contacts service module for myPSNs.

Orchestrates the parallel computation of:
  - Hydrogen bonds (H-bonds)
  - Salt bridges
  - Hydrophobic clusters (HCs)

across all chains and chain-pairs of a molecular dynamics trajectory.

Depends on:
  - mypsns/utils/contacts_utilities.py  for all low-level logic
  - Topology_recognition output folder   (interaction_topology.csv produced
                                          by interactors.py)
"""

from __future__ import annotations

import datetime
import logging
import multiprocessing as mp
import timeit
from logging.handlers import QueueListener
from os import makedirs
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import MDAnalysis as mda
import numpy as np
import pandas as pd

from mypsns.utils.contacts_utilities import (
    generate_combinations,
    get_selection,
    hydrophobics,
    Hbonds,
    init_defaults,
    load_inputs,
    salt_bridges,
    split_frames,
    startContacts,
    build_results,
    get_interaction_energy,
    _SUBFOLDER,
    build_persistency_matrices,
    threshold_analysis,
    load_lj_params,
    energy_to_percentage,
    get_hc_energy,
)
log = logging.getLogger(__name__)

class Contacts:
    """
    Orchestrates contact calculations (H-bonds, salt bridges, HCs) over an
    MD trajectory.

    Parameters
    ----------
    topology : str
        Path to a topology file compatible with MDAnalysis (PDB, GRO, …).
    trj_file : str
        Path to the trajectory file (XTC, DCD, …).
    outpath : str or Path
        Root output folder.  Sub-folders ``Adjacency/`` (results) and the
        expected ``Topology_recognition/`` folder are resolved from here.
    **kwargs
        Any key from ``self.defaults`` can be overridden at construction time.

    Attributes
    ----------
    defaults : dict
        Runtime parameters; see below for the full list.
    inters : pd.DataFrame
        Interaction-topology table loaded from
        ``Topology_recognition/interaction_topology.csv``.
    top : Path
        Path to the corrected topology PDB written by :func:`load_inputs`.

    Default parameters
    ------------------
    Computation
        - ``nc`` (2) – number of worker processes
        - ``ram_thresh`` (80) – RAM-usage ceiling (%) before flushing to disk
        - ``nframes`` (500) – frames per H-bond batch

    Trajectory window
        - ``beg_frm`` (0) – first frame to analyse
        - ``end_frm`` (-1) – last frame (-1 = all)

    Interaction-type selectors
        - ``interactors_type`` (['named_interactors']): PyInteraph - like (name-based) intractors selection, or 'q_interactors': charge-, connectivity-based selection (not fully developed)
        - ``targets`` (['HBonds', 'Salt_bridges', 'HCs'])
        - ``subject_chains`` (list of polymeric chains (Entity = 'P' in inters) in topology)
        - ``surrounding_chains`` (all chains in topology)

    Geometric cutoffs
        - ``d_a_dist`` (5 Å) – H-bond donor–acceptor distance
        - ``angle`` (150°) – H-bond D–H–A angle
        - ``sb_cutoff`` (5 Å) – salt-bridge distance
        - ``hc_cutoff`` (5 Å) – hydrophobic-cluster distance
    """

    def __init__(
        self,
        topology: str,
        trj_file: str,
        outpath: Union[str, Path],
        na_values: Optional[List[str]] = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"],
        **kwargs,
    ) -> None:

        output_root = Path(outpath)
        self.inp  = output_root / "Interactors"
        self.outp = output_root / "Adjacency"
        self.nas = na_values

        makedirs(self.outp, exist_ok=True)

        # ── Logging ──────────────────────────────────────────────────────────
        # Don't call basicConfig: the CLI already configured the root logger.
        # Attach a dedicated per-stage file handler alongside the CLI handlers.
        _fmt = logging.Formatter(
            "%(asctime)s %(processName)s %(levelname)s: %(message)s"
        )
        self._stage_fh = logging.FileHandler(self.outp / "Contacts.log", mode="w")
        self._stage_fh.setLevel(logging.DEBUG)
        self._stage_fh.setFormatter(_fmt)
        logging.getLogger().addHandler(self._stage_fh)

        # Workers can't share a handler => push records through a queue;
        # the listener replays them into the parent's existing handlers.
        self.log_queue = mp.Queue(-1)
        self.listener  = QueueListener(
            self.log_queue, *logging.getLogger().handlers,
            respect_handler_level=True,
        )
        self.listener.start()
 
        self.start = timeit.default_timer()

        # ── Pre-run check ─────────────────────────────────────────────────
        stop, msg, name = startContacts(topology, self.inp)
        log.info(msg)

        if stop:
            log.error(
                "Topology_recognition must complete successfully "
                "before running Contacts."
            )
            return

        # ── Load inputs ───────────────────────────────────────────────────────
        self.inters, self.top, self.trj = load_inputs(
            self.inp, self.outp, name, topology, trj_file, self.nas
        )

        poly_chains = list(np.unique(self.inters[self.inters['Entity'] == 'P']["Chain"].tolist()))
        found_chains = list(np.unique(self.inters["Chain"].tolist()))
        self.defaults = init_defaults(poly_chains, found_chains, self.log_queue)
        self.defaults.update(kwargs)
        self.defaults.update({k: v for k, v in kwargs.items() if v is not None})

        for k in self.defaults.keys():
            print(k, self.defaults[k])
        if not self.defaults.get("subject_chains"):
            self.defaults["subject_chains"] = poly_chains
        if not self.defaults.get("surrounding_chains"):
            self.defaults["surrounding_chains"] = found_chains
        log.info("Subject chains: %s | surrounding chains: %s",
                 self.defaults["subject_chains"], self.defaults["surrounding_chains"])
        # ── Build universe and split frames ───────────────────────────────────
        u = mda.Universe(self.top.resolve(), self.trj.resolve())
        myframes: List[List[int]] = split_frames(u, **self.defaults)
        self.frames = {frame for framelist in myframes for frame in framelist}

    # ──────────────────────────────────────────────────────────────────────────

    def run(
        self
    ) -> None:
        """
        Execute contact calculations over the trajectory.
        """

        log.info(
            f"\nContacts.run starting at "
            f"{datetime.datetime.now().strftime('%I:%M%p on %B %d, %Y')}\n"
            f"Topology          : {self.top}\n"
            f"Interaction table : {self.inp / 'interaction_topology.csv'}\n"
            f"Requested types   : {self.defaults['interactors_type']}\n"
        )

        # ── Build universe and split frames ───────────────────────────────────
        u = mda.Universe(self.top.resolve(), self.trj.resolve())
        myframes = split_frames(u, **self.defaults)

        # ── Chain combinations ────────────────────────────────────────────────
        combos = generate_combinations(
            self.defaults["subject_chains"],
            self.defaults["surrounding_chains"],
        )
        self.defaults["chain_combinations"] = combos
        log.info(
            "Chain combinations (subject | surroundings):\n"
            + "\n".join(f"  {c[0]:^7} | {' '.join(c[1:])}" for c in combos)
        )

        # ── Per interaction-type loop ─────────────────────────────────────────
        for inttype in self.defaults["interactors_type"]:
            if inttype not in self.inters.columns:
                log.warning(
                    f"Interaction type '{inttype}' not found in "
                    f"{self.inp / 'interaction_topology.csv'} – skipping."
                )
                continue

            log.info(f"Processing interaction type: {inttype}")
            myout = Path(self.outp / inttype)
            makedirs(myout, exist_ok=True)
            
            seles = self.build_selections(inttype)
            log.info('Selections built: %s', seles)
            log.info(f'Will run contacts with parameters:\ntopology: {self.top}')
            log.info(f'trajectory: {self.trj.resolve()}')
            boundary_str = [f"{f[0]} to {f[-1]}" for f in myframes]
            log.info(f'frames to be analyzed: {boundary_str}')
            
            log.info(f'SB output folder: {(myout / "Salt_bridges").absolute()}')
            log.info(f'HB output folder: {(myout / "H_bonds").absolute()}')
            log.info(f'HC output folder: {(myout / "HCs").absolute()}')

            if "Salt_bridges" in self.defaults["targets"]:
                sbout = myout / "Salt_bridges"
                makedirs(sbout, exist_ok=True)
                toremove = [f for f in sbout.iterdir() if f.name.endswith('.csv')]
                if len(toremove) > 0:
                    log.info("\n\nATTENTION: Results of previous runs found. Such files will be removed")
                for f in toremove:
                    f.unlink()

                log.info("Computing salt bridges …")
                salt_bridges(
                    self.top.resolve(), self.trj.resolve(), seles, sbout,
                    myframes, **self.defaults,
                )

            if "HBonds" in self.defaults["targets"]:
                hbout = myout / "H_bonds"
                makedirs(hbout, exist_ok=True)
                toremove = [f for f in hbout.iterdir() if f.name.endswith('.csv')]
                if len(toremove) > 0:
                    log.info("\n\nATTENTION: Results of previous runs found. Such files will be removed")
                for f in toremove:
                    f.unlink()

                log.info("Computing H-bonds …")
                Hbonds(
                    self.top.resolve(), self.trj.resolve(), seles, hbout, myframes,
                    **self.defaults,
                )

            if "HCs" in self.defaults["targets"]:
                hcout = myout / "HCs"
                makedirs(hcout, exist_ok=True)
                toremove = [f for f in hcout.iterdir() if f.name.endswith('.csv')]
                if len(toremove) > 0:
                    log.info("\n\nATTENTION: Results of previous runs found. Such files will be removed")
                for f in toremove:
                    f.unlink()

                log.info("Computing hydrophobic clusters …")
                hydrophobics(
                    self.top.resolve(), self.trj.resolve(), seles, hcout, myframes,
                    **self.defaults,
                )

        #self.listener.stop()

    def get_matrices(self, **kwargs) -> None:
        """
        Run the post-contact analyses selected by ``matrix_type``
        ('persistency' and/or 'energy'), then tear the log listener down once.
        Call after :meth:`run`.
        """
        self.defaults.update(kwargs)
        types = self.defaults.get("matrix_type", ["persistency", "energy"])
        try:
            if "persistency" in types:
                self.get_persistence_analysis()
            if "energy" in types:
                self.get_energy_analysis()
        finally:
            if getattr(self, "listener", None) is not None:
                self.listener.stop()
                self.listener = None
    
    def get_persistence_analysis(self,
                        **kwargs
                        ) -> None:
        self.defaults.update(kwargs)
        start = timeit.default_timer()

        log.info(
            f"\nContacts.get_persistence starting at "
            f"{datetime.datetime.now().strftime('%I:%M%p on %B %d, %Y')}\n"
        )

        def_foldernames = [i.split('_')[0] for i in self.defaults['interactors_type']]
        inters_done = [f for f in self.outp.iterdir() if f.is_dir() and f.name in def_foldernames]

        if len(inters_done) == 0:
            log.warning(
                f"\nNeed to determine contacts for at least an interaction type\n"
            )
            exit

        for intfolder in inters_done:
            intfiles: List[Path] = []
            expected: List[str] = [_SUBFOLDER[i] for i in self.defaults['targets']]

            inttypes = [t for t in intfolder.iterdir() if t.is_dir() and t.name in expected]
            log.info(f'Expected following interaction types: {expected}')
            log.info(f'Found following interaction type folders: {inttypes}')
            
            if len(inttypes) == 0 or len(inttypes) != len(expected):
                log.warning(f"No interaction type found for interactors: {intfolder.name}")
                break

            for inttype in inttypes:
                merged = [merged_csv for merged_csv in inttype.iterdir() if merged_csv.is_file() and merged_csv.name.endswith('_merged.csv')]
    
                if len(merged) == 0:
                    log.warning(f"No merged results found for interactors: {intfolder.name}, interaction type: {inttype.name}")
                    break
                intfiles.extend(merged)
            
            log.info(f'The following files will be used for results computation: {[i.resolve() for i in intfiles]}')
            build_results(intfiles, self.inters, self.frames, intfolder)
            
            # ── adjacency matrices (per-type + union-combined) ────────────
            n_frames = len(self.frames)
            final_csvs   = {}
            visited_csvs = {}
            for label, fname in (('SB', 'SB'), ('HC', 'HC'), ('Hbonds', 'Hbonds')):
                f = intfolder / f'{fname}_final.csv'
                v = intfolder / f'{fname}_final_visited_frames.csv'
                if f.exists():
                    final_csvs[label]   = f
                    visited_csvs[label] = v

            if not final_csvs:
                log.warning(f'{intfolder.name}: no *_final.csv → skip adjacency')
                continue

            log.info(f'Building adjacency matrices in {intfolder}')
            _, matrices, self.labels = build_persistency_matrices(
                final_csvs, visited_csvs, self.inters, n_frames, intfolder, na_values=self.nas
            )

            # ── threshold scan on per-type matrices ───────────────────────
            per_type = {k: v for k, v in matrices.items()}
            log.info(f'Running threshold analysis on subsets of {list(per_type)}')
            threshold_analysis(
                per_type,
                outdir=intfolder / 'graph_report',
                step=self.defaults.get('thr_step', 0.1),
                nc=self.defaults.get('nc', 1),
            )

        #self.listener.stop()

    def get_energy_analysis(self, **kwargs) -> None:
        '''
        Per-contact LJ + Coulomb energy matrices from the SB / Hbond merged CSVs.
        Runs after the full pipeline (needs forcefield/ + Adjacency/<inttype>/).
        '''
        self.defaults.update(kwargs)
        log.info("\nContacts.get_energy_analysis starting at "
                 f"{datetime.datetime.now().strftime('%I:%M%p on %B %d, %Y')}\n")

        ff_dir = self.outp.parent / 'forcefield'          # <outpath>/forcefield
        if not (ff_dir / 'atomtypes_se.csv').exists():
            log.error("atomtypes_se.csv not found in %s – run ff_reader first.", ff_dir)
            return

        n_frames  = len(self.frames)
        lj_lookup = load_lj_params(
            ff_dir, comb_rule=self.defaults.get('comb_rule', 2), na_values=self.nas
        )

        def_foldernames = [i.split('_')[0] for i in self.defaults['interactors_type']]
        inters_done = [f for f in self.outp.iterdir()
                       if f.is_dir() and f.name in def_foldernames]
        if not inters_done:
            log.warning("No interaction-type folders found under %s", self.outp)
            return

        for intfolder in inters_done:
            sb_dir = intfolder / _SUBFOLDER['Salt_bridges']   # Salt_bridges/
            hb_dir = intfolder / _SUBFOLDER['HBonds']          # H_bonds/
            hc_dir   = intfolder / _SUBFOLDER['HCs']        # HCs/
            sb_files = sorted(sb_dir.glob('*_merged.csv')) if sb_dir.exists() else []
            hb_files = sorted(hb_dir.glob('*_merged.csv')) if hb_dir.exists() else []
            hc_files = sorted(hc_dir.glob('*_merged.csv')) if hc_dir.exists() else []

            if not sb_files and not hb_files and not hc_files:
                log.warning("%s: no SB/Hbond/HC merged CSVs – skip energy.", intfolder.name)
                continue

            log.info("%s: energy from %d SB + %d HB merged files",
                     intfolder.name, len(sb_files), len(hb_files))
            get_interaction_energy(
                sb_files, hb_files, self.inters, n_frames, intfolder, lj_lookup,
                na_values=self.nas,
            )

            if hc_files:
                inttype = next((t for t in self.defaults['interactors_type']
                                if t.split('_')[0] == intfolder.name), None)
                if inttype is None:
                    log.warning("HC energy: can't resolve interactors_type for %s – skip HC",
                                intfolder.name)
                else:
                    seles = self.build_selections(inttype)
                    hc_label_sel: Dict[str, str] = {}
                    for _ch, roles in seles.items():
                        hc = roles.get('HC')
                        if isinstance(hc, dict):
                            for label, (sel_str, _) in hc.items():
                                if sel_str is not None:
                                    hc_label_sel[label] = sel_str
                    get_hc_energy(hc_files, self.inters, hc_label_sel, self.top, self.trj,
                                  n_frames, intfolder, lj_lookup, na_values=self.nas)

                        # ── energy → strength% transform (+ distribution plot), then elbow ──
            emats: Dict[str, np.ndarray] = {}
            for kind in ('SB', 'Hbonds', 'HC'):
                p = intfolder / f'{kind}_total_energy.dat'
                if p.exists():
                    m = np.loadtxt(p)
                    if m.ndim == 2 and m.size:
                        emats[kind] = m
            if emats:
                report = intfolder / 'energy_graph_report'
                use_abs = self.defaults.get('energy_use_abs', False)
                step    = self.defaults.get('thr_step', 0.1)
                nc      = self.defaults.get('nc', 1)

                # total interaction energy across all types (element-wise sum)
                total_all = sum(emats.values()) if len(emats) > 1 else None
                if total_all is not None:
                    np.savetxt(intfolder / 'all_total_energy.dat', total_all, fmt='%.6e', delimiter=' ')

                # per-type subsets
                pct = energy_to_percentage(emats, outdir=report, use_abs=use_abs, plot=True)
                threshold_analysis(pct, outdir=report, step=step, nc=nc,
                                   label='energy', xlabel=r'strength (%)')

                # all-combined, analysed on its own (NOT folded into the subset combos)
                if total_all is not None:
                    pct_all = energy_to_percentage({'all': total_all}, outdir=report,
                                                   use_abs=use_abs, plot=True)
                    threshold_analysis(pct_all, outdir=report, step=step, nc=nc,
                                       label='energy', xlabel=r'strength (%)')
                        
    # ──────────────────────────────────────────────────────────────────────────
    # Private helpers
    # ──────────────────────────────────────────────────────────────────────────

    # ──────────────────────────────────────────────────────────────────────────
    # Private helpers
    # ──────────────────────────────────────────────────────────────────────────

    def build_selections(
        self,
        inttype: str,
    ) -> Dict[str, Dict[str, Union[Tuple[str, int], Dict[str, Tuple[str, int]]]]]:
        """
        Build MDAnalysis selection strings for every chain and requested
        interaction role (HBond, HAcc, HDon, Acidic, Basic, HC).

        Returns
        -------
        dict
            ``{chain_id: {role: (sel_str, expected_atom_count)}}``, except for
            HC which carries a per-cluster sub-dict
            ``{label: (sel_str, expected_atom_count)}``. Chains/roles with no
            tagged atoms are omitted.
        """
        targets = self.defaults["targets"]
        all_chains = sorted(set(self.defaults["surrounding_chains"]) |
                            set(self.defaults["subject_chains"]))

        role_map: Dict[str, List[str]] = {
            "HBonds":       ["HBond", "HAcc", "HDon"],
            "Salt_bridges": ["Acidic", "Basic"],
            "HCs":          ["HC"],
        }
        roles_needed: List[str] = []
        for t in targets:
            roles_needed.extend(role_map.get(t, []))

        seles: Dict[str, Dict[str, Union[Tuple[str, int],
                                         Dict[str, Tuple[str, int]]]]] = {}
        for chain in all_chains:
            chain_df = self.inters[self.inters["Chain"] == chain]
            if chain_df.empty:
                continue

            chain_seles: Dict[str, Union[Tuple[str, int],
                                         Dict[str, Tuple[str, int]]]] = {}
            for role in roles_needed:
                sel = get_selection(chain_df, inttype, role, chainID=chain)
                if sel is None:
                    log.info(f"  Chain {chain} | {role}: no atoms found")
                    continue
                chain_seles[role] = sel
                log.info(f"  Chain {chain} | {role}: {sel}")

            if chain_seles:
                seles[chain] = chain_seles

        return seles