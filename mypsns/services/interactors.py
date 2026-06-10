"""
interactors.py – Interactors service module for myPSNs.
 
Bridges Topology_recognition output and the Contacts module by:
 
  1. Validating that every residue type has been successfully associated to a
     forcefield entry (``__init__``).
  2. Reconstructing the full whole-topology table with correct bond info,
     chain assignments, polymer classification, and ter chargegroup
     adjustments (``get_complete_topology``).
  3. Assigning interaction-type labels – PIE and/or PyInteraph – to every
     atom in the topology (``add_inttypes``), producing the
     ``interaction_topology.csv`` file consumed by :class:`Contacts`.
 
Depends on:
  - mypsns/utils/interactors_utilities.py  for all heavy lifting
  - Topology_recognition output folder:
      complete_residue_topologies.csv
      Resindex_topology_association.csv
      topologies/<TopName>.csv  (per-residue ff-parameter CSVs)
"""

from __future__ import annotations
 
import logging
import timeit
from os import listdir, makedirs
from pathlib import Path
from typing import List, Optional, Union
 
import MDAnalysis as mda
import numpy as np
import pandas as pd
 
from mypsns.utils.interactors_utilities import (
    get_charge_based_interactors,
    get_named_interactors,
    reconstruct_whole_topology,
)

from mypsns.utils.general_utilities import get_endtime
 
pd.options.mode.chained_assignment = None

log = logging.getLogger(__name__)
# Default na_values shared across every CSV read in this module.
_NA_VALUES: List[str] = [
    "", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN",
    "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A",
    "NULL", "NaN", "n/a", "nan", "null",
]

class Interactors():
    """
    Assign complete bond, chain, and interaction-type information to every
    atom in the topology from Topology_recognition output.
 
    Parameters
    ----------
    topology : str
        Original topology file (PDB, GRO, …), same as passed to
        :class:`Topology_recognition`.
    general_outpath : str or Path
        Root output folder.  The sub-folder ``Topology_recognition/`` must
        already exist and contain the output of
        :meth:`Topology_recognition.get_ff_parameters`.
    na_values : list[str], optional
        Strings recognised as ``pd.NA`` by pandas.  value ``"NA"``
        omitted to preserve sodium atom names.
 
    Attributes
    ----------
    u : mda.Universe
        MDAnalysis Universe built from *topology*.
    out : Path
        ``general_outpath / 'Topology_recognition'`` – where all input CSVs
        live and where output CSVs are written.
    start : float
        ``timeit.default_timer()`` at construction time.
    nas : list[str]
        Active na_values list.
    """

 
    def __init__(
        self,
        topology: str,
        general_outpath: Union[str, Path],
        na_values: Optional[List[str]] = None,
        ) -> None:

        self.u       = mda.Universe(topology)
        self.inp     = Path(general_outpath) / 'Topology_recognition'
        self.out     = Path(general_outpath) / 'Interactors'
        makedirs(self.out, exist_ok=True)
        self.nas     = na_values if na_values is not None else _NA_VALUES
        self.start   = timeit.default_timer()

        log.info('Interactors run starting; reading Topology_recognition from %s', self.out.absolute())
        # validate topology info from Topology_recognition
        self.validate_topology_recognition()

    def validate_topology_recognition(self) -> None:
        """
        Check that ``complete_residue_topologies.csv`` exists and that every
        residue type has non-null ``ff_atoms`` (=>
        successfully associated to a forcefield entry by
        :meth:`Topology_recognition.get_ff_parameters`).
 
        Raises
        ------
        FileNotFoundError
            If the required CSV is absent.
        ValueError
            If one or more residue types could not be associated to the
            forcefield; lists corresponding residue names.
        """
        csv_path = self.inp / 'complete_residue_topologies.csv'
        if not csv_path.exists():
            msg = (
                f'complete_residue_topologies.csv not found in '
                f'{self.out.absolute()}.  Run '
                f'Topology_recognition.get_ff_parameters() successfully first.'
            )
            log.error(msg)
            raise FileNotFoundError(msg)
 
        topologies = pd.read_csv(
            csv_path, header=0, index_col=False, na_values=self.nas
        )
 
        # ff_atoms column that contains per-residue CSV forcefield data.
        if 'ff_atoms' not in topologies.columns:
            log.warning('ff_atoms column not found in complete_residue_topologies.csv - skipping null check')
            return
 
        unmatched = topologies[topologies['ff_atoms'].isna()]['resname'].tolist()
        if unmatched:
            msg = (
                f'The following residue types have no forcefield association '
                f'(ff_atoms is null): {unmatched}.\n'
                f'Re-run Topology_recognition.get_ff_parameters() and ensure '
                f'all residues are present in the forcefield.'
            )
            log.error(msg)
            raise ValueError(msg)
 
        log.info('Validation passed: all %s residue types have forcefield associations', len(topologies))

    def get_complete_topology(self):
        '''
        Function that aggregates the results from get_ff_parameters.
        In particular:
         - Fills bond information
         - assigns chains
        
        '''  
        
        start = timeit.default_timer()
        log.info('Reconstructing complete topology')

        try:
            self.res_topol = pd.read_csv(self.inp / 'complete_residue_topologies.csv', header=0, index_col = False)

        except FileNotFoundError:
            log.error('complete_residue_topologies.csv not found - run Topology_recognition.get_ff_parameters() first')
            return False

        found_tops = np.unique(self.res_topol['ff_atoms'])
        topologies = [self.inp / 'found_topologies' / top for top in found_tops]

        reconstruct_whole_topology(self.u, topologies, self.res_topol, self.inp, self.out, self.nas)

        self.template_dir   = self.out / 'found_topologies'
        self.residue_index  = self.out / 'residue_index.csv'
        self.atom_index_map = self.out / 'atom_index_map.csv'

        script_time, total_time = get_endtime(start, self.start)
        log.info('Topology reconstructed; get_complete_topology took %s (h:m:s), total %s (h:m:s)',
                 script_time, total_time)
        return True

    def add_inttypes(self,
                    na_values: list = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"],
                    **kwargs
                    ):
        '''
        Function that adds info concerning the capabilities of various atoms to form salt bridges, H bonds or hydrophobic interacions.
        2 types: q_interactors (based on charges and connectivity) and named_interactors (PyInteraph - like, exact translation of default selections in PyInteraph 1)
         - q_interactors support following parameters:
              - H bond parameters:
                - max_DA_charge: maximum (less negative) charge that a donor or acceptor may have
                - min_H_charge: minimum charge of an H atom to be able to form an H bond
                - DA_masses: mass of the atoms that can function as donor/acceptor atoms (N, O, F)
              - SB parameters:
                - min_abs_charge: minimum charge of a chargegroup to find and classify atoms as basic/acidic
              - HC parameters:
                - max_abs_charge: maximum charge of atoms to classify them as hydrophobic
                - max_num_heavy_ats: maximum number of heavy atoms to be included in each HC cluster
                - min_bonds_to_charged: minimum number of bonds separating the possibly HC forming atom from a charged (with partial charge > max_abs_charge) atom
                                        to classify it as HC
         - named_interactors: 
             - H-Bonds:
                 - donors HDon = [('atom_name', ['OH', 'OG', 'NE2', 'OG1', 'NE', 'N', 'ND1', 'NZ', 'NH1', 'NH2', 'OW', 'ND2', 'SG', 'OH2', 'NE1'])]
                 - acceptors HAcc = [('atom_name', ['OH', 'OG', 'OD1', 'OD2', 'OG1', 'O', 'ND1', 'NE2', 'OE2', 'OW', 'SG', 'OE1', 'OH2', 'SD'])]
                 - hb_masses = [1]: mass(es) of atoms bound to HDon atoms to be considered as H possibly involved in H-Bonds
             - Hydrophobic: HC = [[('resname', ['ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'PRO', 'TRP', 'MET', 'IIL']),
                                   ('not atom_name', ['CA', 'C', 'O', 'N', 'H', 'H1', 'H2', 'H3', 'O1', 'O2', 'OXT', 'OT1', 'OT2'])]]
             - Salt-bridges: 
                 - Basic = [
                     [('resname',  ['GLY-NH3+', 'NH3+', 'NHT3+', 'PRO-NH2+']),
                      ('atom_name', ['N', 'H1', 'H2', 'H3', 'HT1', 'HT2', 'HT3'])],
                     [('resname',  ['HIS', 'HSD', 'HSE', 'HSP', 'HID', 'HIE', 'HIP']),
                      ('atom_name', ['CG', 'ND1', 'HD1', 'CD2', 'CE1', 'NE2', 'HE2'])],
                     [('resname',  ['LYS']),
                      ('atom_name', ['CE', 'NZ', 'HZ1', 'HZ2', 'HZ3'])],
                     [('resname',  ['ARG']),
                      ('atom_name', ['CD', 'NE', 'CZ', 'NH1', 'HH11', 'HH12',
                                     'NH2', 'HH21', 'HH22'])],
                         ]
                 - Acidic = [
                     [('resname',  ['ASP']),
                      ('atom_name', ['CG', 'OD1', 'OD2'])],
                     [('resname',  ['GLU']),
                      ('atom_name', ['CD', 'OE1', 'OE2'])],
                     [('resname',  ['COO-']),
                      ('atom_name', ['C', 'O1', 'O2', 'O1T', 'O2T', 'OT1', 'OT2'])],
                 ]            
                              
        NW: the list na_values is used from pandas to identify missing values. the entry NA has been removed from such list to keep sodium atoms.
            If any atom inserted in the forcefield has a name comprised in this list, it is necessary to remove it from it.
        '''
        start = timeit.default_timer()
        log.info('Assigning interaction types')


        template_dir = self.out / 'found_topologies'
        variant_files = sorted(
            p for p in template_dir.iterdir()
            if p.is_file() and p.suffix == '.csv'
        )

        for vpath in variant_files:
            log.debug('Processing variant %s', vpath.stem)

            vdf = pd.read_csv(
                vpath, header=0, index_col=False, na_values=na_values
            )

            # Scaffolding for analysis functions (they groupby Resindex
            # and read Entity). Single-residue variant → both fixed values
            # collapse the loops to a single iteration with no neighbours.
            vdf['Resindex'] = 0
            vdf['Entity']   = 'M'

            vdf = get_charge_based_interactors(vdf, **kwargs)

            if 'HC_cluster' in vdf.columns:
                if not any(isinstance(x, str) and x != ''
                           for x in vdf['HC_cluster']):
                    vdf.drop('HC_cluster', axis=1, inplace=True)

            vdf = get_named_interactors(vdf, **kwargs)

            vdf.drop(columns=['Resindex', 'Entity'],
                     inplace=True, errors='ignore')
            vdf.to_csv(vpath, header=True, index=False)

        script_time, total_time = get_endtime(start, self.start)
        log.info('Interaction types added; add_inttypes took %s (h:m:s), total %s (h:m:s)',
                 script_time, total_time)