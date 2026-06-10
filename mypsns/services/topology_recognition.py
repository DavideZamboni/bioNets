from os import path, makedirs
import pandas as pd
from pathlib import Path
import os
import MDAnalysis as mda
import timeit
import logging
import sys
script_path = path.abspath(__file__)
script_dir = path.dirname(script_path)
import warnings

warnings.filterwarnings(
    "ignore",
    message="Invalid atom serials were present, bonds will not be parsed",
    category=UserWarning,
    module="MDAnalysis.topology.PDBParser"
)

warnings.filterwarnings(
    "ignore",
     message="Serial numbers went over 100,000.  Higher serials have been guessed",
     category=UserWarning,
     module="MDAnalysis.topology.PDBParser"
 )

warnings.filterwarnings(
    "ignore",
    message="Unknown element  found for some atoms. These have been given an empty element record. If needed they can be guessed using universe.guess_TopologyAttrs(context='default', to_guess=['elements']).",
    category=UserWarning,
    module="MDAnalysis.topology.PDBParser"
)

sys.path.insert(0, script_dir + '/scripts')

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from mypsns.utils.topology_recognition_utilities import associate_ff_params, get_residue_types, parse_initial_topology_info
from mypsns.utils.general_utilities import get_endtime
pd.options.mode.chained_assignment = None
log = logging.getLogger(__name__)

class Topology_recognition():
    def __init__(
            self, 
            topology: str,
            outpath: Path, 
            ff_path: Path,
            na_values: list = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"]                ):
        '''
         - topology: topology file (pdb, gro, ...)
         - outpath: output folder
         - na_values: specifies strings recognized as pandas.NA values.
                      feature added to avoid problems with sodium

        Obtains useful information from input topology, producing .csv file with all the possible following info:

        - chainID: chain index
        - resindex: internal (unique and 0-based) residue index
        - resid: residue index information (from topology)
        - resname: residue name
        - index: internal (unique and 0-based) atom index
        - id: atom index information (from topology)
        - name: atom name
        - type: atom type (could be element or forcefield atomtype, depending on input topology file)
        - element: atomic element
        - mass: atomic mass 
        - charge: atom's partial charge 

        NW: obtained info depend on info defined in the topology file provided. For more information, visit https://userguide.mdanalysis.org/stable/topology_system.html#format-specific-attributes

        Defines:
         - self.u: reference Universe
         - self.out: output folder, exactly outpath
         - self.nas: na_values list
        '''

        self.u = mda.Universe(topology)

        # Internalize output folder, create log 
        self.out = Path(outpath)

        # Internalize ff folder
        self.ff = ff_path


        # 1 - Obtain information from current topology
        parse_initial_topology_info(self.u, self.out)
    
        if not hasattr(self, 'start'):
            self.start = timeit.default_timer()

        self.nas = na_values
        log.info('Topology_recognition run starting; initial topology info written to %s',
                 (self.out / 'initial_topology_info.csv').absolute())
      
    def get_restypes(self):
        '''
        Produces a .csv file containing a table containing all the possible topologies found for each of the residue names in universe. E.g.:

         resname    atnames             rixs
             ALA    ALA_atoms.csv       ALA_rixs.txt
             ARG    ARG_atoms.csv       ARG_rixs.txt
             ASN    ASN_atoms.csv       ASN_rixs.txt
             ASP    ASP_atoms.csv       ASP_rixs.txt

        where:
         - resname: resname found in the topology file
         - atnames: atom names in the residue(s) indicated in the corresponding index list
         - rixs: list of residue indices containing the residues constituted exactly by the corresponding list of atom names

        Defines
        
         - self.res_topology_folder: self.out / resid_topologies, 2 relevant file types resulting:
            - residuename(_topId)_rsix.txt -> with all residue indices for the respective topology
            - residuename(_topId)_atoms.csv -> csv file with atom names and masses
            NW: (_topId): numer of different topologies (different number/definition of atoms for each residue name)
                          present in the filename only if it is > 0

         - self.restypes: self.out / different_residue_topologies.csv
       '''

        start = timeit.default_timer()
        log.info('Detecting residue types')

        self.res_topology_folder = self.out / "resid_topologies"
        makedirs(self.res_topology_folder, exist_ok = True)

        # Find all different types of residues associated with the same resname:
        restypes = pd.DataFrame.from_dict({'resname':[], 'atnames':[], 'rixs':[], 'masses':[]})
        top_df = pd.read_csv(self.out / 'initial_topology_info.csv', header = 0, index_col = False)
       
        restypes = get_residue_types(top_df, self.res_topology_folder)

        self.restypes = restypes
        self.restypes.to_csv(self.out / 'different_residue_topologies.csv', header = True, index = False)
        proc_time, tot_time = get_endtime(start, self.start)
        log.info('Residue-type detection complete; get_restypes took %s (h:m:s), total %s (h:m:s)',
                 proc_time, tot_time)

    def get_ff_parameters(self, 
                          maxToCheck: int = -1
                         ):
        '''
        Function that associate the unique residue topologies found to specific ff entries

         - verbose: whether to write report file of whole process
         - maxToCheck: max number of residues to check from ff, default (-1) = check all.
                       residue names in forcefield are sorted by similarity to topology residue name

        '''
        # To associate unique restypes to ff entries
        start = timeit.default_timer()

        associate_ff_params(self.ff, self.out, self.nas, maxToCheck)

        proc_time, tot_time = get_endtime(start, self.start)
        log.info('Residue-to-forcefield association complete; get_ff_parameters took %s (h:m:s), total %s (h:m:s)',
                 proc_time, tot_time)