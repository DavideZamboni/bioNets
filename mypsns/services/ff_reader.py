
## Parse atoms and things from Charmm36-jul2022.ff
from os import makedirs
from pathlib import Path
from typing import Optional, List
from pathlib import Path
import logging

## import required functions from ff_reader_utilities
from mypsns.utils.ff_reader_utilities import get_residues, get_ters, get_sigma_eps_mass, merge_data, add_residue
from mypsns.resources import load_files

def ff_reader(
              inpath: str, 
              outpath: str, 
              na_values: Optional[List[str]] = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"],
              nonbonded_name: Optional[str] = "ffnonbonded.itp",
              additional_residues: Optional[bool] = True
                           ) -> None:
    
    """
    Function required to get forcefield information on all the defined residues
    The information that will be associated with each of the atoms defined in the forcefield are:
    - class: molecular class of the defined residue (silicates, aminoacids, lipids, ethers, nucleic acids, ...)
    - resname: residue name (must be unique in the forcefield)
    - atom_name: name of the atom in the topology
    - atom_type: type of atom (defining charge, sigma and eps parameters)
    - charge: atom's patial charge
    - chargegroup: atom's chargegroup within the residue
    - bonds: atom names of the atoms bound to the resiude

    2 principal output files:
    - final_ff_residues.csv:    file with all the information of each residue name found in the forcefield .rtp files
    - final_ff_ters.csv:        file with all the information of each ter entry found in the forcefield .tdb files
                                this has the following nomenclature in atom bonds: 
                                    -AN: AtomName to be removed from normal topology
                                    rAN: AtomName to be replaced from normal topology with the one at the precise row
    """
    path = Path(inpath)
    out = Path(outpath)

    
    outf = out / 'forcefield'

    makedirs(outf, exist_ok=True)

    logging.info(f'Analyzing forcefield: {path.name}\nStarting with residue topology files (.rtp)\n\n')
    get_residues(path, outf)
    logging.info(f'\n\nDone with residue topology files, results can be found in the output folder {out.name}/{outf.name}, with filename residues.csv\n\nMoving to termini database (.tdb) files\n\n')
    get_ters(path, outf)
    logging.info(f'\n\nDone with termini database (.tdb) entries, results can be found in the output folder {out.name}/{outf.name}, with filename ters.csv\n\nRetriving sigma, epsilon and mass constants for each atom type from ffnonbonded.itp file\n\n')
    get_sigma_eps_mass(path, outf, nonbonded_name = nonbonded_name)
    logging.info(f'\n\nConstants from .itp files have been obtained, results can be found in the output folder {out.name}/{outf.name}, with filenames:\n- atomtypes_se.csv for sigma, epsilon and mass data\n- pairtypes.csv for specific pair type interactions\nMerging data\n')
    merge_data(outf, na_values = na_values) 
    logging.info(f'\n\nForcefield table built.\nResults can be found in the output folder: {out.name}/{outf.name}\nFilename: merged_data.csv\n')

    if additional_residues:
        add_residue(load_files(na_values = na_values), outf)

## NW PANDAS DATAFRAMES DONT LIKE ENTRIES LIKE NA (HEME AND SOD ENTRIES)
