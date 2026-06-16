## Parse atoms and things from Charmm36-jul2022.ff
import numpy as np
import pandas as pd
from os import makedirs
from pathlib import Path
import datetime
from typing import Optional, Union, List, Tuple #, Literal, Iterable
from pathlib import Path
import logging
import re
log = logging.getLogger(__name__)

## Helper functions
def get_index(atdict: dict, atom_name: str, resname: str) -> Union[int, None]:
    '''
    given the dictionary of parsed residues, an atom name and its residue name,
    return the index of such atom in the dictionary.

    Used to add bond information after the parsing of the atoms in such residue.
    
    parsed residues is a dictionary with values:
    - class: [...]
    - resname: [...]
    - atom_name: [...]
    - atom_type: [...]
    - charge: [...]
    - chargegroup: [...]
    - bonds: [...]

    all list have the same length   
    '''
    indices = [-i for i in range(-(len(atdict['resname'])-1), 0, 1)]
    indices.append(0)

    for i in indices:
        if (atdict['resname'][i] == resname) and (atdict['atom_name'][i] == atom_name):
            return i

   #log.debug('NOT FOUND: ', atom_name, resname)
    return None

def get_float(datastr: str) -> Union[str, None]:
    '''
    Removes the termini (e.g. \\n) from a string that should be a float
    input: string, 
    output: trimmed string or None if it was not possible to convert the string to float

    '''
    conv = False
    l = len(datastr)
    i = 1

    while not conv and i < l - 1:
        try:
            float(datastr[:-i])
            conv = True
            return datastr[:-i]
        
        except:
            i += 1
    return None

def get_classes(infile: Path):

    separators = ['\n', '\t']
    myclasses = {}
    with open(infile, 'r') as f:
        l = [ll for ll in f.readlines() if len([k for k in ll.split(' ') if k not in ["", "\n"]]) == 2]
        for ll in l:
            for s in separators:
                    if s in ll:
                        ll = ll.replace(s, ' ')
            info = [rns for rns in ll.split(' ') if rns != '']
            if len(info) == 2:
                if info[0] not in myclasses.keys():
                    myclasses[info[0]] = info[1]
    return myclasses

## rtp reader funcions 

def read_atom_entry(residues: dict, ll: List[str], residue_name: str, rtp_file: str) -> dict:
    # always have 4 entries: atom name, atom type, partial charge and chargegroup
    if len(ll) == 4:
        residues["rtp_file"].append(rtp_file)
        residues['resname'].append(residue_name)
        residues['atom_name'].append(ll[0])
        residues['atom_type'].append(ll[1])
        residues['charge'].append(ll[2])
        residues['chargegroup'].append(ll[3])
        # initialize the bond field for each atom
        residues['bonds'].append('')
    return residues

def read_bond_entry(residues: dict, ll: List[str], residue_name: str) -> dict:
    if len(ll) == 2:
        at1 = ll[0]
        at2 = ll[1]
        
        ## If the atom name starts with + or -, it indicates that it is bound to an atom from the precedent (-) or successive (+) residue
        if at1[0] not in ['+', '-']:
            ix = get_index(residues, at1, residue_name)
            if ix is not None:
                residues['bonds'][ix] += f'{at2} '
            else:
                log.warning(f'Warning: in bonding information of residue {residue_name}, the following atom name has not been found: {at2}')

        if at2[0] not in ['+', '-']:
            ix = get_index(residues, at2, residue_name)
            if ix is not None:
                residues['bonds'][ix] += f'{at1} '
            else:
                log.warning(f'Warning: in bonding information of residue {residue_name}, the following atom name has not been found: {at2}')

    return residues

## rtp reader

def read_rtp_file(path: Path, ff: str, residues: dict) -> dict:
    
    # to insert correctly the resname 
    residue_name = None
    
    # Flags to understand if reading bond or atom data
    atoms = False
    bonds = False

    rtp_file = ff.split('.')[0]
    with open(path / ff, 'r') as f:
        mylines = f.readlines()
        # as ; are comments
        lines = [l for l in mylines if (len(l) > 0 and not l.startswith(';'))]
        for l in lines:
            if ';' in l:
                # remove comments
                idx = l.index(';')
                l = l[0 : idx] 
            
            if l.startswith('[ '):
                # then defining something after
                resname = l.split()[1]
                if resname != 'bondedtypes':
                    residue_name = resname
                    atoms = False
                    bonds = False
            
            elif 'bonds' in l:
                atoms = False
                bonds = True
                
            elif 'atoms' in l:
                atoms = True
                bonds = False

            else:
                if residue_name != None:
                    ll = l.split()
                    if len(ll) > 0:
                        if atoms:
                            residues = read_atom_entry(residues, ll, residue_name, rtp_file)
                        elif bonds:
                            residues = read_bond_entry(residues, ll, residue_name)
                            
    return residues

## tdb reader functions

def parse_ter_definition(line: str, replace: bool, add: bool, delete: bool, 
                         toremove: str, replace_count: int, myTername: str
                         )  -> Tuple[
                             bool, bool, bool, str, int, str]:
    if line.lower().startswith("[ replace ]"):
                        replace = True
                        add = False
                        delete = False
                        replace_count += 1
                    
    elif line.lower().startswith("[ add ]"):
        replace = False
        add = True
        delete = False

    elif line.lower().startswith("[ delete ]"):
        replace = False
        add = False
        delete = True
        toremove = ''

    elif line.lower().split()[1] in ["bonds", "angles", "dihedrals", "impropers"]:
        replace = False
        add = False
        delete = False
        toremove = ''
    
    else:
        myTername = line.split()[1]  
        replace_count = 0

    return replace, add, delete, toremove, replace_count, myTername

def read_replace5_entries(ll:List[str], residues: dict[str, List[str]], myfname: str, myTername: str, toremove: str) -> dict:
    for k in residues.keys():
        log.debug(f"{k}, {residues[k]}")
    
    residues["tdb_file"].append(myfname)
    residues['resname'].append(myTername)
    residues['atom_name'].append(ll[1])
    residues['atom_type'].append(ll[2])
    residues['mass'].append(ll[3])
    residues['charge'].append(ll[4])
    residues['chargegroup'].append(-1)
    # initialize the bond field for each atom
    residues['bonds'].append(f"r{ll[0]} ")
    if len(toremove) > 0:
        log.debug(f'len(ll) == 5, len(toremove) > 0, {len(toremove)}, {toremove}\n')
        residues['bonds'][-1] += toremove
    return residues

def read_replace4_entries(ll:List[str], residues: dict, myfname: str, myTername: str, toremove: str) -> dict:
    residues["tdb_file"].append(myfname)
    residues['resname'].append(myTername)
    residues['atom_name'].append(ll[0])
    residues['atom_type'].append(ll[1])
    residues['mass'].append(ll[2])
    residues['charge'].append(ll[3])
    residues['chargegroup'].append(-1)
    # initialize the bond field for each atom
    residues['bonds'].append(f"r{ll[0]} ")

    if len(toremove) > 0:
        log.debug(f'll: {ll}\nlen(ll) == 4, len(toremove) > 0, {len(toremove)}, {toremove}\n')
        residues['bonds'][-1] += toremove
    return residues

def replace_new_atoms(ll:List[str], residues: dict, myTername: str) -> dict:
    # replace newly added atoms. This is done e.g. for COOH termini, where there is C type replacement, terminal O and H addition, terminal O bound to H type replacement
    # Thus, need to modify data stored in the dict
    ix = get_index(residues, ll[0], myTername)
    #log.debug(myTername)
    #log.debug(ll)
     #'log.debug(pd.DataFrame.from_dict(residues).head())'
    if len(ll) == 5:
        residues["atom_name"][ix] = ll[1]
        residues['atom_type'][ix] = ll[2]
        residues['mass'][ix] = ll[3]
        residues['charge'][ix] = ll[4]

    elif len(ll) == 4:
        residues['atom_type'][ix] = ll[1]
        residues['mass'][ix] = ll[2]
        residues['charge'][ix] = ll[3]
    return residues

def read_replace_entry(ll:List[str], replace_count: int, residues:dict, myfname: str, myTername: str, toremove: str) -> dict:
    """
    GROMACS Documentation, Release 2021, pag. 404

    Replace an existing atom by one with a different atom type, atom name, charge, and/or mass.
    This entry can be used to replace an atom that is present both in the input coordinates and in
    the rtp (page 429) database, but also to only rename an atom in the input coordinates such that
    it matches the name in the force field. 

    In the latter case, there should also be a corresponding [ add ] section present that gives 
    instructions to add the same atom, such that the position in the sequence and the bonding is known.
    Such an atom can be present in the input coordinates and kept, or not present and constructed by
    pdb2gmx (page 128). 
    
    For each atom to be replaced on line should be entered with the following fields:
    – name of the atom to be replaced
    – new atom name (optional)
    – new atom type
    – new mass
    – new charge

    [ replace ]
    ll[0]  ll[1]     ll[2]       ll[3]     ll[4]
    C5'    CN9       12.011000  -0.2700    -      # No new name assigned
    HT     H5'       HN9        1.008000   0.0900 # new name assigned

    """
    if replace_count == 1:
        if not 'X' in ll[0]:
            if len(ll) == 5:
                residues = read_replace5_entries(ll, residues, myfname, myTername, toremove)

            elif len(ll) == 4:
                residues = read_replace4_entries(ll, residues, myfname, myTername, toremove)

    elif replace_count > 1:
        # replace newly added atoms. This is done e.g. for COOH termini, where there is C type replacement, terminal O and H addition, terminal O bound to H type replacement
        # Thus, need to modify data stored in the dict
        residues = replace_new_atoms(ll, residues, myTername)
    
    return residues

def get_numfloats(ll:List[str]) -> int:
    num_floats = 0
    for entry in ll:
        try:
            float(entry)
            num_floats += 1
        except:
            True
    return  num_floats

def add_bonds_ctrl_atoms(residues:dict, myTername: str, 
                         toadd: List, atname: str) ->  dict:
    # adding 1st control atom i
    if toadd[2] not in residues['bonds'][-1]:
        residues['bonds'][-1] += f' {toadd[2]}'
        log.debug(f'| {residues["tdb_file"][-1]:^7} | {residues["resname"][-1]:^7} | {residues["atom_name"][-1]:^9} | {residues["atom_type"][-1]:^9} | {residues["charge"][-1]:^6} | {residues["chargegroup"][-1]:^11} | {residues["bonds"][-1]:^11} |\n')
    # Then need to add this atom to the bonds of the 1st control atom i, 
    # also check to add control atom j to bonds of i.
    # If i not found, store info about bond i-j into newly created atom as +j
    ix = get_index(residues, toadd[2], myTername)
    if ix is not None:
        if atname not in residues['bonds'][ix]:
            residues['bonds'][ix] += f' {atname}'

        if toadd[3] not in residues['bonds'][ix]:
            residues['bonds'][ix] += f' {toadd[3]}'
    else:
        if toadd[3] not in residues['bonds'][-1]:
            residues['bonds'][-1] += f' +{toadd[3]}'
        
    # adding 1st control atom i to bond information of 2nd control atom j
    ix_j = get_index(residues, toadd[3], myTername)
    if ix_j is not None:
        if toadd[3] not in residues['bonds'][ix_j]:
            residues['bonds'][ix_j] += f' {toadd[2]}'

    return residues

def add_atoms(ll:List[str], residues:dict, myfname: str, myTername: str, 
              toremove: str, toadd: List) -> dict:
    '''
    In toadd: [number atoms to add, name of the new atom, control atom i, control atom j] (1st line of add definition)
    In ll: 2nd line of add definition: [new atom type, mass, charge, chargegroup]
        '''

    # .rtp-like line, with: atom type, mass, partial charge and chargegroup. The atom name is stored in toadd[1]
    log.debug(f'adding atoms:\ntoadd = {toadd} ([#atoms to add, atomname, 1st control atom, 2nd control atom])\natom_type: {ll[0]}  charge: {ll[2]}\n\n')
    for atcount in range(1, toadd[0] + 1):
        if toadd[0] == 1:
            atname = toadd[1]
        else:
            atname = f"{toadd[1]}{str(atcount)}"

        residues["tdb_file"].append(myfname)
        residues['resname'].append(myTername)
        residues['atom_name'].append(atname)
        residues['atom_type'].append(ll[0])
        residues['mass'].append(ll[1])
        residues['charge'].append(ll[2])
        residues['chargegroup'].append(ll[3])
        # Bound to the first control atom
        residues['bonds'].append(f'')
        if len(toremove) > 0:
            #log.debug(f'num_floats == 3, toadd[0] > 1, len(toremove) > 0, {len(toremove)}, {toremove}\n')
            residues['bonds'][-1] += toremove
        if toadd[2] not in residues['bonds'][-1]:
            residues['bonds'][-1] += f' {toadd[2]} '
        
        residues = add_bonds_ctrl_atoms(residues, myTername, toadd, atname)
    return residues

def read_add_entry(ll:List[str], residues:dict, myfname: str, myTername: str,
                   toremove: str, toadd: List) ->  Tuple[dict, List]:
    """
    Add new atoms.
    For each (group of) added atom(s), a two-line entry is necessary.
    The first line
    contains following entries:

    number of atoms, type of addition, name of the new atom, control atoms

    Type of addition can be int 1 -> 9 (see hydrogen database GROMACS documentation)
     1. one planar hydrogen, e.g. rings or peptide bond (1 H atom n on plane of atoms (i, j, k), at a distance of 0.1 nm from atom i, such that the angles (n-i-j) and (n-i-k) are > 90.
     2. one single hydrogen, e.g. hydroxyl: (1 H atom n on plane of atoms at a distance of 0.1 nm from atom i, such that angle (n-i-j)=109.5 degrees and dihedral (n-i-j-k)=trans
     3. two planar hydrogens, e.g. ethylene -C=CH 2 , or amide -C(=O)NH 2
     4. two or three tetrahedral hydrogens, e.g. -CH 3
     5. one tetrahedral hydrogen, e.g. C 3 CH
     ...
    => can use as bond info only 1st 2 atoms i and j (i bound to generated h tom n, j bound to i)

    After this line, another line follows that specifies the details of the added atom(s), in the same
    way as for replacing atoms:

    atom type  mass  charge  charge group (optional)
    
    Like in the hydrogen database (see rtp (page 429)), when more than one atom is connected to
    an existing one, a number will be appended to the end of the atom name. Note that, like in the
    hydrogen database, the atom name is now on the same line as the control atoms, whereas it was
    at the beginning of the second line prior to GROMACS version 3.3. 
    When the charge group field is left out, the added atom will have the same charge group number 
    as the atom that it is bonded to.

    ex:
    [ add ]
    ll[0]  ll[1]     ll[2]       ll[3]
    1       5         H53'        C5'       In this case
                                           ll[0]: #atoms to add, ll[1]: type of addition, ll[2]: name of the new atom, ll[3]: control atoms
    HN9   1.008000   0.0900       -1       In this case
                                           ll[0]: atom type, ll[1]: mass, ll[2]: charge, ll[3]: charge group (optional)
    
    The two type of lines will be distinguished from the number of entries that can be read as floats
    (2 for 1st line of definition, 3 for the 2nd)
    """
    
    log.debug(f'll in add: {ll}\n')
    
    # get number of floats
    num_floats = get_numfloats(ll)    
    log.debug(f'number of floats: {num_floats}\n')
    
    if num_floats == 2:
        # H database-like line, with: #atoms to add, rule #, atomname, ctrl atoms 
        # 
        toadd = [int(ll[0]), ll[2], ll[3], ll[4]]
        # keep [#atoms to add, atomname, 1st control atom] , 1st control atom is the 1st atom bound to the atoms to add

    elif num_floats == 3:
        residues = add_atoms(ll, residues, myfname, myTername, toremove, toadd)

    return residues, toadd

## tdb reader

def read_tbd_file(path: Path, residues:dict, ff_filename: str) -> dict:

    # Flags to understand whether the line is referring to replace or add entries
    delete = False
    # flags to understand line content
    replace = False
    add = False
    # info on atoms to be added/removed
    toremove = ''
    toadd = []

    myfname = ff_filename.split('.')[0]
    myTername = None

    line = 0
    jj = 1
    replace_count = 0

    with open(path / ff_filename, 'r') as f:
        mylines = f.readlines()
        # as ; are comments
        lines = [l for l in mylines if (len(l) > 0 and not l.startswith(';'))]
        for l in lines:
            line += 1
            log.debug(f'line: {l}\n')

            if ';' in l:
                # remove comments
                idx = l.index(';')
                l = l[0 : idx] 

            # if [ None ]: empty terminus, nothing to add
            if l.startswith('[ ') and not l.startswith("[ None ]"):
                replace, add, delete, toremove, replace_count, myTername = parse_ter_definition(l, replace, add, delete, toremove, replace_count, myTername)
                log.debug(f'''Has set following variables:\nmyresname = {myTername}\nreplace = {replace}\nadd = {add}\ndelete = {delete}\nreplace_count = {replace_count}\ntoremove = {toremove}''')
                
            else:  
                ll = l.split()
                if len(ll) > 0:
                    log.debug(f'resname: {myTername}\nll: {ll}\ntoadd: {toadd}\nreplace_count: {replace_count}\n')
                    log.debug(f'delete: {delete}\nreplace: {replace}\nadd: {add}\n')
                    if delete:
                        if ll[0] != '[':
                            toremove += f"-{ll[0]} "
                    if replace:
                        myTername = f"{myTername}"
                        residues = read_replace_entry(ll, replace_count, residues, myfname, myTername, toremove)

                    elif add:
                        residues, toadd = read_add_entry(ll, residues, myfname, myTername, toremove, toadd)
                    if toremove != '':
                        log.debug(f'toremove: {toremove}\n')
                    if toadd != []:
                        log.debug(f'toadd: {toadd}\n')
        if len(residues['tdb_file']) > 0:
            log.debug(f'Last added:\n|  tdb_file  | resname | atom_name | atom_type | charge | chargegroup |    bonds    |\n')
            log.debug(f'| {residues["tdb_file"][-1]:^8} | {residues["resname"][-1]:^7} | {residues["atom_name"][-1]:^9} | {residues["atom_type"][-1]:^9} | {residues["charge"][-1]:^6} | {residues["chargegroup"][-1]:^11} | {residues["bonds"][-1]:^11} |\n')

    return residues

## adjust tdb output

def get_bonds(g: pd.DataFrame) -> set:
    bond_tokens = set([])
    for b in g["bonds"].fillna(""):
        tokens = re.split(r"\s+", b.strip())
        for t in tokens:
            if t:
                bond_tokens.add(t)
    
    return bond_tokens

def pick_second_or_first(series: pd.Series):
    return series.iloc[1] if len(series) > 1 else series.iloc[0]

def merge_group(g: pd.DataFrame):
    #log.debug(g.head())
    out = {}
    out["tdb_file"] = g["tdb_file"].dropna().iloc[0]
    #out["resname"] = g["resname"].iloc[0] # not necessary as grouping by it
    #out["atom_name"] = pick_second_or_first(g["atom_name"]) # not necessary as grouping by it

    # atom_type
    out["atom_type"] = pick_second_or_first(g["atom_type"])
    
    # mass
    masses = g["mass"].astype(float)
    out["mass"] = pick_second_or_first(masses)
    
    # charge
    charges = g["charge"].astype(float)
    out["charge"] = pick_second_or_first(charges)

    # chargegroup
    out["chargegroup"] = pick_second_or_first(g["chargegroup"])

    # ---- merge bonds ----
    bond_tokens = set()
    for b in g["bonds"].fillna(""):
        tokens = re.split(r"\s+", b.strip())
        for t in tokens:
            if t:
                bond_tokens.add(t)
    # ricostruzione stringa bonds
    out["bonds"] = " ".join(sorted(bond_tokens))
    return pd.Series(out)

def merge_duplicate_atoms(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    merged_df = (
        df.groupby(["resname", "atom_name"], as_index=False)
          .apply(merge_group)
          .reset_index(drop=True)
    )
    #log.debug(f'merged_df:\n{merged_df.head()}')
    return merged_df

## helper of se_data

def merge_se(residues: pd.DataFrame, se_df: pd.DataFrame) -> Tuple[pd.DataFrame, List, List]:
    residues['sigma'] = np.nan
    residues['epsilon'] = np.nan
    residues['mass'] = np.nan

    # Merge with sigma, epsilon and mass data
    ## There are problems with HT and OT, that in ffnonbonded.itp  come with 2 different mass values for:
    #  heavy hydrogens (4.032000 vs 1.008000)
    #  light oxygens   (7.935400 vs 15.999400)

    ht_data = ['HT', 1.008000, 0.040001352445, 0.1924] 
    ot_data = ['OT', 15.999400, 0.315057422683, 0.6363]

    missing_ats = set()
    dup_ats = []

    for i in residues.index:
        myrow = residues.iloc[i]
        if myrow['atom_type'] in ['HT', 'OT']:
            if myrow['atom_type'] == 'HT':
                #log.debug('HT')
                residues.at[i, 'sigma'] = ht_data[2]
                residues.at[i, 'epsilon'] = ht_data[3]
                residues.at[i, 'mass'] = ht_data[1]      
            
            elif myrow['atom_type'] == 'OT':
                #log.debug('OT')
                residues.at[i, 'sigma'] = ot_data[2]
                residues.at[i, 'epsilon'] = ot_data[3]
                residues.at[i, 'mass'] = ot_data[1]      
        else:
            se_data = se_df[se_df['atom_type'] == myrow['atom_type']]

            if se_data.shape[0] == 1:
                residues.at[i, 'sigma'] = float(se_data['sigma'].iloc[0])
                residues.at[i, 'epsilon'] = float(se_data['epsilon'].iloc[0])
                residues.at[i, 'mass'] = float(se_data['mass'].iloc[0])
            else:
                if se_data.shape[0] == 0:
                    missing_ats.add(myrow['atom_type'])
                elif se_data.shape[0] > 1:
                    dup_ats.append(myrow['atom_type'])
                    #dup_ats.add(se_data.to_dict('split'))
                    dup_ats.append(se_data)

    return residues, [i for i in missing_ats], dup_ats

## ffnonbonded.itp parser

def get_sigma_eps_mass(inpath: str, outpath: str, nonbonded_name: str) -> Tuple[Path, Path]:
    '''
    parsing of ffnonbonded.itp

    Just one problem in parsing HT and OT entries
    in ffnonbonded.itp  come with 2 different mass values for:
     - heavy hydrogens (4.032000 vs 1.008000)
     - light oxygens   (7.935400 vs 15.999400) 
    '''

    atypes = {"atom_type":[], "mass":[], "sigma":[], "epsilon":[]}
    ptypes = {"atom1":[], "atom2":[], "sigma":[], "epsilon":[], "func":[]}
    path = Path(inpath)
    out = Path(outpath)
    with open(path / nonbonded_name, 'r') as f:
        mylines = f.readlines()
        # as ; are comments
        lines = [l for l in mylines if (len(l) > 0 and not l.startswith(';'))]
        for l in lines:
            if ';' in l:
                # remove comments
                idx = l.index(';')
                l = l[0 : idx]
            log.debug(f'L: {l}\n')
            if '\t' in l:
                l = l.replace('\t', ' ')

            data = [i for i in l.split(' ') if i not in ['', '\n']]
            log.debug(f'data prefloat: {data}\n')
            if len(data) >= 5:
                data[-1] = get_float(data[-1])
                log.debug(f'data: {data}\n')

                if len(data) == 7:
                    ## Then entry of [ atomtypes ], specifying type, atnum, mass, charge, ptype, sigma, epsilon of each listed atom
                    atypes["atom_type"].append(data[0])
                    atypes["mass"].append(data[2])
                    atypes["sigma"].append(float(data[-2]))
                    atypes["epsilon"].append(float(data[-1]))

                elif len(data) == 5:
                    ## Then entry of [ pairtypes ], specifying type1, type2, func, sigma, epsilon of each listed atom
                    ## func specifies the function type for the interaction between these two atom types, if 1: 12-6 Lennard-Jones potential
                    ptypes["atom1"].append(data[0])
                    ptypes["atom2"].append(data[1])
                    ptypes["func"].append(int(data[2]))
                    ptypes["sigma"].append(float(data[-2]))
                    ptypes["epsilon"].append(float(data[-1]))

                else:
                    log.debug(f'len(data) != 5 and len(data) != 7\n{data}\n{len(data)}\n\n')

    inttypes = set(ptypes["func"])
    if len(inttypes) > 1:
       log.debug(f'Warning: more interaction types found: {" ".join(inttypes)}')
    
    se_df = pd.DataFrame.from_dict(atypes)
    pdf = pd.DataFrame.from_dict(ptypes)
    sig_eps = out / 'atomtypes_se.csv'
    pairtypes = out / 'pairtypes.csv'
    
    se_df.to_csv(sig_eps, header = True, index=False)
    pdf.to_csv(pairtypes, header = True, index=False)
    return sig_eps, pairtypes

## Callers

def get_residues(path: Path, outpath: Path) -> Path:
    '''
    parser of .rtp files from ff folder path
    returns dictionary with entries:
    - rtp_file: filename of origin
    - resname: residue name
    - atom_name: atom name 
    - atom_type: element 
    - charge: partial charge
    - chargegroup: (not necessary) small atom groups with tot charge = 0, 1 or -1
    - bonds: bond info with other atoms
    '''
    log.info(f'Obtaining residue info with get_residues, parameters:\n- forcefield folder: {path.name}\n')
    pathfnames = [f.name for f in path.iterdir() if f.is_file() and f.name.endswith('.rtp')]

    ## PARSE THE RESIDUE TOPOLOGY (.rtp) FILES
    # get list of residue topology files 
    rtps = [f for f in pathfnames if f.endswith('.rtp')]
    log.info(f'The residues included in the following .rtp files will be considered:\n{pathfnames}')
    
    # Initialize dict of results
    residues = {'rtp_file':[], 'resname':[], 'atom_name':[], 'atom_type':[], 'charge':[], 'chargegroup':[], 'bonds':[]}

    # start to iterate over files
    for ff in rtps:
        log.info(f'Analyzing file {ff}')
        residues = read_rtp_file(path, ff, residues)
        
    tmp = pd.DataFrame.from_dict(residues)
    o = outpath / "residues.csv"
    tmp.to_csv(o, index=False, header = True)
    return o

def get_ters(path: Path, 
             outpath: Path,
            ) -> Path:
    f'''
    parsing of .tdb files from ff folder path {path.as_posix()}

    '''
    
    # Initialize dict of results
    residues = {'tdb_file':[], 'resname':[], 'atom_name':[], 'atom_type':[], 'mass':[], 'charge':[], 'chargegroup':[], 'bonds':[]}
    pathfnames = [f.name for f in path.iterdir() if f.is_file()]
    ## PARSE THE TERMINI DATABASE FILES (.tdb)
    tdbs = [f for f in pathfnames if f.endswith('.tdb')]
    tdbs_str = "\n".join([path.as_posix() + t for t in tdbs])
    log.info(f'termini database file parsing will analyze files: {tdbs_str}\n')
    for ff_filename in tdbs:
        residues = read_tbd_file(path, residues, ff_filename)
    
    #log.debug('residues dict keys:')
    #for k in residues.keys():
    #    log.debug(f'{k}, {len(residues[k])}')
    tmp = pd.DataFrame.from_dict(residues)
    #log.debug(f'residues -> tmp df: {tmp.head()}')
    tmp = merge_duplicate_atoms(tmp)

    o = outpath / "ters.csv"
    tmp.to_csv(o, index=False, header = True)

def merge_data(outpath: str, 
               na_values: Optional[List[str]] = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"])-> Path:
    
    out = Path(outpath)

    # retrive residue topologies
    residues = pd.read_csv(out / 'residues.csv', header = 0, index_col = False, na_values = na_values, keep_default_na = False)
    se_df = pd.read_csv(out / 'atomtypes_se.csv', header = 0, index_col = False, na_values = na_values, keep_default_na = False)
    residues = residues.merge(
                              se_df[['atom_type', 'sigma', 'epsilon', 'mass']] ,
                              how="left",
                              on="atom_type"
                              )
    o = out / "final_ff_residues.csv"
    residues.to_csv(o, header = True, index=False)

    # retrive termini
    ters = pd.read_csv(out / 'ters.csv', header = 0, index_col = False, na_values = na_values, keep_default_na = False)    
    ters = ters.merge(               
                      se_df[['atom_type', 'sigma', 'epsilon']] ,
                      how="left",
                      on="atom_type"
                      )

    o = out / "final_ff_ters.csv"
    ters.to_csv(o, header = True, index=False)

def add_residue(data: pd.DataFrame,
                 outpath: str, 
                 na_values: Optional[List[str]] = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"])-> Path:
    
    out = Path(outpath)
    log.debug(f'adding {data.head()}')
    # retrive residue topologies
    residues = pd.read_csv(out / 'final_ff_residues.csv', header = 0, index_col = False, na_values = na_values, keep_default_na = False)
    # concat data
    residues = pd.concat([residues, data])
    residues.to_csv(out / 'final_ff_residues.csv', header = True, index=False)
