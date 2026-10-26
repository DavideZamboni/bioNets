import numpy as np
from os import makedirs
import pandas as pd
from pathlib import Path
import MDAnalysis as mda
from typing import Union, List, Tuple
import difflib
import logging

log = logging.getLogger(__name__)

##

_ATOMS_SUFFIX = '_atoms.csv'

def topology_key(atoms_file: str) -> str:
    '''
    Name of a residue topology from its resid_topologies file:
    "MET_1_atoms.csv" -> "MET_1". The key names the topology in every later
    file (found_topologies/<key>_ff_atoms.csv,
    possible_topologies/<key>_uncertain_atoms_<n>.csv) and is what
    `mypsns accept-topology` takes on the command line.

    Strips the exact suffix. The old str.rstrip("_atoms.csv") stripped a
    character SET instead, so a key ending in any of "_atomscv." lost those
    letters too (e.g. "Na_atoms.csv" -> "N").
    '''
    if atoms_file.endswith(_ATOMS_SUFFIX):
        return atoms_file[:-len(_ATOMS_SUFFIX)]
    return atoms_file

def get_residue_types(top_df: pd.DataFrame, res_topology_folder: Path) -> pd.DataFrame:
    # Iterate over each resname
    restypes_dict = {'resname':[], 'atoms':[], 'rixs':[]}
    for rname, rndf in top_df.groupby(by = 'resname'):
        atsets = []
        resixs = []

        # find all unique atom sets 

        for rix, ridf in rndf.groupby(by = 'resindex'):
            myats = ridf[['name', 'mass']]
            atnames = [an for an in myats['name'].to_list()]

            # Test if atoms with same atom name are present
            found = False
            if len(set(atnames)) != len(atnames):
                log.warning('get_restypes: duplicate atom names for residue %s (resindex %s)', rname, rix)
            for i, ats in enumerate(atsets):
                ot_atnames = [an for an in ats['name'].to_list()]
                if atnames == ot_atnames:
                    resixs[i].append(rix)
                    found = True
                    break

            if not found:
                atsets.append(myats)
                resixs.append([rix])
            
        # add rows to dataframe, writing info
        for i in range(len(atsets)):
            if len(resixs[i]) > 1:
                ixs = ' '.join([str(k) for k in resixs[i]])
            else:
                ixs = str(resixs[i][0])
            basename = f"{rname}_"
            if i > 0:
                basename += f"{i}_"

            ix_p = f"{basename}rixs.txt"
            res_p = f"{basename}atoms.csv"
            with open(res_topology_folder / ix_p, "w") as resixout:
                resixout.write(ixs)
                resixout.close()

            atsets[i].to_csv(res_topology_folder / res_p, header=True, index = False)
            restypes_dict['resname'].append(rname)
            restypes_dict['atoms'].append(res_p)
            restypes_dict['rixs'].append(ix_p)

    restypes = pd.DataFrame.from_dict(restypes_dict)

    return restypes

## Function for initial parsing of provided pdb file
def parse_initial_topology_info(u: mda.Universe, outpath: Path) -> None:
    
    topdata = {'chainID':[], 'resindex':[],  'resid':[], 'resname':[], 'index':[],  'id':[], 'name':[], 'type':[], 'element':[], 'mass':[], 'charge':[]}

    notfound = set()
    for at in u.atoms:
        try:
            topdata['resindex'].append(at.resindex)
        except mda.exceptions.NoDataError:
            notfound.add('resindex')
            topdata['resindex'].append('')

        try:
            topdata['resid'].append(at.resid)
        except mda.exceptions.NoDataError:
            notfound.add('resid')
            topdata['resid'].append('')

        try:
            topdata['index'].append(at.index)
        except mda.exceptions.NoDataError:
            notfound.add('index')
            topdata['index'].append('')

        try:
            topdata['id'].append(at.id)
        except mda.exceptions.NoDataError:
            notfound.add('id')
            topdata['id'].append('')

        try:
            topdata['charge'].append(at.charge)
        except mda.exceptions.NoDataError:
            notfound.add('charge')
            topdata['charge'].append('')

        try:
            topdata['chainID'].append(at.chainID)
        except mda.exceptions.NoDataError:
            notfound.add('chainID')
            topdata['chainID'].append('')

        try:
            topdata['element'].append(at.element)
        except mda.exceptions.NoDataError:
            notfound.add('element')
            topdata['element'].append('')

        try:
            topdata['mass'].append(at.mass)
        except mda.exceptions.NoDataError:
            notfound.add('mass')
            topdata['mass'].append('')

        try:
            topdata['name'].append(at.name)
        except mda.exceptions.NoDataError:
            notfound.add('name')
            topdata['name'].append('')

        try:
            topdata['resname'].append(at.resname)
        except mda.exceptions.NoDataError:
            notfound.add('resname')
            topdata['resname'].append('')

        try:
            topdata['type'].append(at.type)
        except mda.exceptions.NoDataError:
            notfound.add('type')
            topdata['type'].append('')   

    ## Merging obtained info in dataframe
    truecols = [i for i in topdata.keys() if len([j for j in topdata[i] if j != '']) != 0]
    top_df = pd.DataFrame.from_dict(topdata)[truecols]

    ## Sometimes (e.g. from pdb topologies in input) the element is reported as atom type
    ## The two cases can be distinguished as forcefield atom names are longer than the max two characters of the atom symbols
    if 'element' not in truecols:
        if 'type' in truecols:
            maxlen = max([len(k) for k in top_df['type']])
            if maxlen <= 2:
                top_df['element'] = top_df['type']

    top_df.to_csv(outpath / 'initial_topology_info.csv', header = True, index = False)

## Functions used in tolopolgy recognition (assignation of each atom in the topology to a precise atom in the forcefield)

def get_1st_str_letter(mystring: str) -> Union[str, None]:
    '''
    Given a string of alphanumerical characters, returns the 1st letter if present or None
    '''

    for i in mystring:
        try:
            float(i)

        except:
            return i
    
    return None

def check_masses(ff_res: pd.DataFrame, top_atoms_df: pd.DataFrame) -> float:
    #print(ff_res.head())
    #print(top_atoms_df.head())
    ff_masses = ff_res['mass'].to_list()
    top_masses = top_atoms_df['mass'].to_list()
    ff_atnames = ff_res['atom_name'].to_list()
    top_atnames = top_atoms_df['name'].to_list()
    if (len(ff_atnames) != len(top_atnames)) or (len(ff_masses) != len(top_masses)):
        log.debug('check_masses: different number of atoms')
        return -np.inf
    else:
        if ff_masses != top_masses:
            log.debug('check_masses: different masses (top=%s ff=%s)', top_masses, ff_masses)
            return -np.inf

        found = [name for name in top_atnames if name in ff_atnames]
        score_atnames = len(found)

        log.debug('check_masses: common atoms %s -> score %s', found, score_atnames)
        return score_atnames

def check_namemasses(tmp_data: pd.DataFrame, results: dict, top_rn: str, top_atoms_df: pd.DataFrame) -> Tuple[dict, bool]:
    '''
    Function that checks wheter all atoms found in topology are matched to the analyzed ff entry.
    If so, updates the results and sets found to True, passing to the next topologyy entry.

    mytopology: ridxs, top_rn, top_atoms
     
    '''
    log.debug('check_namemasses: starting check for %s', top_rn)

    top_atoms = [k for k in top_atoms_df['name']]
    ff_atnames = [k for k in tmp_data['atom_name']]
    not_found = [k for k in top_atoms if k not in ff_atnames]
    extra = [k for k in ff_atnames if k not in top_atoms]
    log.debug('check_namemasses: atoms not found=%s extra=%s', not_found, extra)

    # check if all atoms have been found and none is in excess after the associations
    if len(not_found) == 0 and len(extra) == 0:
        # Then precisely found
        results[top_rn].append([np.inf, tmp_data])
        log.debug('check_namemasses: forcefield residue accepted as match for %s', top_rn)
        return results, True
    else:
        found = [k for k in top_atoms if k in tmp_data['atom_name'].values]
        ff_not_assigned = tmp_data[(~tmp_data['atom_name'].isin(found))]
        top_not_assigned = top_atoms_df[(top_atoms_df['name'].isin(not_found))]
        score_atnames = check_masses(ff_not_assigned, top_not_assigned)

        results[top_rn].append([score_atnames, tmp_data])

        return results, False

def get_replaced_removed_terAtoms(ter_residue: pd.DataFrame) -> Tuple[set, list, dict]:
    replaced_removed = set()
    to_be_removed = []
    to_be_replaced = {}
    for ix, row in ter_residue.iterrows():
        bonds = [b for b in ter_residue.at[ix, 'bonds'].split(' ') if b != '']
        atname = ter_residue.at[ix, 'atom_name']
        for b in bonds:
            if b.startswith('r'):
                myb = b.lstrip('r')
                to_be_replaced[atname] = myb
                replaced_removed.add(myb)
            elif b.startswith('-'):
                myb = b.lstrip('-')
                to_be_removed.append(myb)
                replaced_removed.add(myb)
    replaced = [to_be_replaced[k] for k in to_be_replaced.keys()]
    log.debug('ter: atoms replaced from current residue: %s (map %s)', replaced, to_be_replaced)

    return replaced_removed, to_be_removed, to_be_replaced

def bond_tokens(bonds) -> List[str]:
    '''The tokens of a bonds cell ('' and NaN give []).'''
    return [t for t in bonds.split() if t] if isinstance(bonds, str) else []

def is_patch_instruction(token: str) -> bool:
    '''
    True for a token of a terminal-patch row (final_ff_ters.csv) that is an
    instruction rather than a bond. ff_reader writes the .tdb sections into the
    bonds column of the patch rows as "rX" ([ replace ]: this atom replaces
    atom X) and "-X" ([ delete ]: the patch deletes atom X, repeated on every
    row of the patch), plus "+X" when an [ add ] control atom is not one of the
    patch's own rows. Inside a patch row none of them is a bond, unlike "-X" /
    "+X" in a residue row, which name atoms of the previous / next residue.
    '''
    return token[0] in '-+r'

def get_ter_ffres(ter_residue: pd.DataFrame, ff_residue: pd.DataFrame,to_be_replaced: dict, to_be_removed: List) -> pd.DataFrame:
    '''
    Apply one terminal patch (a .tdb entry, as ff_reader stores it in
    final_ff_ters.csv) to a forcefield residue: to its atoms AND to every bonds
    list.

     - ter_residue: the patch's rows
     - ff_residue: the residue's rows (final_ff_residues.csv)
     - to_be_replaced: {patch atom name: residue atom it replaces}, from the
       "rX" tokens. A different name renames the atom (e.g. the apr2019 COO-
       replaces O by OT1).
     - to_be_removed: residue atoms the patch deletes (its "-X" tokens that
       the residue has).

    Atoms: deleted atoms are dropped; a replaced atom takes the patch row's
    name, type, mass, charge, chargegroup, sigma and epsilon; the patch's other
    rows are added after the residue's atoms.

    Bonds get the same edits:
     - a deleted atom disappears from every bonds list, and a renamed atom is
       renamed in every bonds list, not only in its own row;
     - a replaced atom keeps its residue bonds, minus "+X" (a patch that
       replaces the atom bonded to the next residue makes a C-terminus), and
       gains the bonds its patch row lists;
     - an added atom gets the bonds its patch row lists, without the patch's
       instruction tokens (is_patch_instruction);
     - bonds within the residue are made two-sided. ff_reader records the
       bond between an [ add ] atom and its control atom in the control atom's
       row, and its substring checks used to drop it from the added row (the
       "N" of an NH3+ H1 is "in" "-HN");
     - a token that names no atom of the patched residue is dropped. Leaving
       it is how a hand-made NHT3+ (NH3+ with its added H1-H3 renamed HT1-HT3
       but N still listing "H1 H2 H3") reached found_topologies. Such a patch
       leaves atoms without bonds, which save_results reports.

    A bonds cell that none of these edits changes is kept byte for byte.

    Before, only the replaced atoms' own rows were edited: renamed and deleted
    names survived in the other rows, and an added atom that its control atom
    did not list by name kept the instruction tokens as bonds (NHT3+ HT1
    "-HN -HN1 -HN2 -HT1 -HT2", which every later stage reads as bonds to the
    previous residue).
    '''
    # a replacement wins over a deletion of the same atom
    removed = set(to_be_removed) - set(to_be_replaced.values())
    tmp = ff_residue[~ff_residue['atom_name'].isin(removed)].copy()

    # resolve every replaced row before renaming any of them, so that a rename
    # onto a name that another replacement frees cannot pick the wrong row
    replaced_ix = {new: tmp.index[tmp['atom_name'] == old][0]
                   for new, old in to_be_replaced.items()}
    renames = {old: new for new, old in to_be_replaced.items() if new != old}
    ter_ix = {an: ix for ix, an in ter_residue['atom_name'].items()}
    log.debug('ter replace: to_be_replaced=%s removed=%s renames=%s',
              to_be_replaced, sorted(removed), renames)

    cols = ['resname', 'atom_name', 'atom_type', 'mass',
            'charge', 'chargegroup', 'sigma', 'epsilon']
    for new, ix in replaced_ix.items():
        tmp.loc[ix, cols] = ter_residue.loc[ter_ix[new], cols].values

    # added atoms: patch rows whose name the patched residue does not have yet
    to_add = ter_residue[~ter_residue['atom_name'].isin(tmp['atom_name'])].copy()
    to_add.rename(columns={'tdb_file': 'rtp_file'}, inplace=True)
    log.debug('ter replace: atoms to be added: %s', to_add['atom_name'].to_list())
    out = pd.concat([tmp, to_add], ignore_index=True)
    names = out['atom_name'].tolist()
    nameset = set(names)
    n_res = tmp.shape[0]
    replaced_rows = {tmp.index.get_loc(ix) for ix in replaced_ix.values()}

    def patch_bonds(atom_name):
        # a patch row names atoms by their final names; map a renamed-away
        # name only when it is not itself an atom of the patched residue
        return [t if t in nameset else renames.get(t, t)
                for t in bond_tokens(ter_residue.at[ter_ix[atom_name], 'bonds'])
                if not is_patch_instruction(t)]

    original = out['bonds'].tolist()
    bonds = []
    for row, (an, cell) in enumerate(zip(names, original)):
        if row < n_res:
            toks = []
            for t in bond_tokens(cell):
                if t[0] == '+' and row in replaced_rows:
                    continue
                if t[0] not in '+-':
                    if t in removed:
                        continue
                    t = renames.get(t, t)
                toks.append(t)
            if row in replaced_rows:
                toks += patch_bonds(an)
        else:
            toks = patch_bonds(an)
        bonds.append(list(dict.fromkeys(toks)))   # drop repeats, keep order

    for row, toks in enumerate(bonds):
        dangling = [t for t in toks if t[0] not in '+-' and t not in nameset]
        if dangling:
            log.debug('ter %s: %s lists %s, not atoms of the patched residue; dropped',
                      ter_residue['resname'].iloc[0], names[row], dangling)
            bonds[row] = [t for t in toks if t not in dangling]

    row_of = {an: row for row, an in enumerate(names)}
    for row, toks in enumerate(bonds):
        for t in list(toks):
            if t[0] in '+-' or t == names[row]:
                continue
            partner = bonds[row_of[t]]
            if names[row] not in partner:
                partner.append(names[row])

    out['bonds'] = [cell if row < n_res and bond_tokens(cell) == toks else ' '.join(toks)
                    for row, (cell, toks) in enumerate(zip(original, bonds))]
    return out

def get_res_w_ter(ter_residue: pd.DataFrame, ff_residue: pd.DataFrame) -> pd.DataFrame:
    '''
    Replaces or substitutes the atoms found in ff_residue according to the information
    stored in ter_residue bonds entries. In particular, if an atom name in bonds starts with:
     - 'r': atom to be replaced with the one here defined
     - '-': atom to be removed
    While the information concerning the atoms to be removed is redundant (present in all bonds of the ter entry),
    only the atom that replaces the original has the atom to be replaced in bond list
    '''
    # check if all ter atom names to be replaced or removed are present in the residue
    # otherwise it cant be the right ter entry

    # key: atom in ter, value: atom in ff_residue
    tn = ter_residue['resname'].to_list()[0]
    rn = ff_residue['resname'].to_list()[0]
    log.debug('ter: trying to apply ter %s to forcefield residue %s', tn, rn)

    replaced_removed, to_be_removed, to_be_replaced = get_replaced_removed_terAtoms(ter_residue)

    # test if all atoms to be replaced or removed are present in the forcefield residue
    atnames = set(ff_residue['atom_name'])

    # Only the rX atoms (replacements) must exist in the ff residue.
    # The -X remove list may over-specify atoms from other ter definitions
    # Filter to_be_removed to only the atoms actually present.
    filtered_remove = [a for a in to_be_removed if a in atnames]


    # All replacement targets (rX) MUST be present
    replace_targets = list(to_be_replaced.values())
    missing_replacements = [a for a in replace_targets if a not in atnames]
    log.debug('ter: replacement targets=%s missing=%s filtered_remove=%s',
              replace_targets, missing_replacements, filtered_remove)

    if missing_replacements:
        log.debug('ter: cannot apply, replacement targets not in residue: %s', missing_replacements)
        tmp = ff_residue.copy()
    else:
        tmp = get_ter_ffres(ter_residue, ff_residue, to_be_replaced, filtered_remove)

    return tmp

def get_minmax_delta(myters: pd.DataFrame):
    min_delta = float('inf')
    max_delta = float('-inf')
    #addter = None
    #remter = None
    # for each ter:
    # - number added atoms: number of defined atoms 
    #                       that do not have any bond starting with 'r'
    # - number of removed atoms: length of set of bonds starting with '-'

    for tername, ter in myters.groupby(by = 'resname'):
        my_add = 0
        removed = set()
        for ix, row in ter.iterrows():
            bonds = myters.at[ix, 'bonds']
            bs = [b for b in bonds.split() if b != '']
            ter_removed = [b for b in bs if b.startswith('-')]
            for r in ter_removed:
                removed.add(r)
            if any(b.startswith('r') for b in bs):
                my_add += 1

        if len(removed) > 3:
            lR = 3
        else:
            lR = len(removed)
        total_delta = my_add - lR
        if total_delta < min_delta:
            min_delta = total_delta
            remter = tername
        if total_delta > max_delta:
            max_delta = total_delta
            addter = tername

    #print(remter, addter)
    return min_delta, max_delta
            

'''
myters = pd.read_csv('/home/davide/storage/lab_dati/psns_out/rec_grk1_ros_edit/charmm36_apr2019_vm/final_ff_ters.csv', header=0, index_col=False)
min_delta, max_delta = get_minmax_delta(myters)
print(min_delta, max_delta)
'''
def check_ters(results: dict, top_atoms_df: pd.DataFrame, 
                ter_params: pd.DataFrame, ff_residue: pd.DataFrame, 
                top_rn: str, top_atoms: List[str]) -> Tuple[dict, bool]:

    #print('entering ter:')
    #print(top_atoms_df.head(50))
    #print(ff_residue.head(50))
    ff_resname = ff_residue['resname'].to_list()[0]

    fname = ff_residue['rtp_file'].to_list()[0]
    myters = ter_params[ter_params['tdb_file'] == fname]

    top_atnums = top_atoms_df.shape[0]
    #if (ff_resname == 'GLN'):
    #    print('GLN in ters')
    #    print(ff_residue.head(50))
    #    print(top_atoms_df.head(50))
    #    print(f'top_atnums: {top_atnums}')

    for tername, ter_residue in myters.groupby(by = 'resname'):
        log.debug('check_ters: applying ter substitutions (%s)', tername)
        tmp = get_res_w_ter(ter_residue, ff_residue)

        if tmp.shape[0] == top_atoms_df.shape[0]:
            results, found = check_namemasses(tmp, results, top_rn, top_atoms_df)
            if found:
                return results, found
    return results, False

def get_results(results: dict, top_residuename: str, top_atoms_df: pd.DataFrame, ff_params: pd.DataFrame, ter_params: pd.DataFrame, ff_resnames: List[str], maxToCheck: int) -> dict:
    log.debug('get_results: starting for topology residue %s', top_residuename)

    if top_residuename not in results:
        results[top_residuename] = []

    top_atnames = top_atoms_df['name'].to_list()

    log.debug('get_results: residue %s atoms = %s', top_residuename, ", ".join(top_atnames))
    found = False
    ## SequenceMatcher: flexible class for comparing pairs of sequences of any type
    #                 - the .ratio() method returns a measure of the sequences’ similarity as a float in the range [0, 1].
    sorted_ff_resnames = sorted(ff_resnames, key=lambda x: difflib.SequenceMatcher(None, x, top_residuename.split('_')[0]).ratio(), reverse=True)
    ### ITERATE OVER ALL POSSIBLE RESIDUE TOPOLOGIES IN FORCEFIELD
    index = 0
    if maxToCheck == -1:
        maxToCheck = len(sorted_ff_resnames) - 1

    while (index < len(sorted_ff_resnames)) and (index < maxToCheck):

        ff_residuename = sorted_ff_resnames[index]
        index += 1

        if not isinstance(ff_residuename, str):
            continue

        ff_res = ff_params[ff_params['resname'] == ff_residuename]

        if ff_res.shape[0] != len(top_atnames):
            continue

        results, found = check_namemasses(
            ff_res, results, top_residuename, top_atoms_df
        )

        if found:
            log.debug('get_results: exact match for %s', top_residuename)
            return results  # EXIT ENTIRE FUNCTION

    index = 0
    while (index < len(sorted_ff_resnames)) and (index < maxToCheck):

        ff_residuename = sorted_ff_resnames[index]
        index += 1
        log.debug('get_results: evaluating forcefield residue %s', ff_residuename)

        if not isinstance(ff_residuename, str):
            continue

        ff_res = ff_params[ff_params['resname'] == ff_residuename]


        results, found = check_ters(
            results, top_atoms_df, ter_params, ff_res, top_residuename, top_atnames
        )

        if found:
            log.debug('get_results: matched %s via TER', top_residuename)
            return results

    if not found:
        log.warning('Could not associate topology residue %s (atoms=%s) to any forcefield entry',
                    top_residuename, top_atnames)

    return results
'''
tests = Path('/home/davide/storage/lab_dati/test_files')
ter_params = pd.read_csv(tests/'final_ff_ters.csv', header=0, index_col=False)

# normal test
ff_params = pd.read_csv(tests/'ff_ala.csv', header=0, index_col=False)
top_atoms_df = pd.read_csv(tests/'ALA_atoms.csv', header=0, index_col=False)
#results, logstr = get_results({}, "ALA", top_atoms_df, ff_params, ter_params, ["ALA"])
#print(results)
#print(logstr)

# ter test
resnames = ['ILE', 'MET', 'LEU']
for rn in resnames:
    ff_params = pd.read_csv(tests/f'ff_{rn.lower()}.csv', header=0, index_col=False)
    top_atoms_df = pd.read_csv(tests/f'{rn.upper()}_1_atoms.csv', header=0, index_col=False)
    results = get_results({}, f"{rn.upper()}", top_atoms_df, ff_params, ter_params, [f"{rn.upper()}"], -1)
    #print(results)
    #print(results)
'''
def collect_results(results: dict) -> Tuple[dict, dict]:

    # results: dictionary:
    #  - keys: residue names found in topology (made unique for topology)
    #           e.g. if ALA at Nter there will be ALA and ALA_1, one referring 
    #                to Nter residue topology and the other to intrachain ALA residue topoligy
    #                                       )
    #  - values: lists of tuples with entries
    #            (score_atnames: int, ff_atoms: pd:DataFrame)

    found_residues = {}
    best_choices = {}
    for rn in results.keys():
        if not results[rn]:
            # No forcefield candidate at all (no ff entry with a comparable atom
            # set). Leave this residue out of both dicts: save_results then records
            # it as unmatched and get_ff_parameters halts the pipeline with a clear
            # message instead of crashing here on max() of an empty list.
            log.warning('No forcefield candidate found for topology residue %s', rn)
            continue
        maxScore = max([r[0] for r in results[rn]])
        topres = sorted(results[rn], key = lambda x:x[0], reverse = True)
        if maxScore == np.inf:
            # then precisely found
            found_residues[rn] = topres[0][1]

        else:
            # store all results with score corresponding to maxScore 
            best_choices[rn] =  [(topres[k][0], topres[k][1]) for k in range(0, len(topres)) if topres[k][0] == maxScore]

    return found_residues, best_choices


def save_results(out : Path, found_residues: dict, best_choices: dict, nas: list) -> Tuple[bool, list]:
    # create directory for results
    found_topologies_path = out / 'found_topologies'
    makedirs(found_topologies_path, exist_ok = True)

    # create directory for top scoring topologies for each missing topology
    possible_topologies_path = out / 'possible_topologies'
    makedirs(possible_topologies_path, exist_ok = True)

    # retrieve df with info about different topologies in the input structure
    different_residue_topologies = pd.read_csv(out / 'different_residue_topologies.csv', header = 0, index_col = False,
                                               na_values = nas, keep_default_na = False)
    
    # initilize result columns
    different_residue_topologies["ff_atoms"] = ""
    different_residue_topologies["possible_atoms"] = ""
    different_residue_topologies["unmatching_masses"] = 0

    '''
    found_residues: dict with:
     - keys: residue names found in topology (made unique for topology)
             e.g. if ALA at Nter there will be ALA and ALA_1, one referring 
                  to Nter residue topology and the other to intrachain ALA residue topoligy
     - values: a ff_residue topology has been found to match all the atom names in the structure_topology
               ff_topology: pd.DataFrame               
    best_choices: dict with:
     - keys: same of found_residues
     - values: no ff_residue topology has been found to perfectly match all the atom names in the structure_topology
               score is 0 if all masses match, else = -1 for each of the mismatching masses
                list of tuples (score: int, ff_topology: pd.DataFrame)
    '''
    complete_match = True
    unmatched: list = []   # residues with no perfect forcefield match (halts pipeline)
    for ix, row in different_residue_topologies.iterrows():
        atoms_file = different_residue_topologies.loc[ix, "atoms"]
        unique_top_resname = topology_key(atoms_file)
        resname = different_residue_topologies.loc[ix, 'resname']
        if unique_top_resname in found_residues.keys():
            ff_topology_name = unique_top_resname + "_ff_atoms.csv"
            log.debug('Perfect atom-name match for residue %s (%s) -> %s',
                      resname, unique_top_resname, ff_topology_name)

            different_residue_topologies.loc[ix, "ff_atoms"] = ff_topology_name
            # an exact match has the structure's atom names, so each atom's
            # structure name is its own; the column makes that explicit and is
            # the one accept-topology fills with differing names
            ff_atoms = with_structure_names(found_residues[unique_top_resname],
                                            found_residues[unique_top_resname]['atom_name'])
            unbonded = unbonded_patch_atoms(ff_atoms)
            if unbonded:
                log.warning('Residue %s (%s): the terminal patch %s leaves %s without a bond to '
                            'any other atom of the residue. Fix that patch\'s rows in the '
                            'forcefield\'s final_ff_ters.csv (bonds column) and re-run '
                            '`mypsns topology --overwrite`.',
                            resname, unique_top_resname,
                            '/'.join(sorted(set(ff_atoms.loc[ff_atoms['atom_name'].isin(unbonded),
                                                             'resname'].astype(str)))),
                            ' '.join(unbonded))
            ff_atoms.to_csv(found_topologies_path / ff_topology_name, header = True, index = False)

        elif unique_top_resname in best_choices.keys():
            complete_match = False
            possile_tops_str = ""
            candidates = best_choices[unique_top_resname]
            # The candidates' score is not a count: it is 0 when the leftover
            # (differently named) atoms match by mass and -inf otherwise (see
            # check_masses). Count the shared atom names directly, taking the
            # best candidate, so the report says how close the match is.
            top_atnames = pd.read_csv(out / 'resid_topologies' / atoms_file, header = 0, index_col = False,
                                      na_values = nas, keep_default_na = False)['name'].astype(str).tolist()
            matching_ats = max(len(set(top_atnames) & set(c[1]['atom_name'].astype(str))) for c in candidates)
            total_ats = len(top_atnames)
            log.warning('No perfect match for residue %s (%s): best candidate shares %s/%s atom names, %s proposed alternative(s)',
                        resname, unique_top_resname,
                        matching_ats, total_ats, len(candidates))
            for count in range(len(best_choices[unique_top_resname])):
                ff_topology_name = unique_top_resname + f"_uncertain_atoms_{count}.csv"
                possile_tops_str += f"{ff_topology_name} "
                best_choices[unique_top_resname][count][1].to_csv(possible_topologies_path / ff_topology_name, header = True, index = False)

            different_residue_topologies.loc[ix, "possible_atoms"] = possile_tops_str
            unmatched.append({'resname': resname, 'topology': unique_top_resname,
                              'matched_atoms': matching_ats, 'total_atoms': total_ats,
                              'n_candidates': len(candidates)})
        else:
            # residue neither perfectly matched nor with any proposed alternative
            complete_match = False
            log.warning('No forcefield candidate at all for residue %s (%s)',
                        resname, unique_top_resname)
            unmatched.append({'resname': resname, 'topology': unique_top_resname,
                              'matched_atoms': None, 'total_atoms': None,
                              'n_candidates': 0})

    if complete_match:
        oname = "complete_residue_topologies.csv"
        stale = "unmatching_residue_topologies.csv"
        columns = different_residue_topologies.columns.to_list()
        columns.remove('possible_atoms')
    else:
        oname = "unmatching_residue_topologies.csv"
        stale = "complete_residue_topologies.csv"
        columns = different_residue_topologies.columns.to_list()

    different_residue_topologies.to_csv(out / oname, header = True, columns=columns, index = False)
    log.debug('save_results: wrote %s (complete_match=%s)', oname, complete_match)

    # Only the outcome of THIS run may be on disk. A complete_residue_topologies.csv
    # left by an earlier run (or an earlier accept-topology) would let the
    # interactors stage run on stale associations after this run halted, and a
    # stale unmatching file would offer accept-topology candidates that no
    # longer exist.
    if (out / stale).exists():
        (out / stale).unlink()
        log.info('Removed %s left over from a previous run', stale)
    return complete_match, unmatched

def associate_ff_params(ff: Path, out: Path, nas: list, maxToCheck: int) -> Tuple[bool, list]:
    log.info('Associating topology residues to forcefield in: %s', ff)

    # 1 - load unique topologies found in input topology file
    try:
        restypes= pd.read_csv(out / 'different_residue_topologies.csv', header = 0, index_col = False,  na_values = nas,
                              keep_default_na = False)
    
    except FileNotFoundError:
        exit

    # 2 - Load interactors table (forcefield parameters)
    ff_params = pd.read_csv(ff / "final_ff_residues.csv", header=0, index_col = False, na_values = nas,
                            keep_default_na = False)
    ter_params = pd.read_csv(ff / "final_ff_ters.csv", header=0, index_col = False, na_values = nas,
                             keep_default_na = False)
    ff_params['bonds'] = ff_params['bonds'].fillna('')
    ter_params['bonds'] = ter_params['bonds'].fillna('')

    log.debug('Starting to associate atoms in topology to atoms in the forcefield')

    poss_names = [poss_rn for poss_rn in set(ff_params['resname'].values) if type(poss_rn) == str]
    log.debug('candidate forcefield residue names: %s', poss_names)

    ### ITERATE OVER RESIDUE NAMES IN TOPOLOGY
    # Start to compute results
    # only perfect matches on atom names will be considered as effectiva matches
    # otherwise the residue defined in the forcfield will be stored
    
    results = {}

    ### ITERATE OVER ALL POSSIBLE RESIDUE TOPOLOGIES (RESIDUES WITH SAME SET OF ATOMS) IN TOPOLOGY
    for at_path, topres_df in restypes.groupby(by = 'atoms'):
        top_residuename = topology_key(at_path)
        top_atoms_df = pd.read_csv(out / "resid_topologies" / at_path, header = 0, index_col = False, na_values = nas,
                                   keep_default_na = False)
        # results: dictionary:
        #  - keys: residue names found in topology
        #  - values: lists of tuples with entries
        #            (ff_rn, score, ff_atoms)
        #   Where:
        #            ff_resname = residue name in forcefield
        #            score: if = 1, perfect match of atomnames
        #                   if < 1, atom names do not match.
        #                       if = 0: all atom masses match
        #                       if < 0: a number of atomic masses = -score do not match
        #            ff_atoms = pandas DataFrame with extracted forcefield information
        #            

        results = get_results(results, top_residuename, top_atoms_df, ff_params, ter_params, poss_names, maxToCheck)
    '''
    found_residues: dict with:
     - keys: residue names found in topology (made unique for topology)
             e.g. if ALA at Nter there will be ALA and ALA_1, one referring 
                  to Nter residue topology and the other to intrachain ALA residue topoligy
    
     - values: a ff_residue topology has been found to match all the atom names in the structure_topology
               tuple with (ff_topology: pd.DataFrame)               
    best_choices: dict with:
     - keys: same of found_residues
     - values: no ff_residue topology has been found to perfectly match all the atom names in the structure_topology
               score is 0 if all masses match, else = -1 for each of the mismatching masses
                list of tuples (score: int, ff_topology: pd.DataFrame)
    '''
    found_residues, best_choices = collect_results(results)

    return save_results(out, found_residues, best_choices, nas)

## Functions used by `mypsns accept-topology`: promote a proposed candidate
## (possible_topologies/) to the residue's forcefield topology (found_topologies/)
## after recognition halted on residues with no exact match.

UNMATCHED_TABLE = 'unmatching_residue_topologies.csv'
COMPLETE_TABLE = 'complete_residue_topologies.csv'

# Columns of every per-residue forcefield CSV (as get_results / check_ters write
# them). A hand-edited replacement must keep all of them: the interactors stage
# reads bonds, charges and chargegroups from it and contacts reads masses and
# the LJ sigma/epsilon used for the energy matrices.
FF_ATOM_COLUMNS = ['rtp_file', 'resname', 'atom_name', 'atom_type', 'charge',
                   'chargegroup', 'bonds', 'sigma', 'epsilon', 'mass']

# Column of found_topologies/<key>_ff_atoms.csv, next to atom_name, with each
# atom's name in the input structure. atom_name (and bonds) keep the
# forcefield's names, so every atom can be traced to both. Recognition only
# writes exact matches (the two names are equal); accept-topology writes the
# names its pairing found. The interactors stage builds the whole topology in
# the structure's names (structure_named), because contacts joins it with the
# structure on (Resindex, atom_name). A file without the column (written
# before it existed, or copied by hand) is used as is.
STRUCTURE_NAME_COLUMN = 'structure_atom_name'

# Max |mass difference| (amu) for a renamed atom pair to count as the same atom.
# The structure's masses can be MDAnalysis guesses (O 15.999 against the
# forcefield's 15.9994), so exact equality is too strict, while the smallest gap
# that matters (H 1.008 against D 2.014) is far above it.
MASS_TOLERANCE = 0.1

def read_str_csv(path: Path, nas: list) -> pd.DataFrame:
    '''
    Read a topology CSV with every column as a string, so that rewriting it
    changes only the cells we edit (no float reformatting). keep_default_na=False
    keeps pandas' default "NA" from turning sodium names into NaN; only the
    configured na_values count as missing.
    '''
    return pd.read_csv(path, header = 0, index_col = False, dtype = str,
                       na_values = nas, keep_default_na = False)

def read_topology_table(path: Path, nas: list) -> pd.DataFrame:
    '''
    Read a *_residue_topologies.csv. The ff_atoms, possible_atoms and
    accepted_from columns come back as '' when empty (or absent:
    accepted_from only exists once accept-topology has run).
    '''
    table = read_str_csv(path, nas)
    for col in ('ff_atoms', 'possible_atoms', 'accepted_from'):
        if col not in table.columns:
            table[col] = ''
        table[col] = table[col].fillna('')
    return table

def candidate_files(row: pd.Series) -> List[str]:
    '''File names (in possible_topologies/) of the candidates proposed for one
    row of unmatching_residue_topologies.csv, in index order.'''
    return [f for f in str(row['possible_atoms']).split() if f]

def candidate_index(candidate_file: str) -> int:
    '''"MET_1_uncertain_atoms_3.csv" -> 3, the number accept-topology takes.'''
    return int(candidate_file[:-len('.csv')].rsplit('_', 1)[1])

def pair_candidate_atoms(top_atoms: pd.DataFrame, cand: pd.DataFrame) -> pd.DataFrame:
    '''
    Pair every atom of a forcefield candidate with an atom of the structure
    residue, the way recognition scored the candidate (check_namemasses and
    check_masses): atoms with the same name pair with each other, and the
    remaining candidate atoms pair with the remaining structure atoms BY
    POSITION, each side in its own file order.

     - top_atoms: resid_topologies/<key>_atoms.csv (columns name, mass)
     - cand: candidate forcefield topology (columns atom_name, mass, ...)

    Returns one row per candidate atom, in candidate order: ff_name, top_name,
    ff_mass, top_mass, renamed (the names differ) and mass_ok
    (|ff_mass - top_mass| <= MASS_TOLERANCE). Raises ValueError when the two
    sides cannot be paired at all (different atom counts, duplicated names).
    '''
    top_names = top_atoms['name'].astype(str).tolist()
    ff_names = cand['atom_name'].astype(str).tolist()
    if len(top_names) != len(ff_names):
        raise ValueError(f'the candidate has {len(ff_names)} atoms, '
                         f'the structure residue {len(top_names)}')
    for side, names in (('structure residue', top_names), ('candidate', ff_names)):
        duplicated = sorted({n for n in names if names.count(n) > 1})
        if duplicated:
            raise ValueError(f'duplicated atom names in the {side}: {duplicated}')

    top_mass = dict(zip(top_names, top_atoms['mass'].astype(float)))
    shared = set(top_names) & set(ff_names)
    # same counts and no duplicates -> both leftover lists have the same length
    top_left = iter([n for n in top_names if n not in shared])
    rows = []
    for ff_name, ff_mass in zip(ff_names, cand['mass'].astype(float)):
        top_name = ff_name if ff_name in shared else next(top_left)
        rows.append({'ff_name': ff_name, 'top_name': top_name,
                     'ff_mass': ff_mass, 'top_mass': top_mass[top_name]})
    pairs = pd.DataFrame(rows)
    pairs['renamed'] = pairs['ff_name'] != pairs['top_name']
    pairs['mass_ok'] = (pairs['ff_mass'] - pairs['top_mass']).abs() <= MASS_TOLERANCE
    return pairs

def rename_candidate_atoms(cand: pd.DataFrame, renames: dict) -> pd.DataFrame:
    '''
    Give the candidate's atoms the structure's names ({ff_name: structure_name}),
    in atom_name AND in every bonds string. The renamed topology is what the
    interactors stage stores, and contacts joins it with the structure on
    (Resindex, atom_name): an atom name the structure does not have makes that
    join fail.

    Bond tokens prefixed with "-" or "+" name an atom of the previous / next
    residue, whose names are not ours to change, so they are left alone. A
    bonds string with no renamed atom is kept byte for byte.
    '''
    renamed = cand.copy()
    renamed['atom_name'] = renamed['atom_name'].map(lambda n: renames.get(n, n))

    def _rename_bonds(bonds):
        if not isinstance(bonds, str):
            return bonds
        tokens = bonds.split()
        if not any(t in renames for t in tokens):
            return bonds
        return ' '.join(renames.get(t, t) for t in tokens)

    renamed['bonds'] = renamed['bonds'].map(_rename_bonds)
    return renamed

def with_structure_names(ff_atoms: pd.DataFrame, structure_names) -> pd.DataFrame:
    '''
    Copy of a forcefield topology with STRUCTURE_NAME_COLUMN set to
    structure_names (one per row, in row order) and placed right after
    atom_name. A column already there (e.g. in a CSV copied from
    found_topologies/) is replaced.
    '''
    out = ff_atoms.drop(columns = [STRUCTURE_NAME_COLUMN], errors = 'ignore')
    out.insert(out.columns.get_loc('atom_name') + 1, STRUCTURE_NAME_COLUMN,
               list(structure_names))
    return out

def structure_named(ff_atoms: pd.DataFrame) -> pd.DataFrame:
    '''
    A found topology as the later stages use it: atom_name, and every bonds
    entry naming those atoms, in the structure's names (rename_candidate_atoms
    with the STRUCTURE_NAME_COLUMN pairs), without that column. A topology
    without the column is returned unchanged.
    '''
    if STRUCTURE_NAME_COLUMN not in ff_atoms.columns:
        return ff_atoms
    renames = {ff: st for ff, st in zip(ff_atoms['atom_name'], ff_atoms[STRUCTURE_NAME_COLUMN])
               if isinstance(st, str) and st and st != ff}
    out = rename_candidate_atoms(ff_atoms, renames) if renames else ff_atoms.copy()
    return out.drop(columns = [STRUCTURE_NAME_COLUMN])

def unbonded_patch_atoms(ff_atoms: pd.DataFrame) -> List[str]:
    '''
    Names of the terminal-patch atoms (chargegroup -1, ff_reader's mark for
    patch rows) of a multi-atom residue that have no bond to another atom of
    the residue. Every atom a .tdb patch adds or replaces is bonded to the
    residue, so a non-empty result means the patch's bonds are wrong in
    final_ff_ters.csv (e.g. a hand-made NHT3+ whose N still lists H1-H3).
    '''
    if ff_atoms.shape[0] <= 1 or 'chargegroup' not in ff_atoms.columns:
        return []
    names = ff_atoms['atom_name'].astype(str).tolist()
    nameset = set(names)
    bonded = {an: False for an in names}
    for an, cell in zip(names, ff_atoms['bonds']):
        for t in bond_tokens(cell):
            if t in nameset and t != an:
                bonded[an] = bonded[t] = True
    patch = pd.to_numeric(ff_atoms['chargegroup'], errors = 'coerce') == -1
    return [an for an, is_patch in zip(names, patch) if is_patch and not bonded[an]]

def dangling_bonds(ff_atoms: pd.DataFrame) -> dict:
    '''
    {atom_name: [bond partners that are not atoms of this residue]}, ignoring
    "-" / "+" tokens (atoms of the previous / next residue).

    After a hand edit, a non-empty result usually means atom_name was changed
    but the bonds column was not. It is only a warning. Recognition's terminal
    patching drops such names (get_ter_ffres); topologies it wrote before
    that can still have them (e.g. an NHT3+ N-terminus whose N lists
    "H1 H2 H3" while its atoms are HT1-HT3).
    '''
    names = set(ff_atoms['atom_name'].astype(str))
    dangling = {}
    for atom, bonds in zip(ff_atoms['atom_name'].astype(str), ff_atoms['bonds']):
        if not isinstance(bonds, str):
            continue
        unknown = [b for b in bonds.split() if b[0] not in '-+' and b not in names]
        if unknown:
            dangling[atom] = unknown
    return dangling