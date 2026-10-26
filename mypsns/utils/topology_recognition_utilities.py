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

def get_ter_ffres(ter_residue: pd.DataFrame, ff_residue: pd.DataFrame,to_be_replaced: dict, to_be_removed: List) -> pd.DataFrame:
    #print('All atoms to be replaced or removed are present in forcefield residue')
    tmp = ff_residue.copy()
    tmp = tmp[~tmp['atom_name'].isin(to_be_removed)]  # now correctly filtered
    rep = []
    for k in to_be_replaced.keys():
        rep.append(k)
        if k != to_be_replaced[k]:
            to_be_removed.append(to_be_replaced[k])

    log.debug('ter replace: to_be_replaced=%s to_be_removed=%s', to_be_replaced, to_be_removed)
    atnames = tmp['atom_name'].to_list()

    bonds_to_edit = {}
    for an in to_be_replaced.keys():
        myat = tmp[tmp['atom_name'] == to_be_replaced[an]]
        terat = ter_residue[ter_residue['atom_name'] == an]
        cols = [
            'resname','atom_name','atom_type','mass',
            'charge','chargegroup','sigma','epsilon'
        ]
        # substitute columns of replaced atoms (no bonds, need to be adjusted)
        tmp.loc[myat.index, cols] = terat[cols].values

        at_b_entry = tmp.at[myat.index[0], 'bonds']
        ter_b_entry = ter_residue.at[terat.index[0], 'bonds']
        resbonds = [b for b in at_b_entry.split(' ') if b != '']
        terbonds = [b for b in ter_b_entry.split(' ') if b != '']

        # also removing e.g. +N when attempting to build a Cter
        truebonds = [b for b in resbonds if b not in to_be_removed and not b.startswith('+')]
        for b in terbonds:
            if (b not in truebonds) and (b[0] not in ['-','r']):
                truebonds.append(b)
                
        for b in truebonds:
            if b in ter_residue['atom_name'].to_list() and b not in atnames:
                bonds_to_edit[b] = an

        mybs = ''
        if len(truebonds) > 1:
            mybs += ' '.join(truebonds) + ' '
        else:
            if len(truebonds)  == 0:
                truebonds.append('')
                
            mybs += f'{truebonds[0]} '
        tmp.at[myat.index[0], 'bonds'] = mybs
    
    atnames = tmp['atom_name'].to_list()
    # now need to add atoms generated by ter definition to ff_residue
    # before need to adjust bonding info according to replaced atoms
    to_add = ter_residue[~ter_residue['atom_name'].isin(atnames)].copy()
    log.debug('ter replace: atoms to be added: %s', to_add['atom_name'].to_list())

    for ix, row in to_add.iterrows():
        atname = to_add.at[ix, 'atom_name']
        if atname in bonds_to_edit.keys():
            to_add.at[ix, 'bonds'] = bonds_to_edit[atname]

    to_add.rename(columns={'tdb_file': 'rtp_file'}, inplace=True)
    tmp = pd.concat([tmp, to_add], ignore_index=True)

    return tmp

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
    different_residue_topologies = pd.read_csv(out / 'different_residue_topologies.csv', header = 0, index_col = False)
    
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
            found_residues[unique_top_resname].to_csv(found_topologies_path / ff_topology_name, header = True, index = False)

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
        restypes= pd.read_csv(out / 'different_residue_topologies.csv', header = 0, index_col = False,  na_values = nas)
    
    except FileNotFoundError:
        exit

    # 2 - Load interactors table (forcefield parameters)
    ff_params = pd.read_csv(ff / "final_ff_residues.csv", header=0, index_col = False, na_values = nas)
    ter_params = pd.read_csv(ff / "final_ff_ters.csv", header=0, index_col = False, na_values = nas)
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
        top_atoms_df = pd.read_csv(out / "resid_topologies" / at_path, header = 0, index_col = False, na_values = nas)
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

def dangling_bonds(ff_atoms: pd.DataFrame) -> dict:
    '''
    {atom_name: [bond partners that are not atoms of this residue]}, ignoring
    "-" / "+" tokens (atoms of the previous / next residue).

    After a hand edit, a non-empty result usually means atom_name was changed
    but the bonds column was not. It is only a warning: recognition's own
    terminal patching can leave such names too (e.g. an NHT3+ N-terminus whose
    N lists "H1 H2 H3" while its atoms are HT1-HT3).
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