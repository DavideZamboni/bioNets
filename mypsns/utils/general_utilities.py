from typing import Optional, Tuple, Union
import timeit
import signal
from contextlib import contextmanager

def get_endtime(start: float, tstart: Optional[float] = None) -> Union[Tuple[str, str], str]:
    """Function that returns string hh:mm:ssss of execution time of the process. Tuple of times if 2 are given, 1st process start time, 2nd whole start time"""

    end = timeit.default_timer()
    script_time = end - start
    # output running time in a nice format.
    mins, secs = divmod(script_time, 60)
    hours, mins = divmod(mins, 60)

    my_t = f"{hours}:{mins}:{secs:4f}"
    
    if tstart is not None:
        total_time = end - tstart
        # output running time in a nice format.
        tmins, tsecs = divmod(total_time, 60)
        thours, tmins = divmod(mins, 60)
        tot_t = f"{thours}:{tmins}:{tsecs:4f}"
        return my_t, tot_t
    return my_t


def timeout(seconds: int):
    def handler(signum, frame):
        raise TimeoutError(f"Timed out after {seconds}s")

    old = signal.signal(signal.SIGALRM, handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)

def fix_h_naming(pdb_path, out_path):
    # Parse all lines, track atom names per residue
    lines = open(pdb_path).readlines()
    
    # Group atom names by (chain, resid)
    from collections import defaultdict
    res_atoms = defaultdict(set)
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            chain = line[21]
            resid = line[22:27].strip()
            name = line[12:16].strip()
            res_atoms[(chain, resid)].add(name)
    
    out = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            chain = line[21]
            resid = line[22:27].strip()
            name = line[12:16].strip()
            
            # Check if name ends with '3' and the '1' variant is missing
            if name.endswith('3') and len(name) >= 3:
                name1 = name[:-1] + '1'
                atoms = res_atoms[(chain, resid)]
                if name1 not in atoms:
                    # Replace: right-justify in columns 13-16 (PDB format)
                    new_name = name1
                    # PDB atom name field is cols 13-16 (0-indexed: 12:16)
                    if len(new_name) < 4:
                        formatted = f' {new_name:<3s}'
                    else:
                        formatted = f'{new_name:<4s}'
                    line = line[:12] + formatted + line[16:]
        out.append(line)
    
    with open(out_path, 'w') as f:
        f.writelines(out)
'''op = "/home/davide/storage/lab_dati/structures"
inp = op + "/old"
name = "/wt_3ca.pdb"
fix_h_naming(inp + name, op + name)'''