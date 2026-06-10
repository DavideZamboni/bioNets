from pathlib import Path
from typing import List
from pathlib import Path
import logging
from collections import deque

log = logging.getLogger(__name__)

def guess_ff(
        inpath: Path,
        nonbonded_name: str = "ffnonbonded.itp",
    ) -> bool:

    rtps = [f for f in inpath.iterdir() if f.is_file() and f.name.endswith('.rtp')]
    if len(rtps)  < 1:
        log.debug(f"No residue topology (*.rtp) files found in {inpath.resolve()}")
        return False
    
    tdbs = [f for f in inpath.iterdir() if f.is_file() and f.name.endswith('.tdb')]
    if len(tdbs)  < 1:
        log.debug(f"No termini database (*.tdb) files found in {inpath.resolve()}")
        return False

    nonbonded = [f for f in inpath.iterdir() if f.is_file() and f.name.endswith(nonbonded_name)]
    if len(nonbonded)  < 1:
        log.debug(f"File {nonbonded} is not present in {inpath.resolve()}")
        return False

    return True

def guess_ff_path(
    parent_dir:     str,
    recursive:      bool = False,
    max_depth:      int  = 10,
    nonbonded_name: str  = "ffnonbonded.itp",
) -> List[Path]:
    '''
    Log (and return) every directory at/under `parent_dir` whose contents match
    a GROMACS-style forcefield layout (passes `guess_ff`: ≥1 *.rtp, ≥1 *.tdb,
    and a file ending with `nonbonded_name`).

    - `parent_dir` itself is always checked.
    - With `recursive=True`, descend BFS to a maximum nesting of `max_depth`
      levels below `parent_dir` (so `max_depth=1` checks parent + direct
      children, `max_depth=10` is up to ten levels below).
    '''
    inpath = Path(parent_dir).expanduser()
    found: List[Path] = []

    if not inpath.exists() or not inpath.is_dir():
        log.warning("guess_ff_path: %s does not exist or is not a directory", inpath)
        return found

    # check the parent itself
    if guess_ff(inpath, nonbonded_name=nonbonded_name):
        log.info("guess_ff_path: forcefield layout at %s", inpath.resolve())
        found.append(inpath)

    if not recursive:
        max_depth = 1

    # BFS over subdirectories, depth-limited (depth = levels below `inpath`)
    queue: deque[tuple[Path, int]] = deque([(inpath, 0)])
    while queue:
        d, depth = queue.popleft()
        if depth >= max_depth:
            continue
        try:
            children = sorted(c for c in d.iterdir() if c.is_dir())
        except (PermissionError, OSError) as e:
            log.debug("guess_ff_path: cannot list %s (%s)", d, e)
            continue
        for child in children:
            if guess_ff(child, nonbonded_name=nonbonded_name):
                log.info("guess_ff_path: forcefield layout at %s", child.resolve())
                found.append(child)
            queue.append((child, depth + 1))

    if not found:
        log.info("guess_ff_path: no forcefield layout found under %s", inpath.resolve())
    return found