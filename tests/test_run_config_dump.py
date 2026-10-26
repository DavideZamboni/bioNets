"""<outpath>/run_config.yaml must keep describing the run after a single stage
re-dumps it. Each stage replaces only its own section(s); before, a stage run
from flags (e.g. `mypsns network -t … -o …`) rewrote every other section with
defaults and dropped the forcefield path, so the file claimed SB 4.5 Å for a
contacts run made at 3.5 Å and `mypsns validate` rejected it."""
import yaml

from mypsns.commands._helpers import dump_run_config
from mypsns.config import build_pipeline_config


def _cfg(tmp_path, **overrides):
    base = {'topology': str(tmp_path / 'system.gro'), 'outpath': str(tmp_path)}
    base.update(overrides)
    return build_pipeline_config(overrides=base, require_forcefield=False)


def _read(tmp_path):
    return yaml.safe_load((tmp_path / 'run_config.yaml').read_text())


def test_stage_dump_keeps_other_sections_and_paths(tmp_path):
    ff = str(tmp_path / 'charmm.ff')
    dump_run_config(_cfg(tmp_path, forcefield=ff,
                         contacts={'sb_cutoff': 3.5, 'hc_cutoff': 5.5}), tmp_path)
    # later: network from flags only (no forcefield, contacts left at defaults)
    dump_run_config(_cfg(tmp_path, network={'min_persistency': 20.0}), tmp_path,
                    sections=['network'])
    out = _read(tmp_path)
    assert out['forcefield'] == ff
    assert out['contacts']['sb_cutoff'] == 3.5 and out['contacts']['hc_cutoff'] == 5.5
    assert out['network']['min_persistency'] == 20.0
    # the merged file is still a valid config
    build_pipeline_config(config_path=tmp_path / 'run_config.yaml')


def test_paths_given_to_the_stage_win(tmp_path):
    dump_run_config(_cfg(tmp_path, trajectory=str(tmp_path / 'a.xtc')), tmp_path)
    dump_run_config(_cfg(tmp_path, trajectory=str(tmp_path / 'b.xtc')), tmp_path,
                    sections=['contacts'])
    assert _read(tmp_path)['trajectory'] == str(tmp_path / 'b.xtc')


def test_full_dump_without_sections(tmp_path):
    dump_run_config(_cfg(tmp_path, contacts={'sb_cutoff': 3.5}), tmp_path)
    dump_run_config(_cfg(tmp_path), tmp_path)        # sections=None: whole config
    assert _read(tmp_path)['contacts']['sb_cutoff'] == 4.5


def test_run_keeps_the_network_section(tmp_path):
    # `mypsns run` owns its four stages only: re-running it after `network`
    # must not reset what `network` recorded
    run_sections = ['ff_reader', 'topology_recognition', 'interactors', 'contacts']
    dump_run_config(_cfg(tmp_path, network={'min_persistency': 20.0}), tmp_path)
    dump_run_config(_cfg(tmp_path, contacts={'sb_cutoff': 3.5}), tmp_path,
                    sections=run_sections)
    out = _read(tmp_path)
    assert out['network']['min_persistency'] == 20.0
    assert out['contacts']['sb_cutoff'] == 3.5


def test_first_stage_dump_writes_a_complete_file(tmp_path):
    dump_run_config(_cfg(tmp_path), tmp_path, sections=['network'])
    out = _read(tmp_path)
    assert {'contacts', 'network', 'interactors'} <= set(out)
