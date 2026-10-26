"""get_charge_based_interactors takes its defaults from InteractorsConfig.

Its own defaults dict lacked min_bonds_to_polar, so a call without that kwarg
raised KeyError (the CLI always passed it, so nothing showed). Each parameter
is also an `interactors` flag and is forwarded by both `interactors` and
`run`.
"""
import inspect

import pytest

from mypsns.commands import interactors_cmd, run_cmd
from mypsns.config import InteractorsConfig
from mypsns.utils.interactors_utilities import (
    CHARGE_BASED_PARAMETERS, add_hbonds, add_hc_atoms, add_sb_atoms,
    get_charge_based_interactors,
)


def hc_partition(df):
    hc = df[df['HC_cluster'] != '']
    return sorted(sorted(g['atom_name']) for _, g in hc.groupby('HC_cluster'))


def test_runs_with_no_kwargs(phe_variant):
    df = get_charge_based_interactors(phe_variant)
    assert 'HC_cluster' in df.columns and hc_partition(df)


def test_parameters_are_config_fields():
    assert set(CHARGE_BASED_PARAMETERS) <= set(InteractorsConfig.model_fields)


@pytest.mark.parametrize('func', [add_hbonds, add_sb_atoms, add_hc_atoms])
def test_every_typing_parameter_is_covered(func):
    params = set(inspect.signature(func).parameters) - {'df'}
    assert params <= set(CHARGE_BASED_PARAMETERS), params - set(CHARGE_BASED_PARAMETERS)


def test_kwargs_override_and_others_are_ignored(phe_variant):
    default = get_charge_based_interactors(phe_variant.copy())
    one_cluster = get_charge_based_interactors(
        phe_variant.copy(), max_num_heavy_ats=7,
        named_rules={'HC': []}, min_residue_charge=0.8)    # the named scheme's: ignored
    assert len(hc_partition(default)) == 3
    assert len(hc_partition(one_cluster)) == 1


@pytest.mark.parametrize('param', list(CHARGE_BASED_PARAMETERS) + ['min_residue_charge'])
def test_cli_flag_and_forwarding(param):
    assert param in inspect.signature(interactors_cmd.interactors).parameters
    for module in (interactors_cmd, run_cmd):
        src = inspect.getsource(module)
        assert f'{param}=cfg.interactors.{param}' in src, f'{module.__name__} does not forward {param}'
