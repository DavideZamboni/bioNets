"""NAMED_INTERACTOR_DEFAULTS has a single source of truth (mypsns/config.py).

It used to be duplicated in utils/interactors_utilities.py and kept in sync by
hand. The utilities module now imports it; these tests make sure nobody
re-introduces a second copy, and that the config default cannot alias it.
"""
import copy

from mypsns import config
from mypsns.config import NAMED_INTERACTOR_DEFAULTS, InteractorsConfig
from mypsns.utils import interactors_utilities


def test_utilities_use_the_config_object():
    # Same object, not merely an equal copy: an equal copy could drift later.
    assert interactors_utilities.NAMED_INTERACTOR_DEFAULTS is config.NAMED_INTERACTOR_DEFAULTS


def test_defaults_cover_every_named_role():
    roles = set(interactors_utilities.NAMED_INTERACTOR_ROLES)
    assert roles <= set(NAMED_INTERACTOR_DEFAULTS)
    assert set(NAMED_INTERACTOR_DEFAULTS) - roles == {"hb_masses"}


def test_config_default_equals_module_defaults():
    assert InteractorsConfig().named_interactors_rules == NAMED_INTERACTOR_DEFAULTS


def test_config_default_is_a_deep_copy():
    # Editing one config's rules in place must not leak into the module
    # constant (and so into every later config / get_named_interactors call).
    before = copy.deepcopy(NAMED_INTERACTOR_DEFAULTS)
    cfg = InteractorsConfig()
    cfg.named_interactors_rules["HC"][0]["residue"].append("GLY")
    cfg.named_interactors_rules["hb_masses"].append(2)
    assert NAMED_INTERACTOR_DEFAULTS == before
    assert InteractorsConfig().named_interactors_rules == before
