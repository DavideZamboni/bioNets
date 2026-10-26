"""`forcefield` is optional for network/compare-variants only; everything else still requires it."""
import pytest

from mypsns.config import ConfigError, build_pipeline_config

BASE = {"topology": "s.pdb", "outpath": "out"}


def test_forcefield_required_by_default():
    with pytest.raises(ConfigError, match="forcefield: Field required"):
        build_pipeline_config(overrides=BASE)


def test_forcefield_optional_for_network_stages():
    cfg = build_pipeline_config(overrides=BASE, require_forcefield=False)
    assert cfg.forcefield is None


def test_forcefield_still_used_when_given():
    cfg = build_pipeline_config(overrides={**BASE, "forcefield": "ff"})
    assert str(cfg.forcefield) == "ff"
