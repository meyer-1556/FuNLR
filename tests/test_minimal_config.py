"""The small template must preserve preset semantics and resolved defaults."""
import yaml

from funlr.config import FunlrConfig, write_example_yaml


def test_minimal_and_full_illumina_resolve_to_same_science():
    minimal = yaml.safe_load(write_example_yaml(minimal=True))
    full = yaml.safe_load(write_example_yaml())
    assert set(minimal) == {"inputs", "databases", "execution", "context"}
    assert FunlrConfig(minimal).scientific_settings() == FunlrConfig(full).scientific_settings()


def test_minimal_profile_switch_reapplies_unwritten_hifi_defaults():
    minimal = yaml.safe_load(write_example_yaml(minimal=True))
    minimal["context"]["profile"] = "HIFI"
    actual = FunlrConfig(minimal).scientific_settings()
    expected = FunlrConfig({"context": {"profile": "HIFI"}}).scientific_settings()
    assert actual == expected
    assert actual != FunlrConfig().scientific_settings()


def test_minimal_hifi_roundtrip_preserves_explicit_threshold(tmp_path):
    minimal = yaml.safe_load(write_example_yaml(profile="HIFI", minimal=True))
    minimal["thresholds"] = {"min_nbd_ali_len": 151}
    path = tmp_path / "settings.yaml"
    path.write_text(yaml.safe_dump(minimal))
    actual = FunlrConfig.from_yaml(path)
    assert actual.scientific_settings()["MIN_NBD_ALI_LEN"] == 151
    assert actual.scientific_settings()["NLR_PROFILE"] == "HIFI"
