"""Profile identity and source evidence remain separate from names and roles."""
import hashlib
import json
from importlib.resources import files

import pytest

from funlr import model_inventory


def hmm(name="sensor", accession=None, nseq=7):
    text = f"HMMER3/f [3.4]\nNAME  {name}\n"
    if accession:
        text += f"ACC   {accession}\n"
    text += "LENG  12\n"
    if nseq is not None:
        text += f"NSEQ  {nseq}\n"
    return (text + "HMM          A C D\n  fake numerical body\n//\n").encode()


def test_inventory_preserves_duplicate_parameters_and_exact_block_identity(tmp_path):
    first, second = hmm("TPR", "PF00001.2"), hmm("TPR__dup2", "PF00001.2")
    path = tmp_path / "test.hmm"
    path.write_bytes(first + second)
    before = (path.read_bytes(), path.stat().st_mtime_ns)
    result = model_inventory.inventory_models(tmp_path)
    assert result["schema_version"] == 1
    library = result["libraries"][0]
    assert library["sha256"] == hashlib.sha256(first + second).hexdigest()
    assert library["profile_count"] == 2
    assert library["identical_except_name_extra_records"] == 1
    assert [r["name"] for r in library["profiles"]] == ["TPR", "TPR__dup2"]
    profile = library["profiles"][0]
    assert (profile["accession"], profile["length"], profile["nseq"]) == ("PF00001.2", 12, 7)
    assert profile["sha256"] == hashlib.sha256(first).hexdigest()
    assert profile["publication_profile_identity"] is None
    assert library["role"] is None and library["source_audit"] is None
    assert library["redistribution_permission"] is None
    assert (path.read_bytes(), path.stat().st_mtime_ns) == before
    assert sorted(p.name for p in tmp_path.iterdir()) == ["test.hmm"]


def test_same_filename_does_not_inherit_known_source_without_matching_bytes(tmp_path, monkeypatch):
    path = tmp_path / "known.hmm"
    path.write_bytes(hmm())
    source = {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
              "role": "test role", "citations": ["https://example.org/source"],
              "source_audit": {"historical_reconstruction": None}}
    monkeypatch.setattr(model_inventory, "_metadata", lambda: {"metadata_version": "test", "libraries": [source]})
    known = model_inventory.inventory_models(tmp_path)["libraries"][0]
    assert known["snapshot_identity"] == "RECOGNIZED"
    assert known["role"] == "test role"
    path.write_bytes(hmm("changed"))
    unknown = model_inventory.inventory_models(tmp_path)["libraries"][0]
    assert unknown["snapshot_identity"] == "UNRECOGNIZED"
    assert unknown["citations"] is None and unknown["role"] is None


def test_absent_optional_profile_headers_are_null(tmp_path):
    (tmp_path / "test.hmm").write_bytes(hmm(nseq=None))
    profile = model_inventory.inventory_models(tmp_path)["libraries"][0]["profiles"][0]
    assert profile["accession"] is None and profile["nseq"] is None


@pytest.mark.parametrize("data,match", [
    (b"", "Empty or unterminated"),
    (hmm().replace(b"//\n", b""), "unterminated"),
    (hmm().replace(b"NAME  sensor\n", b""), "Missing NAME/LENG"),
    (hmm().replace(b"LENG  12", b"LENG  -1"), "Invalid LENG"),
    (hmm().replace(b"NSEQ  7", b"NSEQ  abc"), "Invalid NSEQ"),
    (hmm().replace(b"NAME  sensor", b"NAME  sensor\nNAME  other"), "Invalid NAME"),
    (b"unexpected\n" + hmm(), "Unexpected text"),
    (hmm().replace(b"//\n", b"") + hmm(), "Unterminated"),
])
def test_malformed_model_cannot_produce_a_plausible_inventory(tmp_path, data, match):
    (tmp_path / "bad.hmm").write_bytes(data)
    with pytest.raises(ValueError, match=match):
        model_inventory.inventory_models(tmp_path)


def test_inventory_rejects_symlink_libraries(tmp_path):
    target = tmp_path / "actual.txt"
    target.write_bytes(hmm())
    (tmp_path / "linked.hmm").symlink_to(target)
    with pytest.raises(ValueError, match="non-symlink"):
        model_inventory.inventory_models(tmp_path)


def test_optional_snapshot_check_is_the_existing_fixed_verifier(tmp_path, monkeypatch):
    from funlr import databases
    (tmp_path / "test.hmm").write_bytes(hmm())
    seen = []
    def fail(directory):
        seen.append(directory)
        raise ValueError("fixed snapshot mismatch")
    monkeypatch.setattr(databases, "verify_models", fail)
    with pytest.raises(ValueError, match="fixed snapshot mismatch"):
        model_inventory.inventory_models(tmp_path, verify_snapshot=True)
    assert seen == [tmp_path]


def test_source_metadata_is_packaged_and_keeps_unknown_history_null():
    data = json.loads(files("funlr").joinpath("data/model_sources.json").read_text())
    assert data["schema_version"] == 1
    assert len(data["libraries"]) == 6
    for row in data["libraries"]:
        assert len(row["sha256"]) == 64
        assert row["source_audit"]["per_profile_training_sequences"] is None
        assert row["source_audit"]["redistribution_permission"] is None


def test_snapshot_success_cannot_be_combined_with_different_inventory_bytes(tmp_path, monkeypatch):
    from funlr import databases
    path = tmp_path / "test.hmm"
    original = hmm()
    path.write_bytes(original)
    def mutate_after_check(directory):
        path.write_bytes(hmm("renamed"))
        return {"libraries": [{"file": path.name, "sha256": hashlib.sha256(original).hexdigest(),
                               "bytes": len(original), "profiles": 1}]}
    monkeypatch.setattr(databases, "verify_models", mutate_after_check)
    with pytest.raises(ValueError, match="changed during snapshot inventory"):
        model_inventory.inventory_models(tmp_path, verify_snapshot=True)
