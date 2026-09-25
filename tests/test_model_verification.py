"""The distributed snapshot check is strict, complete, and read-only."""
import hashlib
import json
from importlib.resources import files

import pytest

from funlr.databases import verify_models


@pytest.fixture
def libraries(tmp_path):
    data = {
        "NLR_ASM_combined.hmm": b"HMMER3/f\nNAME  ASM\n//\n",
        "NLR_NBD_combined.hmm": b"HMMER3/f\nNAME  NBD\n//\n",
        "NLR_effectors_combined.hmm": b"HMMER3/f\nNAME  effector\n//\n",
        "NLR_sensors_combined.hmm": b"HMMER3/f\nNAME  sensor\n//\n",
    }
    data["NLR_custom_combined.hmm"] = (
        data["NLR_effectors_combined.hmm"] + data["NLR_sensors_combined.hmm"]
    )
    records = []
    for name, content in sorted(data.items()):
        (tmp_path / name).write_bytes(content)
        records.append({"file": name, "bytes": len(content),
                        "sha256": hashlib.sha256(content).hexdigest(),
                        "profiles": content.count(b"NAME  ")})
    return tmp_path, {"libraries": records}


def test_complete_snapshot_check_does_not_write(libraries):
    directory, snapshot = libraries
    before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in directory.iterdir()}
    result = verify_models(directory, snapshot)
    assert result["status"] == "VERIFIED"
    assert len(result["libraries"]) == 5
    assert result["custom_concatenation_verified"] is True
    assert {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in directory.iterdir()} == before


@pytest.mark.parametrize("filename", ["NLR_ASM_combined.hmm", "NLR_NBD_combined.hmm"])
def test_required_library_cannot_be_omitted(libraries, filename):
    directory, snapshot = libraries
    (directory / filename).unlink()
    with pytest.raises(ValueError, match=f"Missing model library: {filename}"):
        verify_models(directory, snapshot)


def test_same_size_alteration_fails(libraries):
    directory, snapshot = libraries
    target = directory / "NLR_NBD_combined.hmm"
    target.write_bytes(target.read_bytes().replace(b"NBD", b"BAD"))
    with pytest.raises(ValueError, match="snapshot mismatch: NLR_NBD"):
        verify_models(directory, snapshot)


def test_concatenation_order_is_checked_independently(libraries):
    directory, snapshot = libraries
    custom = directory / "NLR_custom_combined.hmm"
    reversed_data = (directory / "NLR_sensors_combined.hmm").read_bytes() + (directory / "NLR_effectors_combined.hmm").read_bytes()
    custom.write_bytes(reversed_data)
    next(r for r in snapshot["libraries"] if r["file"] == custom.name)["sha256"] = hashlib.sha256(reversed_data).hexdigest()
    with pytest.raises(ValueError, match="effector-plus-sensor concatenation"):
        verify_models(directory, snapshot)


def test_snapshot_is_available_as_installed_package_resource():
    snapshot = json.loads(files("funlr").joinpath("data/model_snapshot.json").read_text())
    records = {row["file"]: row for row in snapshot["libraries"]}
    assert len(records) == 5
    assert records["NLR_custom_combined.hmm"]["profiles"] == 127
    assert records["NLR_NBD_combined.hmm"]["sha256"] == "230982ccb1e008aee682d9469eb525a407d007ec104ba5ebb287765b1f276433"
