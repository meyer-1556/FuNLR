"""Fetch the fixed Pfam release verified against the reference input."""
from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
import shutil
import tempfile
from urllib.request import urlopen

PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/Pfam-A.hmm.gz"
PFAM_SHA256 = "4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e"
PFAM_BYTES = 2246909846


def unpack_verified(source, destination, expected_sha256, expected_bytes):
    """Verify before installation; refuse overwrite and clean up failed writes."""
    destination = Path(destination)
    if destination.exists() or destination.is_symlink():
        raise ValueError(f"Destination already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, prefix=".pfam-", delete=False) as out:
        partial = Path(out.name)
        digest = hashlib.sha256()
        size = 0
        try:
            with gzip.open(source, "rb") as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    size += len(block)
                    if size > expected_bytes:
                        raise ValueError("Pfam exceeds the expected uncompressed size")
                    digest.update(block)
                    out.write(block)
            if size != expected_bytes or digest.hexdigest() != expected_sha256:
                raise ValueError("Pfam checksum/size mismatch; database was not installed")
            out.flush()
        except BaseException:
            partial.unlink(missing_ok=True)
            raise
    try:
        # Exclusive creation also handles another process completing meanwhile.
        with destination.open("xb") as target, partial.open("rb") as source_file:
            try:
                shutil.copyfileobj(source_file, target, 1024 * 1024)
            except BaseException:
                destination.unlink(missing_ok=True)
                raise
    finally:
        partial.unlink(missing_ok=True)


def fetch_pfam(directory):
    from .runner import sha256, write_json
    directory = Path(directory).expanduser().resolve()
    target = directory / "Pfam-A.hmm"
    if target.exists() or target.is_symlink():
        if target.is_symlink() or not target.is_file() or target.stat().st_size != PFAM_BYTES or sha256(target) != PFAM_SHA256:
            raise ValueError(f"Existing Pfam file does not match release38.2: {target}")
        return {"status": "VERIFIED_EXISTING", "path": str(target), "sha256": PFAM_SHA256}
    directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".pfam-download-", dir=directory) as working:
        compressed = Path(working) / "Pfam-A.hmm.gz"
        with urlopen(PFAM_URL, timeout=120) as response, compressed.open("wb") as out:
            if not response.geturl().startswith("https://"):
                raise ValueError("Pfam download was redirected away from HTTPS")
            shutil.copyfileobj(response, out, 1024 * 1024)
        unpack_verified(compressed, target, PFAM_SHA256, PFAM_BYTES)
    report = {"status": "DOWNLOADED_AND_VERIFIED", "release": "38.2", "url": PFAM_URL,
              "path": str(target), "bytes": PFAM_BYTES, "sha256": PFAM_SHA256}
    write_json(directory / "Pfam38.2.provenance.json", report)
    return report


def verify_models(directory, snapshot=None):
    """Read-only verification of all five libraries in the fixed reference snapshot."""
    import json
    from importlib.resources import files
    if snapshot is None:
        snapshot = json.loads(files('funlr').joinpath('data/model_snapshot.json').read_text())
    directory = Path(directory).expanduser().resolve()
    if not directory.is_dir():
        raise ValueError(f'Model directory does not exist: {directory}')
    records = []
    errors = []
    for expected in snapshot['libraries']:
        name = expected['file']
        if Path(name).name != name:
            raise ValueError('Model snapshot contains an invalid filename')
        path = directory / name
        if not path.is_file():
            errors.append(f'Missing model library: {name}')
            continue
        digest = hashlib.sha256()
        size = 0
        profiles = 0
        with path.open('rb') as handle:
            for line in handle:
                digest.update(line)
                size += len(line)
                profiles += line.startswith(b'NAME  ')
        actual = {'file': name, 'sha256': digest.hexdigest(), 'bytes': size, 'profiles': profiles}
        records.append(actual)
        if any(actual[key] != expected[key] for key in ('sha256', 'bytes', 'profiles')):
            errors.append(f'Model snapshot mismatch: {name}')
    if errors:
        raise ValueError('; '.join(errors))
    combined = hashlib.sha256()
    for name in ('NLR_effectors_combined.hmm', 'NLR_sensors_combined.hmm'):
        with (directory / name).open('rb') as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                combined.update(block)
    expected_custom = next(row['sha256'] for row in records if row['file'] == 'NLR_custom_combined.hmm')
    if combined.hexdigest() != expected_custom:
        raise ValueError('Custom library is not the verified effector-plus-sensor concatenation')
    return {'status': 'VERIFIED', 'directory': str(directory), 'libraries': records,
            'custom_concatenation_verified': True,
            'scope': 'Byte identity and composition of this snapshot; not a grant of redistribution rights or a rebuild of the original component models.'}
