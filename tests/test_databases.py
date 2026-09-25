import gzip
import hashlib
from pathlib import Path
import pytest
from funlr.databases import unpack_verified


def test_verified_unpack_preserves_exact_bytes(tmp_path):
    data=b'HMMER3/f\nNAME fixture\n//\n'
    source=tmp_path/'source.gz';source.write_bytes(gzip.compress(data))
    dest=tmp_path/'db/Pfam-A.hmm'
    unpack_verified(source,dest,hashlib.sha256(data).hexdigest(),len(data))
    assert dest.read_bytes()==data


def test_corrupt_or_wrong_release_never_installed(tmp_path):
    source=tmp_path/'source.gz';source.write_bytes(gzip.compress(b'wrong data'))
    dest=tmp_path/'db/Pfam-A.hmm'
    with pytest.raises(ValueError,match='mismatch'):
        unpack_verified(source,dest,'0'*64,10)
    assert not dest.exists()
    assert list(dest.parent.iterdir())==[]


def test_existing_database_untouched(tmp_path):
    dest=tmp_path/'Pfam-A.hmm';dest.write_bytes(b'original')
    with pytest.raises(ValueError,match='already exists'):
        unpack_verified(tmp_path/'absent.gz',dest,'0'*64,10)
    assert dest.read_bytes()==b'original'
