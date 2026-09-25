#!/usr/bin/env python3
"""Copy the canonical source example into package resources, or verify it.

Maintainer command; does not download data. Copies only the listed files.
The duplication keeps the browsable source example and installed wheel usable.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
NAMES = ("genome.fa", "proteins.faa", "annotation.gff3", "emapper.annotations", "db/nbd.hmm",
         "db/pfam_mini.hmm", "config.yaml", "README.md", "PROVENANCE.json", "DATA_LICENSES.md")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Check resource equality without writing files")
    args = parser.parse_args()
    source = ROOT / "examples/demo"
    target = ROOT / "src/funlr/data/public_sequences"
    contents = {name: (source / name).read_bytes() for name in NAMES}
    hashes = {name: hashlib.sha256(data).hexdigest() for name, data in contents.items()}
    manifest = json.dumps({"schema_version": 1, "source": "examples/demo", "sha256": hashes}, indent=2, sort_keys=True) + "\n"
    contents["assets.json"] = manifest.encode()
    if args.check:
        missing = [name for name, data in contents.items() if not (target / name).is_file() or (target / name).read_bytes() != data]
        if missing:
            print("Stale/missing packaged example resources: " + ", ".join(missing), file=sys.stderr)
            return 1
    else:
        for name, data in contents.items():
            dest = target / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(data)
    print("Public example resources match source example." if args.check else "Public example resources copied and checksummed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
