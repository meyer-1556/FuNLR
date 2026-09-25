#!/usr/bin/env python3
"""Reconstruct only the verified effector+sensor composition; never overwrite."""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--models-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    metadata=json.loads((Path(__file__).resolve().parents[1]/'docs/MODEL_SNAPSHOT.json').read_text())
    expected={r['file']:r['sha256'] for r in metadata['libraries']}
    pieces=[]
    for name in ('NLR_effectors_combined.hmm','NLR_sensors_combined.hmm'):
        data=(args.models_dir/name).read_bytes()
        if hashlib.sha256(data).hexdigest()!=expected[name]:
            raise SystemExit('Component checksum mismatch: '+name)
        pieces.append(data)
    data=b''.join(pieces)
    actual=hashlib.sha256(data).hexdigest()
    if actual!=expected['NLR_custom_combined.hmm']:
        raise SystemExit('Reconstructed custom library checksum mismatch')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('xb') as stream:stream.write(data)
    print(json.dumps({'output':str(args.output),'sha256':actual,'bytes':len(data),'recipe':'effectors then sensors; verified snapshot'},indent=2))

if __name__=='__main__':main()
