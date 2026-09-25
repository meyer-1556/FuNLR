"""Public example packaging, provenance, and failure boundaries."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

import pytest

from funlr import public_demo

ROOT = Path(__file__).resolve().parents[1]
RESOURCES = Path(public_demo.__file__).parent / 'data/public_sequences'


def test_packaged_assets_equal_canonical_example():
    result = subprocess.run([sys.executable, str(ROOT / 'scripts/sync_public_demo.py'), '--check'], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assets, manifest = public_demo.verified_assets()
    assert public_demo.INPUTS <= set(assets)
    assert all(hashlib.sha256(data).hexdigest() == manifest['sha256'][name] for name, data in assets.items())


def test_resource_corruption_fails_before_demo_output_is_created(tmp_path, monkeypatch):
    fixture = tmp_path / 'corrupt'
    shutil.copytree(RESOURCES, fixture)
    (fixture / 'proteins.faa').write_text('>replaced\nAAA\n')
    with pytest.raises(ValueError, match='checksum mismatch: proteins.faa'):
        public_demo.verified_assets(fixture)
    monkeypatch.setattr(public_demo, 'files', lambda package: tmp_path)
    with pytest.raises(OSError):
        public_demo.run_public_demo(tmp_path / 'output')
    assert not (tmp_path / 'output').exists()


def test_asset_manifest_cannot_write_parent_paths(tmp_path):
    manifest = {'schema_version': 1, 'sha256': {name: '0'*64 for name in sorted(public_demo.INPUTS)}}
    manifest['sha256'] = {'../escape': '0'*64, **manifest['sha256']}
    (tmp_path / 'assets.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='Invalid path'):
        public_demo.verified_assets(tmp_path)


def test_existing_output_is_preserved_without_probing_tools(tmp_path, monkeypatch):
    marker = tmp_path / 'keep';marker.write_text('unchanged')
    monkeypatch.setattr(public_demo, 'verified_assets', lambda: pytest.fail('must reject existing output first'))
    with pytest.raises(ValueError, match='new directory'):
        public_demo.run_public_demo(tmp_path)
    assert marker.read_text() == 'unchanged'


def test_missing_tools_fail_without_creating_output(tmp_path, monkeypatch):
    monkeypatch.setattr(public_demo, 'check_installation', lambda: {'ok': False, 'errors': ['miniprot is unavailable']})
    out = tmp_path / 'new'
    with pytest.raises(ValueError, match='miniprot is unavailable'):
        public_demo.run_public_demo(out)
    assert not out.exists()


def test_public_fixture_claims_have_sequence_and_hmm_evidence():
    assets, _ = public_demo.verified_assets()
    provenance = json.loads(assets['PROVENANCE.json'])
    for name, details in provenance['files'].items():
        assert len(assets[name]) == details['bytes']
        assert hashlib.sha256(assets[name]).hexdigest() == details['sha256']
    records = {}
    for line in assets['proteins.faa'].decode().splitlines():
        if line.startswith('>'):
            key=line[1:].split()[0];records[key]=''
        else: records[key]+=line.strip()
    assert len(records) == 7 and len(set(records.values())) == 6
    assert records['HETE_AAL37299'] == records['HETE_Q8X1P4']
    nbd = [line.split()[1] for line in assets['db/nbd.hmm'].decode().splitlines() if line.startswith('ACC ')]
    assert nbd == ['PF05729.19', 'PF00931.29']
    assert assets['db/pfam_mini.hmm'].decode().count('\nACC ') == 16


def test_public_example_config_resolves_outside_checkout(tmp_path, monkeypatch):
    from funlr.config import FunlrConfig
    monkeypatch.chdir(tmp_path)
    cfg = FunlrConfig.from_yaml(str(ROOT / 'examples/demo/config.yaml'))
    for value in cfg['inputs'].values():
        if value is not None: assert Path(value).is_file()
    for value in cfg['databases'].values():
        if value is not None: assert Path(value).is_file()
    for section, settings in public_demo.DEMO_SETTINGS.items():
        for name, expected in settings.items():
            assert cfg[section][name] == expected


def test_demo_cli_exposes_both_datasets():
    from funlr.cli import parser
    assert parser().parse_args(['demo','--outdir','x']).dataset == 'synthetic'
    assert parser().parse_args(['demo','--dataset','public-sequences','--outdir','x']).dataset == 'public-sequences'


def test_public_checker_accepts_observed_native_refinement_and_rejects_missing_mapping(tmp_path):
    run = tmp_path / 'run'
    shutil.copytree(ROOT / 'tests/fixtures/public_demo_v3', run)
    counts = public_demo.verify_outputs(run)
    assert counts['candidate_count'] == 6
    assert counts['stage6_refine_count'] == 1
    assert counts['stage6_comprehensive_refine_count'] == 3
    summary = run / 'results/stage6/comprehensive/tables/exonerate_summary.tsv'
    lines = summary.read_text().splitlines()
    summary.write_text('\n'.join(lines[:-1]) + '\n')
    with pytest.raises(RuntimeError, match='comprehensive Exonerate mapping mismatch'):
        public_demo.verify_outputs(run)


@pytest.mark.skipif(os.name != 'posix', reason='POSIX scheduler termination and process groups')
@pytest.mark.parametrize('dataset', ['synthetic', 'public-sequences'])
def test_terminated_demo_cli_stops_its_subprocess(tmp_path, dataset):
    """A scheduler SIGTERM must reach cleanup despite detached child groups."""
    ready = tmp_path / 'child.pid'
    child_code = ('import os,time; from pathlib import Path; '
                  f'Path({str(ready)!r}).write_text(str(os.getpid())); time.sleep(60)')
    code = f'''
import os, sys
from pathlib import Path
from funlr import cli, demo, public_demo
from funlr.runner import invoke
root = Path({str(tmp_path)!r})
def wait_demo(outdir):
    invoke([sys.executable, '-c', {child_code!r}], root/'tool.log', os.environ.copy(), root, root/'commands.jsonl')
demo.run_demo = public_demo.run_public_demo = wait_demo
raise SystemExit(cli.main(['demo', '--dataset', {dataset!r}, '--outdir', str(root/'output')]))
'''
    process = subprocess.Popen([sys.executable, '-c', code], cwd=tmp_path,
                               env=dict(os.environ, PYTHONPATH=str(ROOT / 'src')),
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    child_pid = None
    try:
        deadline = time.monotonic() + 10
        while not ready.exists() and process.poll() is None and time.monotonic() < deadline:
            time.sleep(0.02)
        assert ready.exists(), 'Demo did not launch its child'
        child_pid = int(ready.read_text())
        process.terminate()
        _, stderr = process.communicate(timeout=10)
        assert process.returncode == 130, stderr
        with pytest.raises(ProcessLookupError):
            os.kill(child_pid, 0)
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=5)
        if child_pid is not None:
            try:
                os.kill(child_pid, signal.SIGKILL)
            except ProcessLookupError:
                pass


@pytest.mark.skipif(os.environ.get('FUNLR_TEST_REAL_TOOLS') != '1', reason='opt-in real bioinformatics executables')
def test_real_public_example_is_packaged_and_resumable(tmp_path):
    report = public_demo.run_public_demo(tmp_path / 'public demo with spaces')
    assert report['status'] == 'passed'
    assert report['candidate_count'] == 6
    assert report['strict_candidate_count'] == 3
    assert report['spliced_locus_copies'] == 2
    assert report['intron_bp'] == 49
    assert report['stage6_refine_count'] == 1
    assert report['stage6_comprehensive_refine_count'] == 3
    assert report['completed_stages'] == report['resume_reused_stages'] == 8
