"""Installed-style command contracts for read-only evidence utilities."""
import json
from pathlib import Path
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]


def invoke(*args):
    return subprocess.run([sys.executable, '-m', 'funlr', *map(str, args)], capture_output=True, text=True)


def test_qc_config_paths_and_exclusive_output(tmp_path):
    # No HMM libraries or tools are supplied: primary-input QC is independent.
    (tmp_path / 'genome.fa').write_text('>c\nAcgtnN\n')
    (tmp_path / 'proteins.faa').write_text('>p\nMA*\n')
    (tmp_path / 'annotation.gff3').write_text('c\ttest\tmRNA\t1\t6\t.\t+\t.\tID=p;Parent=g\n')
    config = tmp_path / 'settings.yaml'
    config.write_text('inputs:\n  genome: genome.fa\n  proteins: proteins.faa\n  gff3: annotation.gff3\n')
    result = invoke('input-qc', '--config', config)
    assert result.returncode == 0, result.stderr
    report = json.loads(result.stdout)
    assert report['genome']['n_bases'] == 2
    assert report['proteins']['n_with_terminal_stop'] == 1
    assert report['gff3']['n_gene_features'] == 0
    output = tmp_path / 'qc.json'
    result = invoke('input-qc', '--config', config, '--output', output)
    assert result.returncode == 0, result.stderr
    assert not result.stdout
    saved = output.read_bytes()
    result = invoke('input-qc', '--config', config, '--output', output)
    assert result.returncode == 1
    assert output.read_bytes() == saved


def test_qc_requires_primary_files_without_traceback():
    result = invoke('input-qc')
    assert result.returncode == 1
    assert 'Missing input QC files' in result.stderr
    assert 'Traceback' not in result.stderr


def test_reference_comparison_command():
    fixture = REPO / 'examples/benchmark'
    result = invoke('benchmark', '--reference', fixture / 'reference.json', '--predictions', fixture / 'called.tsv', '--input', fixture / 'input.faa')
    assert result.returncode == 0, result.stderr
    report = json.loads(result.stdout)
    assert report['independence'] == 'SYNTHETIC'
    assert report['metrics']['positive_recovery'] == 0.5
    assert report['metrics']['precision'] is None
    assert report['unknown_called_ids'] == ['unresolved']


def test_model_inventory_does_not_require_external_tools():
    result = invoke('model-inventory', '--directory', REPO / 'examples/demo/db')
    assert result.returncode == 0, result.stderr
    report = json.loads(result.stdout)
    assert report['schema_version'] == 1
