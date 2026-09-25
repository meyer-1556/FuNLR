"""Discovery regressions for empty evidence, literal IDs and deterministic ties.

These exercise the active discovery helpers. Existing tests already cover an empty
ALL source and a single explicit-LOW batch; the cases here extend those checks
to gated-only discovery, cohort independence, and process-level hash seeds.
"""
import json
import os
from pathlib import Path
import subprocess
import sys
import textwrap

import funlr
import pandas as pd
import pytest

from funlr.stages._discovery import description_priority
from funlr.stages._settings import DEFAULTS
from funlr.stages._discovery_run import build_union, write_ids
from funlr.stages.tiering import TierPolicy, build_profile_table


def discovery_dirs(tmp_path):
    stage0, stage1 = tmp_path / "stage0", tmp_path / "stage1"
    stage0.mkdir()
    stage1.mkdir()
    return stage0, stage1


def test_missing_annotation_cells_do_not_drop_keyword_columns(tmp_path):
    """A missing annotation in one row must not remove evidence in others."""
    stage0, stage1 = discovery_dirs(tmp_path)
    (stage0 / "eggnog_light.tsv").write_text(
        "protein_id\tDescription\tPreferred_name\n"
        "core\tNACHT domain protein\t\n"
        "repeat\t\tWD40\n"
        "effector\tgasdermin-like protein\t\n"
        "background\tABC transporter\tabcA\n"
        "unannotated\t\t\n"
    )
    description_priority(stage0, stage1, DEFAULTS)
    assert (stage1 / "desc_priority_ids.txt").read_text().splitlines() == [
        "core", "repeat", "effector"
    ]


@pytest.mark.parametrize("mode, expected", [("ALL", []), ("ANY", ["candidate"])])
def test_gated_description_uses_empty_sources_even_when_union_channel_is_off(
    tmp_path, mode, expected
):
    stage0, stage1 = discovery_dirs(tmp_path)
    (stage0 / "proteins_clean.faa").write_text(">candidate\nMMMM\n")
    for name, ids in {
        "nbd_candidate_ids_relaxed.txt": [],
        "nbd_candidate_ids_strict.txt": [],
        "pfam_nbd_candidate_ids.txt": [],
        "pfam_priority_ids.txt": ["candidate"],
        "desc_priority_ids.txt": ["candidate"],
    }.items():
        write_ids(stage1 / name, ids)
    values = {
        **DEFAULTS,
        "USE_PFAM_ENRICHMENT_FOR_UNION": 0,
        "INCLUDE_DESC_KEYWORDS_IN_UNION": 1,
        "DESC_KEYWORD_REQUIRE_STRUCTURAL_HINT": 1,
        "DESC_KEYWORD_STRUCTURAL_HINT_SOURCES": "PFAM_ENRICHED,PFAM_NBD",
        "DESC_KEYWORD_STRUCTURAL_HINT_MODE": mode,
    }
    # The enrichment IDs may support a description without being a direct
    # union channel. ALL still needs support from the requested empty source.
    assert build_union(stage0, stage1, values) == expected
    assert (stage1 / "desc_priority_ids.structural_hint.txt").read_text().splitlines() == expected


def test_explicit_low_confidence_is_independent_of_other_candidates():
    arch = pd.DataFrame([
        dict(protein_id=pid, has_nbd=True, has_nbd_stage1=False,
             has_nbd_pfam=True, has_sensor=True, has_ssfr_sensor=True,
             has_non_ssfr_sensor=False, domains_raw="NACHT;WD40",
             domains_grouped="NACHT;WD40", nbd_confidence=confidence,
             sensor_is_cterm_of_nbd=True)
        for pid, confidence in [("low", "LOW"), ("high", "HIGH")]
    ])
    # Separate scaffolds prevent the extra candidate from changing the
    # intentional genomic-neighbor rescue context of the LOW candidate.
    master = pd.DataFrame([
        dict(protein_id=pid, scaffold=pid, start=40000, end=42700,
             strand="+", protein_length=900)
        for pid in ("low", "high")
    ])
    contigs = pd.DataFrame([
        dict(scaffold=pid, contig_len=200000) for pid in ("low", "high")
    ])
    fields = ["nbd_confidence", "nbd_ok", "tier", "flags", "rescue_priority"]
    alone = build_profile_table(master, arch.iloc[:1], contigs, TierPolicy()).set_index("protein_id")
    together = build_profile_table(master, arch, contigs, TierPolicy()).set_index("protein_id")
    assert alone.loc["low", "nbd_confidence"] == "LOW"
    assert not alone.loc["low", "nbd_ok"]
    pd.testing.assert_series_equal(alone.loc["low", fields], together.loc["low", fields])


def test_enrichment_top_n_candidate_ids_are_stable_across_hash_seeds(tmp_path):
    """An exact enrichment tie must not change the candidate ID universe."""
    code = textwrap.dedent("""
        import json
        from pathlib import Path
        import pandas as pd
        from funlr.stages._discovery import enrich_pfams
        from funlr.stages._settings import DEFAULTS

        stage0, stage1 = Path('stage0'), Path('stage1')
        stage0.mkdir()
        stage1.mkdir()
        master = pd.DataFrame({
            'protein_id': ['fg1', 'fg2', 'only_a', 'only_b', 'bg1', 'bg2'],
            'PFAMs': ['PF_A,PF_B', 'PF_A,PF_B', 'PF_A', 'PF_B', '', ''],
        })
        master.to_csv(stage0 / 'master_table.tsv', sep='\\t', index=False)
        master.iloc[:2].to_csv(stage1 / 'nbd_hits_with_master_strict.tsv', sep='\\t', index=False)
        values = {**DEFAULTS, 'PFAM_ENRICH_MIN_FG_COUNT': 2,
                  'PFAM_ENRICH_MIN_LOG2FC': 0, 'PFAM_ENRICH_TOP_N': 1}
        enrich_pfams(stage0, stage1, values)
        print(json.dumps({
            'top': (stage1 / 'pfam_enriched_top50.txt').read_text().splitlines(),
            'ids': (stage1 / 'pfam_priority_ids.txt').read_text().splitlines(),
        }))
    """)
    package_parent = str(Path(funlr.__file__).resolve().parent.parent)
    outputs = []
    for seed in ("1", "4"):
        run_dir = tmp_path / seed
        run_dir.mkdir()
        env = dict(os.environ, PYTHONHASHSEED=seed, PYTHONPATH=package_parent)
        result = subprocess.run(
            [sys.executable, "-c", code], cwd=run_dir, env=env,
            text=True, capture_output=True, check=True, timeout=30,
        )
        outputs.append(json.loads(result.stdout))
    assert outputs == [{"top": ["PF_A"], "ids": ["fg1", "fg2", "only_a"]}] * 2
