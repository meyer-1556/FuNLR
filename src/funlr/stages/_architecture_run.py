"""Native Pfam, custom-domain and ASM scans for architecture annotation."""
import shutil
from pathlib import Path

import pandas as pd

from funlr.core.errors import StageError
from funlr.core.tools import require_tools
from funlr.parsers.fasta import fasta_iter, write_fasta
from ._settings import enabled, settings
from ._discovery_run import scan, read_ids, write_ids
from ._architecture import build_architecture


def run(ctx):
    values = settings(ctx.config)
    log = ctx.stage_logger(2, "architecture")
    require_tools(ctx.runner, ["hmmscan", "hmmpress"], log)
    stage0 = ctx.paths.results_dir / "stage0"
    stage1, stage2 = ctx.paths.results_dir / "stage1", ctx.paths.stage_dir(2)
    query = stage1 / "union_candidates.faa"
    if not query.is_file():
        raise StageError(f"Stage-1 union FASTA is missing: {query}")
    mode = str(values["NBD_MODE"]).lower()
    names = {"strict": "nbd_candidate_ids_strict.txt", "relaxed": "nbd_candidate_ids_relaxed.txt", "legacy": "nbd_candidate_ids.txt"}
    if mode not in names:
        raise StageError(f"NBD_MODE must be strict, relaxed, or legacy: {mode}")
    ids_path = stage1 / names[mode]
    # A valid empty ID list remains empty; it is not a failed upstream stage.
    if not ids_path.is_file():
        raise StageError(f"Stage-1 NBD ID list is missing: {ids_path}")
    write_ids(stage2 / "stage1_nbd.ids", sorted(set(read_ids(ids_path))))
    hits = stage1 / "nbd_hits.tsv"
    if not hits.is_file():
        raise StageError(f"Stage-1 NBD hit table is missing: {hits}")
    shutil.copyfile(hits, stage2 / "stage1_nbd_hits.tsv")
    # Retain the defined strict interval source even when relaxed membership is selected.
    value = values["MAX_DOM_I_EVAL"]
    values["EVAL_USE"] = values["EVAL_PFAM"] if float(value) == 0 else value
    pfam = ctx.config.get("databases", "pfam")
    scan(ctx, pfam, query, stage2, "pfam", values["EVAL_USE"])
    if str(values["PFAM_ENRICH_BG"]).lower() == "proteome":
        if not scan(ctx, pfam, stage0 / "proteins_clean.faa", stage2, "pfam_bg", values["EVAL_USE"], check=False):
            values["PFAM_ENRICH_BG"] = "candidates"
            (stage2 / "pfam_bg.domtblout").write_text("")
            log.warning("Background Pfam scan failed; retaining reference candidates-background fallback")
    else:
        (stage2 / "pfam_bg.domtblout").write_text("")
        (stage2 / "pfam_bg.hmmscan.txt").write_text("")
    custom = ctx.config.get("databases", "custom_hmms")
    if custom:
        scan(ctx, custom, query, stage2, "custom", values["EVAL_USE"])
        (stage2 / "custom_scan.status").write_text("RUN\n" if query.stat().st_size else "EMPTY_INPUT\n")
    else:
        (stage2 / "custom.domtblout").write_text("")
        (stage2 / "custom.hmmscan.txt").write_text("")
        (stage2 / "custom_scan.status").write_text("DISABLED (no custom database)\n")
    window_faa = stage2 / "asm_window_nterm.faa"
    if enabled(values["ASM_ENABLE"]):
        asm = ctx.config.get("databases", "asm_hmms")
        if not asm:
            raise StageError("ASM_ENABLE requires databases.asm_hmms; set ASM_ENABLE=0 to disable this evidence channel")
        with query.open() as handle, window_faa.open("w") as out:
            write_fasta(((pid, seq[:int(values["ASM_WINDOW_NTERM_AA"])]) for pid, seq in fasta_iter(handle)), out, wrap=60)
        scan(ctx, asm, window_faa, stage2, "asm", values["ASM_MAX_EVALUE"])
        (stage2 / "asm_scan.status").write_text("RUN\n" if query.stat().st_size else "EMPTY_INPUT\n")
    else:
        for name in ("asm_window_nterm.faa", "asm.domtblout", "asm.hmmscan.txt"):
            (stage2 / name).write_text("")
        (stage2 / "asm_scan.status").write_text("DISABLED (config disabled)\n")
    summary = build_architecture(stage0, stage1, stage2, values)
    (stage2 / "run_info.txt").write_text("stage: 2\n" + "".join(f"{key.lower()}: {value}\n" for key, value in values.items()) +
                                         "asm_windows_implemented: N_TERM\nrepeatlike_groups_implemented: SSFR\n" +
                                         f"candidate_count: {len(summary)}\n")
    log.info("Architecture complete: %s candidates, %s NBD calls", len(summary), int(summary["has_nbd"].sum()))
    return sorted(path for path in stage2.iterdir() if path.is_file())
