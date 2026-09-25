"""Native executable orchestration for NBD and enrichment discovery."""
import json
import shutil
from pathlib import Path

import pandas as pd

from funlr.core.errors import StageError
from funlr.core.tools import ensure_pressed, require_tools
from ._settings import enabled, settings
from ._discovery import parse_calls, merge_calls, enrich_pfams, description_priority


def read_ids(path):
    path = Path(path)
    return [line.strip() for line in path.read_text().splitlines() if line.strip()] if path.is_file() else []


def write_ids(path, ids):
    Path(path).write_text("".join(str(pid) + "\n" for pid in ids))


def scan(ctx, database, queries, stage, stem, evalue, *, check=True):
    """Scan a private pressed database and retain both raw HMMER outputs."""
    domtbl, report = stage / f"{stem}.domtblout", stage / f"{stem}.hmmscan.txt"
    if not Path(queries).stat().st_size:
        domtbl.write_text("")
        report.write_text("")
        return True
    pressed = ensure_pressed(database, ctx.paths.pressed_dir, ctx.runner, ctx.logger)
    result = ctx.runner.run([
        ctx.runner.resolve("hmmscan"), "--cpu", str(ctx.cpu_threads),
        "--domtblout", str(domtbl), "--incE", str(evalue), "--incdomE", str(evalue),
        str(pressed), str(queries),
    ], cwd=stage, stdout_path=report, check=check)
    return result.returncode == 0


def build_union(stage0, stage1, values):
    """Apply reference evidence-union rules; ALL includes empty requested sources."""
    names = {"RELAXED_NBD": "nbd_candidate_ids_relaxed.txt", "STRICT_NBD": "nbd_candidate_ids_strict.txt",
             "PFAM_NBD": "pfam_nbd_candidate_ids.txt", "PFAM_ENRICHED": "pfam_priority_ids.txt"}
    evidence = {key: set(read_ids(stage1 / name)) for key, name in names.items()}
    candidate_ids = evidence["RELAXED_NBD"] | evidence["PFAM_NBD"]
    if enabled(values["USE_PFAM_ENRICHMENT_FOR_UNION"]):
        candidate_ids |= evidence["PFAM_ENRICHED"]
    description = set(read_ids(stage1 / "desc_priority_ids.txt"))
    backed_description = set()
    if enabled(values["INCLUDE_DESC_KEYWORDS_IN_UNION"]):
        if enabled(values["DESC_KEYWORD_REQUIRE_STRUCTURAL_HINT"]):
            selected = [key.strip().upper() for key in values["DESC_KEYWORD_STRUCTURAL_HINT_SOURCES"].split(",") if key.strip()]
            unknown = set(selected) - evidence.keys()
            if unknown:
                raise StageError(f"Unknown description structural-hint sources: {sorted(unknown)}")
            sets = [evidence[key] for key in selected]
            if sets:
                hints = set.intersection(*sets) if values["DESC_KEYWORD_STRUCTURAL_HINT_MODE"].upper() == "ALL" else set.union(*sets)
                backed_description = description & hints
            candidate_ids |= backed_description
        else:
            candidate_ids |= description
    write_ids(stage1 / "desc_priority_ids.structural_hint.txt", sorted(backed_description))
    with (stage0 / "proteins_clean.faa").open() as handle:
        proteome = [line[1:].split()[0] for line in handle if line.startswith(">")]
    proteome_set = set(proteome)
    raw = sorted(candidate_ids)
    missing = sorted(candidate_ids - proteome_set)
    filtered = sorted(candidate_ids & proteome_set) if enabled(values["UNION_INTERSECT_WITH_PROTEOME_IDS"]) else raw
    write_ids(stage1 / "proteome_ids.txt", sorted(proteome_set))
    write_ids(stage1 / "union_candidate_ids.raw.txt", raw)
    write_ids(stage1 / "union_candidate_ids.txt", filtered)
    write_ids(stage1 / "union_candidate_ids.missing_from_fasta.txt", missing)
    # The missing-ID audit remains truthful even if intersection is disabled.
    rows = [("proteome_unique_ids", len(proteome_set)), ("union_raw_ids", len(raw)),
            ("union_filtered_ids", len(filtered)), ("union_missing_from_fasta", len(missing))]
    pd.DataFrame(rows, columns=["metric", "value"]).to_csv(stage1 / "union_candidate_accounting.tsv", sep="\t", index=False)
    return filtered


def run(ctx):
    values = settings(ctx.config)
    log = ctx.stage_logger(1, "discover")
    require_tools(ctx.runner, ["hmmscan", "hmmpress", "seqkit"] + (["hmmfetch"] if enabled(values["PFAM_NBD_SCAN"]) else []), log)
    stage0, stage1 = ctx.paths.results_dir / "stage0", ctx.paths.stage_dir(1)
    protein = stage0 / "proteins_clean.faa"
    if not protein.is_file():
        raise StageError(f"Stage-0 protein FASTA missing: {protein}")
    scan(ctx, ctx.config.get("databases", "nbd_hmms"), protein, stage1, "nbd_whole", values["EVAL_NBD"])
    pfam_success = False
    if enabled(values["PFAM_NBD_SCAN"]):
        parts = []
        for model, filename in (("NACHT", "pfam_nacht.hmm"), ("NB-ARC", "pfam_nbarc.hmm")):
            dest = stage1 / filename
            dest.unlink(missing_ok=True)
            result = ctx.runner.run([ctx.runner.resolve("hmmfetch"), "-o", str(dest),
                                     str(ctx.config.get("databases", "pfam")), model], check=False)
            if result.returncode == 0 and dest.is_file() and dest.stat().st_size:
                parts.append(dest)
        if len(parts) == 2:
            combined = stage1 / "pfam_nbd.hmm"
            with combined.open("wb") as handle:
                for part in parts:
                    with part.open("rb") as source:
                        shutil.copyfileobj(source, handle)
            scan(ctx, combined, protein, stage1, "pfam_nbd", values["EVAL_NBD_RELAXED"])
            pfam_success = True
        else:
            log.warning("Could not extract both NACHT and NB-ARC from Pfam; retaining reference fallback to empty Pfam-NBD evidence")
    if not pfam_success:
        (stage1 / "pfam_nbd.domtblout").write_text("")
        (stage1 / "pfam_nbd.hmmscan.txt").write_text("")

    parse_calls(stage0, stage1, values)
    # Additive alias supports NBD_MODE=legacy without treating missing as empty.
    shutil.copyfile(stage1 / "nbd_candidate_ids_strict.txt", stage1 / "nbd_candidate_ids.txt")
    lengths = pd.read_csv(stage0 / "protein_lengths.tsv", sep="\t")
    known = set(lengths["protein_id"].astype(str))
    with (stage1 / "candidate_id_sanity.txt").open("w") as handle:
        for mode in ("strict", "relaxed"):
            ids = set(read_ids(stage1 / f"nbd_candidate_ids_{mode}.txt"))
            missing = sorted(ids - known)
            handle.write(f"{mode}_candidates_total\t{len(ids)}\n{mode}_missing_in_lengths\t{len(missing)}\n{mode}_missing_percent\t{100 * len(missing) / max(1, len(ids)):.2f}\n")
            if missing:
                handle.write(f"{mode}_missing_ids_preview\t{','.join(missing[:50])}\n")
    merge_calls(stage0, stage1, values)
    enrich_pfams(stage0, stage1, values)
    description_priority(stage0, stage1, values)
    ids = build_union(stage0, stage1, values)
    union_faa = stage1 / "union_candidates.faa"
    if ids:
        ctx.runner.run([ctx.runner.resolve("seqkit"), "grep", "-f", str(stage1 / "union_candidate_ids.txt"), str(protein)],
                       cwd=stage1, stdout_path=union_faa, check=True)
    else:
        union_faa.write_text("")
    with union_faa.open() as handle:
        extracted = [line[1:].split()[0] for line in handle if line.startswith(">")]
    if set(extracted) != set(ids) or len(extracted) != len(ids):
        raise StageError(f"Candidate FASTA does not match requested IDs ({len(ids)} requested, {len(extracted)} extracted)")
    with (stage1 / "union_candidate_accounting.tsv").open("a") as handle:
        handle.write(f"union_fasta_seqs\t{len(extracted)}\n")
    sources = {"nbd_strict": "nbd_candidate_ids_strict.txt", "nbd_relaxed": "nbd_candidate_ids_relaxed.txt",
               "pfam_nbd": "pfam_nbd_candidate_ids.txt", "pfam_enriched": "pfam_priority_ids.txt",
               "desc_keywords": "desc_priority_ids.txt", "union_raw": "union_candidate_ids.raw.txt",
               "union_filtered": "union_candidate_ids.txt", "union_missing_from_fasta": "union_candidate_ids.missing_from_fasta.txt"}
    summary = [(key, len(read_ids(stage1 / name)), str(stage1 / name)) for key, name in sources.items()]
    summary.append(("union_fasta_seqs", len(extracted), str(union_faa)))
    pd.DataFrame(summary, columns=["source", "count", "file"]).to_csv(stage1 / "priority_sources_summary.tsv", sep="\t", index=False)
    (stage1 / "run_info.txt").write_text("stage: 1\n" + "".join(f"{key.lower()}: {value}\n" for key, value in values.items()) +
                                         f"pfam_nbd_scan_completed: {pfam_success}\nunion_count: {len(ids)}\n")
    log.info("Strict NBD=%s, relaxed=%s, union=%s", len(read_ids(stage1 / "nbd_candidate_ids_strict.txt")), len(read_ids(stage1 / "nbd_candidate_ids_relaxed.txt")), len(ids))
    return sorted(path for path in stage1.iterdir() if path.is_file())
