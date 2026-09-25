"""Final reports with explicit, deduplicated fusion evidence."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil

import pandas as pd

from funlr.core.errors import StageError
from funlr.stages.architecture_v3 import build_tripartite_architecture
from funlr.stages.rescue_v3 import table, truth, dirs, manifest, EVIDENCE_COLUMNS

BASE_STRICT = {"TIER_1A_HIGH_CONFIDENCE", "TIER_1B_NEEDS_REVIEW", "TIER_2A_HIGH_PRIORITY_RESCUE", "TIER_2B_RESCUE_CANDIDATE", "TIER_3A_INTEGRATED_DECOY", "TIER_3B_ARCHITECTURAL_VARIANT"}
REPORT_COLUMNS = ["protein_id", "gene_id", "transcript_id", "scaffold", "start", "end", "strand", "protein_length", "tier", "rescue_priority", "flags", "nbd_confidence", "domains_raw", "domains_grouped", "NLR_architecture", "Description", "Preferred_name", "GOs", "KEGG_ko", "mp_scaffold", "mp_start", "mp_end", "mp_strand", "exo_scaffold", "exo_start", "exo_end", "exo_strand", "n_features", "comprehensive_mp_scaffold", "comprehensive_mp_start", "comprehensive_mp_end", "comprehensive_mp_strand", "comprehensive_exo_scaffold", "comprehensive_exo_start", "comprehensive_exo_end", "comprehensive_exo_strand", "comprehensive_exo_n_features", "is_fusion_model", "fusion_source", "fusion_raw_id", "fusion_members", "fused_into", "fusion_evidence_rule"]


def mapping_table(path, kind, track):
    frame = table(path)
    if kind == "exo" and (not Path(path).is_file() or not Path(path).stat().st_size):
        from funlr.stages.stage7_finalize import _parse_exonerate_merged_gff
        frame = _parse_exonerate_merged_gff(Path(path).parent.parent / "gff/exonerate_merged.gff3")
        frame = frame.rename(columns={"exo_" + key: key for key in ("scaffold", "start", "end", "strand")})
    prefix = "comprehensive_" if track == "comprehensive" else ""
    wanted = ["protein_id", *[f"{prefix}{kind}_{key}" for key in ("scaffold", "start", "end", "strand")]]
    if frame.empty:
        return pd.DataFrame(columns=wanted)
    if "protein_id" not in frame:
        for key in ("query", "qid", "query_id"):
            if key in frame:
                frame = frame.rename(columns={key: "protein_id"})
                break
    if "protein_id" not in frame:
        raise StageError(f"Mapping summary lacks query identifier: {path}")
    rename = {key: f"{prefix}{kind}_{key}" for key in ("scaffold", "start", "end", "strand")}
    if kind == "exo":
        rename["n_features"] = "comprehensive_exo_n_features" if prefix else "n_features"
        wanted.append(rename["n_features"])
    frame = frame.rename(columns=rename)
    frame = frame.dropna(subset=["protein_id"])
    return frame.reindex(columns=wanted).drop_duplicates("protein_id", keep="first")


def integrate_fusions(out, evidence):
    """Collapse the same members across tracks; never suppress unsupported members."""
    out = out.copy()
    for key in ("fusion_source", "fusion_raw_id", "fusion_members", "fused_into", "fusion_evidence_rule"):
        out[key] = ""
    out["is_fusion_model"] = False
    if evidence.empty:
        return out
    accepted = evidence[evidence["accepted"].map(truth)].copy()
    groups = {}
    for row in accepted.to_dict("records"):
        key = tuple(sorted(str(row["members"]).split(";")))
        groups.setdefault(key, []).append(row)
    appended = []
    for members, rows in sorted(groups.items()):
        # Content identity is independent of track counters and input row order.
        fid = "FUSION_" + hashlib.sha256("\0".join(members).encode()).hexdigest()[:12]
        chosen = max(rows, key=lambda row: float(row.get("score") or 0))
        sources = sorted({row["track"].upper() for row in rows})
        raw_ids = sorted({str(row["fusion_id"]) for row in rows})
        record = {key: "" for key in REPORT_COLUMNS}
        nbd_conf = out.loc[out["protein_id"].isin(members), "nbd_confidence"]
        conf_rank = {"HIGH": 2, "MEDIUM": 1, "LOW": 0}
        record.update(protein_id=fid, gene_id="", transcript_id="", scaffold=chosen["scaffold"], start=chosen["start"], end=chosen["end"], strand=chosen["strand"], protein_length=chosen["query_length"], tier="FUSION_RESCUE", rescue_priority=100, flags="FUSION_SPANNING_ALIGNMENT;REVIEW_REQUIRED", nbd_confidence=max(nbd_conf, key=lambda value: conf_rank.get(str(value), -1), default=""), is_fusion_model=True, fusion_source=";".join(sources), fusion_raw_id=";".join(raw_ids), fusion_members=chosen["members"], fusion_evidence_rule=chosen["acceptance_rule"])
        for source in sources:
            best = max((row for row in rows if row["track"].upper() == source), key=lambda row: float(row.get("score") or 0))
            prefix = "comprehensive_" if source == "COMPREHENSIVE" else ""
            for key in ("scaffold", "start", "end", "strand"):
                record[f"{prefix}exo_{key}"] = best[key]
        for index in out.index[out["protein_id"].isin(members)]:
            existing = out.at[index, "fused_into"]
            out.at[index, "fused_into"] = ";".join(sorted(filter(None, [existing, fid])))
        appended.append(record)
    if appended:
        out = pd.concat([out, pd.DataFrame(appended)], ignore_index=True)
    return out


def strict_sets(out):
    ordinary = out["tier"].isin(BASE_STRICT) | (out["tier"].eq("TIER_2D_PROBABLE_UNKNOWN_SENSOR_REPEATLIKE") & out["nbd_confidence"].fillna("").str.upper().eq("HIGH"))
    ordinary &= out["fused_into"].fillna("").eq("")
    fusion = out["is_fusion_model"].map(truth)
    def present(column):
        return out[column].notna() & out[column].astype(str).ne("")
    masks = {
        "combined": ordinary | fusion,
        "priority": (ordinary & (present("mp_scaffold") | present("exo_scaffold"))) | (fusion & out["fusion_source"].fillna("").str.contains("PRIORITY")),
        "comprehensive": (ordinary & present("comprehensive_mp_scaffold")) | (fusion & out["fusion_source"].fillna("").str.contains("COMPREHENSIVE")),
    }
    result = {}
    for key, mask in masks.items():
        frame = out.loc[mask].copy()
        frame["flags"] = frame["flags"].fillna("").astype(str).str.replace(r"\bNBD_ONLY\b", "NBD_ONLY_UNKNOWN_SENSOR", regex=True)
        result[key] = frame
    return result


def bed_file(frame, path, prefix="", suffix="", original=False):
    with open(path, "w") as handle:
        for row in frame.to_dict("records"):
            chrom, start, end = (row.get(prefix + key) for key in ("scaffold", "start", "end"))
            if pd.isna(chrom) or pd.isna(start) or pd.isna(end) or chrom == "" or start == "" or end == "":
                continue
            name = str(row["protein_id"]) + ("|" + str(row["tier"]) if original else suffix)
            handle.write(f"{chrom}\t{int(float(start))-1}\t{int(float(end))}\t{name}\t0\t{row.get(prefix+'strand', '.')}\n")


def domain_plots(ctx, tiered, stage7, final_dir):
    from importlib.resources import files
    hit_path = ctx.paths.stage_dir(2) / "all_domain_hits.tsv"
    hits = table(hit_path, ["protein_id", "domain", "start", "end"])
    eligible = tiered["tier"].str.match(r"^TIER_[123]")
    if "is_fusion_model" in tiered:
        eligible &= ~tiered["is_fusion_model"].astype(str).str.strip().str.lower().isin({"true", "1", "yes", "t"})
    ids = tiered.loc[eligible, "protein_id"]
    target = hits[hits["protein_id"].isin(ids)]
    for pid in target["protein_id"].unique():
        if (pid in {"", ".", ".."} or "/" in pid or "\\" in pid
                or any(ord(ch) < 32 or ord(ch) == 127 for ch in pid)):
            raise StageError(f"Protein ID cannot be used as a plot filename: {pid!r}")
    status_path = stage7 / "metadata/plot_status.json"
    # These directories contain only products of this stage. Start fresh so an
    # old PNG cannot mask a failed device or survive after its candidate is gone.
    work = stage7 / "plot_work"
    destinations = (stage7 / "domain_plots", final_dir / "domain_plots")
    for directory in (work, *destinations):
        if directory.exists():
            shutil.rmtree(directory)
    if target.empty:
        status_path.write_text(json.dumps({"status": "SKIPPED_NO_DOMAIN_HITS", "proteins": 0}) + "\n")
        return [status_path]
    work.mkdir(exist_ok=True)
    target.to_csv(work / "all_domain_hits.tsv", sep="\t", index=False)
    tiered.to_csv(work / "tiered_candidates.tsv", sep="\t", index=False)
    (work / "candidate_ids.txt").write_text("".join(key + "\n" for key in ids))
    backend = ctx.config.get("reporting", "plot_backend", "matplotlib")
    if backend == "matplotlib":
        # A stage worker owns this cache; do not write into the user's home or
        # rely on a pre-existing machine-wide Matplotlib configuration.
        cache = work / ".matplotlib"
        cache.mkdir()
        os.environ["MPLCONFIGDIR"] = str(cache.resolve())
        from funlr.plots import generate_domain_plots, plotting_provenance
        generate_domain_plots(work / "all_domain_hits.tsv",
                              work / "tiered_candidates.tsv", work, ctx.logger)
        (work / "plot_provenance.json").write_text(
            json.dumps(plotting_provenance(), indent=2, sort_keys=True) + "\n"
        )
    elif backend == "r":
        script = files("funlr.stages").joinpath("domain_plots.R")
        ctx.runner.run([ctx.runner.resolve("Rscript"), "--vanilla", str(script)], cwd=work, stdout_path=stage7 / "logs/domain_plots.stdout.log", stderr_path=stage7 / "logs/domain_plots.stderr.log")
    else:
        raise StageError(f"Unknown domain plotting backend: {backend}")
    expected = [work / name for name in ("nlr_domains_by_tier.png", "sensor_distribution.png", "nbd_distribution.png", "asm_presence.png", "length_vs_priority.png", "order_sanity.png")]
    expected += [work / "individual" / (str(pid) + ".png") for pid in target["protein_id"].unique()]
    for path in expected:
        if not path.is_file():
            raise StageError(f"Domain plotting did not produce expected PNG: {path}")
        with path.open("rb") as handle:
            if handle.read(8) != b"\x89PNG\r\n\x1a\n":
                raise StageError(f"Domain plotting produced an invalid PNG: {path}")
    outputs = []
    for destination in destinations:
        destination.mkdir(exist_ok=True)
        for path in work.rglob("*"):
            if path.is_file() and (path.suffix == ".png" or path.name in {"domain_architecture_summary.tsv", "R_sessionInfo.txt", "plot_provenance.json"}):
                dest = destination / path.relative_to(work)
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(path, dest)
                outputs.append(dest)
    status_path.write_text(json.dumps({"status": "WRITTEN", "backend": backend, "proteins": target["protein_id"].nunique(), "protein_length_source": "Stage3 protein_length", "fusion_models_plotted": False}) + "\n")
    return [*outputs, status_path]


def run_final(ctx):
    if ctx.dry_run:
        return []
    stage7, final = ctx.paths.stage_dir(7), ctx.paths.final_dir
    dirs(stage7)
    (stage7 / "beds").mkdir(exist_ok=True)
    tiered_path = ctx.paths.stage_dir(3) / "tables/tiered_candidates.tsv"
    tiered = table(tiered_path)
    if "protein_id" not in tiered:
        raise StageError(f"Tiered candidate table missing: {tiered_path}")
    out = tiered.copy()
    for track in ("priority", "comprehensive"):
        tail = "comprehensive" if track == "comprehensive" else ""
        for kind, stage, filename in (("mp", 5, "miniprot_summary.tsv"), ("exo", 6, "exonerate_summary.tsv")):
            out = out.merge(mapping_table(ctx.paths.stage_dir(stage) / tail / "tables" / filename, kind, track), on="protein_id", how="left")
    for key in REPORT_COLUMNS:
        if key not in out:
            out[key] = pd.NA
    architecture = build_tripartite_architecture(ctx.paths.stage_dir(2) / "all_domain_hits.tsv")
    out["NLR_architecture"] = out["protein_id"].map(architecture).fillna("")
    evidence = []
    for track in ("priority", "comprehensive"):
        path = ctx.paths.stage_dir(6) / ("comprehensive" if track == "comprehensive" else "") / "tables/fusion_evidence.tsv"
        frame = table(path, EVIDENCE_COLUMNS)
        frame["track"] = track
        evidence.append(frame)
    evidence = pd.concat(evidence, ignore_index=True)
    evidence.to_csv(stage7 / "tables/fusion_evidence.tsv", sep="\t", index=False)
    out = integrate_fusions(out, evidence).reindex(columns=REPORT_COLUMNS)
    strict = strict_sets(out)
    reports = {"nlr_final_report.tsv": out, "nlr_strict_candidates.tsv": strict["combined"], "nlr_strict_candidates.priority.tsv": strict["priority"], "nlr_strict_candidates.comprehensive.tsv": strict["comprehensive"]}
    for name, frame in reports.items():
        for directory in (stage7 / "tables", final):
            frame.to_csv(directory / name, sep="\t", index=False)
    tiers = out["tier"].value_counts().rename_axis("tier").reset_index(name="count")
    qc = [("n_total_candidates", len(tiered)), *[("tier_count::" + key, int(value)) for key, value in tiered["tier"].value_counts().items()], *[("n_strict_candidates_" + key, len(value)) for key, value in strict.items()], ("n_priority_miniprot_hits", int(out["mp_scaffold"].fillna("").ne("").sum())), ("n_comprehensive_miniprot_hits", int(out["comprehensive_mp_scaffold"].fillna("").ne("").sum())), ("n_fusion_models_total", int(out["is_fusion_model"].map(truth).sum())), ("n_fragments_fused", int(out["fused_into"].fillna("").ne("").sum())), ("n_fusion_alignments_accepted", int(evidence["accepted"].map(truth).sum()))]
    for directory in (stage7 / "summaries", final):
        tiers.to_csv(directory / "tier_summary.tsv", sep="\t", index=False)
        pd.DataFrame(qc, columns=["metric", "value"]).to_csv(directory / "qc_summary.tsv", sep="\t", index=False)
    for name, prefix, suffix in (("nlr_candidates.bed", "", ""), ("nlr_miniprot_loci.bed", "mp_", "|MINIPROT"), ("nlr_exonerate_loci.bed", "exo_", "|EXONERATE"), ("nlr_comprehensive_miniprot.bed", "comprehensive_mp_", "|COMP_MINIPROT"), ("nlr_comprehensive_exonerate.bed", "comprehensive_exo_", "|COMP_EXONERATE"), ("nlr_fusions_exonerate.bed", "", "|FUSION")):
        frame = out[out["is_fusion_model"].map(truth)] if "fusions" in name else out
        for directory in (stage7 / "beds", final):
            bed_file(frame, directory / name, prefix, suffix, original=name == "nlr_candidates.bed")
    extra = []
    if ctx.config.get("reporting", "plots", True):
        extra += domain_plots(ctx, tiered, stage7, final)
    if ctx.config.get("reporting", "integration", True):
        from funlr.stages.annotation_export import export_annotations
        extra += export_annotations(ctx, out, final / "integration")
    params = {"stage": 7, "fusion_policy": "member-spanning alignments are review hypotheses, deduplicated by members", "strict_counts": {key: len(value) for key, value in strict.items()}, "reporting": ctx.config.data.get("reporting", {})}
    outputs = manifest(stage7, "stage7_manifest.tsv", params)
    return list(dict.fromkeys([*outputs, *extra, *[path for path in final.rglob("*") if path.is_file()]]))
