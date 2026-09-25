"""Portable reference rescue tracks and explicitly audited fusion hypotheses.

Fusion acceptance is an engineering guard, not a biological validation claim:
one alignment must span every query member and overlap every source locus on
the expected scaffold/strand. A local match to one fragment cannot replace it.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import json
from pathlib import Path
import re
import shlex

import pandas as pd

from funlr.core.errors import StageError
from funlr.parsers.fasta import fasta_iter, write_fasta
from funlr.parsers.gff import read_exonerate_features

FEATURE_COLUMNS = ["query", "scaffold", "start", "end", "strand", "feature", "source", "score", "phase", "attrs_raw", "mrna_id", "parent_id", "identity", "positive", "rank"]
LOCUS_COLUMNS = ["query", "scaffold", "start", "end", "strand", "n_features"]
FUSION_COLUMNS = ["fusion_id", "members", "scaffold", "strand", "fusion_start", "fusion_end", "max_gap", "rationale", "member_lengths", "member_starts", "member_ends"]
EVIDENCE_COLUMNS = ["fusion_id", "members", "accepted", "reason", "query_length", "aligned_query_start", "aligned_query_end", "query_span_fraction", "scaffold", "start", "end", "strand", "score", "acceptance_rule"]
ACCEPTANCE_RULE = "single_alignment_spans_all_members_and_overlaps_all_source_loci"


def table(path, columns=()):
    path = Path(path)
    if not path.is_file() or not path.stat().st_size:
        return pd.DataFrame(columns=columns)
    return pd.read_csv(path, sep="\t")


def records(path):
    if not Path(path).is_file():
        return {}
    with open(path) as handle:
        result = {}
        for key, seq in fasta_iter(handle):
            result.setdefault(key, seq)
        return result


def record_list(path):
    if not Path(path).is_file():
        return []
    with open(path) as handle:
        return list(fasta_iter(handle))


def write_records(path, values):
    with open(path, "w") as handle:
        write_fasta(iter(values), handle)


def write_ids(path, values):
    Path(path).write_text("".join(str(v) + "\n" for v in values))


def read_ids(path):
    return list(dict.fromkeys(x.strip() for x in Path(path).read_text().splitlines() if x.strip())) if Path(path).is_file() else []


def truth(value):
    return str(value).strip().lower() in {"true", "1", "yes"}


def flag_tokens(value):
    return set(re.split(r"[;,|]", str(value))) if pd.notna(value) else set()


def dirs(root):
    for name in ("queries", "tables", "gff", "summaries", "logs", "metadata", "qc"):
        (root / name).mkdir(parents=True, exist_ok=True)


def manifest(root, name, params):
    info = root / "metadata/run_info.json"
    info.write_text(json.dumps(params, indent=2) + "\n")
    out = root / name
    paths = sorted(p for p in root.rglob("*") if p.is_file() and p != out and "comprehensive" not in p.relative_to(root).parts)
    pd.DataFrame([(str(p.relative_to(root)), str(p)) for p in paths], columns=["label", "path"]).to_csv(out, sep="\t", index=False)
    return [*paths, out]


def miniprot_features(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        fields = line.split("\t")
        if line.startswith("#") or len(fields) != 9:
            continue
        chrom, source, feature, start, end, score, strand, phase, attrs = fields
        attr = dict(re.findall(r"([A-Za-z0-9_]+)=([^;]+)", attrs))
        query = attr.get("Target", "").split()
        if not query:
            continue
        try:
            start, end = int(start), int(end)
        except ValueError:
            continue
        stats = []
        for key in ("Identity", "Positive", "Rank"):
            fallback = re.search(rf"{key.lower()}:([0-9.]+)", attr.get("Note", ""))
            stats.append(attr.get(key, fallback.group(1) if fallback else ""))
        rows.append([query[0], chrom, start, end, strand, feature, source, score, phase, attrs,
                     attr.get("ID", "") if feature == "mRNA" else "",
                     attr.get("Parent", "") if feature == "CDS" else "", *stats])
    return pd.DataFrame(rows, columns=FEATURE_COLUMNS)


def loci(features, id_col="query"):
    columns = [id_col, *LOCUS_COLUMNS[1:]]
    if features.empty:
        return pd.DataFrame(columns=columns)
    # Deliberately preserve the defined first-scaffold / min-max aggregation.
    return features.groupby(id_col, sort=True, as_index=False).agg(scaffold=("scaffold", "first"), start=("start", "min"), end=("end", "max"), strand=("strand", "first"), n_features=("feature", "size"))[columns]


def refinement_reasons(query_ids, features, tiered):
    result = {}
    tier_map = tiered.set_index("protein_id").to_dict("index") if "protein_id" in tiered else {}
    for query in sorted(set(query_ids)):
        group = features[features["query"] == query]
        mrna = group[group["feature"] == "mRNA"]
        if mrna.empty:
            result[query] = "no_miniprot_hit" if features.empty else "no_mrna"
            continue
        # Preserve reference metrics: last mRNA and all CDS per query.
        row = mrna.iloc[-1]
        count = int((group["feature"] == "CDS").sum())
        span = int(row["end"] - row["start"] + 1)
        try:
            score = float(row["score"])
        except (ValueError, TypeError):
            score = 0
        reasons = []
        if count >= 12:
            reasons.append(f"MANY_CDS({count})")
        if span >= 15000:
            reasons.append(f"LONG_SPAN({span})")
        if score <= 20:
            reasons.append(f"LOW_SCORE({score})")
        scaffolds = group["scaffold"].nunique()
        if scaffolds > 1:
            reasons.append(f"MULTI_SCAFFOLD({scaffolds})")
        if reasons:
            meta = tier_map.get(query, {})
            flags = str(meta.get("flags", ""))
            if str(meta.get("nbd_confidence", "")).upper() == "LOW":
                reasons.append("LOW_CONF")
            for flag, explanation in (("ORIENTATION_INVALID", "ORIENTATION_INVALID"), ("FP_DOMAIN_OVERLAPS_NBD", "FP_OVERLAP")):
                if flag in flags:
                    reasons.append(explanation)
            if meta.get("tier") in {"TIER_2B_RESCUE_CANDIDATE", "TIER_2C_LOW_PRIORITY_FRAGMENT", "TIER_4G_SUSPECTED_PSEUDOGENE"}:
                reasons.append("TIER_FLAGGED")
            result[query] = ";".join(reasons)
    return result


def run_miniprot(ctx):
    if ctx.dry_run:
        return []
    cfg = ctx.config
    stage5 = ctx.paths.stage_dir(5)
    stage4 = ctx.paths.stage_dir(4)
    include_tier3 = bool(cfg.get("rescue", "include_tier3_in_rescue", False) or cfg.get("rescue", "include_tier3", False))
    tiered = table(ctx.paths.stage_dir(3) / "tables/tiered_candidates.tsv")
    proteins = records(ctx.paths.stage_dir(0) / "proteins_clean.faa")
    preferred = stage4 / "rescue_priority.faa"
    if preferred.is_file():
        priority_raw = record_list(preferred)
        priority = records(preferred)  # An explicitly empty export remains empty.
    else:
        files = sorted(stage4.glob("TIER_2*.faa")) + sorted(stage4.glob("TIER_1B*.faa"))
        if include_tier3:
            files += sorted(stage4.glob("TIER_3*.faa"))
        if not files:  # Previous-release stage exports remain readable.
            files = [stage4 / name for name in ("tier2A_high_priority_rescue.faa", "tier2B_rescue_candidate.faa", "tier1B_needs_review.faa")]
            if include_tier3:
                files.append(stage4 / "tier3_architectural_variant.faa")
        priority = {}
        priority_raw = []
        for path in files:
            priority_raw.extend(record_list(path))
            for key, seq in records(path).items():
                priority.setdefault(key, seq)
    if cfg.get("rescue", "filter_low_nbd", False):
        rescue = table(ctx.paths.stage_dir(3) / "tables/rescue_priority.tsv")
        if "protein_id" not in rescue:
            raise StageError("Low-NBD filtering requires the Stage 3 rescue priority table")
        mask = rescue["nbd_ok"].map(truth) if "nbd_ok" in rescue else rescue["nbd_confidence"].astype(str).str.upper().isin(["HIGH", "MEDIUM"])
        keep = set(rescue.loc[mask, "protein_id"])
        priority = {key: seq for key, seq in priority.items() if key in keep}
        priority_raw = [(key, seq) for key, seq in priority_raw if key in keep]
    comprehensive_ids = tiered.loc[tiered["tier"].str.match(r"^TIER_[123]"), "protein_id"].tolist() if "tier" in tiered else []
    missing = set(comprehensive_ids) - proteins.keys()
    if missing:
        raise StageError("Comprehensive proteins absent from Stage 0: " + ", ".join(sorted(missing)[:10]))
    comprehensive = {key: proteins[key] for key in comprehensive_ids}
    outputs = []
    extra = shlex.split(cfg.get("rescue", "miniprot_extra", "") or "")
    for name, root, queries in (("priority", stage5, priority), ("comprehensive", stage5 / "comprehensive", comprehensive)):
        dirs(root)
        raw = root / ("queries/rescue_queries.raw.faa" if name == "priority" else "queries/comprehensive_queries.faa")
        query_fasta = root / ("queries/rescue_queries.faa" if name == "priority" else "queries/comprehensive_queries.dedup.faa")
        raw_records = priority_raw if name == "priority" else list(queries.items())
        write_records(raw, raw_records)
        write_records(query_fasta, queries.items())
        write_ids(root / f"queries/miniprot_{name}_ids.txt", queries)
        gff, log = root / "gff/miniprot.gff3", root / "logs/miniprot.log"
        if queries:
            ctx.runner.run([ctx.runner.resolve("miniprot"), "-t", ctx.cpu_threads, "--gff", *extra, cfg.get("inputs", "genome"), query_fasta], stdout_path=gff, stderr_path=log)
        else:
            gff.write_text("##gff-version 3\n")
            log.write_text("No queries; miniprot was not run.\n")
        features = miniprot_features(gff)
        summary = loci(features)
        features.to_csv(root / "tables/miniprot_features.tsv", sep="\t", index=False)
        summary.to_csv(root / "tables/miniprot_summary.tsv", sep="\t", index=False)
        pd.DataFrame(Counter(features["feature"]).most_common(), columns=["feature", "count"]).to_csv(root / "tables/miniprot_hit_counts.tsv", sep="\t", index=False)
        reasons = refinement_reasons(queries, features, tiered)
        write_ids(root / "queries/exonerate_refine_ids.txt", reasons)
        pd.DataFrame(reasons.items(), columns=["query", "reason"]).to_csv(root / "qc/exonerate_refine_reasons.tsv", sep="\t", index=False)
        # An isolated CDS row is not a complete miniprot hit. Match the defined
        # requested-versus-observed audit; keep the summary and refine rules.
        observed = set(features.loc[features["feature"].eq("mRNA"), "query"])
        write_ids(root / "summaries/miniprot_missing_hits.txt", sorted(set(queries) - observed))
        write_ids(root / "summaries/miniprot_extra_hits.txt", sorted(observed - set(queries)))
        pd.DataFrame([("raw_concat", str(raw), len(raw_records)), ("dedup", str(query_fasta), len(queries))], columns=["label", "file", "seqs"]).to_csv(root / "summaries/query_counts.tsv", sep="\t", index=False)
        outputs += manifest(root, "stage5_manifest.tsv" if name == "priority" else "comprehensive_manifest.tsv", {"stage": 5, "track": name, "query_count": len(queries), "refine_count": len(reasons), "miniprot_extra": extra})
    return outputs


def fusion_candidates(tiered, proteins, settings, track):
    """reference connected-component heuristic, with strand-correct member order."""
    rank = {"LOW": 0, "MEDIUM": 1, "HIGH": 2}
    min_conf = rank[settings.get("min_nbd_conf", "MEDIUM")]
    excludes = settings.get("exclude_flags", "LRR,FP_DOMAIN_OVERLAPS_NBD")
    excludes = set(excludes.split(",")) if isinstance(excludes, str) else set(excludes)
    groups = defaultdict(list)
    features = {}
    for row in tiered.to_dict("records"):
        if flag_tokens(row.get("flags")) & excludes:
            continue
        if any(pd.isna(row.get(k)) for k in ("scaffold", "start", "end", "strand")):
            continue
        row = dict(row)
        row["start"], row["end"] = int(row["start"]), int(row["end"])
        row["confident_nbd"] = truth(row.get("has_nbd")) and rank.get(str(row.get("nbd_confidence")).upper(), 0) >= min_conf
        features[row["protein_id"]] = row
        groups[(row["scaffold"], row["strand"])].append(row)
    adjacency = defaultdict(set)
    for (_, strand), group in groups.items():
        group.sort(key=lambda x: x["start"], reverse=strand == "-")
        for i, a in enumerate(group):
            if not a["confident_nbd"] or truth(a.get("has_sensor")):
                continue
            for b in group[i + 1:]:
                if not truth(b.get("has_sensor")) or b["confident_nbd"]:
                    continue
                gap = b["start"] - a["end"] if strand != "-" else a["start"] - b["end"]
                flags = flag_tokens(a.get("flags")) | flag_tokens(b.get("flags"))
                allowed = settings.get("max_gap", 50000)
                if "HARD_END" in flags:
                    allowed = max(allowed, settings.get("gap_boost_hard_end", 100000))
                if "REPEAT_NEARBY" in flags:
                    allowed = max(allowed, settings.get("gap_boost_repeat", 50000))
                if gap < -settings.get("max_overlap", 200) or max(0, gap) > allowed:
                    continue
                adjacency[a["protein_id"]].add(b["protein_id"])
                adjacency[b["protein_id"]].add(a["protein_id"])
    visited, result, query_records = set(), [], {}
    for seed in sorted(adjacency):
        if seed in visited:
            continue
        pending, component = [seed], []
        while pending:
            key = pending.pop()
            if key in visited:
                continue
            visited.add(key)
            component.append(features[key])
            pending.extend(sorted(adjacency[key] - visited))
        if not 2 <= len(component) <= settings.get("max_members", 3):
            continue
        component.sort(key=lambda x: (x["start"], x["end"]), reverse=component[0]["strand"] == "-")
        if not component[0]["confident_nbd"] or truth(component[0].get("has_sensor")) or not truth(component[-1].get("has_sensor")):
            continue
        members = [row["protein_id"] for row in component]
        missing = set(members) - proteins.keys()
        if missing:
            raise StageError("Fusion member proteins missing: " + ", ".join(sorted(missing)))
        prefix = settings.get(track + "_prefix", track.upper())
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", prefix):
            raise StageError("Fusion namespace prefixes must be filename-safe identifiers")
        padding = int(settings.get("id_pad_width", 6))
        if not 1 <= padding <= 20:
            raise StageError("Fusion ID padding must be between 1 and 20")
        fid = f"{prefix}__FUSION_{len(result) + 1:0{padding}d}"
        gaps = [b["start"] - a["end"] if a["strand"] != "-" else a["start"] - b["end"] for a, b in zip(component, component[1:])]
        result.append(dict(zip(FUSION_COLUMNS, [fid, ";".join(members), component[0]["scaffold"], component[0]["strand"], min(x["start"] for x in component), max(x["end"] for x in component), max([0, *gaps]), "NBD+sensor_split", ";".join(str(len(proteins[m])) for m in members), ";".join(str(x["start"]) for x in component), ";".join(str(x["end"]) for x in component)])))
        query_records[fid] = "".join(proteins[m] for m in members)
    return pd.DataFrame(result, columns=FUSION_COLUMNS), query_records


def fusion_evidence(fusion, raw_path):
    """Require one local alignment to support every member on its source locus.

    No empirical biological accuracy threshold is implied. This guard rejects
    fragment-only matches; accepted hypotheses still explicitly require review.
    """
    lengths = list(map(int, str(fusion["member_lengths"]).split(";")))
    starts = list(map(int, str(fusion["member_starts"]).split(";")))
    ends = list(map(int, str(fusion["member_ends"]).split(";")))
    offsets = [0]
    for length in lengths:
        offsets.append(offsets[-1] + length)
    result = {key: "" for key in EVIDENCE_COLUMNS}
    result.update(fusion_id=fusion["fusion_id"], members=fusion["members"], accepted=False, reason="NO_ALIGNMENT", query_length=sum(lengths), acceptance_rule=ACCEPTANCE_RULE)
    best = None
    for line in Path(raw_path).read_text().splitlines():
        if not line.startswith("vulgar:"):
            continue
        fields = line.split()
        if len(fields) < 10 or fields[1] != fusion["fusion_id"]:
            continue
        try:
            qs, qe, ts, te, score = int(fields[2]), int(fields[3]), int(fields[6]), int(fields[7]), int(fields[9])
        except ValueError:
            continue
        qlo, qhi = sorted((qs, qe))
        tlo, thi = sorted((ts, te))
        operations = fields[10:]
        if not operations or len(operations) % 3:
            continue
        qpos, tpos = qs, ts
        matched_query, matched_target = [], []
        try:
            for index in range(0, len(operations), 3):
                op, qa, ta = operations[index:index + 3]
                qa, ta = int(qa), int(ta)
                qnext = qpos + qa if qe >= qs else qpos - qa
                tnext = tpos + ta if te >= ts else tpos - ta
                if op in {"M", "C", "S"} and qa > 0 and ta > 0:
                    matched_query.append(tuple(sorted((qpos, qnext))))
                    matched_target.append(tuple(sorted((tpos, tnext))))
                qpos, tpos = qnext, tnext
        except ValueError:
            continue
        if qpos != qe or tpos != te:
            continue
        query_spans = all(any(max(lo, a) < min(hi, b) for lo, hi in matched_query) for a, b in zip(offsets, offsets[1:]))
        locus_spans = all(any(max(lo, s - 1) < min(hi, e) for lo, hi in matched_target) for s, e in zip(starts, ends))
        same_locus = fields[5] == fusion["scaffold"] and fields[8] == fusion["strand"]
        accepted = query_spans and locus_spans and same_locus
        # Prefer accepted evidence, then score, without combining alignments.
        priority = (int(accepted), score)
        if best is not None and priority <= best:
            continue
        best = priority
        result.update(accepted=accepted, reason="MEMBER_SPANNING_ALIGNMENT_REVIEW_REQUIRED" if accepted else "FRAGMENT_ONLY_OR_WRONG_LOCUS", aligned_query_start=qlo, aligned_query_end=qhi, query_span_fraction=(qhi-qlo)/sum(lengths), scaffold=fields[5], start=tlo + 1, end=thi, strand=fields[8], score=score)
    return result


def run_exonerate(ctx):
    if ctx.dry_run:
        return []
    from funlr.stages.stage6_exonerate import merge_gffs
    cfg = ctx.config
    tiered = table(ctx.paths.stage_dir(3) / "tables/tiered_candidates.tsv")
    proteins = records(ctx.paths.stage_dir(0) / "proteins_clean.faa")
    settings = dict(cfg.data.get("fusion", {}))
    profile = cfg.get("context", "profile", "ILLUMINA")
    settings["max_gap"] = settings.get("max_gap_hifi", 20000) if profile == "HIFI" else settings.get("max_gap_illumina", 50000)
    outputs = []
    for name in ("priority", "comprehensive"):
        source = ctx.paths.stage_dir(5) / ("comprehensive" if name == "comprehensive" else "")
        root = ctx.paths.stage_dir(6) / ("comprehensive" if name == "comprehensive" else "")
        dirs(root)
        per_query = root / "gff/per_query"
        per_query.mkdir(exist_ok=True)
        requested = read_ids(source / "queries/exonerate_refine_ids.txt")
        cap = int(cfg.get("rescue", "exon_maxn", 200))
        selected = requested[:cap] if cap else requested
        source_fa = source / ("queries/comprehensive_queries.dedup.faa" if name == "comprehensive" else "queries/rescue_queries.faa")
        source_records = records(source_fa)
        missing = set(selected) - source_records.keys()
        if missing:
            raise StageError("Exonerate refine IDs missing from rescue FASTA: " + ", ".join(sorted(missing)))
        original = {key: source_records[key] for key in selected}
        write_ids(root / "queries/exonerate_refine_ids.work.txt", selected)
        write_records(root / "queries/exonerate_queries.faa", original.items())
        write_ids(root / "queries/exonerate_queries.ids", original)
        if settings.get("enabled", True) and not tiered.empty:
            fusion, fusion_queries = fusion_candidates(tiered, proteins, settings, name)
        else:
            fusion, fusion_queries = pd.DataFrame(columns=FUSION_COLUMNS), {}
        fusion.to_csv(root / "tables/fusion_manifest.tsv", sep="\t", index=False)
        write_records(root / "queries/fusion_queries.faa", fusion_queries.items())
        write_ids(root / "queries/fusion_refine_ids.txt", fusion_queries)
        all_queries = {**original, **fusion_queries}
        write_records(root / "queries/exonerate_all.faa", all_queries.items())
        write_ids(root / "queries/exonerate_all.ids", sorted(all_queries))
        write_ids(root / "queries/exonerate_run.ids", sorted(all_queries))
        pairs, evidence = [], []
        fusion_map = fusion.set_index("fusion_id", drop=False).to_dict("index")
        for key in sorted(all_queries):
            query_fasta = root / "queries" / f"{key}.faa"
            write_records(query_fasta, [(key, all_queries[key])])
            raw = per_query / f"{key}.exonerate.gff"
            log = root / "logs" / f"{key}.exonerate.log"
            ctx.runner.run([ctx.runner.resolve("exonerate"), "--model", "protein2genome", "--query", query_fasta, "--target", cfg.get("inputs", "genome"), "--showtargetgff", "yes", "--maxintron", cfg.get("rescue", "exon_max_intron", 5000), "--percent", cfg.get("rescue", "exon_min_pct", 70), "--score", cfg.get("rescue", "exon_min_score", 150)], stdout_path=raw, stderr_path=log)
            pairs.append((key, raw))
            if key in fusion_map:
                evidence.append(fusion_evidence(fusion_map[key], raw))
            query_fasta.unlink()
        merged = root / "gff/exonerate_merged.gff3"
        merge_gffs(pairs, merged)
        features = read_exonerate_features(merged)
        summary = loci(features, "protein_id")
        summary.to_csv(root / "tables/exonerate_summary.tsv", sep="\t", index=False)
        pd.DataFrame(Counter(features["feature"]).most_common(), columns=["feature", "count"]).to_csv(root / "summaries/exonerate_counts.tsv", sep="\t", index=False)
        pd.DataFrame(evidence, columns=EVIDENCE_COLUMNS).to_csv(root / "tables/fusion_evidence.tsv", sep="\t", index=False)
        observed = set(summary["protein_id"])
        absent = sorted(set(all_queries) - observed)
        write_ids(root / "summaries/exonerate_missing_ids.txt", absent)
        write_ids(root / "summaries/exonerate_extra_ids.txt", sorted(observed - set(all_queries)))
        pd.DataFrame([(key, "NO_EXONERATE_FEATURES") for key in absent], columns=["protein_id", "reason"]).to_csv(root / "qc/missing_ids.tsv", sep="\t", index=False)
        pd.DataFrame([(key, "EXONERATE_HIT") for key in sorted(observed)], columns=["protein_id", "status"]).to_csv(root / "qc/processed_ids.tsv", sep="\t", index=False)
        outputs += manifest(root, "stage6_manifest.tsv", {"stage": 6, "track": name, "fusion_settings": settings, "requested_original": len(requested), "selected_original": len(original), "fusion_queries": len(fusion_queries), "accepted_fusion_evidence": sum(row["accepted"] for row in evidence), "acceptance_rule": ACCEPTANCE_RULE})
    return outputs
