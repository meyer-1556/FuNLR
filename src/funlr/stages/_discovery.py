"""Stage 1 NBD parsing, enrichment and description-evidence transforms.

Source: supplied 01_discover_candidates.sh. External programs are orchestrated
by stage1_discover.run; this module only reads/writes stage data. Empty output
schemas and deterministic enrichment ties are deliberate portability fixes.
"""
from pathlib import Path

def parse_calls(stage0, stage1, settings):
    import pandas as pd, sys, os

    DOMTBL = str(stage1 / 'nbd_whole.domtblout')
    EVAL_STRICT = float(settings["EVAL_NBD"])
    EVAL_RELAXED = float(settings["EVAL_NBD_RELAXED"])
    MIN_ALI_STRICT = int(settings["MIN_NBD_ALI_LEN"])
    MIN_ALI_RELAXED = int(settings["MIN_NBD_ALI_LEN_RELAXED"])

    print(f"STRICT thresholds: eval <= {EVAL_STRICT}, ali_len >= {MIN_ALI_STRICT}", file=sys.stderr)
    print(f"RELAXED thresholds: eval <= {EVAL_RELAXED}, ali_len >= {MIN_ALI_RELAXED}", file=sys.stderr)

    rows=[]
    total=0
    strict_kept=0
    relaxed_kept=0
    filtered_eval_strict=0
    filtered_len_strict=0
    filtered_eval_relaxed=0
    filtered_len_relaxed=0
    badline=0

    with open(DOMTBL) as f:
        for line in f:
            if line.startswith("#"):
                continue
            total += 1
            x=line.split()
            # domtblout has fixed first 22+ columns
            if len(x) < 22:
                badline += 1
                continue

            # hmmscan domtblout:
            #   query name (protein) = x[0]
            #   target name (HMM)    = x[3]
            hmm  = x[0]   # target name (HMM)
            prot = x[3]   # query name (protein)

            try:
                i_eval = float(x[12])  # domain i-Evalue
                ali_start = int(x[17])
                ali_end   = int(x[18])
            except Exception:
                badline += 1
                continue

            alilen = abs(ali_end - ali_start) + 1
        
            # Determine which filters pass
            passes_eval_strict = i_eval <= EVAL_STRICT
            passes_len_strict = alilen >= MIN_ALI_STRICT
            passes_eval_relaxed = i_eval <= EVAL_RELAXED
            passes_len_relaxed = alilen >= MIN_ALI_RELAXED
        
            passes_strict = passes_eval_strict and passes_len_strict
            passes_relaxed = passes_eval_relaxed and passes_len_relaxed
        
            # Count filtering reasons
            if not passes_eval_strict:
                filtered_eval_strict += 1
            if not passes_len_strict:
                filtered_len_strict += 1
            if not passes_eval_relaxed:
                filtered_eval_relaxed += 1
            if not passes_len_relaxed:
                filtered_len_relaxed += 1
            
            if passes_strict:
                strict_kept += 1
            if passes_relaxed:
                relaxed_kept += 1
        
            # Always keep all hits for analysis, with flags
            rows.append({
                "protein_id": prot,
                "nbd_hmm": hmm,
                "dom_i_evalue": i_eval,
                "ali_start": ali_start,
                "ali_end": ali_end,
                "ali_len": alilen,
                "passes_eval_strict": passes_eval_strict,
                "passes_len_strict": passes_len_strict,
                "passes_strict": passes_strict,
                "passes_eval_relaxed": passes_eval_relaxed,
                "passes_len_relaxed": passes_len_relaxed,
                "passes_relaxed": passes_relaxed
            })

    df=pd.DataFrame(rows, columns=[
        "protein_id", "nbd_hmm", "dom_i_evalue", "ali_start", "ali_end", "ali_len",
        "passes_eval_strict", "passes_len_strict", "passes_strict", "passes_eval_relaxed",
        "passes_len_relaxed", "passes_relaxed"])
    for flag in [c for c in df.columns if c.startswith("passes_")]:
        df[flag] = df[flag].astype(bool)
    df.to_csv(str(stage1 / 'nbd_hits_detailed.tsv'), sep="\t", index=False)

    # Extract candidate IDs for strict and relaxed criteria
    strict_cand = df[df["passes_strict"]][["protein_id"]].drop_duplicates()
    relaxed_cand = df[df["passes_relaxed"]][["protein_id"]].drop_duplicates()

    strict_cand.to_csv(str(stage1 / 'nbd_candidate_ids_strict.txt'), index=False, header=False)
    relaxed_cand.to_csv(str(stage1 / 'nbd_candidate_ids_relaxed.txt'), index=False, header=False)

    # Also write the filtered (strict) hits for backward compatibility
    df_strict = df[df["passes_strict"]].drop(columns=[
        "passes_eval_strict", "passes_len_strict", "passes_strict",
        "passes_eval_relaxed", "passes_len_relaxed", "passes_relaxed"
    ])
    df_strict.to_csv(str(stage1 / 'nbd_hits.tsv'), sep="\t", index=False)

    # Parse Pfam NBD hits if available
    pfam_domtbl = str(stage1 / 'pfam_nbd.domtblout')
    pfam_rows = []
    if os.path.exists(pfam_domtbl) and os.path.getsize(pfam_domtbl) > 0:
        print("Parsing Pfam NBD hits...", file=sys.stderr)
        with open(pfam_domtbl) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                x=line.split()
                if len(x) < 22:
                    continue
                hmm = x[0]
                prot = x[3]
                # Only keep if hmm is NACHT or NB-ARC
                if hmm in ["NACHT", "PF05729", "NB-ARC", "PF00931"]:
                    try:
                        i_eval = float(x[12])
                        ali_start = int(x[17])
                        ali_end   = int(x[18])
                    except Exception:
                        continue
                    alilen = abs(ali_end - ali_start) + 1
                    if i_eval <= EVAL_RELAXED and alilen >= MIN_ALI_RELAXED:
                        pfam_rows.append({
                            "protein_id": prot,
                            "pfam_nbd_hmm": hmm,
                            "dom_i_evalue": i_eval,
                            "ali_start": ali_start,
                            "ali_end": ali_end,
                            "ali_len": alilen
                        })
    
        if pfam_rows:
            pfam_df = pd.DataFrame(pfam_rows)
            pfam_df.to_csv(str(stage1 / 'pfam_nbd_hits.tsv'), sep="\t", index=False)
            pfam_cand = pfam_df[["protein_id"]].drop_duplicates()
            pfam_cand.to_csv(str(stage1 / 'pfam_nbd_candidate_ids.txt'), index=False, header=False)
            print(f"Pfam NBD candidates: {pfam_cand.shape[0]}", file=sys.stderr)
        else:
            print("No Pfam NBD hits found.", file=sys.stderr)

    if not pfam_rows:
        pd.DataFrame(columns=["protein_id", "pfam_nbd_hmm", "dom_i_evalue", "ali_start", "ali_end", "ali_len"]).to_csv(stage1 / "pfam_nbd_hits.tsv", sep="\t", index=False)
        (stage1 / "pfam_nbd_candidate_ids.txt").write_text("")

    with open(str(stage1 / 'nbd_parse_report.txt'),"w") as out:
        out.write(f"total_domtbl_lines\t{total}\n")
        out.write(f"bad_short_lines\t{badline}\n")
        out.write(f"strict_kept\t{strict_kept}\n")
        out.write(f"relaxed_kept\t{relaxed_kept}\n")
        out.write(f"filtered_high_i_eval_strict\t{filtered_eval_strict}\n")
        out.write(f"filtered_short_align_strict\t{filtered_len_strict}\n")
        out.write(f"filtered_high_i_eval_relaxed\t{filtered_eval_relaxed}\n")
        out.write(f"filtered_short_align_relaxed\t{filtered_len_relaxed}\n")
        out.write(f"unique_strict_candidate_proteins\t{strict_cand.shape[0]}\n")
        out.write(f"unique_relaxed_candidate_proteins\t{relaxed_cand.shape[0]}\n")
    
        # Add distribution information
        if not df.empty:
            eval_stats = df["dom_i_evalue"].describe()
            len_stats = df["ali_len"].describe()
            out.write(f"eval_mean\t{eval_stats['mean']:.2e}\n")
            out.write(f"eval_median\t{eval_stats['50%']:.2e}\n")
            out.write(f"len_mean\t{len_stats['mean']:.1f}\n")
            out.write(f"len_median\t{len_stats['50%']:.1f}\n")
        
            # Borderline bucket: eval <= 1e-3 and len >= 80
            borderline = df[(df["dom_i_evalue"] <= 1e-3) & (df["ali_len"] >= 80)]
            out.write(f"borderline_hits_1e3_80bp\t{len(borderline)}\n")

    print(f"STRICT NBD candidate proteins: {strict_cand.shape[0]}", file=sys.stderr)
    print(f"RELAXED NBD candidate proteins: {relaxed_cand.shape[0]}", file=sys.stderr)
    print(f"Wrote: {stage1}/nbd_parse_report.txt", file=sys.stderr)


def merge_calls(stage0, stage1, settings):
    import pandas as pd, os, sys

    hits = pd.read_csv(str(stage1 / 'nbd_hits_detailed.tsv'), sep="\t")
    for flag in ("passes_strict", "passes_relaxed"):
        hits[flag] = hits[flag].astype(bool)
    master_path = str(stage0 / 'master_table.tsv')
    out_path_strict = str(stage1 / 'nbd_hits_with_master_strict.tsv')
    out_path_relaxed = str(stage1 / 'nbd_hits_with_master_relaxed.tsv')

    if not os.path.exists(master_path) or os.path.getsize(master_path) == 0:
        # Write both strict and relaxed subsets
        hits[hits["passes_strict"]].to_csv(out_path_strict, sep="\t", index=False)
        hits[hits["passes_relaxed"]].to_csv(out_path_relaxed, sep="\t", index=False)
        print("master_table missing/empty; wrote hits only.", file=sys.stderr)
    else:
        master = pd.read_csv(master_path, sep="\t")
    
        # Strict hits
        strict_hits = hits[hits["passes_strict"]].drop_duplicates(subset=["protein_id"])
        merged_strict = strict_hits.merge(master, on="protein_id", how="left")
        merged_strict.to_csv(out_path_strict, sep="\t", index=False)
    
        # Relaxed hits
        relaxed_hits = hits[hits["passes_relaxed"]].drop_duplicates(subset=["protein_id"])
        merged_relaxed = relaxed_hits.merge(master, on="protein_id", how="left")
        merged_relaxed.to_csv(out_path_relaxed, sep="\t", index=False)
    
        print(f"Wrote strict: {out_path_strict} rows: {merged_strict.shape[0]}", file=sys.stderr)
        print(f"Wrote relaxed: {out_path_relaxed} rows: {merged_relaxed.shape[0]}", file=sys.stderr)


def enrich_pfams(stage0, stage1, settings):
    import math
    import os
    import re
    import sys
    from collections import Counter

    import pandas as pd

    STAGE0_DIR = str(stage0)
    STAGE1_DIR = str(stage1)
    MIN_FG_COUNT = int(settings.get("PFAM_ENRICH_MIN_FG_COUNT", "2"))
    MIN_LOG2FC = float(settings.get("PFAM_ENRICH_MIN_LOG2FC", "1.0"))
    TOP_N = int(settings.get("PFAM_ENRICH_TOP_N", "50"))

    MERGED = f"{STAGE1_DIR}/nbd_hits_with_master_strict.tsv"
    MASTER = f"{STAGE0_DIR}/master_table.tsv"

    OUT_ENRICH = f"{STAGE1_DIR}/pfam_enriched_in_nbd.tsv"
    OUT_TOP    = f"{STAGE1_DIR}/pfam_enriched_top50.txt"
    OUT_IDS    = f"{STAGE1_DIR}/pfam_priority_ids.txt"

    if not os.path.exists(MASTER) or os.path.getsize(MASTER) == 0:
        open(OUT_ENRICH, "w").close()
        open(OUT_TOP, "w").close()
        open(OUT_IDS, "w").close()
        raise SystemExit(f"MASTER missing/empty: {MASTER}")

    master = pd.read_csv(MASTER, sep="\t")
    if "protein_id" not in master.columns or "PFAMs" not in master.columns:
        raise SystemExit("Missing protein_id or PFAMs in master_table.tsv")

    if os.path.exists(MERGED) and os.path.getsize(MERGED) > 0:
        hits = pd.read_csv(MERGED, sep="\t").drop_duplicates(subset=["protein_id"])
    else:
        hits = pd.DataFrame(columns=["protein_id", "PFAMs"])

    def split_pfams(v):
        if not isinstance(v, str):
            return []
        v = v.strip()
        if v in {"", "-", "--"}:
            return []
        return [p.strip() for p in re.split(r"[,\|;]\s*", v) if p.strip() and p.strip() not in {"-", "--"}]

    bg_count = Counter()
    for v in master["PFAMs"].fillna("").astype(str):
        ps = split_pfams(v)
        if ps:
            bg_count.update(sorted(set(ps)))

    fg_count = Counter()
    for v in hits.get("PFAMs", pd.Series([], dtype=str)).fillna("").astype(str):
        ps = split_pfams(v)
        if ps:
            fg_count.update(sorted(set(ps)))

    N_bg = master.shape[0]
    N_fg = hits.shape[0]

    rows = []
    for pf, fg_n in fg_count.items():
        bg_n = bg_count.get(pf, 0)
        if fg_n < MIN_FG_COUNT:
            continue
        fg_rate = (fg_n + 1) / (N_fg + 2)
        bg_rate = (bg_n + 1) / (N_bg + 2)
        log2fc = math.log2(fg_rate / bg_rate)
        if log2fc < MIN_LOG2FC:
            continue
        rows.append((pf, fg_n, bg_n, fg_rate, bg_rate, log2fc))

    enr = pd.DataFrame(
        rows,
        columns=["PFAM", "NBD_hits_with_PFAM", "Proteome_with_PFAM", "NBD_rate", "Proteome_rate", "log2FC"]
    )
    if not enr.empty:
        enr = enr.sort_values(["log2FC", "NBD_hits_with_PFAM", "PFAM"], ascending=[False, False, True])

    enr.to_csv(OUT_ENRICH, sep="\t", index=False)

    top = enr.head(TOP_N)["PFAM"].tolist() if not enr.empty else []
    with open(OUT_TOP, "w") as out:
        out.write("\n".join(top) + ("\n" if top else ""))

    topset = set(top)

    def has_top(v):
        return any(p in topset for p in split_pfams(v))

    prio_ids = master[master["PFAMs"].fillna("").astype(str).apply(has_top)][["protein_id"]].dropna().drop_duplicates()
    prio_ids.to_csv(OUT_IDS, index=False, header=False)

    print(f"Foreground STRICT NBD proteins: {N_fg}", file=sys.stderr)
    print(f"PFAM enrichment thresholds: min_fg_count={MIN_FG_COUNT}, min_log2fc={MIN_LOG2FC}, top_n={TOP_N}", file=sys.stderr)
    print(f"Wrote: {OUT_ENRICH} (n={enr.shape[0]})", file=sys.stderr)
    print(f"Wrote: {OUT_TOP} (n={len(top)})", file=sys.stderr)
    print(f"Wrote: {OUT_IDS} (n={prio_ids.shape[0]})", file=sys.stderr)


def description_priority(stage0, stage1, settings):
    import os
    import re
    import sys

    import pandas as pd

    EGGNOG_LIGHT=str(stage0 / 'eggnog_light.tsv')
    OUT=str(stage1 / 'desc_priority_ids.txt')
    PATTERN = settings.get("DESC_PRIORITY_KEYWORDS_REGEX", "").strip()

    if (not os.path.exists(EGGNOG_LIGHT)) or os.path.getsize(EGGNOG_LIGHT) == 0:
        open(OUT, "w").close()
        print("eggNOG_LIGHT missing/empty; wrote empty:", OUT, file=sys.stderr)
        return

    eggnog = pd.read_csv(EGGNOG_LIGHT, sep="\t")
    if "protein_id" not in eggnog.columns:
        open(OUT, "w").close()
        print("eggNOG_LIGHT lacks protein_id column; wrote empty:", OUT, file=sys.stderr)
        return

    if not PATTERN:
        open(OUT, "w").close()
        print("DESC_PRIORITY_KEYWORDS_REGEX is empty; wrote empty:", OUT, file=sys.stderr)
        return

    cols = [col for col in eggnog.columns if col != "protein_id" and eggnog[col].dtype == object]
    pat = re.compile(PATTERN, re.IGNORECASE)

    def row_hit(r):
        for c in cols:
            v = r.get(c, "")
            if isinstance(v, str) and pat.search(v):
                return True
        return False

    if cols and not eggnog.empty:
        prio = eggnog[eggnog.apply(row_hit, axis=1)][["protein_id"]].dropna().drop_duplicates()
    else:
        prio = eggnog[["protein_id"]].head(0)

    prio.to_csv(OUT, index=False, header=False)
    print(f"Description priority proteins: {prio.shape[0]}", file=sys.stderr)
