"""Domain architecture and orientation rules as a package function.

Source: supplied 02_domain_architecture.sh, stage 2.4 Python. Domain vocabularies,
interval rules, thresholds, source priority, and confidence behavior are retained.
Fixes: per-protein effector state, empty schemas, deterministic enrichment ties.
The source's SSFR-only repeat-like flag remains unchanged.
"""
from pathlib import Path
from ._settings import enabled

def build_architecture(stage0, stage1, stage2, settings):
    import os, re, sys
    import pandas as pd
    FASTA = str(stage1 / "union_candidates.faa")                       # union_candidates.faa
    # Helper to strip taxonomic suffixes (e.g., "__class__Agaricomycetes")
    def domain_base(d):
        if "__" in d:
            return d.split("__")[0]
        return d

    # Effector domain groups (names from combined effector HMM)
    EFFECTOR_GROUPS = {
        "DUF676", "Goodbye", "HeLo", "HeLo-like", "Helo_like_N", "HET", "het-s",
        "HET-s_218-289", "NACHT_sigma", "NAD1", "Patatin", "PNP_UDP_1", "PP",
        "RelA_SpoT", "sesA", "SesA", "sesB", "sigma"
    }

    # Effector placement heuristics
    # Allow a small boundary buffer because domtblout coordinates can vary slightly
    EFFECTOR_NTERM_BUFFER_AA = 10

    # Optional: if you want to allow effectors that slightly overlap the NBD start
    # set to >0 (I recommend 0 or 5). Keep conservative by default.
    EFFECTOR_ALLOW_OVERLAP_AA = 0

    EVAL  = float(settings["EVAL_USE"])         # same threshold used for hmmscan + parsing
    MIN_ALI = int(settings["MIN_DOM_ALI_LEN"])  # same minimum alignment length as parser

    PROTEOME_FASTA = str(stage0 / "proteins_clean.faa")
    PFAM_BG_DOMTBL = str(stage2 / 'pfam_bg.domtblout')

    STAGE1_NBD_IDS = str(stage2 / 'stage1_nbd.ids')
    STAGE1_NBD_HITS = str(stage2 / 'stage1_nbd_hits.tsv')

    OUT_ALL = str(stage2 / 'all_domain_hits.tsv')
    OUT_SUM = str(stage2 / 'architecture_summary.tsv')
    OUT_REP = str(stage2 / 'parse_filter_report.tsv')
    OUT_ASM_SUM = str(stage2 / 'asm_hits_summary.tsv')

    # --------------------------------------------------------------------------
    # NBD confidence thresholds (from nlr_config.sh, passed via environment)
    # --------------------------------------------------------------------------
    NBD_CONF_EV_HIGH  = float(settings.get("NBD_CONF_EV_HIGH", "1e-6"))
    NBD_CONF_COV_HIGH = float(settings.get("NBD_CONF_COV_HIGH", "0.45"))
    NBD_CONF_EV_MED   = float(settings.get("NBD_CONF_EV_MED", "1e-3"))
    NBD_CONF_COV_MED  = float(settings.get("NBD_CONF_COV_MED", "0.25"))

    def fasta_ids(path):
        ids=[]
        with open(path) as f:
            for line in f:
                if line.startswith(">"):
                    ids.append(line[1:].split()[0])  # Only first token
        return ids

    all_prots = fasta_ids(FASTA)
    all_set = set(all_prots)

    # Stage1 NBD IDs
    stage1_nbd = set()
    if os.path.exists(STAGE1_NBD_IDS) and os.path.getsize(STAGE1_NBD_IDS) > 0:
        with open(STAGE1_NBD_IDS) as f:
            stage1_nbd = {l.strip() for l in f if l.strip()}

    # Stage1 NBD align intervals (optional, best union interval per protein)
    stage1_nbd_iv = {}
    if os.path.exists(STAGE1_NBD_HITS) and os.path.getsize(STAGE1_NBD_HITS) > 0:
        try:
            h = pd.read_csv(STAGE1_NBD_HITS, sep="\t")
            # Expect: protein_id, ... ali_start, ali_end
            for pid, g in h.groupby("protein_id"):
                if pid not in all_set: 
                    continue
                a = int(g["ali_start"].min())
                b = int(g["ali_end"].max())
                if a > b: a, b = b, a
                stage1_nbd_iv[pid] = (a, b)
        except Exception:
            stage1_nbd_iv = {}

    def parse_domtbl(path, source, eval_thresh, min_ali, min_bitscore=None):
        rows=[]
        total=0; kept=0; bade=0; short=0; badline=0; lowbits=0
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            return pd.DataFrame(columns=["protein_id","domain","acc","hmm_len","dom_i_evalue","bitscore","start","end","ali_len","source"]), {
                "source": source, "total_domtbl_lines": 0, "kept_after_filters": 0,
                "filtered_high_i_eval": 0, "filtered_short_align": 0, "bad_short_lines": 0,
                "filtered_low_bitscore": 0
            }
        with open(path) as f:
            for line in f:
                if line.startswith("#"): 
                    continue
                total += 1
                x=line.split()
                if len(x) < 22:
                    badline += 1
                    continue
                dom_name = x[0]
                dom_acc  = x[1]
                hmm_len  = int(x[2])          # target length (HMM)
                prot = x[3].split()[0]        # normalize query id to first token (matches fasta_ids)
            
                # For ASM, we need to check if the protein is in our set
                # (it will be, because we created the window FASTA from all_prots)
                if source == "asm" and prot not in all_set:
                    # Try to extract just the protein ID (before any space)
                    prot_base = prot.split()[0] if " " in prot else prot
                    if prot_base in all_set:
                        prot = prot_base
                    else:
                        continue
            
                try:
                    i_eval = float(x[12])  # domain i-Evalue
                    bitscore = float(x[13])  # domain score (bits)
                    ali_start = int(x[17])
                    ali_end   = int(x[18])
                except (IndexError, ValueError):
                    badline += 1
                    continue
                alilen = abs(ali_end - ali_start) + 1
                if i_eval > eval_thresh:
                    bade += 1
                    continue
                if alilen < min_ali:
                    short += 1
                    continue
                if min_bitscore is not None and bitscore < min_bitscore:
                    lowbits += 1
                    continue
                kept += 1
                s, e = ali_start, ali_end
                if s > e: s, e = e, s
                rows.append((prot, dom_name, dom_acc, hmm_len, i_eval, bitscore, s, e, alilen, source))
        df = pd.DataFrame(rows, columns=["protein_id","domain","acc","hmm_len","dom_i_evalue","bitscore","start","end","ali_len","source"])
        rep = {"source": source, "total_domtbl_lines": total, "kept_after_filters": kept,
               "filtered_high_i_eval": bade, "filtered_short_align": short, "bad_short_lines": badline,
               "filtered_low_bitscore": lowbits}
        return df, rep

    pf, pf_rep = parse_domtbl(str(stage2 / 'pfam.domtblout'), "pfam", EVAL, MIN_ALI)
    cu, cu_rep = parse_domtbl(str(stage2 / 'custom.domtblout'), "custom", EVAL, MIN_ALI)

    # Parse ASM with its own thresholds
    asm_eval = float(settings["ASM_MAX_EVALUE"])
    asm_min_ali = int(settings["ASM_MIN_ALI_LEN"])
    asm_min_bits = float(settings["ASM_MIN_BITSCORE"])
    asm, asm_rep = parse_domtbl(str(stage2 / 'asm.domtblout'), "asm", asm_eval, asm_min_ali, asm_min_bits)

    # ------------------------------------------------------------------------------
    # Parse Stage1 NBD hits and add as domain entries (FIX B1)
    # ------------------------------------------------------------------------------
    stage1_rows = []
    stage1_nbd_hits_path = STAGE1_NBD_HITS
    if os.path.exists(stage1_nbd_hits_path) and os.path.getsize(stage1_nbd_hits_path) > 0:
        try:
            stage1_hits = pd.read_csv(stage1_nbd_hits_path, sep="\t")
            for _, row in stage1_hits.iterrows():
                pid = row["protein_id"]
                hmm = row["nbd_hmm"]
                # Map Sordariales HMM names to standard domain names
                if "NACHT" in hmm:
                    domain = "NACHT"
                elif "NB-ARC" in hmm:
                    domain = "NB-ARC"
                else:
                    domain = hmm   # fallback (rare)
                acc = ""                     # no Pfam accession
                hmm_len = 0                 # not available
                evalue = row["dom_i_evalue"]
                bitscore = 0.0             # not available
                start = int(row["ali_start"])
                end   = int(row["ali_end"])
                ali_len = int(row["ali_len"])
                source = "stage1_nbd"
                stage1_rows.append([pid, domain, acc, hmm_len, evalue, bitscore, start, end, ali_len, source])
        except Exception as e:
            print(f"Warning: could not parse Stage1 NBD hits: {e}", file=sys.stderr)

    stage1_df = pd.DataFrame(stage1_rows, columns=[
        "protein_id","domain","acc","hmm_len","dom_i_evalue","bitscore","start","end","ali_len","source"
    ])

    # ------------------------------------------------------------------------------
    # Avoid FutureWarning by only concatenating non-empty frames
    # ------------------------------------------------------------------------------
    frames = [f for f in [pf, cu, stage1_df] if not f.empty]
    if frames:
        df = pd.concat(frames, ignore_index=True)
    else:
        df = pd.DataFrame(columns=["protein_id","domain","acc","hmm_len","dom_i_evalue","bitscore","start","end","ali_len","source"])
    
    df.to_csv(OUT_ALL, sep="\t", index=False)

    # --- Domain vocab ---
    def is_lrr(dom):
        return bool(re.search(r"\bLRR\b|Leucine[_ -]?rich|LRR_", dom, re.I))

    def group_domain(d):
        d_str = str(d).upper()
    
        # NBD domains (fungal NLR core - conservative)
        if re.search(r"\bNACHT\b", d_str): return "NACHT"
        if re.search(r"\bNB-?ARC\b", d_str): return "NB-ARC"
    
        # STAND-like (separate tracking, not for NBD calls)
        if re.search(r"\bSTAND\b", d_str): return "STAND_LIKE"
        if re.search(r"\bNBD\b", d_str) and not re.search(r"NB-?ARC", d_str): return "NBD_LIKE"
    
        # SSFR (Sequence-Similar Fold Repeats) - fungal NLR sensors
        if d_str.startswith("ANK"): return "ANK"
        if d_str.startswith("TPR"): return "TPR"
        if re.search(r"WD40|WD_40|WD.?REPEAT", d_str): return "WD40"
        if d_str.startswith("HEAT") or "ARM" in d_str: return "HEAT"
        if re.search(r"KELCH", d_str): return "KELCH"
    
        # LRR domains (EXCLUDE in fungi)
        if is_lrr(d): return "LRR"
    
        # --------------------------------------------------------------------------
        # Correct kinase grouping logic + Add kinase-like synonym
        # --------------------------------------------------------------------------
        if re.search(r"PROTEIN_KINASE|PKINASE|KINASE[-_ ]?LIKE", d_str):
            return "KINASE"
    
        # Non-SSFR sensors (valid fungal NLR sensors)
        if re.search(r"C2H2|ZINC.?FINGER", d_str): return "C2H2_ZF"
        if re.search(r"ZZ", d_str): return "ZZ"
        if re.search(r"SPRY|B30\.2", d_str): return "SPRY"
        if re.search(r"RING", d_str): return "RING"
    
        # Other domains for categorization (keep original names for debugging)
        if re.search(r"SAM", d_str): return "SAM"
        if re.search(r"SMC", d_str): return "SMC"
        if re.search(r"ABC", d_str): return "ABC"
        if re.search(r"SEPTIN", d_str): return "SEPTIN"
        if re.search(r"DYNAMIN", d_str): return "DYNAMIN"
        if re.search(r"HELICASE", d_str): return "HELICASE"
    
        # Return original domain name for unknown domains (not "OTHER")
        return d

    NBD_GROUPS = {"NACHT", "NB-ARC"}  # Conservative: only true fungal NLR NBDs
    STAND_LIKE_GROUPS = {"STAND_LIKE", "NBD_LIKE"}  # Track separately
    SSFR_GROUPS = {"WD40", "ANK", "TPR", "HEAT", "KELCH"}
    NON_SSFR_SENSOR_GROUPS = {"C2H2_ZF", "ZZ", "SPRY", "RING", "KINASE"}
    SENSOR_GROUPS = SSFR_GROUPS.union(NON_SSFR_SENSOR_GROUPS)
    LRR_GROUP = {"LRR"}

    # Housekeeping / false-positive domain patterns
    FP_DOM_RE = re.compile(r"(\bSMC\b|ABC_tran|ABC_transporter|Septin|Dynamin|Helicase_C|\bDEAD\b|DEXDc|ResIII|Restriction|Toprim|FtsK_SpoIIIE|CDC6|ORC[0-9]*|MCM[0-9]*|RAD[0-9]*|\bDNA\b.*\bREPAIR\b)", re.I)

    def interval_overlap(a1, a2, b1, b2):
        return max(0, min(a2, b2) - max(a1, b1) + 1)

    # Build per-protein sorted hits
    hit_map = {}
    if not df.empty:
        df["domain_group"] = df["domain"].astype(str).apply(group_domain)
        for pid, g in df.groupby("protein_id"):
            hit_map[pid] = g.sort_values(["start","end","dom_i_evalue"], ascending=[True, True, True]).copy()

    # ASM hits (if available) - only consider hits in N-terminal window
    asm_hits = {}
    if ("1" if enabled(settings["ASM_ENABLE"]) else "0") == "1" and not asm.empty:
        # ASM queries will have the same protein IDs as the original
        for pid, g in asm.groupby("protein_id"):
            if pid not in all_set:
                # Try to clean the protein ID (remove any window annotation)
                pid_clean = pid.split()[0] if " " in pid else pid
                if pid_clean in all_set:
                    pid = pid_clean
                else:
                    continue
            # Take best hit per protein
            best = g.loc[g["dom_i_evalue"].idxmin()]
            asm_hits[pid] = {
                "asm_domain": best["domain"],
                "asm_evalue": float(best["dom_i_evalue"]),
                "asm_bitscore": float(best["bitscore"]),
                "asm_start": int(best["start"]),
                "asm_end": int(best["end"]),
                "asm_window": "N_TERM"  # All ASM hits are in N-terminal window
            }

    summ=[]
    for pid in all_prots:
        g = hit_map.get(pid, None)

        doms_raw=[]
        doms_grp=[]
    
        # NBD interval sources
        nbd_iv_stage1 = stage1_nbd_iv.get(pid, None)
        nbd_iv_pfam = None
    
        has_nbd_pfam = False
        has_stand_like = False
        has_ssfr = False
        has_non_ssfr_sensor = False
        has_lrr = False
        flag_fp_domains = False
        flag_fp_overlaps_nbd = False
        integrated_domain_cterm = False
    
        order_invalid = False
        sensor_is_cterm = False
        has_effector_any = False
        has_effector_nterm = False
        effector_is_nterm_of_nbd = False

        if g is not None and not g.empty:
            doms_raw = g["domain"].astype(str).tolist()
            doms_grp = g["domain_group"].astype(str).tolist()

            has_lrr = any(d == "LRR" for d in doms_grp)
        
            # NBD from pfam hits (conservative)
            nbd_hits = g[g["domain_group"].isin(NBD_GROUPS)]
            if not nbd_hits.empty:
                has_nbd_pfam = True
                a = int(nbd_hits["start"].min())
                b = int(nbd_hits["end"].max())
                nbd_iv_pfam = (a, b) if a <= b else (b, a)
        
            # STAND-like domains (track separately)
            has_stand_like = any(d in STAND_LIKE_GROUPS for d in doms_grp)

            # Choose NBD interval priority: stage1 > pfam
            nbd_final = nbd_iv_stage1 if nbd_iv_stage1 is not None else nbd_iv_pfam

            # Sensors
            has_ssfr = any(d in SSFR_GROUPS for d in doms_grp)
            has_non_ssfr_sensor = any(d in NON_SSFR_SENSOR_GROUPS for d in doms_grp)
        
            # FP flags
            flag_fp_domains = any(FP_DOM_RE.search(d) for d in doms_raw)
        
            # Domain order and FP overlap
            if nbd_final is not None:
                nbd_start, nbd_end = nbd_final
            
                # ------------------------------------------------------------------
                # More nuanced order_invalid logic
                # Check if ANY sensor domain is clearly C-terminal to the NBD
                # Buffer of 10aa to allow for domain boundaries
                # ------------------------------------------------------------------
                if has_ssfr or has_non_ssfr_sensor:
                    sensor_hits = g[g["domain_group"].isin(SENSOR_GROUPS)]
                    if not sensor_hits.empty:
                        buffer = 10
                        # Find sensors that start at least buffer aa after the NBD ends
                        c_term_sensors = sensor_hits[sensor_hits["start"] >= (nbd_end + buffer)]
                        sensor_is_cterm = not c_term_sensors.empty
                        # Order is only invalid if NO sensor is C-terminal
                        # (some sensors may be N-terminal or overlapping, which is okay)
                        order_invalid = not sensor_is_cterm
            
                # Check for integrated kinase at C-terminus (not overlapping NBD)
                kinase_hits = g[g["domain_group"] == "KINASE"]
                if not kinase_hits.empty:
                    kinase_min_start = int(kinase_hits["start"].min())
                    if kinase_min_start > nbd_end:
                        integrated_domain_cterm = True
            
                # FP overlap with NBD
                if flag_fp_domains:
                    for _, r in g.iterrows():
                        if FP_DOM_RE.search(r["domain"]):
                            ov = interval_overlap(nbd_start, nbd_end, int(r["start"]), int(r["end"]))
                            if ov >= 20:  # meaningful overlap
                                flag_fp_overlaps_nbd = True
                                break

            # Effector detection (track ANY vs N-term-of-NBD)
            has_effector_any = False
            has_effector_nterm = False
            effector_is_nterm_of_nbd = False  # keep old column name for backwards compatibility

            eff_hits = []
            if g is not None and not g.empty:
                for _, row in g.iterrows():
                    if domain_base(row["domain"]) in EFFECTOR_GROUPS:
                        eff_hits.append((int(row["start"]), int(row["end"]), str(row["domain"])))

            has_effector_any = (len(eff_hits) > 0)

            if has_effector_any and nbd_final is not None:
                nbd_start, nbd_end = nbd_final
                # "N-term effector" if ANY effector ends before (NBD start - buffer),
                # optionally allowing a tiny overlap margin.
                cutoff = (nbd_start - EFFECTOR_NTERM_BUFFER_AA) + EFFECTOR_ALLOW_OVERLAP_AA
                has_effector_nterm = any(e_end <= cutoff for (e_start, e_end, _) in eff_hits)
                effector_is_nterm_of_nbd = has_effector_nterm

        # Stage1 membership
        has_nbd_stage1 = pid in stage1_nbd
    
        # Final NBD call (conservative: only stage1 or conservative Pfam NBDs)
        has_nbd = has_nbd_stage1 or has_nbd_pfam
    
        # Sensor present?
        has_sensor = has_ssfr or has_non_ssfr_sensor
    
        # ASM information
        asm_present = pid in asm_hits
        asm_info = asm_hits.get(pid, {})

        # ----------------------------------------------------------------------
        # NBD confidence: HIGH (stage1) or based on Pfam NBD e-value + coverage
        # ----------------------------------------------------------------------
        nbd_confidence = "LOW"
        if pid in stage1_nbd:
            nbd_confidence = "HIGH"
        elif has_nbd_pfam and not nbd_hits.empty:
            # Get the best Pfam NBD hit (lowest i-evalue)
            best = nbd_hits.loc[nbd_hits["dom_i_evalue"].idxmin()]
            best_eval = best["dom_i_evalue"]
            best_ali_len = best["ali_len"]
            best_hmm_len = best["hmm_len"]
            coverage = best_ali_len / best_hmm_len if best_hmm_len > 0 else 0.0

            # Check thresholds
            if best_eval <= NBD_CONF_EV_HIGH and coverage >= NBD_CONF_COV_HIGH:
                nbd_confidence = "HIGH"
            elif best_eval <= NBD_CONF_EV_MED and coverage >= NBD_CONF_COV_MED:
                nbd_confidence = "MEDIUM"
            # else stays LOW

        # ----------------------------------------------------------------------
        # has_any_repeat_region (based on SSFR groups detected)
        # ----------------------------------------------------------------------
        has_any_repeat_region = any(d in SSFR_GROUPS for d in doms_grp)

        summ.append({
            "protein_id": pid,
            "domains_raw": ";".join(doms_raw),
            "domains_grouped": ";".join(doms_grp),
            "has_nbd_stage1": bool(has_nbd_stage1),
            "has_nbd_pfam": bool(has_nbd_pfam),
            "has_nbd": bool(has_nbd),
            "has_stand_like": bool(has_stand_like),
            "nbd_interval_stage1": "" if nbd_iv_stage1 is None else f"{nbd_iv_stage1[0]}-{nbd_iv_stage1[1]}",
            "nbd_interval_pfam": "" if nbd_iv_pfam is None else f"{nbd_iv_pfam[0]}-{nbd_iv_pfam[1]}",
            "has_ssfr_sensor": bool(has_ssfr),
            "has_non_ssfr_sensor": bool(has_non_ssfr_sensor),
            "has_sensor": bool(has_sensor),
            "has_lrr": bool(has_lrr),
            "sensor_is_cterm_of_nbd": bool(sensor_is_cterm) if has_sensor and has_nbd else False,
            "order_invalid": bool(order_invalid) if has_sensor and has_nbd else False,
            "flag_fp_domains": bool(flag_fp_domains),
            "flag_fp_overlaps_nbd": bool(flag_fp_overlaps_nbd),
            "integrated_domain_cterm": bool(integrated_domain_cterm),
            "asm_present": bool(asm_present),
            "asm_domain": asm_info.get("asm_domain", ""),
            "asm_evalue": asm_info.get("asm_evalue", ""),
            "asm_bitscore": asm_info.get("asm_bitscore", ""),
            "asm_start": asm_info.get("asm_start", ""),
            "asm_end": asm_info.get("asm_end", ""),
            "asm_window": asm_info.get("asm_window", ""),
            "nbd_confidence": nbd_confidence,
            "has_any_repeat_region": bool(has_any_repeat_region),
            "has_effector": bool(has_effector_any),
            "has_effector_nterm": bool(has_effector_nterm),
            "effector_is_nterm_of_nbd": bool(effector_is_nterm_of_nbd),
        })

    out = pd.DataFrame(summ)

    # Sort columns logically
    cols_order = [
        "protein_id", "has_nbd", "has_nbd_stage1", "has_nbd_pfam", "has_stand_like",
        "has_sensor", "has_ssfr_sensor", "has_non_ssfr_sensor", "has_lrr",
        "has_effector", "has_effector_nterm", "effector_is_nterm_of_nbd",
        "sensor_is_cterm_of_nbd", "order_invalid",
        "flag_fp_domains", "flag_fp_overlaps_nbd", "integrated_domain_cterm",
        "asm_present", "asm_domain", "asm_evalue", "asm_bitscore", "asm_start", "asm_end", "asm_window",
        "nbd_confidence", "has_any_repeat_region",
        "nbd_interval_stage1", "nbd_interval_pfam",
        "domains_raw", "domains_grouped"
    ]

    # ------------------------------------------------------------------------------
    # Stage 2 PFAM enrichment (keep Stage 1 enrichment unchanged)
    # Foreground = STRICT stage1 NBD proteins (or relaxed if NBD_MODE=relaxed)
    # Background = candidates (union set) OR proteome (optional pfam_bg.domtblout)
    # Outputs:
    #   pfam_enriched_in_nbd_stage2.tsv
    #   pfam_enriched_top50_stage2.txt
    #   pfam_priority_ids_stage2.txt
    # ------------------------------------------------------------------------------
    OUT_ENRICH = str(stage2 / 'pfam_enriched_in_nbd_stage2.tsv')
    OUT_TOP    = str(stage2 / 'pfam_enriched_top50_stage2.txt')
    OUT_IDS    = str(stage2 / 'pfam_priority_ids_stage2.txt')

    PFAM_BG_MODE = settings.get("PFAM_ENRICH_BG", "proteome").strip().lower()
    if PFAM_BG_MODE not in ("proteome", "candidates"):
        PFAM_BG_MODE = "proteome"

    MIN_FG_COUNT = int(settings.get("PFAM_ENRICH_MIN_FG_COUNT", "2"))
    MIN_LOG2FC = float(settings.get("PFAM_ENRICH_MIN_LOG2FC", "1.0"))
    TOP_N = int(settings.get("PFAM_ENRICH_TOP_N", "50"))

    def extract_pfam_id(domain, acc):
        # Prefer PF accession like PFxxxxx
        a = str(acc) if acc is not None else ""
        m = re.search(r"(PF\d{5})", a)
        if m:
            return m.group(1)
        # Fall back: sometimes domain name includes PFxxxxx
        d = str(domain) if domain is not None else ""
        m = re.search(r"(PF\d{5})", d)
        if m:
            return m.group(1)
        return None

    def per_protein_pfam_sets(df_hits):
        # returns: dict pid -> set(PFxxxxx)
        mp = {}
        if df_hits is None or df_hits.empty:
            return mp
        sub = df_hits[df_hits["source"] == "pfam"].copy()
        if sub.empty:
            return mp
        # Ensure expected cols exist
        for c in ["protein_id", "domain", "acc"]:
            if c not in sub.columns:
                return mp
        for pid, g in sub.groupby("protein_id"):
            s = set()
            for _, r in g.iterrows():
                pfid = extract_pfam_id(r.get("domain"), r.get("acc"))
                if pfid:
                    s.add(pfid)
            if s:
                mp[pid] = s
        return mp

    # Foreground proteins: stage1 strict IDs you already wrote to stage1_nbd.ids
    fg_ids = set()
    if os.path.exists(STAGE1_NBD_IDS) and os.path.getsize(STAGE1_NBD_IDS) > 0:
        with open(STAGE1_NBD_IDS) as f:
            fg_ids = {l.strip() for l in f if l.strip()}

    # Candidate-set PFAM calls (union candidates scan)
    cand_pfam = per_protein_pfam_sets(df)

    # Background PFAM calls:
    #   - candidates: background universe = union candidates (all_set)
    #   - proteome:  background universe = full proteome FASTA (all_bg_prots)
    bg_pfam = cand_pfam
    bg_ids  = set(all_set)

    all_bg_prots = []
    if PFAM_BG_MODE == "proteome":
        try:
            all_bg_prots = fasta_ids(PROTEOME_FASTA)
        except Exception:
            all_bg_prots = []

        # If proteome FASTA not readable, force candidates mode
        if not all_bg_prots:
            PFAM_BG_MODE = "candidates"
        else:
            # Parse background domtblout (pfam vs proteome)
            bg_hits, _ = parse_domtbl(PFAM_BG_DOMTBL, "pfam", EVAL, MIN_ALI)
            if not bg_hits.empty:
                bg_hits["source"] = "pfam"
                bg_pfam = per_protein_pfam_sets(bg_hits)
                bg_ids  = set(all_bg_prots)   # denominator = full proteome
            else:
                print(f"[PFAM_ENRICH_STAGE2] WARNING: pfam_bg.domtblout empty; falling back to candidates background", file=sys.stderr)
                PFAM_BG_MODE = "candidates"
                bg_pfam = cand_pfam
                bg_ids  = set(all_set)

    # If we fell back to candidates, ensure bg_ids/bg_pfam are consistent
    if PFAM_BG_MODE != "proteome":
        bg_pfam = cand_pfam
        bg_ids  = set(all_set)

    # Foreground pfam sets restricted to final fg_ids
    fg_pfam = {pid: s for pid, s in cand_pfam.items() if pid in fg_ids}

    from collections import Counter
    bg_count = Counter()
    fg_count = Counter()

    N_fg = len(fg_ids)
    N_bg = len(bg_ids)

    for pid in sorted(bg_ids):
        s = bg_pfam.get(pid, set())
        if s:
            bg_count.update(sorted(s))

    for pid in sorted(fg_ids):
        s = fg_pfam.get(pid, set())
        if s:
            fg_count.update(sorted(s))

    rows = []
    if N_fg > 0 and N_bg > 0:
        import math
        for pfid, fg_n in fg_count.items():
            bg_n = bg_count.get(pfid, 0)
            if fg_n < MIN_FG_COUNT:
                continue
            fg_rate = (fg_n + 1) / (N_fg + 2)
            bg_rate = (bg_n + 1) / (N_bg + 2)
            log2fc = math.log2(fg_rate / bg_rate)
            if log2fc < MIN_LOG2FC:
                continue
            rows.append((pfid, fg_n, bg_n, fg_rate, bg_rate, log2fc))

    enr = pd.DataFrame(rows, columns=[
        "PFAM", "NBD_fg_with_PFAM", "BG_with_PFAM", "NBD_rate", "BG_rate", "log2FC"
    ]).sort_values(["log2FC", "NBD_fg_with_PFAM", "PFAM"], ascending=[False, False, True])

    enr.to_csv(OUT_ENRICH, sep="\t", index=False)

    top = enr.head(TOP_N)["PFAM"].tolist() if not enr.empty else []
    with open(OUT_TOP, "w") as f:
        f.write("\n".join(top) + ("\n" if top else ""))

    topset = set(top)
    prio_ids = []
    # Prioritize within the *union* candidates (these are what Stage 2/3 act on)
    for pid in all_prots:
        s = cand_pfam.get(pid, set())
        if s and (s & topset):
            prio_ids.append(pid)

    with open(OUT_IDS, "w") as f:
        f.write("\n".join(prio_ids) + ("\n" if prio_ids else ""))

    print(f"[PFAM_ENRICH_STAGE2] bg_mode={PFAM_BG_MODE}", file=sys.stderr)
    print(f"[PFAM_ENRICH_STAGE2] min_fg_count={MIN_FG_COUNT}", file=sys.stderr)
    print(f"[PFAM_ENRICH_STAGE2] min_log2fc={MIN_LOG2FC}", file=sys.stderr)
    print(f"[PFAM_ENRICH_STAGE2] top_n={TOP_N}", file=sys.stderr)
    print(f"[PFAM_ENRICH_STAGE2] N_fg={N_fg} (stage1 NBD ids)", file=sys.stderr)
    print(f"[PFAM_ENRICH_STAGE2] N_bg={N_bg}", file=sys.stderr)
    print(f"[PFAM_ENRICH_STAGE2] wrote: {OUT_ENRICH} rows={enr.shape[0]}", file=sys.stderr)
    print(f"[PFAM_ENRICH_STAGE2] wrote: {OUT_TOP} n={len(top)}", file=sys.stderr)
    print(f"[PFAM_ENRICH_STAGE2] wrote: {OUT_IDS} n={len(prio_ids)}", file=sys.stderr)

    if out.empty:
        out = out.reindex(columns=cols_order)

    out = out[[c for c in cols_order if c in out.columns] + 
              [c for c in out.columns if c not in cols_order]]
    out.to_csv(OUT_SUM, sep="\t", index=False)

    # Write parse reports
    rep_df = pd.DataFrame([pf_rep, cu_rep])
    if enabled(settings["ASM_ENABLE"]):
        rep_df = pd.concat([rep_df, pd.DataFrame([asm_rep])], ignore_index=True)
    rep_df.to_csv(OUT_REP, sep="\t", index=False)

    # Write ASM summary if available
    if True:  # Always write a headed summary, including zero ASM hits.
        asm_summary = out.loc[out["asm_present"].astype(bool)][["protein_id", "asm_domain", "asm_evalue", "asm_bitscore", "asm_window"]]
        asm_summary.to_csv(OUT_ASM_SUM, sep="\t", index=False)

    print("Domain architecture summary:", file=sys.stderr)
    print(f"Proteins analyzed: {out.shape[0]}", file=sys.stderr)
    print(f"Has NBD (conservative): {int(out['has_nbd'].sum())}", file=sys.stderr)
    print(f"Has STAND-like domains: {int(out['has_stand_like'].sum())}", file=sys.stderr)
    print(f"Has SSFR sensor: {int(out['has_ssfr_sensor'].sum())}", file=sys.stderr)
    print(f"Has non-SSFR sensor: {int(out['has_non_ssfr_sensor'].sum())}", file=sys.stderr)
    print(f"Has LRR (exclude): {int(out['has_lrr'].sum())}", file=sys.stderr)
    print(f"Order invalid (no C-terminal sensor): {int(out['order_invalid'].sum())}", file=sys.stderr)
    print(f"FP domains present: {int(out['flag_fp_domains'].sum())}", file=sys.stderr)
    print(f"FP domains overlap NBD: {int(out['flag_fp_overlaps_nbd'].sum())}", file=sys.stderr)
    print(f"Integrated domain C-terminal: {int(out['integrated_domain_cterm'].sum())}", file=sys.stderr)
    if ("1" if enabled(settings["ASM_ENABLE"]) else "0") == "1":
        print(f"ASM motifs detected (N-terminal window): {int(out['asm_present'].sum())}", file=sys.stderr)

    print(f"Wrote: {OUT_ALL}", file=sys.stderr)
    print(f"Wrote: {OUT_SUM}", file=sys.stderr)
    print(f"Wrote: {OUT_REP}", file=sys.stderr)
    if ("1" if enabled(settings["ASM_ENABLE"]) else "0") == "1" and not asm.empty:
        print(f"Wrote: {OUT_ASM_SUM}", file=sys.stderr)
    return out
