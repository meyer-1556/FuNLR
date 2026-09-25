"""reference dominant-span tripartite reporting labels; slot order is conceptual."""
import re
import sys
from pathlib import Path
import pandas as pd

EFFECTOR_LIKE = {
    "Goodbye", "HeLo", "HeLo-like", "HET", "HET-s", "TIR",
    "Patatin", "PNP_UDP", "RelA_SpoT", "Ses", "CHAT",
    "Crinkler", "SAM", "C2", "Peptidase_S8", "CARD/Pyrin", "CC",
    # Fold ASM-like motifs into the effector-side slot for tripartite output
    "HRAM", "PP", "sigma", "PUASM", "Basidio_ASM", "BASS", "Lineage_ASM"
}

NBD_LIKE = {
    "NACHT", "NB-ARC", "AAA"
}

SENSOR_LIKE = {
    "WD40", "ANK", "TPR", "HEAT", "LRR", "SPRY",
    "Kinase", "Zinc_finger", "ZZ", "HMA/WRKY/LIM"
}

EFFECTOR_TIEBREAK = {
    "Goodbye": 0, "HeLo": 1, "HeLo-like": 2, "HET": 3, "HET-s": 4, "TIR": 5,
    "Patatin": 6, "PNP_UDP": 7, "RelA_SpoT": 8, "Ses": 9, "CHAT": 10,
    "Crinkler": 11, "SAM": 12, "C2": 13, "Peptidase_S8": 14, "CARD/Pyrin": 15, "CC": 16,
    "HRAM": 17, "PP": 18, "sigma": 19, "PUASM": 20, "Basidio_ASM": 21, "BASS": 22, "Lineage_ASM": 23
}
NBD_TIEBREAK = {
    "NACHT": 0, "NB-ARC": 1, "AAA": 2
}
SENSOR_TIEBREAK = {
    "WD40": 0, "ANK": 1, "TPR": 2, "HEAT": 3, "LRR": 4,
    "SPRY": 5, "Kinase": 6, "Zinc_finger": 7, "ZZ": 8, "HMA/WRKY/LIM": 9
}

def simplify_domain_token(tok):
    t = re.sub(r"__.*$", "", str(tok).strip())
    if not t or t.lower() in {"na", "none", "nan", "-", "--"}:
        return None

    # ---------- Effector-like ----------
    if re.search(r"Goodbye", t, re.IGNORECASE):
        return "Goodbye"
    if re.search(r"HeLo-like|HELL", t, re.IGNORECASE):
        return "HeLo-like"
    if re.search(r"\bHeLo\b", t, re.IGNORECASE):
        return "HeLo"
    if re.search(r"HET-s", t, re.IGNORECASE):
        return "HET-s"
    if re.search(r"\bHET\b", t, re.IGNORECASE):
        return "HET"
    if re.search(r"\bTIR\b", t, re.IGNORECASE):
        return "TIR"
    if re.search(r"Patatin", t, re.IGNORECASE):
        return "Patatin"
    if re.search(r"PNP_UDP|PUP", t, re.IGNORECASE):
        return "PNP_UDP"
    if re.search(r"RelA_SpoT", t, re.IGNORECASE):
        return "RelA_SpoT"
    if re.search(r"Ses[AB]|ses[ab]", t, re.IGNORECASE):
        return "Ses"
    if re.search(r"CHAT", t, re.IGNORECASE):
        return "CHAT"
    if re.search(r"Crinkler", t, re.IGNORECASE):
        return "Crinkler"
    if re.search(r"\bSAM\b", t, re.IGNORECASE):
        return "SAM"
    if re.search(r"C2[ -]?domain|PF00168|\bC2\b", t, re.IGNORECASE):
        return "C2"
    if re.search(r"Peptidase_S8", t, re.IGNORECASE):
        return "Peptidase_S8"
    if re.search(r"CARD|Pyrin", t, re.IGNORECASE):
        return "CARD/Pyrin"
    if re.search(r"CC[_-]?|Coiled[_-]?coil", t, re.IGNORECASE):
        return "CC"

    # ---------- NBD-like ----------
    if re.search(r"NACHT", t, re.IGNORECASE):
        return "NACHT"
    if re.search(r"NB-ARC", t, re.IGNORECASE):
        return "NB-ARC"
    if re.search(r"AAA_16|AAA_22|\bAAA\b", t, re.IGNORECASE):
        return "AAA"

    # ---------- Sensor-like ----------
    if re.search(r"WD40|WD_40|WD-repeat", t, re.IGNORECASE):
        return "WD40"
    if re.search(r"ANK|Ank|ANKYRIN", t, re.IGNORECASE):
        return "ANK"
    if re.search(r"\bTPR\b", t, re.IGNORECASE):
        return "TPR"
    if re.search(r"HEAT", t, re.IGNORECASE):
        return "HEAT"
    if re.search(r"LRR|Leucine[ -]?rich", t, re.IGNORECASE):
        return "LRR"
    if re.search(r"SPRY", t, re.IGNORECASE):
        return "SPRY"
    if re.search(r"PKinase|Kinase", t, re.IGNORECASE):
        return "Kinase"
    if re.search(r"C2H2|Zinc[_ -]?finger|ZINC_FINGER", t, re.IGNORECASE):
        return "Zinc_finger"
    if re.search(r"ZZ[_-]?|ZZ-type", t, re.IGNORECASE):
        return "ZZ"
    if re.search(r"HMA|WRKY|LIM", t, re.IGNORECASE):
        return "HMA/WRKY/LIM"

    # ---------- ASM-like; collapse into effector-side slot ----------
    if re.search(r"HRAM", t, re.IGNORECASE):
        return "HRAM"
    if re.search(r"PP[ -]?motif|NLR07|NLR39", t, re.IGNORECASE):
        return "PP"
    if re.search(r"\bsigma\b", t, re.IGNORECASE):
        return "sigma"
    if re.search(r"PUASM|NLR32", t, re.IGNORECASE):
        return "PUASM"
    if re.search(r"NLR05|NLR08|NLR22|NLR29|NLR44", t, re.IGNORECASE):
        return "Basidio_ASM"
    if re.search(r"BASS", t, re.IGNORECASE):
        return "BASS"
    if re.search(r"NLR17|NLR19|NLR34", t, re.IGNORECASE):
        return "Lineage_ASM"

    return None

def simplified_category(label):
    if label in EFFECTOR_LIKE:
        return "effector"
    if label in NBD_LIKE:
        return "nbd"
    if label in SENSOR_LIKE:
        return "sensor"
    return None

def merged_span(intervals):
    spans = []
    for s, e in intervals:
        try:
            s = int(s)
            e = int(e)
        except Exception:
            continue
        if e < s:
            s, e = e, s
        spans.append((s, e))

    if not spans:
        return 0

    spans.sort()
    total = 0
    cur_s, cur_e = spans[0]
    for s, e in spans[1:]:
        if s <= cur_e:
            if e > cur_e:
                cur_e = e
        else:
            total += (cur_e - cur_s + 1)
            cur_s, cur_e = s, e
    total += (cur_e - cur_s + 1)
    return total

def build_tripartite_architecture(domain_hits_path):
    p = Path(domain_hits_path)
    if not p.exists() or p.stat().st_size == 0:
        print(f"No Stage 2 domain-hit table found at {p}; NLR_architecture will be blank.", file=sys.stderr)
        return {}

    hits = pd.read_csv(p, sep="\t")
    required = {"protein_id", "domain", "start", "end"}
    missing = required - set(hits.columns)
    if missing:
        raise SystemExit(f"ERROR: {p} missing required columns for architecture calculation: {sorted(missing)}")

    hits = hits[["protein_id", "domain", "start", "end"]].copy()
    hits["protein_id"] = hits["protein_id"].astype(str).str.strip()
    hits["domain_simple"] = hits["domain"].apply(simplify_domain_token)
    hits["slot"] = hits["domain_simple"].apply(simplified_category)
    hits["start"] = pd.to_numeric(hits["start"], errors="coerce")
    hits["end"] = pd.to_numeric(hits["end"], errors="coerce")

    hits = hits[
        hits["protein_id"].notna() &
        (hits["protein_id"] != "") &
        hits["domain_simple"].notna() &
        hits["slot"].notna() &
        hits["start"].notna() &
        hits["end"].notna()
    ].copy()

    architecture_map = {}

    for pid, g in hits.groupby("protein_id", sort=False):
        slot_best = {}

        for slot, sg in g.groupby("slot", sort=False):
            if slot == "effector":
                tie_map = EFFECTOR_TIEBREAK
            elif slot == "nbd":
                tie_map = NBD_TIEBREAK
            else:
                tie_map = SENSOR_TIEBREAK

            ranked = []
            for label, lg in sg.groupby("domain_simple", sort=False):
                span = merged_span(zip(lg["start"], lg["end"]))
                start_min = int(lg["start"].min())
                tie_rank = tie_map.get(label, 999)
                ranked.append((span, tie_rank, start_min, label))

            ranked.sort(key=lambda x: (-x[0], x[1], x[2], x[3]))
            slot_best[slot] = ranked[0][3]

        ordered = [
            slot_best.get("effector", ""),
            slot_best.get("nbd", ""),
            slot_best.get("sensor", "")
        ]
        architecture_map[pid] = ";".join([x for x in ordered if x])

    return architecture_map

