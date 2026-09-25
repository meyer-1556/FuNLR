#!/usr/bin/env python3
"""Generate the tiny synthetic FuNLR integration fixture (deterministic).

Usage: python3 tests/fixtures/make_synthetic_fixture.py tests/fixtures/synthetic

Layout produced under the target dir:
  genome.fa, proteins.faa, annotation.gff3, emapper.annotations,
  db/nbd.hmm, db/pfam.hmm, db/custom.hmm        (dummy HMM text; shim ignores content)
  canned/nbd.domtblout, canned/pfam.domtblout, canned/custom.domtblout (empty),
  canned/miniprot.gff3, canned/exonerate/<pid>.gff

Expected tiers (asserted by tests/integration/test_smoke.py):
  FUN_000001-T1 1A   FUN_000002-T1 2A   FUN_000003-T1 4B   FUN_000004-T1 4A
  FUN_000005-T1 1B   FUN_000006-T1 3    FUN_000007-T1 4A   FUN_000008-T1 4C
  FUN_000009-T1 2B   FUN_000010-T1 4C   FUN_000011-T1 1B   FUN_000012-T1 3
"""
import random
import sys
from pathlib import Path

random.seed(42)
OUT = Path(sys.argv[1])
(OUT / "db").mkdir(parents=True, exist_ok=True)
(OUT / "canned" / "exonerate").mkdir(parents=True, exist_ok=True)

# scaffold_3 / scaffold_4 each host exactly one gene so that the
# repeat-proximity scan (+-10 genes AND +-10 kb) finds no repeat neighbour
# and contig-end flags stay off for the interior tier-3 / tier-2B cases.
SCAF = {"scaffold_1": 60000, "scaffold_2": 30000, "scaffold_3": 50000, "scaffold_4": 50000}

def rand_dna(n):
    return "".join(random.choice("ACGT") for _ in range(n))

def rand_aa(n):
    return "".join(random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n))

def write_fasta(path, records, wrap=60):
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n")
            for i in range(0, len(seq), wrap):
                fh.write(seq[i:i + wrap] + "\n")

write_fasta(OUT / "genome.fa", [(k, rand_dna(v)) for k, v in SCAF.items()])

# (id, aa length, gff scaffold, start, end, strand, eggnog description, preferred, gos)
PROTS = [
    # canonical NLR mid-scaffold -> 1A
    ("FUN_000001-T1", 900, "scaffold_1", 20000, 24000, "+", "NACHT and WD40 domain-containing protein", "-", "-"),
    # NBD-only near contig end -> 2A
    ("FUN_000002-T1", 400, "scaffold_1", 2000, 4000, "+", "NACHT domain protein", "-", "-"),
    # housekeeping STAND false positive -> 4B
    ("FUN_000003-T1", 700, "scaffold_1", 30000, 33000, "-", "Component of the origin recognition complex (ORC) that binds origins of replication", "ORC4", "GO:0006260"),
    # Kelch-only repeat via eggnog priority -> 4A
    ("FUN_000004-T1", 800, "scaffold_1", 40000, 44000, "+", "Kelch motif", "-", "-"),
    # NLR architecture at scaffold edge -> 1B
    ("FUN_000005-T1", 850, "scaffold_2", 28500, 29900, "-", "NACHT and WD40 domain-containing protein", "-", "-"),
    # NBD + other domains, no repeat, isolated interior scaffold -> 3
    ("FUN_000006-T1", 600, "scaffold_3", 20000, 23000, "+", "Nacht domain-containing protein", "-", "-"),
    # short WD40-only -> 4A (+SHORT_PROTEIN)
    ("FUN_000007-T1", 300, "scaffold_2", 10000, 12000, "+", "WD40 repeat-like protein", "-", "-"),
    # TPR-only via eggnog; TPR_1/TPR_2 do NOT group (legacy \bTPR\b quirk) -> 4C
    ("FUN_000008-T1", 750, "scaffold_2", 15000, 18000, "-", "TPR-like protein", "-", "-"),
    # strict NBD-only, isolated interior scaffold -> 2B
    ("FUN_000009-T1", 500, "scaffold_4", 20000, 22500, "+", "hypothetical protein", "-", "-"),
    # unrelated; enters union only via the legacy HET-keyword quirk ("synthetase") -> 4C
    ("FUN_000010-T1", 450, "scaffold_1", 10000, 14000, "+", "glutamine synthetase", "GLN1", "GO:0006542"),
    # reversed order (Ank before NACHT) + contig end -> 1B with DOMAIN_ORDER_WEIRD
    ("FUN_000011-T1", 950, "scaffold_1", 55000, 58500, "-", "NACHT and ankyrin repeat protein", "-", "-"),
    # NBD + WD40 + Septin: FP domain that is NOT a non-NLR domain, interior -> 3
    ("FUN_000012-T1", 880, "scaffold_1", 46000, 49000, "+", "NACHT and WD40 domain-containing protein", "-", "-"),
]

write_fasta(OUT / "proteins.faa", [(p[0], rand_aa(p[1])) for p in PROTS])

with open(OUT / "annotation.gff3", "w") as fh:
    fh.write("##gff-version 3\n")
    for pid, _ln, scaf, s, e, strand, *_ in PROTS:
        gid = pid.split("-")[0]
        fh.write(f"{scaf}\tfixture\tgene\t{s}\t{e}\t.\t{strand}\t.\tID={gid}\n")
        fh.write(f"{scaf}\tfixture\tmRNA\t{s}\t{e}\t.\t{strand}\t.\tID={pid};Parent={gid}\n")

with open(OUT / "emapper.annotations", "w") as fh:
    fh.write("#query\tseed_ortholog\tevalue\tscore\teggnog_ogs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\n")
    for pid, _ln, _s, _st, _e, _str, desc, pref, gos in PROTS:
        fh.write(f"{pid}\t-\t0\t0\t-\t-\t-\t{desc}\t{pref}\t{gos}\t-\t-\n")

for name in ("nbd", "pfam", "custom"):
    (OUT / "db" / f"{name}.hmm").write_text(f"HMMER3/f [3.3.2]\nNAME  {name}_dummy\n//\n")

# ---- canned domtblout ----------------------------------------------------
# domtblout whitespace columns (only indices used by the parser matter):
#  0 target name, 1 accession, 3 query name, 12 domain i-Evalue,
# 17 ali start, 18 ali end  -> we pad with placeholders.
def dom_row(target, query, i_eval, ali_s, ali_e):
    cols = ["-"] * 22
    cols[0] = target
    cols[1] = f"PF00000.{random.randint(1, 9)}"
    cols[3] = query
    cols[12] = str(i_eval)
    cols[17] = str(ali_s)
    cols[18] = str(ali_e)
    return " ".join(cols) + " description\n"

nbd_rows = [
    dom_row("NACHT", "FUN_000001-T1", "1e-40", 200, 380),
    dom_row("NACHT", "FUN_000002-T1", "1e-30", 50, 230),
    dom_row("NACHT", "FUN_000003-T1", "5e-05", 300, 460),
    dom_row("NACHT", "FUN_000005-T1", "1e-38", 210, 390),
    dom_row("NACHT", "FUN_000006-T1", "1e-25", 100, 280),
    dom_row("NACHT", "FUN_000009-T1", "1e-20", 120, 300),
    dom_row("NACHT", "FUN_000011-T1", "1e-35", 600, 780),
    dom_row("NACHT", "FUN_000012-T1", "1e-33", 200, 380),
    # below-threshold decoys (i-Evalue too weak / alignment too short)
    dom_row("NACHT", "FUN_000010-T1", "0.5", 10, 160),
    dom_row("NACHT", "FUN_000001-T1", "1e-50", 200, 250),  # too short (51 < 150)
]
(OUT / "canned" / "nbd.domtblout").write_text("# canned\n" + "".join(nbd_rows))

pfam_rows = [
    # FUN_000001: NACHT + WD40 correct order
    dom_row("NACHT", "FUN_000001-T1", "1e-40", 200, 380),
    dom_row("WD40", "FUN_000001-T1", "1e-10", 500, 540),
    dom_row("WD40", "FUN_000001-T1", "1e-09", 600, 640),
    # FUN_000002: NACHT only
    dom_row("NACHT", "FUN_000002-T1", "1e-30", 50, 230),
    # FUN_000003: housekeeping signature + NACHT
    dom_row("AAA_16", "FUN_000003-T1", "1e-20", 100, 250),
    dom_row("ORC4_C", "FUN_000003-T1", "1e-15", 400, 550),
    dom_row("NACHT", "FUN_000003-T1", "5e-05", 300, 460),
    # FUN_000004: Kelch only
    dom_row("Kelch_1", "FUN_000004-T1", "1e-12", 100, 150),
    dom_row("Kelch_2", "FUN_000004-T1", "1e-11", 300, 350),
    # FUN_000005: NACHT + WD40 correct order
    dom_row("NACHT", "FUN_000005-T1", "1e-38", 210, 390),
    dom_row("WD40", "FUN_000005-T1", "1e-08", 700, 740),
    # FUN_000006: AAA + NACHT, no repeat
    dom_row("AAA_16", "FUN_000006-T1", "1e-22", 80, 200),
    dom_row("NACHT", "FUN_000006-T1", "1e-25", 300, 480),
    # FUN_000007: WD40 only
    dom_row("WD40", "FUN_000007-T1", "1e-07", 50, 100),
    # FUN_000008: TPR only (ungrouped by legacy \bTPR\b regex)
    dom_row("TPR_1", "FUN_000008-T1", "1e-09", 100, 140),
    dom_row("TPR_2", "FUN_000008-T1", "1e-08", 400, 440),
    # FUN_000009: strict NBD only
    dom_row("NACHT", "FUN_000009-T1", "1e-20", 120, 300),
    # FUN_000011: Ank BEFORE NACHT (reversed)
    dom_row("Ank_2", "FUN_000011-T1", "1e-12", 100, 140),
    dom_row("NACHT", "FUN_000011-T1", "1e-35", 600, 780),
    # FUN_000012: NACHT + WD40 + Septin (FP domain, NOT a non-NLR domain)
    dom_row("NACHT", "FUN_000012-T1", "1e-33", 200, 380),
    dom_row("WD40", "FUN_000012-T1", "1e-07", 700, 740),
    dom_row("Septin", "FUN_000012-T1", "1e-18", 500, 650),
]
(OUT / "canned" / "pfam.domtblout").write_text("# canned\n" + "".join(pfam_rows))
(OUT / "canned" / "custom.domtblout").write_text("")

# canned miniprot gff: resolves FUN_000002-T1 and FUN_000005-T1; FUN_000009-T1
# (and FUN_000011-T1) stay unresolved so they land on the exonerate refine list.
def mp_gff(pid, scaf, s, e, strand):
    return (
        f"{scaf}\tminiprot\tmRNA\t{s}\t{e}\t1000\t{strand}\t.\tID=mp_{pid};Target={pid} 1 400\n"
        f"{scaf}\tminiprot\tCDS\t{s}\t{e}\t1000\t{strand}\t0\tParent=mp_{pid};Target={pid} 1 400\n"
    )

(OUT / "canned" / "miniprot.gff3").write_text(
    "##gff-version 3\n"
    + mp_gff("FUN_000002-T1", "scaffold_1", 1500, 4500, "+")
    + mp_gff("FUN_000005-T1", "scaffold_2", 28300, 29950, "-")
)

# attribute style matches real exonerate --showtargetgff (`sequence <pid>`),
# which the legacy stage-6.3 parser regex keys on
(OUT / "canned" / "exonerate" / "FUN_000009-T1.gff").write_text(
    "scaffold_4\texonerate\tgene\t19800\t22600\t600\t+\t.\tgene_id 0 ; sequence FUN_000009-T1\n"
    "scaffold_4\texonerate\texon\t19800\t22600\t600\t+\t.\tsequence FUN_000009-T1\n"
)

print(f"fixture written to {OUT}")
