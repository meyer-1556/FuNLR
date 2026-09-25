#!/usr/bin/env python3
"""Deliberate executable test doubles. These do not perform biological analysis.

Installed under nine tool names by test_pipeline.py in a temporary PATH. Never
used by the installed package or real example commands.
"""

import json
import csv
import os
from pathlib import Path
import re
import sys
import time
import struct
import zlib


def fasta(path):
    records = []
    name, sequence = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(sequence)))
            name, sequence = line[1:].split()[0], []
        else:
            sequence.append(line.strip())
    if name is not None:
        records.append((name, "".join(sequence)))
    return records


def emit_fasta(records):
    for name, sequence in records:
        print(f">{name}")
        for start in range(0, len(sequence), 60):
            print(sequence[start : start + 60])


def option(args, flag):
    return args[args.index(flag) + 1]


def domain(name, pid, start, end, evalue="1e-30"):
    """HMMER3 domtblout fixed fields, followed by a synthetic description."""
    return " ".join(map(str, [
        name, "PF00000.1", 201, pid, "-", 800, "1e-30", 200, 0,
        1, 1, "1e-30", evalue, 200, 0, 1, 201,
        start, end, start, end, 0.99, "SYNTHETIC TEST HIT",
    ]))


def main():
    tool = Path(sys.argv[0]).name
    args = sys.argv[1:]
    if os.environ.get("FUNLR_TEST_DOUBLES") != "1":
        raise SystemExit("Refusing to run: synthetic test double, not a bioinformatics tool.")
    if args in (["--version"], ["-h"], ["--help"], ["version"]):
        print(f"{tool} SYNTHETIC-TEST-DOUBLE 1.0")
        return
    if tool == "Rscript" and "-e" in args:
        print("LIBRARIES\nSYNTHETIC_LIBRARY\nPACKAGES\ndplyr\tSYNTHETIC\nggplot2\tSYNTHETIC")
        return
    with open(os.environ["FUNLR_TEST_TOOL_LOG"], "a") as log:
        log.write(json.dumps({"tool": tool, "args": args, "pid": os.getpid()}) + "\n")
    if os.environ.get("FUNLR_TEST_WAIT_TOOL") == tool:
        time.sleep(30)  # Parent test interrupts this isolated CLI/process group.
    if os.environ.get("FUNLR_TEST_FAIL_TOOL") == tool:
        print(f"SYNTHETIC_REQUESTED_FAILURE: {tool}", file=sys.stderr)
        raise SystemExit(42)

    if tool == "hmmpress":
        for suffix in ("h3m", "h3i", "h3f", "h3p"):
            Path(args[-1] + "." + suffix).write_text("SYNTHETIC TEST INDEX\n")
        print("Pressed synthetic fixture database")
    elif tool == "hmmfetch":
        Path(option(args, "-o")).write_text(f"HMMER3/f SYNTHETIC TEST\nNAME  {args[-1]}\nLENG  201\n//\n")
    elif tool == "samtools":
        assert args[0] == "faidx", args
        path = args[1]
        assert len(args) == 2, "This double implements indexing only"
        with open(path + ".fai", "w") as out:
            for name, sequence in fasta(path):
                out.write(f"{name}\t{len(sequence)}\t0\t60\t61\n")
    elif tool == "seqkit":
        assert args[0] == "grep", args
        if "-f" in args:
            ids = set(Path(option(args, "-f")).read_text().splitlines())
        elif "-p" in args:
            pattern = option(args, "-p")
            # Exact ID matching is SeqKit's default. Regex needs -r/--use-regexp.
            if "-r" in args or "--use-regexp" in args:
                ids = {name for name, _ in fasta(args[-1]) if re.search(pattern, name)}
            else:
                ids = {pattern}
        else:
            raise AssertionError(args)
        emit_fasta((name, seq) for name, seq in fasta(args[-1]) if name in ids)
    elif tool == "hmmscan":
        out = Path(option(args, "--domtblout"))
        ids = {name for name, _ in fasta(args[-1])}
        rows = []
        if out.name == "nbd_whole.domtblout":
            for pid in ("canonical", "review", "rescue_edge", "rescue_internal", "variant"):
                if pid in ids:
                    rows.append(domain("NACHT", pid, 100, 279, "1e-5"))
            if "housekeeping" in ids:
                rows.append(domain("NACHT", "housekeeping", 100, 219, "1e-3"))
            if "bad_eval" in ids:
                rows.append(domain("NACHT", "bad_eval", 100, 279, "1.01e-3"))
            if "short_hit" in ids:
                rows.append(domain("NACHT", "short_hit", 100, 218, "1e-30"))
        elif out.name == "pfam_nbd.domtblout":
            for pid in ("canonical", "review", "rescue_edge", "rescue_internal", "variant"):
                if pid in ids:
                    rows.append(domain("NACHT", pid, 100, 300))
            if "housekeeping" in ids:
                rows.append(domain("NACHT", "housekeeping", 100, 219, "1e-4"))
        elif out.name in {"pfam.domtblout", "pfam_bg.domtblout"}:
            for pid in ("canonical", "review", "rescue_edge", "rescue_internal", "variant"):
                if pid in ids:
                    rows.append(domain("NACHT", pid, 100, 300))
            if "housekeeping" in ids:
                rows.append(domain("NACHT", "housekeeping", 100, 219, "1e-4"))
            for pid in ("canonical", "review", "repeat_only"):
                if pid in ids:
                    rows.append(domain("WD40", pid, 400, 419, "1e-3"))
            if "variant" in ids:
                rows.append(domain("HET", "variant", 20, 80))
            if "housekeeping" in ids:
                rows.append(domain("ABC_tran", "housekeeping", 400, 450))
            if "canonical" in ids:
                # Either should be removed before housekeeping flag/tier assignment.
                rows.append(domain("ABC_tran", "canonical", 500, 519, "1.01e-3"))
                rows.append(domain("Helicase_C", "canonical", 550, 568))
        out.write_text("# SYNTHETIC HMMER3 domtblout\n" + "".join(r + "\n" for r in rows))
        print("SYNTHETIC hmmscan output; not a similarity search")
    elif tool == "miniprot":
        assert "--gff" in args, args
        print("##gff-version 3")
        if os.environ.get("FUNLR_TEST_MINIPROT_NO_HITS") == "1":
            return
        for pid, _ in fasta(args[-1]):
            if pid == "rescue_internal":
                continue
            chrom = "chr_" + pid
            print(f"{chrom}\tminiprot\tmRNA\t450\t3000\t200\t+\t.\tID=MP_{pid};Target={pid} 1 800")
            print(f"{chrom}\tminiprot\tCDS\t450\t3000\t200\t+\t0\tParent=MP_{pid};Target={pid} 1 800")
    elif tool == "exonerate":
        pid = fasta(option(args, "--query"))[0][0]
        # Includes non-GFF prose, as exonerate's default output commonly does.
        print("Command line: SYNTHETIC TEST DOUBLE")
        print("# --- START OF GFF DUMP ---")
        print(f"chr_{pid}\texonerate:protein2genome\tgene\t19951\t22500\t200\t+\t.\tgene_id 0 ; sequence {pid} ;")
        print(f"chr_{pid}\texonerate:protein2genome\tcds\t19951\t22500\t.\t+\t.\t")
        print("# --- END OF GFF DUMP ---")
    elif tool == "gffread":
        ids = []
        for line in Path(args[0]).read_text().splitlines():
            fields = line.split("\t")
            if len(fields) == 9 and fields[2] == "mRNA":
                attrs = dict(field.split("=", 1) for field in fields[8].split(";") if "=" in field)
                ids.append(attrs["ID"])
        Path(option(args, "-y")).write_text("".join(f">{pid}\n{'M' + 'A' * 799}\n" for pid in ids))
    elif tool == "Rscript":
        # Plot orchestration is exercised without claiming that R was run.
        Path("domain_architecture_summary.tsv").write_text("status\nSYNTHETIC_PLOT_OUTPUT\n")
        Path("R_sessionInfo.txt").write_text("SYNTHETIC R TEST DOUBLE\n")
        def chunk(kind, data):
            return struct.pack('!I', len(data)) + kind + data + struct.pack('!I', zlib.crc32(kind + data))
        png = (b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('!2I5B', 1, 1, 8, 2, 0, 0, 0)) +
               chunk(b'IDAT', zlib.compress(b'\x00\xff\x00\xff')) + chunk(b'IEND', b''))
        names = ['nlr_domains_by_tier', 'sensor_distribution', 'nbd_distribution',
                 'asm_presence', 'length_vs_priority', 'order_sanity']
        for name in names:
            Path(name + '.png').write_bytes(png)
        Path('individual').mkdir(exist_ok=True)
        with Path('all_domain_hits.tsv').open() as handle:
            ids = {row['protein_id'] for row in csv.DictReader(handle, delimiter='\t')}
        for pid in ids:
            (Path('individual') / (pid + '.png')).write_bytes(png)
    else:
        raise AssertionError(f"Unexpected executable name: {tool}")


if __name__ == "__main__":
    main()
