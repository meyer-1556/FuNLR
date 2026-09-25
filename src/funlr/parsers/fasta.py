"""FASTA helpers ported verbatim-in-spirit from the legacy pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Iterator, TextIO


def fasta_iter(fp: TextIO) -> Iterator[tuple[str, str]]:
    """Yield (id, sequence); ID is the first whitespace-delimited header token."""
    name = None
    seq: list[str] = []
    for line in fp:
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def fasta_ids(path: str | Path) -> list[str]:
    ids: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def write_fasta(records: Iterator[tuple[str, str]], out: TextIO, wrap: int = 60) -> int:
    """Write (id, seq) pairs wrapped at `wrap` columns; returns record count."""
    n = 0
    for rid, seq in records:
        out.write(f">{rid}\n")
        for i in range(0, len(seq), wrap):
            out.write(seq[i : i + wrap] + "\n")
        n += 1
    return n


def dedup_keep_first(in_path: str | Path, out_path: str | Path) -> tuple[int, int]:
    """Deduplicate a FASTA by header ID, keeping first occurrence (stable).

    Returns (total_records, kept_records). Ported from legacy stage 5.
    """
    seen: set[str] = set()
    kept = 0
    total = 0
    with open(in_path) as f, open(out_path, "w") as o:
        cur_id = None
        buf: list[str] = []

        def flush() -> None:
            nonlocal total, kept
            if cur_id is None:
                return
            total += 1
            if cur_id not in seen:
                seen.add(cur_id)
                kept += 1
                o.writelines(buf)

        for line in f:
            if line.startswith(">"):
                flush()
                cur_id = line[1:].split()[0]
                buf = [line]
            else:
                buf.append(line)
        flush()
    return total, kept
