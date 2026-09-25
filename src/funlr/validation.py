"""Fail early when the original workflow's identifier assumptions do not hold."""
import re
from pathlib import Path

INPUT_FLAGS = {
    "genome": "GENOME_FILE", "proteins": "PROTEINS_FILE", "gff3": "GFF_FILE",
    "eggnog": "EGGNOG_FILE", "nbd_hmms": "NBD_HMMS", "pfam": "PFAM_DB",
    "custom_hmms": "CUSTOM_HMMS", "annotations": "ANNOTATIONS_FILE",
    "asm_hmms": "ASM_HMMS", "effector_hmms": "COMBINED_EFFECTOR_HMMS", "sensor_hmms": "COMBINED_SENSOR_HMMS",
}


def fasta_lengths(path, *, proteins=False):
    lengths = {}
    name = None
    with Path(path).open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"Empty FASTA identifier in {path}:{line_number}")
                name = fields[0]
                if name in lengths:
                    raise ValueError(f"Duplicate FASTA identifier {name!r} in {path}")
                if proteins and (not re.fullmatch(r"[A-Za-z0-9_.:|+\-]+", name) or name in (".", "..") or name.startswith("-")):
                    raise ValueError(f"Unsupported protein ID {name!r}; use simple IDs without whitespace, slashes, or shell metacharacters")
                lengths[name] = 0
            else:
                if name is None:
                    raise ValueError(f"Expected uncompressed FASTA in {path}:{line_number}")
                if re.search(r"\s", line) or not re.fullmatch(r"[A-Za-z*]+", line):
                    raise ValueError(f"Invalid unaligned FASTA sequence in {path}:{line_number}")
                lengths[name] += len(line)
    if not lengths or any(length == 0 for length in lengths.values()):
        raise ValueError(f"FASTA must contain nonempty sequences: {path}")
    return lengths


def validate_inputs(args):
    inputs = {}
    warnings = []
    for flag, var in INPUT_FLAGS.items():
        value = getattr(args, flag, None)
        if not value:
            if flag in ("eggnog", "custom_hmms", "annotations", "effector_hmms", "sensor_hmms") or (flag == "asm_hmms" and not getattr(args, "scientific", {}).get("ASM_ENABLE", 1)):
                notes = {
                    "eggnog": "annotation-description and annotated Pfam enrichment evidence unavailable",
                    "custom_hmms": "additional custom-domain evidence unavailable; candidate classification may differ",
                    "annotations": "annotation-table join unavailable; this input does not determine candidate calls",
                    "effector_hmms": "separate effector-library provenance unavailable; architecture scanning uses custom_hmms",
                    "sensor_hmms": "separate sensor-library provenance unavailable; architecture scanning uses custom_hmms",
                    "asm_hmms": "ASM scanning is disabled in the resolved configuration",
                }
                warnings.append(f"{flag} omitted; {notes[flag]}")
                continue
            raise ValueError(f"Missing required --{flag.replace('_', '-')}")
        path = Path(value).expanduser().resolve()
        if any(c in str(path) for c in "\n\r\x00"):
            raise ValueError("Input paths cannot contain newlines or NUL characters")
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty --{flag.replace('_', '-')}: {path}")
        inputs[var] = path

    contigs = fasta_lengths(inputs["GENOME_FILE"])
    proteins = fasta_lengths(inputs["PROTEINS_FILE"], proteins=True)
    transcripts = set()
    with inputs["GFF_FILE"].open() as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                raise ValueError("GFF3 with embedded FASTA is unsupported; supply annotation-only GFF3")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GFF3 line {line_number} must have 9 tab-separated fields")
            scaffold, _, kind, start, end, _, strand, _, attrs = fields
            try:
                start, end = int(start), int(end)
            except ValueError:
                raise ValueError(f"Non-integer GFF3 coordinates at line {line_number}") from None
            if scaffold not in contigs or not (1 <= start <= end <= contigs[scaffold]):
                raise ValueError(f"GFF3 coordinates outside genome at line {line_number}: {scaffold}:{start}-{end}")
            if kind in ("mRNA", "transcript"):
                attributes = dict(item.split("=", 1) for item in attrs.split(";") if "=" in item)
                tid = attributes.get("ID")
                if not tid or tid in transcripts or strand not in ("+", "-"):
                    raise ValueError(f"GFF3 transcript needs a unique ID and +/- strand at line {line_number}")
                transcripts.add(tid)
    from .stages._inputs import reconcile_coordinates
    coordinates, mapping_counts = reconcile_coordinates(inputs["GFF_FILE"], inputs["PROTEINS_FILE"])
    mapped = set(coordinates["protein_id"]) & set(proteins)
    missing = set(proteins) - mapped
    if missing:
        raise ValueError("Protein IDs have no GFF3 transcript, CDS protein_id or gene-ID mapping; unmatched: " + ", ".join(sorted(missing)[:8]))
    if not transcripts:
        raise ValueError("GFF3 contains no mRNA/transcript features")
    unmatched_rows = int(mapping_counts.get("rows_with_protein_id_not_in_fasta", 0))
    if unmatched_rows:
        warnings.append(f"{unmatched_rows} GFF3 coordinate rows have no protein sequence")
    return inputs, warnings, {"contigs": len(contigs), "proteins": len(proteins), "transcripts": len(transcripts)}
