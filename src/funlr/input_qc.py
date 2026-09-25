"""Descriptive input measurements that do not alter candidate selection."""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote

SCHEMA_VERSION = 1


def _percent(numerator, denominator):
    return 100.0 * numerator / denominator if denominator else None


def _fasta_metrics(path, *, genome):
    """Read once, retaining lengths and counters instead of whole sequences.

    Ordinary input validation owns alphabet, duplicate-ID and coordinate
    checks. These measurements preserve sequence case and stop characters.
    """
    path = Path(path)
    digest = hashlib.sha256()
    lengths = []
    bases = Counter()
    current_length = 0
    current_stops = 0
    first = last = ""
    have_record = False
    internal_stop_records = terminal_stop_records = starts_m_records = 0

    def finish():
        nonlocal internal_stop_records, terminal_stop_records, starts_m_records
        if not have_record:
            return
        lengths.append(current_length)
        terminal = int(last == "*")
        terminal_stop_records += terminal
        internal_stop_records += int(current_stops > terminal)
        starts_m_records += int(first.upper() == "M")

    with path.open("rb") as handle:
        for line_number, raw in enumerate(handle, 1):
            digest.update(raw)
            line = raw.decode("utf-8").strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                if not line[1:].split():
                    raise ValueError(f"Empty FASTA identifier at line {line_number}: {path}")
                have_record = True
                current_length = current_stops = 0
                first = last = ""
                continue
            if not have_record:
                raise ValueError(f"Sequence before FASTA header at line {line_number}: {path}")
            if not first:
                first = line[0]
            last = line[-1]
            current_length += len(line)
            current_stops += line.count("*")
            if genome:
                bases.update(line)
    finish()
    total = sum(lengths)
    result = {
        "source": {"path": str(path.resolve()), "sha256": digest.hexdigest()},
        "n_sequences": len(lengths),
        "total_sequence_characters": total,
        "min_length": min(lengths) if lengths else None,
        "max_length": max(lengths) if lengths else None,
        "mean_length": total / len(lengths) if lengths else None,
    }
    if genome:
        cumulative = 0
        n50 = None
        for length in sorted(lengths, reverse=True):
            cumulative += length
            if total and cumulative * 2 >= total:
                n50 = length
                break
        acgt = sum(bases[base] for base in "ACGTacgt")
        lower_acgt = sum(bases[base] for base in "acgt")
        n_bases = bases["N"] + bases["n"]
        result.update(
            n50_bp=n50,
            n_bases=n_bases,
            pct_n_bases=_percent(n_bases, total),
            acgt_bases=acgt,
            lowercase_acgt_bases=lower_acgt,
            pct_lowercase_acgt=_percent(lower_acgt, acgt),
            other_sequence_characters=total - acgt - n_bases,
        )
    else:
        result.update(
            n_with_internal_stop=internal_stop_records,
            pct_with_internal_stop=_percent(internal_stop_records, len(lengths)),
            n_with_terminal_stop=terminal_stop_records,
            n_starting_with_m=starts_m_records,
        )
    return result


def _gff_metrics(path, *, genome_lengths=None):
    """Count declared features and transcript-parent relationships.

    Parent groups are not assumed to be independent biological loci. Feature
    coordinates, CDS phases and translated gene models are not evaluated here.
    """
    path = Path(path)
    digest = hashlib.sha256()
    features = Counter()
    gene_ids, transcript_ids = set(), set()
    parent_transcripts = defaultdict(set)
    without_parent = multiple_parents = 0
    with path.open("rb") as handle:
        for line_number, raw in enumerate(handle, 1):
            digest.update(raw)
            line = raw.decode("utf-8").rstrip("\r\n")
            if line.startswith("##FASTA"):
                raise ValueError("Supply annotation-only GFF3; embedded FASTA is unsupported")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"GFF3 line {line_number} needs 9 fields: {path}")
            if genome_lengths is not None:
                try:
                    start, end = int(fields[3]), int(fields[4])
                except ValueError:
                    raise ValueError(f"Non-integer GFF3 coordinates at line {line_number}: {path}") from None
                if fields[0] not in genome_lengths or not (1 <= start <= end <= genome_lengths[fields[0]]):
                    raise ValueError(f"GFF3 coordinates outside genome at line {line_number}: {path}")
            feature = fields[2]
            features[feature] += 1
            attributes = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            identifier = unquote(attributes.get("ID", ""))
            if feature == "gene" and identifier:
                gene_ids.add(identifier)
            if feature not in {"mRNA", "transcript"}:
                continue
            if genome_lengths is not None and (not identifier or identifier in transcript_ids or fields[6] not in {"+", "-"}):
                raise ValueError(f"GFF3 transcript needs a unique ID and +/- strand at line {line_number}: {path}")
            if identifier:
                transcript_ids.add(identifier)
            # Split before decoding: %2C can belong to one identifier.
            parents = {unquote(value) for value in attributes.get("Parent", "").split(",") if value}
            without_parent += int(not parents)
            multiple_parents += int(len(parents) > 1)
            for parent in parents:
                if identifier:
                    parent_transcripts[parent].add(identifier)
    return {
        "source": {"path": str(path.resolve()), "sha256": digest.hexdigest()},
        "n_features": sum(features.values()),
        "feature_counts": dict(sorted(features.items())),
        "n_gene_features": features["gene"],
        "n_distinct_gene_ids": len(gene_ids),
        "n_transcript_features": features["mRNA"] + features["transcript"],
        "n_distinct_transcript_ids": len(transcript_ids),
        "n_transcript_features_without_parent": without_parent,
        "n_transcript_features_with_multiple_parents": multiple_parents,
        "n_parent_ids_referenced_by_transcripts": len(parent_transcripts),
        "n_parent_ids_with_multiple_transcripts": sum(len(ids) > 1 for ids in parent_transcripts.values()),
        "n_declared_genes_with_multiple_transcripts": sum(
            len(parent_transcripts[gene]) > 1 for gene in gene_ids
        ),
    }


def collect_input_qc(genome, proteins, gff3, *, validated=False):
    """Return schema-1 measurements of the three primary inputs, read only.

    Standalone use validates FASTA syntax and GFF coordinates. Set ``validated``
    only after normal pipeline preflight. No models or external tools are
    required. Notes provide interpretation, not thresholds or a quality grade;
    these statistics never change tier rules.
    """
    genome, proteins, gff3 = (Path(path).expanduser().resolve() for path in (genome, proteins, gff3))
    genome_lengths = None
    if not validated:
        from .validation import fasta_lengths
        genome_lengths = fasta_lengths(genome)
        fasta_lengths(proteins, proteins=True)
    genome_report = _fasta_metrics(genome, genome=True)
    proteins_report = _fasta_metrics(proteins, genome=False)
    gff_report = _gff_metrics(gff3, genome_lengths=genome_lengths)
    notes = []
    if genome_report["lowercase_acgt_bases"]:
        notes.append("Lowercase A/C/G/T bases are present. Case information alone does not establish how the assembly was masked.")
    if genome_report["n_bases"]:
        notes.append("N/n bases represent unknown nucleotides; their presence does not establish why those bases are unknown.")
    if proteins_report["n_with_internal_stop"]:
        notes.append("Some input proteins contain internal '*' characters. Review their annotations and translation history; this observation alone does not establish a pseudogene.")
    if not gff_report["n_gene_features"]:
        notes.append("No gene feature rows are declared. Transcript-parent counts remain available, but are not an inferred count of biological loci.")
    if gff_report["n_transcript_features_without_parent"]:
        notes.append("Some transcript features have no Parent attribute, so their parent groups cannot be counted.")
    return {
        "schema_version": SCHEMA_VERSION,
        "purpose": "Descriptive primary-input measurements; no quality grade, biological validation or changes to candidate selection.",
        "validation_scope": "FASTA syntax and GFF3 structure/coordinates; complete protein-to-annotation reconciliation belongs to normal run preflight.",
        "genome": genome_report,
        "proteins": proteins_report,
        "gff3": gff_report,
        "definitions": {
            "sequence_lengths": "Counts of sequence characters, including N/ambiguity symbols and protein '*' characters; FASTA headers and line wrapping are excluded.",
            "pct_n_bases": "100 * N/n bases / all genome sequence characters.",
            "pct_lowercase_acgt": "100 * lowercase a/c/g/t bases / all uppercase or lowercase A/C/G/T bases.",
            "internal_stop": "At least one '*' before the final sequence character; one terminal '*' is counted separately.",
            "starting_with_m": "First sequence character is M or m. This is not an ORF-completeness assessment.",
            "n50_bp": "Contig length at which descending cumulative sequence length reaches half the assembly length; includes ambiguous bases.",
            "annotation_counts": "Declared GFF3 features, distinct IDs and explicit transcript-parent groups. Counts do not collapse alleles, haplotypes, isoforms or duplicated annotations into biological loci.",
            "null": "A measurement with no available denominator or sequence length is null rather than zero.",
        },
        "notes": notes,
    }


def write_input_qc(report, output):
    """Write one new JSON artifact in an existing directory, without overwrite."""
    path = Path(output).expanduser()
    if path.is_symlink() or any(parent.is_symlink() for parent in path.parents):
        raise ValueError("Input-QC output and its parents cannot be symlinks")
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    with path.open("x") as handle:
        handle.write(encoded)
    return path
