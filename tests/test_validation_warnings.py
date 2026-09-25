"""Optional-input diagnostics distinguish analysis, export and provenance roles."""
from types import SimpleNamespace
from funlr.validation import validate_inputs


def test_optional_input_messages_describe_actual_roles(tmp_path):
    genome = tmp_path / 'genome.fa'
    proteins = tmp_path / 'proteins.faa'
    gff = tmp_path / 'genes.gff3'
    model = tmp_path / 'placeholder.hmm'
    genome.write_text('>c\nATGGCTTAA\n')
    proteins.write_text('>p\nMA\n')
    gff.write_text('c\ttest\tmRNA\t1\t9\t.\t+\t.\tID=p;Parent=g\n')
    model.write_text('placeholder for input-file validation only\n')
    args = SimpleNamespace(genome=genome, proteins=proteins, gff3=gff,
                           pfam=model, nbd_hmms=model, scientific={'ASM_ENABLE': 0})
    _, messages, _ = validate_inputs(args)
    by_input = {line.split(' omitted;', 1)[0]: line for line in messages}
    assert 'does not determine candidate calls' in by_input['annotations']
    for role in ('effector_hmms', 'sensor_hmms'):
        assert 'provenance unavailable' in by_input[role]
        assert 'scanning uses custom_hmms' in by_input[role]
    assert 'classification may differ' in by_input['custom_hmms']
    assert 'enrichment evidence unavailable' in by_input['eggnog']
    assert 'disabled in the resolved configuration' in by_input['asm_hmms']
