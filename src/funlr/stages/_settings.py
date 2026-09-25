"""Settings consumed by the discovery stages.

The public configuration is the authoritative validated interface. These
defaults also support pure data-transformation helpers.
"""

DEFAULTS = {
    "EVAL_NBD": "1e-5", "MIN_NBD_ALI_LEN": 180,
    "EVAL_NBD_RELAXED": "1e-3", "MIN_NBD_ALI_LEN_RELAXED": 120,
    "EVAL_PFAM": "1e-3", "MIN_DOM_ALI_LEN": 20, "MAX_DOM_I_EVAL": "0",
    "PFAM_NBD_SCAN": 1, "DISCOVERY_MODE": "BALANCED",
    "INCLUDE_DESC_KEYWORDS_IN_UNION": 0, "DESC_KEYWORD_REQUIRE_STRUCTURAL_HINT": 1,
    "DESC_KEYWORD_STRUCTURAL_HINT_MODE": "ANY",
    "DESC_KEYWORD_STRUCTURAL_HINT_SOURCES": "RELAXED_NBD,PFAM_NBD,PFAM_ENRICHED",
    "DESC_PRIORITY_KEYWORDS_REGEX": r"NACHT|NB-ARC|STAND|NOD-like|NOD-like receptor|NLR|P-loop NTPase|AAA\+|HET|het-e|incompatibility|vegetative incompatibility|HeLo|HELL|Goodbye|gasdermin|MLKL|necrosis|cell death|amyloid|prion|amyloid signaling|HRAM|PP-motif|PUASM|FASS|ASM|WD40|WD-repeat|ANK|ankyrin|TPR|tetratricopeptide|HEAT|kelch|SPRY|B30\.2|zinc finger|C2H2|ZZ-type|protein kinase|PKinase",
    "UNION_INTERSECT_WITH_PROTEOME_IDS": 1, "USE_PFAM_ENRICHMENT_FOR_UNION": 1,
    "PFAM_ENRICH_MIN_FG_COUNT": 2, "PFAM_ENRICH_MIN_LOG2FC": 1.0,
    "PFAM_ENRICH_TOP_N": 50, "PFAM_ENRICH_BG": "candidates", "NBD_MODE": "strict",
    "ASM_ENABLE": 1, "ASM_WINDOW_NTERM_AA": 200, "ASM_MAX_EVALUE": "1e-3",
    "ASM_MIN_BITSCORE": 20, "ASM_MIN_ALI_LEN": 12,
    "NBD_CONF_EV_HIGH": "1e-6", "NBD_CONF_COV_HIGH": 0.45,
    "NBD_CONF_EV_MED": "1e-3", "NBD_CONF_COV_MED": 0.25,
}


def settings(config):
    values = DEFAULTS.copy()
    values.update(config.scientific_settings())
    return values


def enabled(value):
    return str(value).lower() in {"1", "true"}
