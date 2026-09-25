"""Original scientific defaults, separated from execution paths."""
import copy
import yaml
import math
import re
from pathlib import Path
from .updated_defaults import UPDATED, PARAMETERS as NEW_PARAMETERS, FLAGS, ENUMS, STRINGS, NUMERIC, EVALUE, HIFI, DISCOVERY_PRESETS, UNUSED_LEGACY

DEFAULTS = {
    "EVAL_NBD": "1e-4", "EVAL_PFAM": "1e-3",
    "MIN_NBD_ALI_LEN": 150, "MIN_DOM_ALI_LEN": 20, "MAX_DOM_I_EVAL": "0",
    "CONTIG_END_BP": 10000, "MIN_RESCUE_FLANK_BP": 5000,
    "REPEAT_PROX_BP": 10000, "REPEAT_PROX_NEIGHBORS": 10,
    "NON_NLR_KEYWORDS_REGEX": "CDC6|ORC[0-9]*|ANAPC|APC|MCM[0-9]*|SMC|FMO|monooxygenase|helicase|dynein|kinesin|chaperone|HSP|ribosome|proteasome|ubiquitin|SNAP|clathrin|coatom|COPI|COPII|Sec[0-9]+|Rad[0-9]+|DNA repair|replication|cell cycle",
    "NON_NLR_DOMAIN_REGEX": "KAP_NTPase|ATPase_2|Helicase_C|DEAD|DEXDc|ResIII|Toprim|Pkinase|FAD|NAD|ABC_tran|SMC",
    "NON_NLR_GO_REGEX": "DNA replication|cell cycle|chromosome|helicase|repair|ubiquitin|proteasome|vesicle|Golgi|ER",
    "EXON_MAX_INTRON": 5000, "EXON_MIN_SCORE": 150, "EXON_MIN_PCT": 70,
    "EXON_MAXN": 200, "INCLUDE_TIER3": 0,
}

DEFAULTS.update(UPDATED)


def _validate_scientific(overrides=None):
    result = DEFAULTS.copy()
    result.update(overrides or {})
    for key, value in result.items():
        if key.endswith("_REGEX"):
            if not isinstance(value, str):
                raise ValueError(f"{key} must be a regular-expression string")
            try:
                re.compile(value)
            except re.error as exc:
                raise ValueError(f"Invalid {key}: {exc}") from exc
        elif key in FLAGS:
            if value not in (0, 1) or not isinstance(value, (int, bool)):
                raise ValueError(f"{key} must be 0/1 or boolean")
            result[key] = int(value)
        elif key in STRINGS:
            if not isinstance(value, str) or any(c in value for c in "\x00\n\r"):
                raise ValueError(f"{key} must be a single-line string")
            if key in ENUMS and value not in ENUMS[key]:
                raise ValueError(f"{key} must be one of {sorted(ENUMS[key])}")
        elif key in NUMERIC:
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"{key} must be a finite number")
            if key != "PFAM_ENRICH_MIN_LOG2FC" and value < 0:
                raise ValueError(f"{key} must be nonnegative")
            if key.startswith("NBD_CONF_COV_") and value > 1:
                raise ValueError(f"{key} must be between 0 and 1")
        elif key in EVALUE:
            try:
                number = float(value)
            except (ValueError, TypeError):
                raise ValueError(f"{key} must be numeric") from None
            if isinstance(value, bool) or not math.isfinite(number) or number < 0 or (key != "MAX_DOM_I_EVAL" and number == 0):
                raise ValueError(f"Invalid E-value for {key}: {value}")
            if key == "MAX_DOM_I_EVAL" and number == 0:
                result[key] = "0"
        elif not isinstance(value, int) or isinstance(value, bool) or value < 0:
            raise ValueError(f"{key} must be a nonnegative integer")
    for key in ("MIN_NBD_ALI_LEN", "MIN_DOM_ALI_LEN", "EXON_MAX_INTRON"):
        if result[key] == 0:
            raise ValueError(f"{key} must be positive")
    if result["INCLUDE_TIER3"] not in (0, 1) or result["EXON_MIN_PCT"] > 100:
        raise ValueError("INCLUDE_TIER3 must be 0/1 and EXON_MIN_PCT must be <=100")
    for key in UNUSED_LEGACY:
        if result[key] != DEFAULTS[key]:
            raise ValueError(f"{key} was not implemented by the analysis workflow; only its recorded default {DEFAULTS[key]!r} is supported")
    retained_retry = 10000 if result["NLR_PROFILE"] == "HIFI" else 0
    if result["EXON_MAX_INTRON_RETRY"] != retained_retry:
        raise ValueError("EXON_MAX_INTRON_RETRY was not implemented by the analysis workflow; retain the profile default")
    return result


# One source of scientific defaults, exposed through readable configuration sections.
PARAMETERS = {}
for section, keys in {
    "thresholds": ("EVAL_NBD", "EVAL_PFAM", "MIN_NBD_ALI_LEN", "MIN_DOM_ALI_LEN", "MAX_DOM_I_EVAL"),
    "context": ("CONTIG_END_BP", "MIN_RESCUE_FLANK_BP", "REPEAT_PROX_BP", "REPEAT_PROX_NEIGHBORS"),
    "false_positive": ("NON_NLR_KEYWORDS_REGEX", "NON_NLR_DOMAIN_REGEX", "NON_NLR_GO_REGEX"),
    "rescue": ("EXON_MAX_INTRON", "EXON_MIN_SCORE", "EXON_MIN_PCT", "EXON_MAXN", "INCLUDE_TIER3"),
}.items():
    PARAMETERS.update({key: (section, key.lower()) for key in keys})

PARAMETERS.update(NEW_PARAMETERS)

NESTED_DEFAULTS = {
    "sample": {"sample_id": "sample", "species_id": None, "assembly_version": "unspecified"},
    "inputs": {"genome": None, "proteins": None, "gff3": None, "eggnog": None, "annotations": None},
    "databases": {"pfam": None, "nbd_hmms": None, "custom_hmms": None, "asm_hmms": None, "effector_hmms": None, "sensor_hmms": None},
    "execution": {"output_dir": None, "cpu_threads": 1, "ram_gb": 8, "resume": False, "dry_run": False},
    "tools": {key: key for key in ("samtools", "hmmpress", "hmmscan", "seqkit", "miniprot", "exonerate", "hmmfetch", "gffread", "Rscript")},
}
for legacy, (section, name) in PARAMETERS.items():
    NESTED_DEFAULTS.setdefault(section, {})[name] = DEFAULTS[legacy]


class FunlrConfig:
    """Validated defaults < YAML/JSON < explicit CLI overrides.

    The stage API uses the nested configuration object. The earlier flat JSON
    scientific configuration remains accepted without duplicated default values.
    """

    def __init__(self, data=None):
        self.data = copy.deepcopy(NESTED_DEFAULTS)
        self._explicit = set()
        if data is not None:
            if not isinstance(data, dict):
                raise ValueError("Configuration must be a YAML/JSON mapping")
            if data and all(key in DEFAULTS for key in data):
                for legacy, value in data.items():
                    section, name = PARAMETERS[legacy]
                    self.data[section][name] = value
                    self._explicit.add((section, name))
            else:
                self._overlay(data)
        self._validate()

    def _overlay(self, data):
        for section, values in data.items():
            if section not in NESTED_DEFAULTS:
                raise ValueError(f"Unknown configuration key: {section}")
            if not isinstance(values, dict):
                raise ValueError(f"Configuration section {section} must be a mapping")
            for name, value in values.items():
                if name not in NESTED_DEFAULTS[section]:
                    raise ValueError(f"Unknown configuration key: {section}.{name}")
                self.data[section][name] = value
                self._explicit.add((section, name))

    def scientific_settings(self):
        values = {legacy: self.data[section][name] for legacy, (section, name) in PARAMETERS.items()}
        if isinstance(values["INCLUDE_TIER3"], bool):
            values["INCLUDE_TIER3"] = int(values["INCLUDE_TIER3"])
        return _validate_scientific(values)

    def _validate(self):
        preset = HIFI if self.data['context']['profile'] == 'HIFI' else {}
        discovery = self.data['discovery']['discovery_mode']
        preset = {**preset, **DISCOVERY_PRESETS.get(discovery, {})}
        for legacy in set(HIFI) | {k for values in DISCOVERY_PRESETS.values() for k in values}:
            section, name = PARAMETERS[legacy]
            if (section, name) not in self._explicit:
                self.data[section][name] = preset.get(legacy, DEFAULTS[legacy])
        for old, new in [('CONTIG_END_BP', 'HARD_END_BP'), ('INCLUDE_TIER3', 'INCLUDE_TIER3_IN_RESCUE')]:
            old_section, old_name = PARAMETERS[old]
            new_section, new_name = PARAMETERS[new]
            old_explicit = (old_section, old_name) in self._explicit
            new_explicit = (new_section, new_name) in self._explicit
            old_value = self.data[old_section][old_name]
            new_value = self.data[new_section][new_name]
            if old_explicit and new_explicit and old_value != new_value:
                raise ValueError(f"Conflicting legacy/new settings: {old} and {new}")
            if old_explicit and not new_explicit:
                self.data[new_section][new_name] = old_value
            else:
                self.data[old_section][old_name] = new_value
        for legacy, value in self.scientific_settings().items():
            section, name = PARAMETERS[legacy]
            self.data[section][name] = value
        for section in ("sample", "inputs", "databases", "tools"):
            for name, value in self.data[section].items():
                if value is not None and (not isinstance(value, str) or any(c in value for c in "\x00\n\r")):
                    raise ValueError(f"{section}.{name} must be a single-line string or null")
                if section == "tools" and not value:
                    raise ValueError(f"tools.{name} requires an executable name or path")
        execution = self.data["execution"]
        for name in ("cpu_threads", "ram_gb"):
            if type(execution[name]) is not int or execution[name] < 1:
                raise ValueError(f"execution.{name} must be a positive integer")
        for name in ("resume", "dry_run"):
            if type(execution[name]) is not bool:
                raise ValueError(f"execution.{name} must be a boolean")
        output = execution["output_dir"]
        if output is not None and (not isinstance(output, str) or any(c in output for c in "\x00\n\r")):
            raise ValueError("execution.output_dir must be a path string or null")

    @classmethod
    def from_yaml(cls, path):
        path = Path(path).expanduser().resolve()
        try:
            data = yaml.safe_load(path.read_text())
        except yaml.YAMLError as exc:
            raise ValueError(f"Invalid YAML/JSON configuration: {exc}") from exc
        config = cls({} if data is None else data)
        # Relative input/output paths belong to the config file's directory.
        for section in ("inputs", "databases"):
            for name, value in config.data[section].items():
                if value:
                    candidate = Path(value).expanduser()
                    config.data[section][name] = str((path.parent / candidate).resolve())
        value = config.get("execution", "output_dir")
        if value:
            config.data["execution"]["output_dir"] = str((path.parent / Path(value).expanduser()).resolve())
        for name, value in config.data["tools"].items():
            if "/" in value:
                config.data["tools"][name] = str((path.parent / Path(value).expanduser()).resolve())
        return config

    def apply_overrides(self, overrides):
        self._overlay({section: {key: value for key, value in values.items() if value is not None}
                       for section, values in overrides.items()})
        self._validate()

    def get(self, section, key, default=None):
        return self.data.get(section, {}).get(key, default)

    def __getitem__(self, section):
        return self.data[section]

    def section(self, section):
        return self.data[section]

    def to_dict(self):
        return copy.deepcopy(self.data)


def read_config(path=None):
    return (FunlrConfig.from_yaml(path) if path else FunlrConfig()).scientific_settings()


def write_example_yaml(profile="ILLUMINA", minimal=False):
    config = FunlrConfig({"context": {"profile": profile}}).to_dict()
    comment = "# FuNLR: paths are relative to this file; CLI flags override settings.\n"
    if minimal:
        config = {key: config[key] for key in ("inputs", "databases", "execution")}
        config["context"] = {"profile": profile}
        comment += "# Omitted scientific defaults are resolved by this FuNLR version and recorded in each run.\n"
    return comment + yaml.safe_dump(config, sort_keys=False)
