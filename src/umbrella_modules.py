"""Selection, identities and cache checks for the optional umbrella modules.

Modules (MAJIQ, DaPars2, QAPA) are opt-in tools layered on the umbrella's
staged reads and shared HISAT2 alignment. Each keeps small per-run caches
that later comparisons read, so regrouping never needs the reads again.

Identities, all computed, never typed by hand:

  module_reference_id  tool, pinned software, module settings and the base
                       reference_id; a new version or catalogue invalidates
                       only that tool (the base reference_id never changes)
  cache_id             a run's processing identity (umbrella_manifest.
                       RUN_IDENTITY) plus the module_reference_id; labels,
                       comparisons and the selection list stay out
  analysis_id          the selected runs, their biological replicates, the
                       module_reference_id and the analysis settings

Modes (`umbrella_modules_mode`):
  cached_only  (default) every selected run must already have its caches;
               nothing is staged, aligned or quantified for a module
  ingest       missing caches are made; already processed runs still need
               `umbrella_allow_restage: true` for their reads
"""

import hashlib
import json

MODULES = ("majiq", "dapars2", "qapa")
MODES = ("cached_only", "ingest")

# pinned software, part of each module_reference_id
SOFTWARE = {
    "majiq": "majiq v3 (licensed source, see majiq_source)",
    "dapars2": "DaPars2 3UTR/DaPars2@fb81c6cce1a1bbd8093cba327040f16e4ccc7cb6",
    "qapa": "QAPA morrislab/qapa@v1.4.1 + salmon 1.10.3",
}
DAPARS2_COMMIT = "fb81c6cce1a1bbd8093cba327040f16e4ccc7cb6"


def _digest(payload, length=16):
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()[:length]


def parse_selection(config):
    """(selected modules, mode) from config, with actionable errors."""
    raw = config.get("umbrella_modules", []) or []
    if isinstance(raw, str):
        raw = [item.strip() for item in raw.split(",") if item.strip()]
    if not isinstance(raw, list):
        raise ValueError("umbrella_modules must be a list, e.g. [majiq, dapars2, qapa]")
    unknown = sorted(set(raw) - set(MODULES))
    if unknown:
        raise ValueError("unknown umbrella_modules {}; choose from {}".format(
            ", ".join(unknown), ", ".join(MODULES)))
    if len(raw) != len(set(raw)):
        raise ValueError("umbrella_modules lists a module twice")
    mode = str(config.get("umbrella_modules_mode", "cached_only"))
    if mode not in MODES:
        raise ValueError("umbrella_modules_mode must be cached_only or ingest, not {}".format(mode))
    restage = str(config.get("umbrella_allow_restage", False)).lower() in ("true", "t", "1", "yes")
    if mode == "cached_only" and restage:
        raise ValueError("umbrella_modules_mode: cached_only never stages reads; drop "
                         "umbrella_allow_restage or use umbrella_modules_mode: ingest")
    return [module for module in MODULES if module in raw], mode


def module_settings(config, module):
    """The config that changes a module's per-run caches (not its comparisons)."""
    if module == "majiq":
        return {"strandness": config.get("majiq_strandness", "manifest"),
                "source": str(config.get("majiq_source", config.get("majiq_bin_folder", "")))}
    if module == "dapars2":
        return {"coverage": "bedtools genomecov -split -bg, UTR windows, primary mapped reads"}
    return {"library": config.get("qapa_library_type", "A"),
            "polya_sites": str(config.get("qapa_polya_sites", "")),
            "decoys": bool(config.get("qapa_decoys", False))}


def module_reference_id(reference_id, module, config):
    return _digest({"reference_id": reference_id, "module": module,
                    "software": SOFTWARE[module], "settings": module_settings(config, module)})


def cache_id(run, module_reference):
    from src.umbrella_manifest import RUN_IDENTITY
    return _digest({"run": {field: getattr(run, field) for field in RUN_IDENTITY},
                    "module_reference_id": module_reference})


def analysis_id(preflight, module_reference, settings):
    return _digest({"runs": preflight["runs"], "collapse": preflight["collapse"],
                    "module_reference_id": module_reference, "settings": settings})


# file names each module keeps per run (relative to its cache directory)
CACHE_FILES = {
    "majiq": ("run.sj", "cache.json"),
    "dapars2": ("utr_coverage.bedgraph.gz", "depth.json", "cache.json"),
    "qapa": ("quant.sf.gz", "lib_format_counts.json", "cache.json"),
}


def missing_caches(runs, cache_dir_of, modules):
    """[(run_id, module, path)] for caches that do not exist (cached_only check)."""
    import os
    missing = []
    for module in modules:
        for run in runs:
            for name in CACHE_FILES[module]:
                path = "{}/{}".format(cache_dir_of(run, module), name)
                if not os.path.exists(path):
                    missing.append((run.run_id, module, path))
                    break
    return missing


def write_status(path, tool, status, reasons=(), **extra):
    """Compact status JSON: ok, unsupported (design) or no_annotation; failures raise."""
    record = dict(tool=tool, status=status, reasons=list(reasons), **extra)
    with open(path, "w") as stream:
        json.dump(record, stream, sort_keys=True, indent=2)
        stream.write("\n")


# normalized result columns shared by every module (A minus B throughout)
NORMALIZED = ("tool", "tool_version", "reference_id", "module_reference_id", "analysis_id",
              "project_id", "comparison_id", "feature_id", "gene_id", "gene_name", "chrom",
              "start", "end", "coordinates", "strand", "event_class", "effect", "effect_definition",
              "n_a", "n_b", "coverage", "native_statistic", "statistic_type", "p_value",
              "q_value", "status")
