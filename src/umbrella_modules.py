"""Selection, identities and cache checks for the optional umbrella modules.

Modules (MAJIQ, DaPars2, QAPA, vast-tools) are opt-in tools layered on the
umbrella's staged reads and shared HISAT2 alignment. Each keeps small per-run caches
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

MODULES = ("majiq", "dapars2", "qapa", "vasttools", "vasttools_turbo", "rmats_conservative", "psisigma",
           "ciriquant", "circexplorer2", "dcc", "whippet_circ")
# shared producers that are never selected on their own (their own IDs and caches)
INTERNAL_MODULES = ("circrna_star",)
MODES = ("cached_only", "ingest")
SCOPES = ("all", "selected")
# modules that run only on the comparisons of umbrella_module_selection unless told otherwise
DEFAULT_SCOPE = {"ciriquant": "selected", "circexplorer2": "selected", "dcc": "selected",
                 "whippet_circ": "selected"}
CIRCRNA_STAR = ("STAR 2.7.11b chimeric-only (--outSAMtype None, --chimOutType Junctions, legacy format); "
                "DCC joint-mapping settings; profile v1")

# pinned software, part of each module_reference_id
SOFTWARE = {
    "rmats_conservative": "rMATS-turbo 4.3.0, no novel splice sites; profile v1",
    "psisigma": ("PSI-Sigma wososa/PSI-Sigma@bfcb1423f1ec61ac48b6f92250505ec5855f4a53 (v2.3); PDL 2.106, "
                 "PDL-GSL 2.103, PDL-Stats 0.856, Statistics-Multtest 0.14; HISAT2 adapter v1"),
    "majiq": "majiq v3 (licensed source, see majiq_source)",
    "dapars2": "DaPars2 3UTR/DaPars2@fb81c6cce1a1bbd8093cba327040f16e4ccc7cb6",
    "qapa": "QAPA morrislab/qapa@v1.4.1 + salmon 1.10.3",
    "vasttools": "vast-tools vastgroup/vast-tools@5598e1f1b0067611dc9e5abb285b4eacaf543ee3 (2.5.1) + bowtie 1.3.1",
    "vasttools_turbo": ("vast-tools vastgroup/vast-tools@5598e1f1b0067611dc9e5abb285b4eacaf543ee3 (2.5.1) + "
                        "bowtie 1.3.1; chunked align + vast-tools merge, profile v1"),
    "ciriquant": ("CIRIquant 1.1.3 (bioconda) + hisat2 2.2.1; native counting with the load_fai and grouper "
                  "Python 3 fixes (src/run_ciriquant.py); pseudo-circular index built once per catalogue with "
                  "CIRIquant's generate_index/build_index; counts export v1"),
    "circrna_star": CIRCRNA_STAR,
    "circexplorer2": "CIRCexplorer2 2.3.8 parse -t STAR + annotate; on " + CIRCRNA_STAR,
    "dcc": "DCC 0.5.0 (dieterich-lab/DCC@afd5cdd6) detection, CircSkip off; on " + CIRCRNA_STAR,
    "whippet_circ": ("Whippet.jl whippet-quant --circ (whippet_bin_folder; version recorded per cache) on the "
                     "umbrella Whippet index; keeps psi, jnc, gene and isoform TPM; profile v1"),
}
DAPARS2_COMMIT = "fb81c6cce1a1bbd8093cba327040f16e4ccc7cb6"
# vast-tools master at 3 Oct 2026 (VERSION 2.5.1; tagged releases stop at v2.5.1, 2020)
VASTTOOLS_COMMIT = "5598e1f1b0067611dc9e5abb285b4eacaf543ee3"
# VASTDB library per vast-tools assembly code (vast-tools install.R; current release)
VASTDB_RELEASE = "23.06.20"
VASTDB_LIBRARIES = {"hg38": "hs2", "hg19": "hsa", "mm10": "mm2", "mm9": "mmu", "dm6": "dme"}


def vastdb_library(config):
    """(species, library file) for vasttools_species; unknown assemblies are refused."""
    species = str(config.get("vasttools_species", "hg38"))
    if species not in VASTDB_LIBRARIES:
        raise ValueError("vasttools_species {} has no VASTDB library here; choose from {}".format(
            species, ", ".join(sorted(VASTDB_LIBRARIES))))
    return species, "vastdb.{}.{}.tar.gz".format(VASTDB_LIBRARIES[species], VASTDB_RELEASE)


def _digest(payload, length=16):
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()[:length]


def file_sha256(path, _cache={}):
    """sha256 of a file's bytes, or a marker naming the missing path."""
    import os
    path = str(path)
    if not os.path.isfile(path):
        return "missing:" + path
    key = (path, os.path.getmtime(path), os.path.getsize(path))
    if key not in _cache:
        digest = hashlib.sha256()
        with open(path, "rb") as stream:
            for block in iter(lambda: stream.read(1 << 20), b""):
                digest.update(block)
        _cache[key] = digest.hexdigest()
    return _cache[key]


def circrna_catalogue_settings(config):
    """What the frozen circRNA catalogue is made from: sources (content hashes), the
    liftover chain for sources on another assembly, and the build filters."""
    assembly = str(config.get("circrna_assembly", "hg38"))
    sources = []
    for entry in config.get("circrna_catalogue_sources", []) or []:
        if isinstance(entry, str):
            entry = {"path": entry}
        if not isinstance(entry, dict) or "path" not in entry:
            raise ValueError("circrna_catalogue_sources entries need a path (and optionally assembly, name)")
        source_assembly = str(entry.get("assembly", assembly))
        role = str(entry.get("role", "quantify"))
        if role not in ("quantify", "annotate"):
            raise ValueError("circrna_catalogue_sources role must be quantify or annotate, not {}".format(role))
        fmt = str(entry.get("format", "bed6"))
        if fmt not in ("bed6", "circatlas", "circbase", "circpedia"):
            raise ValueError("circrna_catalogue_sources format must be bed6, circatlas, circbase or circpedia")
        sources.append({"name": str(entry.get("name", "")), "assembly": source_assembly, "role": role,
                        "format": fmt, "tissue": str(entry.get("tissue", "")), "sha256": file_sha256(entry["path"])})
    if sources and not any(source["role"] == "quantify" for source in sources):
        raise ValueError("circrna_catalogue_sources needs at least one source with role quantify")
    lifted = any(source["assembly"] != assembly for source in sources)
    chain = config.get("circrna_liftover_chain", "")
    return {"assembly": assembly, "sources": sources,
            "liftover_chain": file_sha256(chain) if lifted and chain else ("" if not lifted else "missing"),
            "max_span": int(config.get("circrna_max_span", 0)),
            "max_pseudo_bases": int(config.get("circrna_max_pseudo_bases", 2_000_000_000)),
            "build": "src/circrna_catalogue.py v2 (format converters; BSJ boundaries lifted one by one; "
                     "annotate sources add names only)"}


def circrna_catalogue_id(reference_id, config):
    return _digest({"reference_id": reference_id, "catalogue": circrna_catalogue_settings(config)})


def parse_scope(config):
    """{module: 'all' | 'selected'} for every module (umbrella_module_scope overrides defaults)."""
    raw = config.get("umbrella_module_scope", {}) or {}
    if not isinstance(raw, dict):
        raise ValueError("umbrella_module_scope must be a mapping, e.g. {ciriquant: selected}")
    unknown = sorted(set(raw) - set(MODULES))
    if unknown:
        raise ValueError("umbrella_module_scope names unknown modules {}; choose from {}".format(
            ", ".join(unknown), ", ".join(MODULES)))
    bad = sorted(module for module, value in raw.items() if value not in SCOPES)
    if bad:
        raise ValueError("umbrella_module_scope values must be all or selected (check {})".format(", ".join(bad)))
    scope = {module: DEFAULT_SCOPE.get(module, "all") for module in MODULES}
    scope.update(raw)
    return scope


SELECTION_COLUMNS = ("comparison_id", "tool", "role", "library_evidence", "decided_by", "approved_by", "reason")


def load_module_selection(path, comparison_ids, scope):
    """Approved per-comparison selection for 'selected'-scope modules.

    Returns {tool: {comparison_id: row}}. Unknown comparisons or tools, all-scope tools,
    duplicates and rows without approved_by are refused; the selection is provenance,
    not identity (it enters no cache, reference or analysis ID)."""
    import csv
    selection = {}
    if not path:
        return selection
    with open(path, newline="") as stream:
        reader = csv.DictReader((line for line in stream if not line.startswith("#")), delimiter="\t")
        missing = set(SELECTION_COLUMNS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError("umbrella_module_selection {} lacks columns {}".format(path, ", ".join(sorted(missing))))
        for number, row in enumerate(reader, start=2):
            row = {key: (value or "").strip() for key, value in row.items() if key}
            where = "umbrella_module_selection line {}".format(number)
            if row["comparison_id"] not in comparison_ids:
                raise ValueError("{}: comparison {} is not in umbrella_comparisons".format(where, row["comparison_id"]))
            if row["tool"] not in scope:
                raise ValueError("{}: unknown tool {}".format(where, row["tool"]))
            if scope[row["tool"]] != "selected":
                raise ValueError("{}: {} runs on every comparison (scope all); set umbrella_module_scope: "
                                 "{{{}: selected}} to select comparisons for it".format(where, row["tool"], row["tool"]))
            if not row["approved_by"]:
                raise ValueError("{}: approved_by is empty; a person approves each row before a launch".format(where))
            if row["comparison_id"] in selection.get(row["tool"], {}):
                raise ValueError("{}: {} is selected twice for {}".format(where, row["comparison_id"], row["tool"]))
            selection.setdefault(row["tool"], {})[row["comparison_id"]] = row
    return selection


def parse_selection(config):
    """(selected modules, mode) from config, with actionable errors."""
    raw = config.get("umbrella_modules", []) or []
    if isinstance(raw, str):
        raw = [item.strip() for item in raw.split(",") if item.strip()]
    if not isinstance(raw, list):
        raise ValueError("umbrella_modules must be a list, e.g. [majiq, dapars2, qapa, vasttools]")
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


def flag(value):
    """A config boolean, also when given as a string on the command line."""
    return str(value).strip().lower() in ("true", "t", "1", "yes")


def module_settings(config, module):
    """The config that changes a module's per-run caches (not its comparisons)."""
    if module == "ciriquant":
        return {"catalogue": circrna_catalogue_settings(config), "anchor": int(config.get("ciriquant_anchor", 5)),
                "library_type": "from umbrella strandedness (auto: inferred junction strandedness)",
                "counts": "every catalogue circRNA; BSJ and FSJ in read pairs"}
    if module == "circrna_star":
        return {"sjdb_overhang": int(config.get("circrna_star_sjdb_overhang", 100)),
                "annotation": "umbrella annotation_gtf", "genome": "umbrella genome_fasta"}
    if module == "circexplorer2":
        return dict(module_settings(config, "circrna_star"), refflat="umbrella annotation_gtf, gene_name or gene_id")
    if module == "dcc":
        return dict(module_settings(config, "circrna_star"))
    if module == "whippet_circ":
        # the core Whippet quant's flags plus --circ, so linear events stay comparable to it
        return {"flags": str(config.get("whippet_flags", "") or ""), "circ": True,
                "whippet": str(config.get("whippet_bin_folder", "")), "index": "umbrella Whippet index"}
    if module == 'psisigma':
        # irmode 0 only: intron retention then depends on the annotation alone
        return {'canonical_only': flag(config.get('psisigma_canonical_only', True)), 'irmode': 0,
                'annotation': 'umbrella annotation_gtf with ME_DB microexons, primary chromosomes',
                'adapter': 'NH:i:1 as MAPQ 255; spliced records grouped by samtools collate'}
    if module == 'rmats_conservative':
        return {'novel_splice_sites': False, 'library_type': 'fr-unstranded',
                'variable_read_length': True, 'profile': 1}
    if module == "majiq":
        return {"strandness": config.get("majiq_strandness", "manifest"),
                "source": str(config.get("majiq_source", config.get("majiq_bin_folder", "")))}
    if module == "vasttools":
        species = str(config.get("vasttools_species", "hg38"))
        library = "vastdb.{}.{}.tar.gz".format(VASTDB_LIBRARIES.get(species, "?"), VASTDB_RELEASE)
        return {"species": species, "vastdb": library, "ir_version": int(config.get("vasttools_ir_version", 2)),
                "expr": bool(config.get("vasttools_expr", False))}
    if module == "vasttools_turbo":
        # chunks change the per-run counts slightly (rounding), so they are identity;
        # threads per chunk do not
        return dict(module_settings(config, "vasttools"), chunks=int(config.get("vasttools_turbo_chunks", 8)),
                    profile=1)
    # the APA annotation ("" = the umbrella annotation_gtf, already in reference_id)
    annotation = {"annotation": str(config.get("apa_annotation_gtf", "") or ""),
                  "transcripts": str(config.get("apa_transcripts", "basic")),
                  "utr_rules": "no CDS or incomplete 3' end (cds_end_NF, mRNA_end_NF) dropped"}
    if annotation["transcripts"] == "megasearch":
        annotation["apa_set"] = apa_set_settings(config)
    if module == "dapars2":
        return dict(annotation, coverage="bedtools genomecov -split -bg, UTR windows, primary mapped reads")
    return dict(annotation, library=config.get("qapa_library_type", "A"),
                polya_sites=str(config.get("qapa_polya_sites", "")),
                gencode_polya=str(config.get("qapa_gencode_polya", "")),
                polyasite=str(config.get("qapa_polyasite", "")),
                decoys=bool(config.get("qapa_decoys", False)))


APA_SET_RULE = ("protein-coding, CDS end in last exon, no *_NF tag, not readthrough, 3' end "
                "Tier-1 or within slop of PolyASite (>= 2 protocols) or GENCODE polyA_site")


def apa_set_settings(config):
    """What the MegaSearch APA set is made from (apa_transcripts: megasearch)."""
    return {"rule": APA_SET_RULE, "slop": int(config.get("apa_set_slop", 50)),
            "gencode_polya": str(config.get("qapa_gencode_polya", "")),
            "polyasite": str(config.get("qapa_polyasite", ""))}


def apa_set_id(reference_id, config):
    return _digest({"reference_id": reference_id, "annotation": str(config.get("apa_annotation_gtf", "") or ""),
                    "settings": apa_set_settings(config)})


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
    "rmats_conservative": ("run.rmats", "cache.json"),
    "psisigma": ("run.SJ.tab.gz", "ir_positions.tsv.gz", "qc.json", "cache.json"),
    "majiq": ("run.sj", "cache.json"),
    "dapars2": ("utr_coverage.bedgraph.gz", "depth.json", "cache.json"),
    "qapa": ("quant.sf.gz", "lib_format_counts.json", "cache.json"),
    "vasttools": ("to_combine.tar", "cache.json"),
    "vasttools_turbo": ("to_combine.tar", "cache.json"),
    "ciriquant": ("counts.tsv.gz", "cache.json"),
    "circrna_star": ("Chimeric.out.junction.gz", "SJ.out.tab.gz", "Log.final.out", "cache.json"),
    "circexplorer2": ("calls.tsv.gz", "circularRNA_known.txt.gz", "back_spliced_junction.bed.gz", "cache.json"),
    # DCC runs per comparison on circrna_star caches; it keeps nothing per run
    "dcc": (),
    "whippet_circ": ("quant.psi.gz", "quant.jnc.gz", "quant.gene.tpm.gz", "quant.isoform.tpm.gz", "cache.json"),
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
    """Compact status JSON; failures raise.

    ok, unsupported (design), native_only (unrecognised native columns), empty
    (the tool reported no features) or no_tests (features, but none tested)."""
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
