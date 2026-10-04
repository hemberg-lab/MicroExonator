"""Opt-in packet execution history and immutable retained-product guards.

No existing Snakemake invocation is enrolled implicitly. The ledger separates
packet identity, execution selection, attempts, and validated output completion.
Completion means validated files, not successful differential inference.
"""
import fcntl
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

DISCOVERY = {"microexonator", "whippet"}
COMPLEMENT = {"hisat2", "salmon", "rmats", "coverage", "suppa2", "majiq", "dapars2", "qapa", "multiqc"}


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    """Hash bytes and relative names, including every member of an index directory."""
    path = Path(path)
    if not path.exists():
        raise ValueError("missing retained/input path: {}".format(path))
    result = hashlib.sha256()
    if path.is_dir():
        files = sorted(p for p in path.rglob("*") if p.is_file())
        if not files:
            raise ValueError("empty retained/index directory: {}".format(path))
        for member in files:
            result.update(str(member.relative_to(path)).encode() + b"\0")
            result.update(digest(member).encode() + b"\0")
    else:
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                result.update(chunk)
    return result.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.{}".format(os.getpid()))
    with temporary.open("w") as stream:
        json.dump(value, stream, sort_keys=True, indent=2)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


class PacketLedger:
    def __init__(self, root, identity):
        self.root = Path(root).resolve()
        self.identity = identity
        self.path = self.root / "packet.json"
        self.lock = None

    def __enter__(self):
        self.root.mkdir(parents=True, exist_ok=True)
        self.lock = (self.root / ".lock").open("a+")
        try:
            fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            self.lock.close()
            self.lock = None
            raise ValueError("packet is locked by another execution")
        try:
            if self.path.exists():
                self.state = json.loads(self.path.read_text())
                if self.state.get("schema") != 1:
                    raise ValueError("unsupported packet ledger schema")
                if self.state["identity"] != self.identity:
                    changed = sorted(key for key in set(self.identity) | set(self.state["identity"])
                                     if self.identity.get(key) != self.state["identity"].get(key))
                    raise ValueError("incompatible packet identity ({}); use an explicit packet revision".format(
                        ", ".join(changed)))
                self.verify()
            else:
                self.state = {"schema": 1, "identity": self.identity, "created": timestamp(), "executions": {}}
                self.save()
            return self
        except BaseException:
            self.__exit__(None, None, None)
            raise

    def __exit__(self, *unused):
        if self.lock:
            fcntl.flock(self.lock, fcntl.LOCK_UN)
            self.lock.close()
            self.lock = None

    def save(self):
        atomic_json(self.path, self.state)

    def verify(self):
        for execution in self.state["executions"].values():
            for files in list(execution.get("outputs", {}).values()) + list(execution.get("artifacts", {}).values()):
                for path, expected in files.items():
                    if digest(path) != expected:
                        raise ValueError("changed retained output: {}".format(path))

    def completed_tools(self):
        return sorted({tool for execution in self.state["executions"].values()
                       if execution["status"] == "complete"
                       for tool in execution.get("completed_tools", [])})

    def begin(self, execution_id, mode, libraries, comparisons, tools):
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", execution_id):
            raise ValueError("unsafe execution ID")
        allowed = DISCOVERY if mode == "discovery" else COMPLEMENT if mode == "complement" else set()
        if not tools or set(tools) - allowed or not libraries:
            raise ValueError("invalid execution profile or empty library selection")
        self.verify()
        if mode == "discovery":
            discovered = {library for item in self.state["executions"].values()
                          if item["mode"] == "discovery" and item["status"] == "complete"
                          for library in item["libraries"]}
            if set(libraries) & discovered:
                raise ValueError("completed discovery producers cannot be repeated; use an explicit packet revision")
        if mode == "complement":
            discovered = {library for execution in self.state["executions"].values()
                          if execution["mode"] == "discovery" and execution["status"] == "complete"
                          for library in execution["libraries"]}
            if set(libraries) - discovered:
                raise ValueError("complement requires completed discovery for every selected library")
        selection = {"mode": mode, "libraries": sorted(set(libraries)),
                     "comparisons": sorted(set(comparisons)), "requested_tools": sorted(set(tools))}
        existing = self.state["executions"].get(execution_id)
        if existing:
            if any(existing[key] != value for key, value in selection.items()):
                raise ValueError("execution selection changed; use a new execution ID")
            if existing["status"] == "complete":
                raise ValueError("execution already complete; completed producers cannot be rerun")
            execution = existing
            # A process killed without an exception handler leaves a running attempt.
            if execution["attempts"][-1]["status"] == "running":
                execution["attempts"][-1].update(status="interrupted", ended=timestamp())
        else:
            execution = dict(selection, execution_id=execution_id, completed_tools=[], outputs={}, attempts=[])
            self.state["executions"][execution_id] = execution
        execution["status"] = "running"
        execution["attempts"].append({"started": timestamp(), "status": "running"})
        self.save()
        return execution

    def fail(self, execution, reason):
        execution["status"] = "failed"
        execution["attempts"][-1].update(status="failed", ended=timestamp(), reason=str(reason))
        self.save()

    def complete(self, execution, outputs, artifacts=None):
        if set(outputs) - set(execution["requested_tools"]):
            raise ValueError("completion includes a tool that was not requested")
        if set(outputs) != set(execution["requested_tools"]) or any(not paths for paths in outputs.values()):
            raise ValueError("incomplete retained tool inventory")
        self.verify()
        def inventory(groups):
            result = {}
            for tool, paths in groups.items():
                result[tool] = {}
                for path in paths:
                    path = Path(path).resolve()
                    if not path.is_file() or path.stat().st_size == 0:
                        raise ValueError("missing/empty retained product: {}".format(path))
                    try:
                        if path.suffix == ".gz":
                            with gzip.open(path, "rb") as stream:
                                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                                    pass
                        elif path.suffix == ".json":
                            json.loads(path.read_text())
                    except (OSError, EOFError, ValueError) as error:
                        raise ValueError("invalid retained product {}: {}".format(path, error))
                    result[tool][str(path)] = digest(path)
            return result
        inventories = inventory(outputs)
        artifact_inventory = inventory(artifacts or {})
        # Do not alter permissions on a directory, index, or external symlink target.
        for files in list(inventories.values()) + list(artifact_inventory.values()):
            for path in files:
                file = Path(path)
                file.chmod(file.stat().st_mode & ~0o222)
        execution.update(outputs=inventories, artifacts=artifact_inventory, completed_tools=sorted(outputs), status="complete",
                         completion_semantics="validated output inventory; inference status is tool-specific")
        execution["attempts"][-1].update(status="complete", ended=timestamp())
        self.save()
        self.write_cumulative()

    def write_cumulative(self):
        """Keep selection and execution identity attached to every retained product."""
        rows = []
        for execution_id, execution in sorted(self.state["executions"].items()):
            if execution["status"] != "complete":
                continue
            for tool, files in sorted(execution["outputs"].items()):
                for path, checksum in sorted(files.items()):
                    rows.append({"execution_id": execution_id, "mode": execution["mode"], "tool": tool,
                                 "libraries": execution["libraries"], "comparisons": execution["comparisons"],
                                 "path": path, "sha256": checksum})
        atomic_json(self.root / "cumulative_inventory.json", rows)
        reports = []
        for execution_id, execution in sorted(self.state["executions"].items()):
            if execution["status"] != "complete":
                continue
            for category, files in execution.get("artifacts", {}).items():
                if category == "reference":
                    continue
                for path, checksum in files.items():
                    reports.append({"execution_id": execution_id, "mode": execution["mode"],
                                    "libraries": execution["libraries"], "comparisons": execution["comparisons"],
                                    "category": category, "path": path, "sha256": checksum})
        atomic_json(self.root / "cumulative_reports.json", reports)


def execution_view(manifest_path, comparisons_path, libraries=None, comparisons=None):
    """Selection is a view, never a mutation or relabelling of the master packet."""
    from src.umbrella_manifest import UmbrellaManifest, load_umbrella_manifest
    from src.comparison_preflight import load_comparisons, preflight
    manifest = load_umbrella_manifest(manifest_path)
    if libraries is not None and not libraries:
        raise ValueError("empty library selection")
    definitions = load_comparisons(comparisons_path) if comparisons_path else []
    comparison_ids = set(comparisons if comparisons is not None else
                         [item["comparison_id"] for item in definitions])
    known_comparisons = {item["comparison_id"] for item in definitions}
    if comparison_ids - known_comparisons:
        raise ValueError("unknown comparison selection")
    definitions = [item for item in definitions if item["comparison_id"] in comparison_ids]
    if libraries is None and comparisons is not None:
        if not definitions:
            raise ValueError("empty comparison selection requires explicit libraries")
        libraries = set()
        for definition in definitions:
            participants = preflight(manifest, definition)["runs"]
            libraries.update(participants["a"] + participants["b"])
    records = manifest.included_runs()
    selected = set(libraries or [record.run_id for record in records])
    known = {record.run_id for record in records} | set(manifest.by_member)
    if selected - known:
        raise ValueError("unknown library selection: {}".format(sorted(selected - known)))
    units = [record for record in records if record.run_id in selected
             or any(member[0] in selected for member in record.members)]
    member_ids = {member[0] for record in units for member in record.members}
    member_ids.update(record.run_id for record in units)
    with open(manifest_path) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fieldnames = reader.fieldnames
        rows = [row for row in reader if row["run_id"] in member_ids]
    view = UmbrellaManifest(manifest_path, units)
    # Validate both sides against the selected view; never silently drop a side.
    for definition in definitions:
        preflight(view, definition)
    return {"libraries": sorted(record.run_id for record in units),
            "rows": rows, "fieldnames": fieldnames, "comparisons": definitions}


DISCOVERY_PRODUCERS = {
    "umbrella_legacy_fastq", "umbrella_legacy_bridge", "umbrella_microexonator_shard",
    "umbrella_microexonator_delta", "umbrella_microexonator_whippet_delta",
    "umbrella_whippet_quant", "umbrella_whippet_shard", "umbrella_whippet_delta",
    "correct_quant", "get_PSI_sparse_quants_se", "coverage_to_PSI_report",
    "ME_SJ_coverage", "ME_reads", "bowtie_to_genome", "read_length_report",
    "ambiguous_positions", "consecutive_runs", "detection_filter", "detection_filter_se",
    "Get_ME_from_annotation", "merge_ME_centric", "merge_tags", "Micro_Exon_Tags",
    "Splice_Junction_Library", "junction_genomic_copies", "generate_fasta_from_bed12",
}


def validate_dry_run(text, mode):
    rules = set(re.findall(r"(?:^|\n)(?:localrule|rule|checkpoint) ([A-Za-z0-9_]+):", text))
    if mode == "complement":
        forbidden = rules & DISCOVERY_PRODUCERS
        forbidden.update(rule for rule in rules if rule.startswith(("Round1_", "Round2_")))
        if forbidden:
            raise ValueError("complement schedules a discovery producer: {}".format(sorted(forbidden)))
    if "umbrella_leafcutter" in rules:
        raise ValueError("LeafCutter is not admitted in the two-pass execution profile")
    return sorted(rules)


def resolve_inputs(value, root):
    """Resolve existing input paths, leaving IDs, flags and scratch variables alone."""
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            if key.endswith("index_prefix") and isinstance(item, str):
                result[key] = str((root / item).resolve())
            else:
                result[key] = resolve_inputs(item, root)
        return result
    if isinstance(value, list):
        return [resolve_inputs(item, root) for item in value]
    if isinstance(value, str) and value and "$" not in value:
        candidate = Path(value) if Path(value).is_absolute() else root / value
        if candidate.exists():
            return str(candidate.resolve())
    return value


def packet_identity(config, source_root, snakefile):
    """Freeze configuration and file contents, not mtime-based reference caches."""
    ignored = {"umbrella_execution", "umbrella_optional", "umbrella_modules", "umbrella_modules_mode",
               "umbrella_allow_restage", "working_directory", "umbrella_stage_tmpdir"}
    producer = {key: value for key, value in config.items() if key not in ignored}
    identity = {"producer_config": hashlib.sha256(json.dumps(producer, sort_keys=True).encode()).hexdigest()}
    def visit(value, key):
        if isinstance(value, dict):
            for name, item in value.items():
                visit(item, key + "." + name)
        elif isinstance(value, list):
            for index, item in enumerate(value):
                visit(item, key + "." + str(index))
        elif isinstance(value, str) and Path(value).is_absolute() and Path(value).exists():
            # Executable and external software directories are identities too.
            identity[key] = digest(value)
            if key.endswith("whippet_index"):
                companion = Path(value + ".exons.tab.gz")
                if not companion.exists():
                    raise ValueError("configured Whippet index is missing its exon sidecar: {}".format(companion))
                identity[key + ".exons"] = digest(companion)
        elif isinstance(value, str) and key.endswith("index_prefix"):
            prefix = Path(value)
            members = sorted(prefix.parent.glob(prefix.name + ".*"))
            if not members:
                raise ValueError("configured index prefix has no files: {}".format(prefix))
            identity[key] = {str(path): digest(path) for path in members if path.is_file()}
    visit(producer, "config")
    identity["seeded_reference_assets"] = {
        str(path.relative_to(source_root)): digest(path) for path in reference_files(source_root)}
    identity["workflow"] = digest(source_root / snakefile)
    for folder in ("src", "rules", "envs"):
        if (source_root / folder).exists():
            files = sorted(path for path in (source_root / folder).rglob("*")
                           if path.is_file() and "__pycache__" not in path.parts and path.suffix != ".pyc")
            identity[folder] = hashlib.sha256(json.dumps(
                [(str(path.relative_to(source_root)), digest(path)) for path in files]).encode()).hexdigest()
    return identity


def tool_for(rule, path):
    for prefix, tool in (("umbrella_microexonator", "microexonator"), ("umbrella_whippet", "whippet"),
                         ("umbrella_salmon", "salmon"), ("umbrella_rmats", "rmats"),
                         ("umbrella_coverage", "coverage"), ("umbrella_suppa", "suppa2"),
                         ("umbrella_majiq", "majiq"), ("umbrella_dapars2", "dapars2"),
                         ("umbrella_qapa", "qapa"), ("umbrella_multiqc", "multiqc")):
        if rule.startswith(prefix):
            return tool
    if rule == "umbrella_deseq2":
        return "salmon" if "tximport" in path else "hisat2"
    if rule in ("umbrella_hisat2", "umbrella_junctions", "umbrella_junction_shard",
                "umbrella_featurecounts", "umbrella_featurecounts_shard", "umbrella_capture_shard"):
        return "hisat2"
    if rule in DISCOVERY_PRODUCERS or rule.startswith(("Round1_", "Round2_")):
        return "microexonator"
    return None


def summary_rows(text):
    lines = text.splitlines()
    start = next((i for i, line in enumerate(lines) if line.startswith("output_file\t")), None)
    if start is None:
        raise ValueError("Snakemake output summary is missing; cannot validate retained inventory")
    return list(csv.DictReader(lines[start:], delimiter="\t"))


def reference_files(source):
    """One inventory controls both reference identity and reuse eligibility."""
    for folder in ("data", "Round1", "Round2", "umbrella/reference", "umbrella/modules/reference", "apa_reference"):
        origin = source / folder
        if not origin.exists():
            continue
        for file in origin.rglob("*"):
            if not file.is_file() or file.name.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
                continue
            if folder == "Round1" and file.name != "ME_TAGs.fa":
                continue
            if folder == "Round2" and not (file.name.startswith("ME_canonical_SJ_tags") or file.name == "TOTAL.ME_centric.txt"):
                continue
            if file.name == "checksum_cache.json":
                continue
            yield file


def seed_references(source, work):
    """Link only inventoried assets; scheduled writes through them are forbidden."""
    for file in reference_files(source):
        destination = work / file.relative_to(source)
        if not destination.exists():
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.symlink_to(file.resolve())


def validate_arguments(arguments):
    import argparse
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument("--profile")
    parser.add_argument("-j", "--jobs", "--cores", type=int)
    parser.add_argument("--use-conda", action="store_true")
    parser.add_argument("--conda-frontend", choices=("conda", "mamba"))
    parser.add_argument("--conda-prefix")
    parser.add_argument("-k", "--keep-going", action="store_true")
    parser.add_argument("--rerun-incomplete", action="store_true")
    parser.add_argument("--resources", nargs="+")
    parser.add_argument("--latency-wait", type=int)
    parser.add_argument("--force-use-threads", action="store_true")
    try:
        unused, unknown = parser.parse_known_args(arguments)
    except SystemExit as error:
        raise ValueError("invalid guarded launch options") from error
    if unknown:
        raise ValueError("unsafe/force/target override options: {}".format(unknown))
    return unused


def validated_profile(profile, source_root):
    import yaml
    folder = Path(profile)
    if not folder.is_absolute():
        folder = source_root / folder
    if not (folder / "config.yaml").exists():
        from snakemake import get_profile_file
        found = get_profile_file(profile, "config.yaml")
        if not found:
            raise ValueError("cannot locate profile {}; provide an absolute profile directory".format(profile))
        folder = Path(found).resolve().parent
    values = yaml.safe_load((folder / "config.yaml").read_text()) or {}
    safe = {"jobs", "cores", "local-cores", "cluster", "cluster-status", "cluster-cancel", "cluster-config",
            "cluster-sync", "jobscript", "jobname", "use-conda", "conda-frontend", "conda-prefix",
            "keep-going", "rerun-incomplete", "restart-times", "latency-wait", "resources",
            "default-resources", "set-resources", "set-threads", "max-jobs-per-second",
            "max-status-checks-per-second", "printshellcmds", "scheduler", "scheduler-ilp-solver"}
    unknown = {key for key in values if key.replace("_", "-") not in safe}
    if unknown:
        raise ValueError("profile has unchecked execution/config override settings: {}".format(sorted(unknown)))
    return folder.resolve()


def run_execution(config_path, source_root, snakemake, arguments, snakefile="MicroExonator.smk", dry_run=False):
    """Guard a phase-local quant_umbrella launch. Existing direct launches are unchanged."""
    import yaml
    config_path, source_root = Path(config_path).resolve(), Path(source_root).resolve()
    config = resolve_inputs(yaml.safe_load(config_path.read_text()), source_root)
    settings = config.get("umbrella_execution", {})
    if settings.get("enabled") is not True:
        raise ValueError("guarded launcher requires umbrella_execution.enabled: true")
    packet_id = settings.get("packet_id", "")
    execution_id = settings.get("execution_id", "")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", packet_id):
        raise ValueError("unsafe/missing packet ID")
    mode = settings.get("mode")
    if mode not in ("discovery", "complement"):
        raise ValueError("execution mode must be discovery or complement")
    launch = validate_arguments(arguments)
    arguments = list(arguments)
    profile = validated_profile(launch.profile, source_root) if launch.profile else None
    if profile:
        if "--profile" in arguments:
            position = arguments.index("--profile")
            arguments[position + 1] = str(profile)
        else:
            arguments = ["--profile=" + str(profile) if item.startswith("--profile=") else item
                         for item in arguments]
    view = execution_view(config["umbrella_manifest"], config.get("umbrella_comparisons"),
                          settings.get("libraries"), settings.get("comparisons"))
    tools = sorted(settings.get("tools", DISCOVERY if mode == "discovery" else COMPLEMENT))
    if mode == "discovery" and set(tools) != DISCOVERY:
        raise ValueError("discovery profile requires MicroExonator and Whippet")
    if mode == "complement" and set(tools) - COMPLEMENT:
        raise ValueError("invalid complement profile")
    identity = packet_identity(config, source_root, snakefile)
    packet_root = source_root / "umbrella/executions" / packet_id
    with PacketLedger(packet_root, identity) as ledger:
        execution = ledger.begin(execution_id, mode, view["libraries"],
                                 [item["comparison_id"] for item in view["comparisons"]], tools)
        work = packet_root / execution_id / "work"
        work.mkdir(parents=True, exist_ok=True)
        (work / "logs").mkdir(exist_ok=True)
        if (source_root / "cluster.PBS.json").is_file():
            shutil.copy2(source_root / "cluster.PBS.json", work / "cluster.PBS.json")
            execution["attempts"][-1]["pbs_config_sha256"] = digest(work / "cluster.PBS.json")
        if profile:
            execution["attempts"][-1]["profile_sha256"] = digest(profile / "config.yaml")
        for folder in ("src", "rules", "envs"):
            if (source_root / folder).exists() and not (work / folder).exists():
                (work / folder).symlink_to(source_root / folder, target_is_directory=True)
        seed_references(source_root, work)
        for previous in ledger.state["executions"].values():
            if previous["status"] == "complete":
                seed_references(packet_root / previous["execution_id"] / "work", work)
        generated = dict(config)
        generated.pop("umbrella_execution", None)
        generated["working_directory"] = str(work) + "/"
        generated["umbrella_allow_restage"] = mode == "complement"
        # MegaSearch microexon inference is Whippet-independent. Ordinary
        # Whippet delta stays enabled in discovery; the hybrid route stays off.
        generated["delta_method"] = "microexonator"
        generated["umbrella_optional"] = {tool: tool in tools for tool in
                                           ("microexonator", "whippet", "salmon", "hisat2", "rmats", "coverage",
                                            "suppa2", "multiqc", "leafcutter", "microexonator_whippet_delta")}
        generated["umbrella_modules"] = sorted(set(tools) & {"majiq", "dapars2", "qapa"})
        generated["umbrella_modules_mode"] = "ingest"
        manifest_path = work / "manifest.tsv"
        with manifest_path.open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=view["fieldnames"], delimiter="\t")
            writer.writeheader()
            for row in view["rows"]:
                row = dict(row)
                for field in ("source_1", "source_2"):
                    if row[field] and row["source_type"] != "sra":
                        row[field] = str((source_root / row[field]).resolve())
                writer.writerow(row)
        comparison_path = work / "comparisons.yaml"
        comparison_path.write_text(yaml.safe_dump({"comparisons": view["comparisons"]}))
        generated["umbrella_manifest"] = str(manifest_path)
        if view["comparisons"]:
            generated["umbrella_comparisons"] = str(comparison_path)
        else:
            generated.pop("umbrella_comparisons", None)
        (work / "config.yaml").write_text(yaml.safe_dump(generated, sort_keys=True))
        atomic_json(packet_root / execution_id / "selection.json", view)
        base = list(snakemake) + ["-s", str(source_root / snakefile), "quant_umbrella"] + list(arguments)
        child_env = dict(os.environ, XDG_CACHE_HOME=str(work / ".cache"))
        # An inherited default profile must not bypass the explicit profile audit.
        child_env.pop("SNAKEMAKE_PROFILE", None)
        def invoke(extra, log_name):
            result = subprocess.run(base + extra, cwd=work, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, env=child_env)
            (packet_root / execution_id / log_name).write_text(result.stdout)
            if result.returncode:
                raise ValueError("Snakemake {} failed; see {}\n{}".format(
                    log_name, packet_root / execution_id / log_name, result.stdout[-2500:]))
            return result.stdout
        try:
            version_log = invoke(["--version"], "snakemake-version.log")
            match = re.search(r"(?m)^(7\.\d+(?:\.\d+)?[^\n]*)$", version_log)
            if not match:
                raise ValueError("guarded execution currently supports validated Snakemake 7 only")
            version = match.group(1)
            execution["versions"] = {"snakemake": version, "workflow_sha256": identity["workflow"],
                                     "reference_versions": config.get("umbrella_reference", {}).get("versions", {})}
            execution["attempts"][-1]["command"] = base
            ledger.save()
            invoke(["--dry-run"], "dry-run.log")
            planned = summary_rows(invoke(["--summary"], "plan-summary.log"))
            scheduled = sorted({row["rule"] for row in planned if row.get("plan") == "update pending"})
            validate_dry_run("\n".join("rule {}:".format(rule) for rule in scheduled), mode)
            # A job must never write through a reused reference symlink.
            for row in planned:
                if row.get("plan") == "update pending":
                    name = row["output_file"]
                    path = work / name
                    try:
                        path.resolve().relative_to(work)
                    except ValueError:
                        raise ValueError("scheduled write to reused reference {}; use an explicit reference revision".format(name))
            ledger.verify()
            if dry_run:
                ledger.fail(execution, "dry-run only; no sequencing jobs launched")
                return {"scheduled_rules": scheduled, "workspace": str(work), "delta_method": generated["delta_method"]}
            invoke([], "run.log")
            if packet_identity(config, source_root, snakefile) != identity:
                raise ValueError("packet identity changed during execution; results cannot be certified")
            rows = summary_rows(invoke(["--summary"], "summary.log"))
            invalid = [row["output_file"] for row in rows
                       if row.get("status") == "missing" or row.get("plan") == "update pending"]
            if invalid:
                raise ValueError("required output inventory is missing or still pending: {}".format(invalid))
            outputs = {tool: [] for tool in tools}
            artifacts = {"synthesis": [], "qc": [], "reference": []}
            for row in rows:
                name, rule = row.get("output_file", ""), row.get("rule", "")
                path = work / name
                tool = tool_for(rule, name)
                if not path.is_file() or path.is_symlink():
                    continue
                if name.startswith(("umbrella/reference/", "umbrella/modules/reference/", "data/", "Round1/", "Round2/ME_canonical_SJ_tags", "Round2/TOTAL.ME_centric")):
                    if path.stat().st_size:
                        artifacts["reference"].append(path)
                elif rule in ("umbrella_synthesis", "umbrella_module_inventory"):
                    artifacts["synthesis"].append(path)
                elif rule.startswith("umbrella_qc") or rule in ("umbrella_fastqc", "umbrella_multiqc"):
                    if path.stat().st_size:
                        artifacts["qc"].append(path)
                    if tool in outputs:
                        outputs[tool].append(path)
                elif tool in outputs:
                    outputs[tool].append(path)
            # Directory outputs (indices, MultiQC and rMATS preparation) need their
            # actual files in the retained inventory, not just a directory timestamp.
            for row in rows:
                name, rule = row.get("output_file", ""), row.get("rule", "")
                directory = work / name
                if not directory.is_dir() or directory.is_symlink():
                    continue
                tool = tool_for(rule, name)
                for path in directory.rglob("*"):
                    if not path.is_file() or path.is_symlink() or path.name == ".snakemake_timestamp" or not path.stat().st_size:
                        continue
                    if name.startswith(("umbrella/reference/", "umbrella/modules/reference/")):
                        artifacts["reference"].append(path)
                    elif tool in outputs:
                        outputs[tool].append(path)
                        if tool == "multiqc":
                            artifacts["qc"].append(path)
            ledger.complete(execution, outputs, artifacts)
            execution["scheduled_rules"] = scheduled
            execution["workspace"] = str(work)
            ledger.save()
            return execution
        except BaseException as error:
            if execution["status"] != "complete":
                ledger.fail(execution, error)
            raise


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--configfile", default="config.yaml")
    parser.add_argument("--source-root", default=".")
    parser.add_argument("--snakemake", default="snakemake")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    arguments = args.arguments[1:] if args.arguments[:1] == ["--"] else args.arguments
    result = run_execution(args.configfile, args.source_root, [args.snakemake], arguments, dry_run=args.dry_run)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    # Executing as src/umbrella_execution.py also needs the repository package root.
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    main()
