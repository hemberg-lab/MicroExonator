"""Validated, immutable intake metadata for the opt-in umbrella workflow."""

import csv
import hashlib
import json
import re
from dataclasses import dataclass, replace
from pathlib import Path


REQUIRED = ("sample_id", "run_id", "biological_replicate_id", "project_id",
            "batch_id", "group", "source_type", "source_1", "source_2",
            "layout", "strandedness", "reference_id", "include")
ALLOWED = {
    "source_type": {"sra", "fastq", "cram", "bam"},
    "layout": {"SE", "PE"},
    "strandedness": {"unstranded", "firststrand", "secondstrand", "auto"},
    "include": {"true", "false"},
}
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
SRA_ID = re.compile(r"^(?:SRR|ERR|DRR)[0-9]+$")
# reference_id "auto": the workflow computes it from config umbrella_reference
AUTO_REFERENCE_ID = "auto"
# The fields that change what is computed for a run. Everything else in a
# manifest row (group, sample_id, biological_replicate_id, include, extra
# metadata columns) only says how runs are selected and combined later, so
# editing it must never restage, realign or requantify a run.
RUN_IDENTITY = ("run_id", "source_type", "source_1", "source_2", "layout",
                "strandedness", "reference_id", "project_id", "batch_id")


@dataclass(frozen=True)
class UmbrellaRun:
    sample_id: str
    run_id: str
    biological_replicate_id: str
    project_id: str
    batch_id: str
    group: str
    source_type: str
    source_1: str
    source_2: str
    layout: str
    strandedness: str
    reference_id: str
    include: bool
    metadata: dict

    @property
    def work_dir(self):
        return "umbrella/work/{}/{}/{}/{}".format(
            self.reference_id, self.project_id, self.batch_id, self.run_id)


class UmbrellaManifest:
    def __init__(self, path, runs):
        self.path = Path(path)
        self.runs = tuple(runs)
        self.by_run = {run.run_id: run for run in runs}

    def with_reference_id(self, reference_id):
        """The manifest with every reference_id "auto" replaced by `reference_id`."""
        return UmbrellaManifest(self.path, [
            replace(run, reference_id=reference_id) if run.reference_id == AUTO_REFERENCE_ID else run
            for run in self.runs])

    def needs_reference_id(self):
        return any(run.reference_id == AUTO_REFERENCE_ID for run in self.runs)

    def included_runs(self):
        return [run for run in self.runs if run.include]

    def groups(self):
        return sorted({run.group for run in self.included_runs()})

    def batches(self):
        return sorted({run.batch_id for run in self.included_runs()})

    def runs_for(self, project_id, group, batch_id):
        return [run for run in self.included_runs() if
                (run.project_id, run.group, run.batch_id) ==
                (project_id, group, batch_id)]

    @staticmethod
    def _digest(runs):
        """Hash the processing identity of the runs that own an output.

        Only RUN_IDENTITY counts: relabelling a run or adding metadata columns
        leaves every per-run and per-shard output untouched.
        """
        payload = [{field: getattr(run, field) for field in RUN_IDENTITY}
                   for run in sorted(runs, key=lambda run: run.run_id)]
        return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

    def run_sha256(self, run_id):
        return self._digest([self._run(run_id)])

    def shard_sha256(self, reference_id, project_id, group, batch_id):
        runs = self.runs_for(project_id, group, batch_id)
        if not runs or any(run.reference_id != reference_id for run in runs):
            raise ValueError("empty or mismatched umbrella shard")
        return self._digest(runs)

    def _run(self, run_id):
        run = self.by_run[run_id]
        if not run.include:
            raise ValueError("run_id {} is excluded".format(run_id))
        return run

    def native_reads(self, run_id):
        run = self._run(run_id)
        reads = (run.work_dir + "/R1.fastq.gz",)
        if run.layout == "PE":
            reads += (run.work_dir + "/R2.fastq.gz",)
        return reads

    def legacy_fastq(self, run_id):
        run = self._run(run_id)
        if run.layout == "SE":
            return self.native_reads(run_id)[0]
        return run.work_dir + "/legacy.fastq.gz"


def link_legacy_fastq(staged, destination):
    """Give MicroExonator its FASTQ without tying it to the temporary staged file.

    A staged symlink (a local .gz source) is followed, and the legacy path
    links to the original source, which the workflow never deletes. A
    workflow-owned staged file is hard-linked, so it survives the staged path
    being removed; across file systems it is copied.
    """
    import os
    import shutil
    staged, destination = Path(staged), Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if staged.is_symlink():
        destination.symlink_to(staged.resolve())
        return
    try:
        os.link(staged, destination)
    except OSError:
        shutil.copyfile(staged, destination)


def load_umbrella_manifest(path):
    path = Path(path).expanduser().resolve()
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        columns = reader.fieldnames or []
        missing = set(REQUIRED) - set(columns)
        if missing:
            raise ValueError("manifest missing required columns: {}".format(", ".join(sorted(missing))))
        if len(columns) != len(set(columns)):
            raise ValueError("manifest contains duplicate column names")
        runs = []
        seen_runs = set()
        samples = {}
        for line_number, raw in enumerate(reader, start=2):
            if None in raw or any(value is None for value in raw.values()):
                raise ValueError("manifest line {} has wrong number of fields".format(line_number))
            values = {key: value.strip() for key, value in raw.items()}
            for key in REQUIRED:
                if key != "source_2" and not values[key]:
                    raise ValueError("manifest line {}: {} is empty".format(line_number, key))
            for key, allowed in ALLOWED.items():
                if values[key] not in allowed:
                    raise ValueError("manifest line {}: invalid {}: {}".format(
                        line_number, key, values[key]))
            for key in ("sample_id", "run_id", "biological_replicate_id",
                        "project_id", "batch_id", "group", "reference_id"):
                if not SAFE_ID.fullmatch(values[key]):
                    raise ValueError("manifest line {}: unsafe {}: {}".format(
                        line_number, key, values[key]))
            if values["run_id"] in seen_runs:
                raise ValueError("duplicate run_id: {}".format(values["run_id"]))
            seen_runs.add(values["run_id"])
            sample_key = values["sample_id"]
            identity = tuple(values[key] for key in (
                "biological_replicate_id", "project_id", "group",
                "layout", "strandedness", "reference_id"))
            if sample_key in samples and samples[sample_key] != identity:
                raise ValueError("sample_id {} has conflicting technical-run metadata".format(sample_key))
            samples[sample_key] = identity
            if values["include"] == "false" and not (values.get("exclusion_reason") or values.get("reason")):
                raise ValueError("manifest line {}: excluded run requires exclusion_reason or reason".format(line_number))
            if values["layout"] == "SE" and values["source_2"]:
                raise ValueError("manifest line {}: SE cannot have source_2".format(line_number))
            if values["source_type"] == "fastq":
                if values["layout"] == "PE" and not values["source_2"]:
                    raise ValueError("manifest line {}: PE FASTQ requires source_2".format(line_number))
            elif values["source_2"]:
                raise ValueError("manifest line {}: source_2 is only for local PE FASTQ".format(line_number))
            if values["source_type"] == "sra":
                if not SRA_ID.fullmatch(values["source_1"]):
                    raise ValueError("manifest line {}: invalid SRA run accession".format(line_number))
            else:
                original_values = dict(values)
                for key in ("source_1", "source_2"):
                    if values[key]:
                        values[key] = str((path.parent / values[key]).resolve())
            if values["source_type"] == "sra":
                original_values = dict(values)
            runs.append(UmbrellaRun(*(values[key] for key in REQUIRED[:-1]),
                                    values["include"] == "true", original_values))
    return UmbrellaManifest(path, runs)
