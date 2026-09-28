"""Contract tests for canonical umbrella sample metadata."""

import pathlib
from pathlib import Path
import tempfile
import unittest

from src.umbrella_manifest import load_umbrella_manifest


COLUMNS = "sample_id run_id biological_replicate_id project_id batch_id group source_type source_1 source_2 layout strandedness reference_id include exclusion_reason phenotype".split()


def row(**changes):
    record = dict(sample_id="sample", run_id="run1", biological_replicate_id="bio1",
                  project_id="project", batch_id="batch", group="case", source_type="fastq",
                  source_1="reads/R1.fastq.gz", source_2="reads/R2.fastq.gz", layout="PE",
                  strandedness="unstranded", reference_id="ref", include="true",
                  exclusion_reason="", phenotype="neural")
    record.update(changes)
    return record


class ManifestTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = pathlib.Path(self.temp.name) / "manifest.tsv"

    def write(self, rows, columns=COLUMNS):
        self.path.write_text("\t".join(columns) + "\n" + "".join(
            "\t".join(str(record.get(column, "")) for column in columns) + "\n"
            for record in rows))
        return self.path

    def test_paths_metadata_and_queries(self):
        manifest = load_umbrella_manifest(self.write([
            row(), row(run_id="run2", source_1="reads/X1.fastq.gz", source_2="reads/X2.fastq.gz"),
            row(run_id="excluded", include="false", exclusion_reason="low quality"),
        ]))
        self.assertEqual([run.run_id for run in manifest.included_runs()], ["run1", "run2"])
        self.assertEqual(manifest.groups(), ["case"])
        self.assertEqual(manifest.batches(), ["batch"])
        self.assertEqual(len(manifest.runs_for("project", "case", "batch")), 2)
        self.assertEqual(manifest.included_runs()[0].metadata["phenotype"], "neural")
        self.assertEqual(manifest.included_runs()[0].source_1,
                         str((self.path.parent / "reads/R1.fastq.gz").resolve()))
        self.assertEqual(manifest.native_reads("run1"), (
            "umbrella/work/ref/project/batch/run1/R1.fastq.gz",
            "umbrella/work/ref/project/batch/run1/R2.fastq.gz"))
        self.assertEqual(manifest.legacy_fastq("run1"),
                         "umbrella/work/ref/project/batch/run1/legacy.fastq.gz")

    def test_technical_run_can_arrive_in_another_batch(self):
        manifest = load_umbrella_manifest(self.write([
            row(), row(run_id="run2", batch_id="later", source_1="reads/X1.fastq.gz",
                       source_2="reads/X2.fastq.gz")]))
        self.assertEqual(manifest.batches(), ["batch", "later"])
        self.assertEqual(len(manifest.runs_for("project", "case", "later")), 1)

    def test_run_and_shard_identities_ignore_appended_batches(self):
        first = load_umbrella_manifest(self.write([row()]))
        run_hash = first.run_sha256("run1")
        shard_hash = first.shard_sha256("ref", "project", "case", "batch")
        extended = load_umbrella_manifest(self.write([
            row(), row(run_id="run2", sample_id="sample2", biological_replicate_id="bio2",
                       batch_id="later", source_1="reads/X1.fastq.gz", source_2="reads/X2.fastq.gz"),
        ]))
        self.assertEqual(extended.run_sha256("run1"), run_hash)
        self.assertEqual(extended.shard_sha256("ref", "project", "case", "batch"), shard_hash)
        self.assertNotEqual(extended.shard_sha256("ref", "project", "case", "later"), shard_hash)
        changed = load_umbrella_manifest(self.write([
            row(source_1="reads/changed.fastq.gz")]))
        self.assertNotEqual(changed.run_sha256("run1"), run_hash)
        self.assertNotEqual(changed.shard_sha256("ref", "project", "case", "batch"), shard_hash)

    def test_sra_accession_is_not_resolved_as_path(self):
        manifest = load_umbrella_manifest(self.write([
            row(source_type="sra", source_1="SRR12345", source_2="")]))
        self.assertEqual(manifest.included_runs()[0].source_1, "SRR12345")

    def test_required_columns_and_values(self):
        with self.assertRaisesRegex(ValueError, "reference_id"):
            load_umbrella_manifest(self.write([row()], [c for c in COLUMNS if c != "reference_id"]))
        for changes in (dict(layout="other"), dict(source_type="other"),
                        dict(strandedness="other"), dict(include="maybe"),
                        dict(source_1=""), dict(include="false", exclusion_reason=""),
                        dict(layout="SE", source_2="reads/R2.fastq.gz"),
                        dict(layout="PE", source_2="")):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                load_umbrella_manifest(self.write([row(**changes)]))

    def test_duplicate_run_and_conflicting_sample_metadata(self):
        with self.assertRaisesRegex(ValueError, "duplicate run_id"):
            load_umbrella_manifest(self.write([row(), row()]))
        for changes in (dict(group="control"), dict(biological_replicate_id="bio2"),
                        dict(project_id="other"), dict(reference_id="other")):
            with self.subTest(changes=changes), self.assertRaisesRegex(ValueError, "sample_id"):
                load_umbrella_manifest(self.write([row(), row(run_id="run2", **changes)]))


if __name__ == "__main__":
    unittest.main()


class LegacyFastqLinkTests(unittest.TestCase):
    def test_legacy_fastq_outlives_the_temporary_staged_file(self):
        import os
        from src.umbrella_manifest import link_legacy_fastq
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.fastq.gz"
            source.write_bytes(b"reads")
            # local .gz source: staged as a symlink, legacy path points at the source
            staged = root / "staged_link.fastq.gz"
            staged.symlink_to(source)
            link_legacy_fastq(staged, root / "FASTQ" / "a.fastq.gz")
            staged.unlink()
            self.assertEqual((root / "FASTQ" / "a.fastq.gz").read_bytes(), b"reads")
            # workflow-owned staged file (download): hard link keeps the data
            owned = root / "owned.fastq.gz"
            owned.write_bytes(b"downloaded")
            link_legacy_fastq(owned, root / "FASTQ" / "b.fastq.gz")
            owned.unlink()
            self.assertEqual((root / "FASTQ" / "b.fastq.gz").read_bytes(), b"downloaded")
            self.assertFalse(os.path.islink(root / "FASTQ" / "b.fastq.gz"))
