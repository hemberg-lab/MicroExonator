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

    def test_run_identity_ignores_labels_and_metadata(self):
        base = load_umbrella_manifest(self.write([row()])).run_sha256("run1")
        for change in (dict(group="ctrl"), dict(sample_id="other", biological_replicate_id="bio9"),
                       dict(phenotype="glia")):
            self.assertEqual(load_umbrella_manifest(self.write([row(**change)])).run_sha256("run1"),
                             base, change)
        for change in (dict(batch_id="batch2"), dict(strandedness="firststrand"),
                       dict(source_1="reads/other_R1.fastq.gz")):
            self.assertNotEqual(load_umbrella_manifest(self.write([row(**change)])).run_sha256("run1"),
                                base, change)
        shard = load_umbrella_manifest(self.write([row()])).shard_sha256("ref", "project", "case", "batch")
        self.assertEqual(load_umbrella_manifest(self.write([row(phenotype="glia")])).shard_sha256(
            "ref", "project", "case", "batch"), shard)

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

    def test_reference_id_auto_is_replaced_by_the_computed_id(self):
        manifest = load_umbrella_manifest(self.write([
            row(reference_id="auto"),
            row(run_id="run2", reference_id="auto", source_1="reads/X1.fastq.gz", source_2="reads/X2.fastq.gz")]))
        self.assertTrue(manifest.needs_reference_id())
        resolved = manifest.with_reference_id("0123456789abcdef")
        self.assertFalse(resolved.needs_reference_id())
        self.assertEqual({run.reference_id for run in resolved.included_runs()}, {"0123456789abcdef"})
        self.assertEqual(resolved.native_reads("run1")[0],
                         "umbrella/work/0123456789abcdef/project/batch/run1/R1.fastq.gz")
        self.assertEqual(resolved.included_runs()[0].metadata["phenotype"], "neural")
        # the ID is part of each run's identity, so a new reference reruns its outputs
        self.assertNotEqual(resolved.run_sha256("run1"),
                            manifest.with_reference_id("fedcba9876543210").run_sha256("run1"))

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


LIBRARY_COLUMNS = COLUMNS + ["library_id", "lane"]


class LibraryTests(unittest.TestCase):
    """Runs sharing library_id become one processing unit (technical-replicate audit cases)."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = pathlib.Path(self.temp.name) / "manifest.tsv"

    def load(self, rows):
        self.path.write_text("\t".join(LIBRARY_COLUMNS) + "\n" + "".join(
            "\t".join(str(record.get(column, "")) for column in LIBRARY_COLUMNS) + "\n" for record in rows))
        return load_umbrella_manifest(self.path)

    @staticmethod
    def sra(run, library="", **changes):
        return row(run_id=run, sample_id=changes.pop("sample_id", "GSM_" + run), source_type="sra",
                   source_1=run, source_2="", library_id=library, **changes)

    def test_lanes_of_one_experiment_become_one_unit(self):
        # PRJNA552045: four runs of SRX6385921
        manifest = self.load([self.sra(run, "SRX6385921") for run in
                              ("SRR9623442", "SRR9623439", "SRR9623441", "SRR9623440")]
                             + [self.sra("SRR9", sample_id="GSM2", biological_replicate_id="bio2")])
        self.assertEqual(sorted(manifest.by_run), ["SRR9", "SRX6385921"])
        unit = manifest.by_run["SRX6385921"]
        self.assertEqual([member[0] for member in unit.members],
                         ["SRR9623439", "SRR9623440", "SRR9623441", "SRR9623442"])
        self.assertEqual(unit.source_1, "SRR9623439;SRR9623440;SRR9623441;SRR9623442")
        self.assertEqual(unit.metadata["library_members"], unit.source_1)
        self.assertIs(manifest.by_member["SRR9623441"], unit)
        self.assertEqual(manifest.native_reads("SRX6385921"),
                         ("umbrella/work/ref/project/batch/SRX6385921/R1.fastq.gz",
                          "umbrella/work/ref/project/batch/SRX6385921/R2.fastq.gz"))
        self.assertEqual(manifest.by_run["SRR9"].members, ())

    def test_lanes_with_their_own_gsm_keep_every_sample_id(self):
        # MYT1L: the two lanes of one culture have different GSMs
        manifest = self.load([self.sra("SRR14130233", "iN_31_dcre_d7", sample_id="GSM5222893", lane="1"),
                              self.sra("SRR14130234", "iN_31_dcre_d7", sample_id="GSM5222894", lane="2")])
        unit = manifest.by_run["iN_31_dcre_d7"]
        self.assertEqual(unit.sample_id, "GSM5222893;GSM5222894")
        self.assertEqual(unit.metadata["lane"], "1;2")
        self.assertEqual(unit.metadata["phenotype"], "neural")

    def test_hashes_unchanged_without_library_id_and_follow_members(self):
        plain = self.load([self.sra("SRR1"), self.sra("SRR2", sample_id="GSM2", biological_replicate_id="b2")])
        with_empty = load_umbrella_manifest(self.path)
        self.assertEqual(plain.run_sha256("SRR1"), with_empty.run_sha256("SRR1"))
        two = self.load([self.sra("SRR1", "LIB"), self.sra("SRR2", "LIB")]).run_sha256("LIB")
        three = self.load([self.sra("SRR1", "LIB"), self.sra("SRR2", "LIB"), self.sra("SRR3", "LIB")])
        self.assertNotEqual(three.run_sha256("LIB"), two)

    def test_excluded_lanes_are_dropped_from_their_library(self):
        manifest = self.load([self.sra("SRR1", "LIB"),
                              self.sra("SRR2", "LIB", include="false", exclusion_reason="low quality")])
        self.assertEqual([member[0] for member in manifest.by_run["LIB"].members], ["SRR1"])
        self.assertTrue(manifest.by_run["LIB"].include)

    def test_inconsistent_libraries_are_refused(self):
        for rows, message in (
                ([self.sra("SRR1", "LIB"), self.sra("SRR2", "LIB", group="ctrl")], "disagree on group"),
                ([self.sra("SRR1", "LIB"), self.sra("SRR2", "LIB", biological_replicate_id="b2")],
                 "disagree on biological_replicate_id"),
                ([self.sra("SRR1", "LIB"), self.sra("SRR2", "LIB", layout="SE")], "disagree on layout"),
                ([self.sra("SRR1", "SRR3"), self.sra("SRR2", "SRR3"), self.sra("SRR3", sample_id="GSM3")],
                 "also another row's run_id"),
                ([self.sra("SRR1", "bad id")], "unsafe library_id")):
            with self.assertRaisesRegex(ValueError, message):
                self.load(rows)

    def test_selections_match_a_library_through_its_members(self):
        from src.comparison_preflight import preflight
        manifest = self.load(
            [self.sra("SRR1", "L1", sample_id="GSM1a", group="case", biological_replicate_id="c1"),
             self.sra("SRR2", "L1", sample_id="GSM1b", group="case", biological_replicate_id="c1"),
             self.sra("SRR3", sample_id="GSM3", group="case", biological_replicate_id="c3"),
             self.sra("SRR4", sample_id="GSM4", group="ctrl", biological_replicate_id="k4"),
             self.sra("SRR5", sample_id="GSM5", group="ctrl", biological_replicate_id="k5")])
        result = preflight(manifest, {"comparison_id": "c", "project_id": "project",
                                      "a": {"label": "picked", "samples": ["GSM1b", "GSM3"]},
                                      "b": {"label": "ctrl", "runs": ["SRR4", "SRR5"]}})
        self.assertEqual(result["runs"]["a"], ["L1", "SRR3"])
        self.assertEqual(result["library_members"], {"L1": ["SRR1", "SRR2"]})
        self.assertEqual(result["technical_run_replicates"], [])
        plain = preflight(manifest, {"comparison_id": "c2", "project_id": "project",
                                     "a": {"label": "one", "runs": ["SRR3", "SRR1"]},
                                     "b": {"label": "ctrl", "groups": ["ctrl"]}})
        self.assertEqual(plain["runs"]["a"], ["L1", "SRR3"])
        self.assertNotIn("library_members", preflight(manifest, {
            "comparison_id": "c4", "project_id": "project",
            "a": {"label": "k", "runs": ["SRR4", "SRR5"]}, "b": {"label": "c", "runs": ["SRR3"]}}))
