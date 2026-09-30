"""Existing shards stay immutable; the one-time params migration keeps results."""

import gzip
import json
import tempfile
import unittest
from pathlib import Path

from src.umbrella_identity import check_existing_shards, migrate_identity_params
from src.umbrella_manifest import load_umbrella_manifest


HEADER = ("sample_id\trun_id\tbiological_replicate_id\tproject_id\tbatch_id\tgroup\t"
          "source_type\tsource_1\tsource_2\tlayout\tstrandedness\treference_id\tinclude\t"
          "exclusion_reason\tsex\n")


def manifest(root, rows):
    path = root / "manifest.tsv"
    path.write_text(HEADER + "".join(
        "{0}\t{0}\t{0}\tproj\t{1}\t{2}\tsra\tSRR{3}\t\tSE\tunstranded\tref\t{4}\t{5}\t{6}\n".format(
            run, batch, group, number, include, "" if include == "true" else "qc", sex)
        for number, (run, batch, group, include, sex) in enumerate(rows, start=1)))
    return load_umbrella_manifest(path)


def shard(root, group, runs, batch="b1"):
    path = root / "junctions" / "ref" / "proj" / group / (batch + ".junctions.tsv.gz")
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as stream:
        stream.write("chrom\tstart\tend\tstrand\tlabel\tmax_anchor\tmulti_total\t"
                     "short_anchor_total\t" + "\t".join(runs) + "\n")


class GuardTests(unittest.TestCase):
    def test_unchanged_shards_and_metadata_edits_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            shard(root, "ctrl", ["a", "b"])
            rows = [("a", "b1", "ctrl", "true", "F"), ("b", "b1", "ctrl", "true", "M"),
                    ("c", "b2", "ctrl", "true", "F")]      # a new batch is fine
            self.assertEqual(check_existing_shards(manifest(root, rows), root), [])

    def test_moved_added_and_removed_runs_are_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            shard(root, "ctrl", ["a", "b", "x"])
            problems = check_existing_shards(manifest(root, [
                ("a", "b1", "ctrl", "true", "F"),
                ("b", "b1", "case", "true", "M"),      # moved
                ("x", "b1", "ctrl", "false", "M"),     # excluded
                ("n", "b1", "ctrl", "true", "F")]),    # added to a written shard
                root)
        text = "\n".join(problems)
        self.assertIn("run b is stored in shard ref/proj/ctrl/b1", text)
        self.assertIn("already exists without run(s) n", text)
        self.assertIn("holds run(s) x that the manifest no longer includes", text)
        self.assertNotIn("holds run(s) b", text)


class MigrationTests(unittest.TestCase):
    def test_drops_identity_params_once(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "metadata"
            (metadata / "sub").mkdir(parents=True)
            staged = {"rule": "umbrella_validate_reads", "params": ["abc"], "starttime": 1.0}
            other = {"rule": "umbrella_leafcutter", "params": ["x"]}
            (metadata / "sub" / "one").write_text(json.dumps(staged))
            (metadata / "two").write_text(json.dumps(other))
            marker = root / "umbrella" / ".run_identity_v2"
            self.assertEqual(migrate_identity_params(str(metadata), marker), 1)
            self.assertEqual(json.loads((metadata / "sub" / "one").read_text()),
                             {"rule": "umbrella_validate_reads", "starttime": 1.0})
            self.assertEqual(json.loads((metadata / "two").read_text()), other)
            self.assertTrue(marker.exists())
            (metadata / "two").write_text(json.dumps(dict(staged, rule="umbrella_qc_shard")))
            self.assertEqual(migrate_identity_params(str(metadata), marker), 0)


if __name__ == "__main__":
    unittest.main()
