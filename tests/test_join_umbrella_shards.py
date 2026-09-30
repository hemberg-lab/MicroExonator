import gzip
import tempfile
import unittest
from pathlib import Path

from src.join_umbrella_shards import collapse_technical_runs, join
from src.shard_guard import write_immutable_bundle


def write_shard(root, name, text, reference="ref", guard_name=None):
    path = root / name
    guard = root / (guard_name or name.replace(".tsv.gz", "_checksums.json"))
    write_immutable_bundle({str(path): gzip.compress(text.encode(), mtime=0)}, str(guard),
                           {"reference_id": reference})
    return str(path)


JUNCTION_HEADER = "chrom\tstart\tend\tstrand\tlabel\tmax_anchor\tmulti_total\tshort_anchor_total"


class JoinTests(unittest.TestCase):
    def test_junction_shards_join_across_batches(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            one = write_shard(root, "b1.junctions.tsv.gz",
                              JUNCTION_HEADER + "\tr1\n"
                              "chr1\t201\t300\t+\tindex_intron\t40\t1\t0\t6\n")
            two = write_shard(root, "b2.junctions.tsv.gz",
                              JUNCTION_HEADER + "\tr2\n"
                              "chr1\t201\t300\t+\tindex_intron\t25\t2\t1\t3\n"
                              "chr1\t50\t90\t+\tnovel\t9\t0\t0\t2\n")
            text, digests = join([two, one], "junctions", "ref")
            self.assertEqual(text.splitlines(), [
                JUNCTION_HEADER + "\tr1\tr2",
                "chr1\t50\t90\t+\tnovel\t9\t0\t0\t0\t2",
                "chr1\t201\t300\t+\tindex_intron\t40\t3\t1\t6\t3",
            ])
            self.assertEqual(len(digests), 2)
            # joining is deterministic whatever the input order
            self.assertEqual(join([one, two], "junctions", "ref")[0], text)

    def test_refuses_edited_shards_other_references_and_repeated_runs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            one = write_shard(root, "b1.featurecounts.tsv.gz", "gene_id\tlength\tr1\ngA\t10\t4\n")
            other = write_shard(root, "b2.featurecounts.tsv.gz", "gene_id\tlength\tr2\ngA\t10\t4\n",
                                reference="other")
            with self.assertRaisesRegex(ValueError, "belongs to reference"):
                join([one, other], "featurecounts", "ref")
            again = write_shard(root, "b3.featurecounts.tsv.gz", "gene_id\tlength\tr1\ngA\t10\t1\n")
            with self.assertRaisesRegex(ValueError, "more than one shard"):
                join([one, again], "featurecounts", "ref")
            Path(one).chmod(0o644)
            Path(one).write_bytes(gzip.compress(b"gene_id\tlength\tr1\ngA\t10\t5\n", mtime=0))
            with self.assertRaisesRegex(ValueError, "checksum guard"):
                join([one], "featurecounts", "ref")

    def test_select_keeps_only_the_chosen_runs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            counts = write_shard(root, "b1.featurecounts.tsv.gz",
                                 "gene_id\tlength\tr1\tr2\tr3\ngA\t10\t4\t5\t6\n")
            text, _ = join([counts], "featurecounts", "ref", runs=["r1", "r3"])
            self.assertEqual(text.splitlines(), ["gene_id\tlength\tr1\tr3", "gA\t10\t4\t6"])
            fields = ("ME_coverages", "excluding_covs", "PSI", "CI_Lo", "CI_Hi")
            header = "ME\t" + "\t".join("{}.{}".format(run, f) for run in ("r1.x", "r2") for f in fields)
            me = write_shard(root, "b1.microexonator.tsv.gz", header + "\n"
                             "me1\t1\t2\t0.3\t0.1\t0.5\t\t\t\t\t\n"
                             "me2\t\t\t\t\t\t3\t4\t0.4\t0.2\t0.6\n")
            text, _ = join([me], "microexonator", "ref", runs=["r1.x"])
            # run IDs may contain dots; rows empty in every kept run are dropped
            self.assertEqual(text.splitlines(), [
                "ME\t" + "\t".join("r1.x." + f for f in fields), "me1\t1\t2\t0.3\t0.1\t0.5"])
            with self.assertRaisesRegex(ValueError, "none of the selected runs"):
                join([counts], "featurecounts", "ref", runs=["r9"])

    def test_salmon_matrices_share_one_guard(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "b1.salmon_counts.tsv.gz"
            write_immutable_bundle({str(path): gzip.compress(b"gene_id\tr1\ngA\t4\n", mtime=0)},
                                   str(root / "b1.salmon_checksums.json"), {"reference_id": "ref"})
            self.assertEqual(join([str(path)], "salmon", "ref")[0], "gene_id\tr1\ngA\t4\n")

    def test_technical_runs_are_summed_into_one_replicate(self):
        text = "gene_id\tlength\tr1\tr2\tr3\ngA\t10\t4\t1\t7\n"
        self.assertEqual(collapse_technical_runs(text, "featurecounts",
                                                 {"r1": "rep1", "r2": "rep1", "r3": "rep2"}),
                         "gene_id\tlength\trep1\trep2\ngA\t10\t5\t7\n")
        with self.assertRaisesRegex(ValueError, "only count"):
            collapse_technical_runs("ME\tr1.PSI\n", "microexonator", {})


if __name__ == "__main__":
    unittest.main()
