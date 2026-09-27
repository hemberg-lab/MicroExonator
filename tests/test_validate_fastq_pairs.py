"""Strict FASTQ structure and mate pairing tests."""

import gzip
import pathlib
import subprocess
import sys
import tempfile
import unittest


SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "src/validate_fastq_pairs.py"


def write_reads(path, names, complete=True):
    with gzip.open(path, "wt") as out:
        for name in names:
            out.write("@%s\nACGT\n+\nFFFF\n" % name)
        if not complete:
            out.write("@incomplete\nACGT\n+")


class FastqPairTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = pathlib.Path(self.temp.name)
        self.r1 = self.root / "R1.fastq.gz"
        self.r2 = self.root / "R2.fastq.gz"
        self.marker = self.root / "valid"

    def validate(self, paired=True):
        args = [sys.executable, str(SCRIPT), str(self.r1)]
        if paired:
            args.append(str(self.r2))
        return subprocess.run(args + ["--marker", str(self.marker), "--run-id", "run1"],
                              capture_output=True, text=True)

    def test_se_and_pe_complete_without_source_mutation(self):
        write_reads(self.r1, ["a/1", "b/1"])
        before = self.r1.read_bytes()
        self.assertEqual(self.validate(False).returncode, 0)
        self.assertEqual(self.r1.read_bytes(), before)
        self.marker.unlink()
        write_reads(self.r2, ["a/2", "b/2"])
        before2 = self.r2.read_bytes()
        self.assertEqual(self.validate().returncode, 0)
        self.assertEqual(self.r2.read_bytes(), before2)
        self.assertTrue(self.marker.exists())

    def test_illumina_mates_and_mismatch(self):
        write_reads(self.r1, ["a 1:N:0:ATCG", "b 1:N:0:ATCG"])
        write_reads(self.r2, ["a 2:N:0:ATCG", "b 2:N:0:ATCG"])
        self.assertEqual(self.validate().returncode, 0)
        self.marker.unlink()
        write_reads(self.r2, ["b 2:N:0:ATCG", "a 2:N:0:ATCG"])
        result = self.validate()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("run1", result.stderr)
        self.assertIn("record 1", result.stderr)
        self.assertFalse(self.marker.exists())

    def test_unequal_counts_and_incomplete_record(self):
        write_reads(self.r1, ["a", "b"])
        write_reads(self.r2, ["a"])
        self.assertIn("record 2", self.validate().stderr)
        write_reads(self.r1, ["a"], complete=False)
        self.assertIn("record 2", self.validate(False).stderr)
        self.assertFalse(self.marker.exists())

    def test_swapped_mate_designations_are_rejected(self):
        write_reads(self.r1, ["a/2"])
        write_reads(self.r2, ["a/1"])
        result = self.validate()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("mate designation mismatch", result.stderr)

    def test_empty_read_id_is_rejected_in_se_and_pe(self):
        write_reads(self.r1, [""])
        write_reads(self.r2, [""])
        for paired in (False, True):
            with self.subTest(paired=paired):
                result = self.validate(paired)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("run1 record 1", result.stderr)
                self.assertIn("empty read ID", result.stderr)


if __name__ == "__main__":
    unittest.main()
