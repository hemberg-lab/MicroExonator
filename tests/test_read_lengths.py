"""Tests for read trimming before tag alignment and the read-length report."""

import pathlib
import subprocess
import sys
import tempfile
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

import read_length_report  # noqa: E402


class ReadTrimmingTests(unittest.TestCase):
    def test_reads_and_qualities_longer_than_max_read_len_are_trimmed(self):
        with tempfile.TemporaryDirectory() as directory:
            lengths = pathlib.Path(directory) / "S1.tsv"
            fastq = "@long\n" + "A" * 12 + "\n+\n" + "I" * 12 + "\n@short\nACGT\n+\nIIII\n"
            result = subprocess.run(
                ["awk", "-v", "L=10", "-v", "out={}".format(lengths), "-f", str(REPOSITORY / "src" / "trim_reads.awk")],
                input=fastq, capture_output=True, text=True, check=True,
            )
            self.assertEqual(result.stdout, "@long\n" + "A" * 10 + "\n+\n" + "I" * 10 + "\n@short\nACGT\n+\nIIII\n")
            self.assertEqual(sorted(lengths.read_text().splitlines()), sorted(["length\treads", "12\t1", "4\t1"]))


class ReadLengthReportTests(unittest.TestCase):
    def write(self, directory, name, hist):
        path = pathlib.Path(directory) / (name + ".tsv")
        path.write_text("length\treads\n" + "".join("{}\t{}\n".format(l, n) for l, n in hist.items()))
        return path

    def test_bulk_sample_within_max_read_len_has_no_note(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.write(directory, "S1", {100: 90, 98: 10})
            self.assertEqual(read_length_report.summarise(path, 150), (100, 100, 100, 0.0, ""))

    def test_long_and_short_libraries_are_noted(self):
        with tempfile.TemporaryDirectory() as directory:
            long_reads = self.write(directory, "S2", {250: 80, 150: 20})
            short_reads = self.write(directory, "S3", {28: 100})
            self.assertEqual(read_length_report.summarise(long_reads, 150)[3:], (80.0, "most reads longer than max_read_len and trimmed"))
            self.assertEqual(read_length_report.summarise(short_reads, 150)[4], "median read shorter than 40 nt")


if __name__ == "__main__":
    unittest.main()
