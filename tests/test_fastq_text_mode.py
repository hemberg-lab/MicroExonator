"""Regression tests for Python 3 FASTQ handling in discovery scripts."""

import gzip
import importlib.util
import pathlib
import subprocess
import sys
import tempfile
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
BIOPYTHON_AVAILABLE = importlib.util.find_spec("Bio") is not None


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "requires a pipeline Biopython environment")
class FastqTextModeTests(unittest.TestCase):
    def test_round1_row_microexon_reads_from_gzipped_fastq(self):
        with tempfile.TemporaryDirectory() as directory:
            workdir = pathlib.Path(directory)
            genome = workdir / "genome.fa"
            genome.write_text(">chr1\n" + "A" * 100 + "AGCGT" + "A" * 95 + "\n")
            alignments = workdir / "sample.sam.pre_processed"
            alignments.write_text(
                "\t".join(
                    (
                        "read1", "0", "chr1:100+120|tx1|100_100", "1", "2M1I1M",
                        "ACGT", "IIII", "0", "4", "C", "2", "1", "1", "C",
                    )
                ) + "\n"
            )
            reads = workdir / "reads.fastq.gz"
            with gzip.open(reads, "wt") as handle:
                handle.write("@read1\nACGT\n+\nIIII\n@read2\nTTTT\n+\nIIII\n")

            result = subprocess.run(
                [sys.executable, str(REPOSITORY / "src" / "row_ME2.py"),
                 str(genome), str(alignments), str(reads)],
                cwd=workdir,
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("chr1_+_102_103", result.stdout)
            output = workdir / "sample.sam.row_ME.fastq"
            self.assertEqual(output.read_text(), "@read1\nACGT\n+\nIIII\n")

    def test_round2_selected_reads_from_gzipped_fastq(self):
        with tempfile.TemporaryDirectory() as directory:
            workdir = pathlib.Path(directory)
            alignments = workdir / "sample.sam.pre_processed"
            alignments.write_text("read1\t0\ttag\t1\t4M\tACGT\tIIII\n")
            reads = workdir / "reads.fastq.gz"
            with gzip.open(reads, "wt") as handle:
                handle.write("@read1\nACGT\n+\nIIII\n@read2\nTTTT\n+\nIIII\n")

            result = subprocess.run(
                [sys.executable, str(REPOSITORY / "src" / "round2_ME_reads_fastq2.py"),
                 str(alignments), str(reads)],
                cwd=workdir,
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                pathlib.Path(str(alignments) + ".fastq").read_text(),
                "@read1\nACGT\n+\nIIII\n",
            )


if __name__ == "__main__":
    unittest.main()
