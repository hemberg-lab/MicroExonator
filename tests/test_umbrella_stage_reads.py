"""Small local tests for Conda-backed umbrella read staging."""

import bz2
import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))
from umbrella_stage_reads import stage_run  # noqa: E402


HEADER = ("sample_id", "run_id", "biological_replicate_id", "project_id",
          "batch_id", "group", "source_type", "source_1", "source_2",
          "layout", "strandedness", "reference_id", "include")


def manifest_for(directory, source_type, source_1, source_2="", layout="PE"):
    manifest = directory / "umbrella.tsv"
    row = ("sample", "run", "rep", "project", "batch", "group", source_type,
           source_1, source_2, layout, "unstranded", "ref", "true")
    manifest.write_text("\t".join(HEADER) + "\n" + "\t".join(row) + "\n")
    return manifest


class StageReadsTests(unittest.TestCase):
    def test_local_gzip_mates_remain_separate_symlinks(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for name in ("mate1.fastq.gz", "mate2.fastq.gz"):
                with gzip.open(directory / name, "wb") as stream:
                    stream.write(b"@read\nAC\n+\nII\n")
            manifest = manifest_for(directory, "fastq", "mate1.fastq.gz", "mate2.fastq.gz")
            r1, r2 = directory / "work/R1.fastq.gz", directory / "work/R2.fastq.gz"
            stage_run(manifest, "run", r1, r2)
            self.assertEqual(r1.resolve(), (directory / "mate1.fastq.gz").resolve())
            self.assertEqual(r2.resolve(), (directory / "mate2.fastq.gz").resolve())
            self.assertTrue(r1.is_symlink() and r2.is_symlink())

    def test_single_bzip2_is_converted_to_gzip(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            with bz2.open(directory / "reads.fastq.bz2", "wb") as stream:
                stream.write(b"@read\nAC\n+\nII\n")
            manifest = manifest_for(directory, "fastq", "reads.fastq.bz2", layout="SE")
            r1 = directory / "work/R1.fastq.gz"
            stage_run(manifest, "run", r1)
            self.assertFalse(r1.is_symlink())
            with gzip.open(r1, "rb") as stream:
                self.assertEqual(stream.read(), b"@read\nAC\n+\nII\n")

    def test_sra_invocation_produces_two_mates(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            manifest = manifest_for(directory, "sra", "SRR000001")
            r1, r2 = directory / "work/R1.fastq.gz", directory / "work/R2.fastq.gz"

            def fake_run(command, check):
                self.assertEqual(command[:3], ["fasterq-dump", "--split-files", "-O"])
                self.assertEqual(command[-1], "SRR000001")
                work = Path(command[3])
                (work / "SRR000001_1.fastq").write_bytes(b"mate1")
                (work / "SRR000001_2.fastq").write_bytes(b"mate2")

            with patch("umbrella_stage_reads.subprocess.run", side_effect=fake_run):
                stage_run(manifest, "run", r1, r2)
            with gzip.open(r1, "rb") as stream:
                self.assertEqual(stream.read(), b"mate1")
            with gzip.open(r2, "rb") as stream:
                self.assertEqual(stream.read(), b"mate2")

    def test_alignment_conversion_uses_samtools(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "reads.bam").touch()
            manifest = manifest_for(directory, "bam", "reads.bam")
            r1, r2 = directory / "work/R1.fastq.gz", directory / "work/R2.fastq.gz"
            calls = []

            def fake_run(command, check):
                calls.append(command)
                if command[1] == "sort":
                    Path(command[command.index("-o") + 1]).touch()
                else:
                    Path(command[command.index("-1") + 1]).write_bytes(b"mate1")
                    Path(command[command.index("-2") + 1]).write_bytes(b"mate2")

            with patch("umbrella_stage_reads.subprocess.run", side_effect=fake_run):
                stage_run(manifest, "run", r1, r2)
            self.assertEqual([call[:2] for call in calls],
                             [["samtools", "sort"], ["samtools", "fastq"]])
            self.assertEqual(calls[1][calls[1].index("-F") + 1], "0x900")
            self.assertTrue(r1.is_file() and r2.is_file())


if __name__ == "__main__":
    unittest.main()
