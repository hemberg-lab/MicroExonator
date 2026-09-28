"""Small local tests for Conda-backed umbrella read staging."""

import bz2
import gzip
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))
from umbrella_stage_reads import compress_fastq, stage_run  # noqa: E402


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

    def fake_sra(self, calls, files):
        """subprocess.run stand-in: prefetch makes <acc>/<acc>.sra, fasterq-dump writes `files`."""
        def fake_run(command, check):
            calls.append(command)
            if command[0] == "prefetch":
                folder = Path(command[command.index("-O") + 1]) / command[-1]
                folder.mkdir()
                (folder / (command[-1] + ".sra")).write_bytes(b"sra")
            else:
                work = Path(command[command.index("-O") + 1])
                for name, content in files.items():
                    (work / name).write_bytes(content)
        return fake_run

    def test_sra_run_is_prefetched_then_dumped_with_threads(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            manifest = manifest_for(directory, "sra", "SRR000001")
            r1, r2 = directory / "work/R1.fastq.gz", directory / "work/R2.fastq.gz"
            calls = []
            files = {"SRR000001_1.fastq": b"mate1", "SRR000001_2.fastq": b"mate2",
                     "SRR000001.fastq": b"unpaired"}
            with patch("umbrella_stage_reads.subprocess.run", side_effect=self.fake_sra(calls, files)):
                stage_run(manifest, "run", r1, r2, threads=6)
            self.assertEqual([call[0] for call in calls], ["prefetch", "fasterq-dump"])
            dump = calls[1]
            self.assertIn("--split-3", dump)
            self.assertEqual(dump[dump.index("--threads") + 1], "6")
            self.assertTrue(dump[-1].endswith("SRR000001/SRR000001.sra"))
            # fasterq-dump's own scratch (fasterq.tmp.*) goes in the run's work
            # directory, never the working directory of the whole workflow
            self.assertTrue(dump[dump.index("--temp") + 1].startswith(str(r1.parent) + "/stage_run_"))
            with gzip.open(r1, "rb") as stream:
                self.assertEqual(stream.read(), b"mate1")
            with gzip.open(r2, "rb") as stream:
                self.assertEqual(stream.read(), b"mate2")
            # only the two gzipped mates remain: no .sra, no uncompressed FASTQ, no scratch
            self.assertEqual(sorted(path.name for path in r1.parent.iterdir()),
                             ["R1.fastq.gz", "R2.fastq.gz"])

    def test_single_end_sra_and_scratch_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            manifest = manifest_for(directory, "sra", "SRR000002", layout="SE")
            r1, scratch = directory / "work/R1.fastq.gz", directory / "scratch"
            calls = []
            with patch("umbrella_stage_reads.subprocess.run",
                       side_effect=self.fake_sra(calls, {"SRR000002.fastq": b"single"})):
                stage_run(manifest, "run", r1, threads=2, tmpdir=scratch)
            self.assertTrue(calls[1][calls[1].index("-O") + 1].startswith(str(scratch)))
            with gzip.open(r1, "rb") as stream:
                self.assertEqual(stream.read(), b"single")
            self.assertEqual(list(scratch.iterdir()), [])

    def test_paired_manifest_row_for_a_single_end_accession_fails_clearly(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            manifest = manifest_for(directory, "sra", "SRR000003")
            with patch("umbrella_stage_reads.subprocess.run",
                       side_effect=self.fake_sra([], {"SRR000003.fastq": b"single"})):
                with self.assertRaisesRegex(FileNotFoundError, "paired-end in SRA"):
                    stage_run(manifest, "run", directory / "work/R1.fastq.gz", directory / "work/R2.fastq.gz")

    def test_pigz_is_used_when_available_and_gives_valid_gzip(self):
        pigz = shutil.which("pigz") or shutil.which("gzip")   # gzip -c stands in for pigz -c
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fake = directory / "pigz"
            fake.write_text("#!/bin/sh\nexec {} -c\n".format(pigz))
            fake.chmod(0o755)
            source = directory / "reads.fastq"
            source.write_bytes(b"@read\nAC\n+\nII\n")
            with patch("umbrella_stage_reads.shutil.which", return_value=str(fake)):
                compress_fastq(source, directory / "reads.fastq.gz", threads=4)
            with gzip.open(directory / "reads.fastq.gz", "rb") as stream:
                self.assertEqual(stream.read(), b"@read\nAC\n+\nII\n")
            self.assertFalse(source.exists())

    def test_paired_sra_reaches_microexonator_with_mate_suffixes(self):
        # the mate fix: MicroExonator's single FASTQ must name mates <id>_1 and <id>_2
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            manifest = manifest_for(directory, "sra", "SRR000004")
            r1, r2 = directory / "work/R1.fastq.gz", directory / "work/R2.fastq.gz"
            files = {"SRR000004_1.fastq": b"@SRR000004.1 1 length=2\nAC\n+\nII\n",
                     "SRR000004_2.fastq": b"@SRR000004.1 1 length=2\nGT\n+\nII\n"}
            with patch("umbrella_stage_reads.subprocess.run", side_effect=self.fake_sra([], files)):
                stage_run(manifest, "run", r1, r2, threads=2)
            legacy = directory / "legacy.fastq.gz"
            subprocess.run(["bash", str(REPOSITORY / "src" / "concat_mates.sh"),
                            str(r1), str(r2), str(legacy), "2"], check=True)
            with gzip.open(legacy, "rt") as stream:
                lines = stream.read().splitlines()
            self.assertEqual([lines[0].split()[0], lines[4].split()[0]],
                             ["@SRR000004.1_1", "@SRR000004.1_2"])
            self.assertEqual([lines[1], lines[5]], ["AC", "GT"])

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
                    self.assertIn("-@", command)
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
