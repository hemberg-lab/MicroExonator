"""Tests for copying local FASTQ files with clean headers and bare '+' lines."""

import bz2
import gzip
import os
import pathlib
import subprocess
import tempfile
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY / "src" / "stage_local_fastq.sh"

CLEAN = b"@SRR1.1 1 length=4\nACGN\n+\nFFFF\n@SRR1.2 2 length=4\nTTTT\n+\nFFFF\n"


def stage(source, destination, env=None):
    subprocess.run(["bash", str(SCRIPT), str(source), str(destination)], check=True, env=env)
    with gzip.open(destination, "rb") as handle:
        return handle.read()


class StageLocalFastqTests(unittest.TestCase):
    def test_clean_gzip_is_copied_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            source = tmp / "in.fastq.gz"
            with gzip.open(source, "wb") as handle:
                handle.write(CLEAN)
            stage(source, tmp / "out.fastq.gz")
            self.assertEqual((tmp / "out.fastq.gz").read_bytes(), source.read_bytes())

    def test_long_clean_gzip_stops_checking_and_is_copied(self):
        # the check reads 10,000 records and closes the pipe; that must not fail the copy
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            source = tmp / "in.fastq.gz"
            with gzip.open(source, "wb") as handle:
                for number in range(50000):
                    handle.write(b"@r%d\nACGTACGTAC\n+\nFFFFFFFFFF\n" % number)
            stage(source, tmp / "out.fastq.gz")
            self.assertEqual((tmp / "out.fastq.gz").read_bytes(), source.read_bytes())

    def test_named_separator_becomes_bare_plus(self):
        # Whippet 1.6 fails on '+<name>' separators when reads contain N
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            source = tmp / "in.fastq.gz"
            with gzip.open(source, "wb") as handle:
                handle.write(CLEAN.replace(b"\n+\n", b"\n+SRR1 again\n"))
            self.assertEqual(stage(source, tmp / "out.fastq.gz"), CLEAN)

    def test_header_bytes_outside_printable_ascii_become_underscores(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            source = tmp / "in.fastq.gz"
            with gzip.open(source, "wb") as handle:
                handle.write("@read\tné 1\nACGT\n+\nFFFF\n".encode("utf-8"))
            # the two UTF-8 bytes of "é" each become "_"
            self.assertEqual(stage(source, tmp / "out.fastq.gz"), b"@read_n__ 1\nACGT\n+\nFFFF\n")

    def test_dirty_record_after_the_first_is_found(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            source = tmp / "in.fastq.gz"
            with gzip.open(source, "wb") as handle:
                handle.write(CLEAN + b"@SRR1.3\nACGT\n+SRR1.3\nFFFF\n")
            self.assertEqual(stage(source, tmp / "out.fastq.gz"), CLEAN + b"@SRR1.3\nACGT\n+\nFFFF\n")

    def test_plain_and_bz2_inputs_are_cleaned_and_compressed(self):
        dirty = CLEAN.replace(b"\n+\n", b"\n+x\n")
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            (tmp / "in.fastq").write_bytes(dirty)
            (tmp / "in.fastq.bz2").write_bytes(bz2.compress(dirty))
            self.assertEqual(stage(tmp / "in.fastq", tmp / "a.fastq.gz"), CLEAN)
            self.assertEqual(stage(tmp / "in.fastq.bz2", tmp / "b.fastq.gz"), CLEAN)

    def test_pigz_is_used_when_available(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            fake = tmp / "bin" / "pigz"
            fake.parent.mkdir()
            fake.write_text("#!/bin/sh\n[ \"$1\" = -p ] || exit 3\ntouch \"$(dirname \"$0\")/used\"\nexec gzip -c\n")
            fake.chmod(0o755)
            (tmp / "in.fastq").write_bytes(CLEAN)
            env = dict(os.environ, PATH=str(fake.parent) + os.pathsep + os.environ["PATH"])
            self.assertEqual(stage(tmp / "in.fastq", tmp / "out.fastq.gz", env=env), CLEAN)
            self.assertTrue((fake.parent / "used").exists())

    def test_unsupported_extension_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            (tmp / "in.bam").write_bytes(b"")
            result = subprocess.run(["bash", str(SCRIPT), str(tmp / "in.bam"), str(tmp / "out.fastq.gz")],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("only fastq", result.stderr)

    def test_local_sample_scripts_use_it(self):
        init = (REPOSITORY / "rules" / "init.smk").read_text()
        self.assertIn('"bash src/stage_local_fastq.sh " + row["path"]', init)
        self.assertNotIn('"cp " + row["path"]', init)


if __name__ == "__main__":
    unittest.main()
