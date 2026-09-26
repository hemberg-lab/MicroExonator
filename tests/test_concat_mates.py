"""Tests for concatenating paired-end mates with distinct read names."""

import gzip
import pathlib
import subprocess
import tempfile
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY / "src" / "concat_mates.sh"


def write_fastq(path, records):
    with gzip.open(path, "wt") as handle:
        for name, seq in records:
            handle.write("@%s\n%s\n+\n%s\n" % (name, seq, "F" * len(seq)))


def read_fastq(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    return [(lines[i], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)]


class ConcatMatesTests(unittest.TestCase):
    def run_script(self, mate1, mate2):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            write_fastq(tmp / "r_1.fastq.gz", mate1)
            write_fastq(tmp / "r_2.fastq.gz", mate2)
            subprocess.run(["bash", str(SCRIPT), str(tmp / "r_1.fastq.gz"), str(tmp / "r_2.fastq.gz"),
                            str(tmp / "r.fastq.gz")], check=True)
            return read_fastq(tmp / "r.fastq.gz")

    def test_mates_get_distinct_names_and_keep_sequences(self):
        records = self.run_script(
            [("SRR1.1 1 length=6", "ACGTAC"), ("SRR1.2 2 length=6", "GGGCCC")],
            [("SRR1.1 1 length=6", "TTTAAA"), ("SRR1.2 2 length=6", "CCCGGG")],
        )
        self.assertEqual([r[0].split()[0] for r in records],
                         ["@SRR1.1_1", "@SRR1.2_1", "@SRR1.1_2", "@SRR1.2_2"])
        self.assertEqual([r[1] for r in records], ["ACGTAC", "GGGCCC", "TTTAAA", "CCCGGG"])
        self.assertEqual([r[2] for r in records], ["FFFFFF"] * 4)
        self.assertEqual(len({r[0].split()[0] for r in records}), 4)

    def test_suffix_survives_names_that_already_end_in_slash_mate(self):
        # ENA-style names can already carry /1 and /2; the added suffix keeps them distinct
        # even for aligners that strip a trailing /1 or /2.
        records = self.run_script([("SRR1.1/1", "ACGT")], [("SRR1.1/2", "TTTT")])
        self.assertEqual([r[0] for r in records], ["@SRR1.1/1_1", "@SRR1.1/2_2"])

    def test_sra_download_script_uses_it(self):
        init = (REPOSITORY / "rules" / "init.smk").read_text()
        self.assertIn("bash src/concat_mates.sh FASTQ/${srr}_1.fastq.gz FASTQ/${srr}_2.fastq.gz", init)
        self.assertNotIn("then cat FASTQ/${srr}_1.fastq.gz", init)


if __name__ == "__main__":
    unittest.main()
