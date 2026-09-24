"""Tests for the Round2 genome blacklist (src/Filter1_round2.py)."""

import pathlib
import subprocess
import sys
import tempfile
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
TAG = "chr1:100+900|tx1|100_100"
TAG_SEQ = "ACGT" * 50


def genome_row(read, nm):
    return "\t".join([read, "0", "chr9", "500", "255", "40M", "*", "0", "0", "A" * 40, "I" * 40,
                      "XA:i:0", "MD:Z:40", "NM:i:{}".format(nm)])


class GenomeBlacklistTests(unittest.TestCase):
    def run_filter(self, use_tags):
        with tempfile.TemporaryDirectory() as directory:
            directory = pathlib.Path(directory)
            read_seq = TAG_SEQ[60:100]
            pre = directory / "S1.sam.pre_processed"
            pre.write_text("".join(
                "\t".join([read, "0", TAG, "60", "40M", seq, "I" * 40]) + "\n"
                for read, seq in [
                    ("better_on_tag", read_seq),
                    ("tie", read_seq),
                    ("unspliced_nowhere", read_seq),
                    ("worse_on_tag", "T" + read_seq[1:]),
                ]
            ))
            sam = directory / "S1.genome.sam"
            sam.write_text("@HD\tVN:1.0\n" + "\n".join([
                genome_row("better_on_tag", 1),
                genome_row("tie", 0),
                genome_row("worse_on_tag", 0),
            ]) + "\n")
            tags = directory / "tags.fa"
            tags.write_text(">{}\n{}\n".format(TAG, TAG_SEQ))
            command = [sys.executable, str(REPOSITORY / "src" / "Filter1_round2.py"), str(pre), str(sam)]
            if use_tags:
                command.append(str(tags))
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            return {line.split("\t")[0] for line in result.stdout.splitlines()}

    def test_read_matching_the_tag_better_than_the_genome_is_kept(self):
        self.assertEqual(self.run_filter(True), {"better_on_tag", "unspliced_nowhere"})

    def test_without_tags_every_genome_hit_is_removed(self):
        self.assertEqual(self.run_filter(False), {"unspliced_nowhere"})


if __name__ == "__main__":
    unittest.main()
