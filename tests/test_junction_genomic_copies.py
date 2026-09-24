"""Tests for flagging junctions that have an unspliced copy in the genome."""

import pathlib
import sys
import tempfile
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

import junction_genomic_copies as jgc  # noqa: E402

ME = "chr1_+_1200_1209"
ME_SEQ = "AATCTATTG"
UP = "".join("ACGT"[i % 4] for i in range(100))
DOWN = "".join("TTGCA"[i % 5] for i in range(100))


class JunctionGenomicCopiesTests(unittest.TestCase):
    def write_inputs(self, directory):
        tags = directory / "tags.fa"
        tags.write_text(
            ">chr1:1000+2000|tx1|100_100\n{}\n".format(UP + DOWN)
            + ">chr1:1000+2000|tx1|100_{}_100\n{}\n".format(ME_SEQ, UP + ME_SEQ + DOWN)
            + ">chr2:5+900|tx2|100_100\n{}\n".format(UP + DOWN)
        )
        centric = directory / "TOTAL.ME_centric.txt"
        centric.write_text("\t".join(
            [ME, "tx1", "0", "chr1:1000+2000", "0", "9", ME_SEQ, "1", "80", "0", "0", ME + "|80|0"]
        ) + "\n")
        return tags, centric

    def test_cores_span_thirty_nt_either_side_of_each_junction(self):
        with tempfile.TemporaryDirectory() as directory:
            cores = jgc.junction_cores(*self.write_inputs(pathlib.Path(directory)))
        self.assertEqual(cores, {
            ME + "|skipping|chr1:1000+2000": UP[70:] + DOWN[:30],
            ME + "|inclusion|chr1:1000+2000": UP[70:] + ME_SEQ + DOWN[:30],
        })

    def test_report_lists_only_cores_with_genome_hits(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = pathlib.Path(directory)
            cores = directory / "cores.fa"
            cores.write_text(">{0}|skipping|chr1:1000+2000\nACGT\n>{0}|inclusion|chr1:1000+2000\nACGT\n".format(ME))
            hits = directory / "hits.txt"
            hits.write_text("{}|skipping|chr1:1000+2000\t-\tchr10\t127948282\tACGT\tIIII\t0\t\n".format(ME))
            lines = jgc.report(cores, hits)
        self.assertEqual(lines, [
            "ME\tjunction\tintron\tgenomic_copies",
            ME + "\tskipping\tchr1:1000+2000\tchr10:127948283:-",
        ])


if __name__ == "__main__":
    unittest.main()
