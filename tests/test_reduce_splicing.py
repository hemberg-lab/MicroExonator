import gzip
import tempfile
import unittest
from pathlib import Path

from src.reduce_splicing import (me_shard_text, rebuild_me_table,
                                 rebuild_whippet_psi, whippet_shard_text)


ME_HEADER = "sample\tME\tME_coverages\texcluding_covs\tPSI\tCI_Lo\tCI_Hi\n"
PSI_HEADER = ("Gene\tNode\tCoord\tStrand\tType\tPsi\tCI_Width\tCI_Lo,Hi\tTotal_Reads\t"
              "Complexity\tEntropy\tInc_Paths\tExc_Paths\tEdges\n")


def gz(path, text):
    path.write_bytes(gzip.compress(text.encode()))
    return str(path)


class MicroExonatorShardTests(unittest.TestCase):
    def test_round_trip_keeps_values_verbatim_and_sparsity(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            a = ME_HEADER + ("r1\tchr1_+_10_20\t3.0\t41.0\t0.068\t0.01\t0.2\n"
                             "r1\tchr1_-_5_8\t0.0\t2.0\t0.0\t0.0\t0.8\n")
            b = ME_HEADER + "r2\tchr1_+_10_20\t10.5\t1.0\t0.913\t0.7\t0.99\n"
            tables = {"r2": gz(root / "b.gz", b), "r1": gz(root / "a.gz", a)}
            shard = root / "shard.tsv"
            shard.write_text(me_shard_text(tables))
            self.assertEqual(shard.read_text().splitlines()[0].split("\t")[:3],
                             ["ME", "r1.ME_coverages", "r1.excluding_covs"])
            self.assertEqual(rebuild_me_table(shard, "r2"), b)
            # rebuilt tables are sorted by microexon ID
            self.assertEqual(sorted(rebuild_me_table(shard, "r1").splitlines()[1:]),
                             sorted(a.splitlines()[1:]))


class WhippetShardTests(unittest.TestCase):
    def psi(self, root, name, psi, reads):
        return gz(root / name, PSI_HEADER +
                  "g1\t2\tchr1:1-10\t+\tCE\t{}\t0.1\t0.4,0.5\t{}\tK1\t0.5\t2-3:1\t2-4:0\t1-2:3\n".format(psi, reads) +
                  "g1\t3\tchr1:20-30\t+\tAA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")

    def test_round_trip_keeps_delta_columns_in_row_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            files = {"r1": self.psi(root, "1.gz", "0.45", "12"),
                     "r2": self.psi(root, "2.gz", "0.9", "30")}
            shard = root / "shard.tsv"
            shard.write_text(whippet_shard_text(files))
            rebuilt = rebuild_whippet_psi(shard, "r2").splitlines()
            self.assertEqual(rebuilt[0], PSI_HEADER.rstrip("\n"))
            self.assertEqual(rebuilt[1].split("\t")[:11],
                             ["g1", "2", "chr1:1-10", "+", "CE", "0.9", "0.1", "0.4,0.5", "30", "K1", "0.5"])
            self.assertEqual(rebuilt[1].split("\t")[11:], ["NA", "NA", "NA"])
            self.assertEqual(len(rebuilt), 3)

    def test_runs_from_different_indexes_are_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            other = gz(root / "o.gz", PSI_HEADER + "g9\t1\tchr2:1-5\t-\tCE\t1\t0\t1,1\t5\tK0\t0\tNA\tNA\tNA\n")
            with self.assertRaisesRegex(ValueError, "not the same index"):
                whippet_shard_text({"r1": self.psi(root, "1.gz", "0.4", "9"), "r2": other})


if __name__ == "__main__":
    unittest.main()
