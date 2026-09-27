import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from src.reduce_salmon import reduce_group


class ReduceSalmonTests(unittest.TestCase):
    def test_script_entrypoint_loads_without_pythonpath_configuration(self):
        script = Path(__file__).resolve().parents[1] / "src/reduce_salmon.py"
        result = subprocess.run([sys.executable, str(script), "--help"],
                                text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        self.assertEqual(result.returncode, 0, result.stdout)

    def test_gene_counts_tpm_and_abundance_weighted_effective_length(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tx2gene = root / "tx2gene.tsv"
            tx2gene.write_text("TXNAME\tGENEID\nt1\tgA\nt2\tgA\nt3\tgB\n")
            quant = root / "quant.sf"
            quant.write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
                             "t1\t100\t80\t3\t9\n"
                             "t2\t200\t160\t1\t2\n"
                             "t3\t300\t250\t5\t4\n")
            prefix = root / "out"
            reduce_group({"run1": str(quant)}, str(tx2gene), str(prefix),
                         "ref", "manifesthash")
            expected = {
                "counts": "gene_id\trun1\ngA\t11\ngB\t4\n",
                "tpm": "gene_id\trun1\ngA\t4\ngB\t5\n",
                "length": "gene_id\trun1\ngA\t100\ngB\t250\n",
            }
            for kind, content in expected.items():
                with gzip.open(str(prefix) + ".salmon_" + kind + ".tsv.gz", "rt") as stream:
                    self.assertEqual(stream.read(), content)

    def test_unexpressed_gene_length_falls_back_to_mean_effective_length(self):
        # tximport uses gene lengths as offsets, so a zero length would give log(0).
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tx2gene = root / "tx2gene.tsv"
            tx2gene.write_text("TXNAME\tGENEID\nt1\tgA\nt2\tgA\n")
            quant = root / "quant.sf"
            quant.write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
                             "t1\t100\t80\t0\t0\n"
                             "t2\t200\t160\t0\t0\n")
            prefix = root / "out"
            reduce_group({"run1": str(quant)}, str(tx2gene), str(prefix),
                         "ref", "manifesthash")
            with gzip.open(str(prefix) + ".salmon_length.tsv.gz", "rt") as stream:
                self.assertEqual(stream.read(), "gene_id\trun1\ngA\t120\n")

    def test_transcript_level_tpm_and_counts_are_kept_for_suppa2(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tx2gene = root / "tx2gene.tsv"
            tx2gene.write_text("TXNAME\tGENEID\nt1\tgA\nt2\tgA\n")
            quants = {}
            for run, (tpm1, reads1) in (("run2", ("2.5", "7")), ("run1", ("3", "9"))):
                quant = root / (run + ".sf")
                quant.write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
                                 "t1\t100\t80\t{}\t{}\n"
                                 "t2\t200\t160\t1\t2\n".format(tpm1, reads1))
                quants[run] = str(quant)
            prefix = root / "out"
            reduce_group(quants, str(tx2gene), str(prefix), "ref", "manifesthash")
            expected = {
                "tx_tpm": "transcript_id\trun1\trun2\nt1\t3\t2.5\nt2\t1\t1\n",
                "tx_counts": "transcript_id\trun1\trun2\nt1\t9\t7\nt2\t2\t2\n",
            }
            for kind, content in expected.items():
                with gzip.open(str(prefix) + ".salmon_" + kind + ".tsv.gz", "rt") as stream:
                    self.assertEqual(stream.read(), content)

    def test_reduction_rejects_unmapped_transcripts_and_changed_quant(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tx2gene = root / "tx2gene.tsv"
            tx2gene.write_text("TXNAME\tGENEID\nt1\tgA\n")
            quant = root / "quant.sf"
            quant.write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
                             "t2\t100\t80\t3\t9\n")
            prefix = root / "out"
            with self.assertRaisesRegex(ValueError, "unmapped transcript"):
                reduce_group({"run1": str(quant)}, str(tx2gene), str(prefix),
                             "ref", "manifesthash")
            quant.write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
                             "t1\t100\t80\t3\t9\n")
            reduce_group({"run1": str(quant)}, str(tx2gene), str(prefix),
                         "ref", "manifesthash")
            quant.write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
                             "t1\t100\t80\t4\t10\n")
            with self.assertRaisesRegex(ValueError, "immutable"):
                reduce_group({"run1": str(quant)}, str(tx2gene), str(prefix),
                             "ref", "manifesthash")


if __name__ == "__main__":
    unittest.main()
