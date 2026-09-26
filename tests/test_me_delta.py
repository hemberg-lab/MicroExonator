"""Tests for the Whippet-free differential inclusion model (src/me_delta.py)."""

import contextlib
import gzip
import io
import os
import pathlib
import sys
import tempfile
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

try:
    import numpy
    import me_delta  # noqa: E402
except ImportError:  # numpy lives in the pipeline env
    me_delta = None


@unittest.skipIf(me_delta is None, "needs numpy")
class MeDeltaModelTests(unittest.TestCase):
    def test_method_of_moments_recovers_beta_parameters(self):
        rng = numpy.random.default_rng(1)
        alpha, beta = me_delta.fit_beta(rng.beta(8.0, 2.0, 200000))
        self.assertAlmostEqual(alpha, 8.0, delta=0.15)
        self.assertAlmostEqual(beta, 2.0, delta=0.05)

    def test_identical_groups_have_no_difference(self):
        values = [(0.4, 200.0), (0.42, 180.0), (0.38, 220.0)]
        result = me_delta.compare(values, values, me_delta.event_rng(1, "chr1_+_10_20"))
        self.assertAlmostEqual(result["DeltaPsi"], 0.0, delta=0.02)
        self.assertLess(result["Probability"], 0.8)

    def test_clear_difference_is_called_with_high_probability(self):
        A = [(0.9, 300.0), (0.88, 250.0), (0.92, 280.0)]
        B = [(0.2, 300.0), (0.25, 250.0), (0.18, 280.0)]
        result = me_delta.compare(A, B, me_delta.event_rng(1, "chr1_+_10_20"))
        self.assertAlmostEqual(result["DeltaPsi"], 0.69, delta=0.03)
        self.assertEqual(result["Probability"], 1.0)
        self.assertEqual((result["Samples_A"], result["Samples_B"]), (3, 3))

    def test_samples_below_min_reads_are_dropped(self):
        result = me_delta.compare([(0.9, 300.0), (0.1, 4.0)], [(0.5, 100.0)], me_delta.event_rng(1, "x"))
        self.assertEqual(result["Samples_A"], 1)
        self.assertEqual(result["Reads_A"], 300.0)
        self.assertIsNone(me_delta.compare([(0.9, 4.0)], [(0.5, 100.0)], me_delta.event_rng(1, "x")))

    def test_results_do_not_depend_on_other_events(self):
        values = ([(0.6, 40.0), (0.7, 30.0)], [(0.4, 50.0)])
        first = me_delta.compare(*values, me_delta.event_rng(123456, "chr2_-_5_8"))
        second = me_delta.compare(*values, me_delta.event_rng(123456, "chr2_-_5_8"))
        self.assertEqual(first, second)


@unittest.skipIf(me_delta is None, "needs numpy")
class MeDeltaScriptTests(unittest.TestCase):
    def write_psi(self, path, rows):
        with gzip.open(path, "wt") as f:
            f.write("sample\tME\tME_coverages\texcluding_covs\tPSI\tCI_Lo\tCI_Hi\n")
            for ME, inclusion, exclusion in rows:
                total = inclusion + exclusion
                psi = "NA" if total == 0 else str(inclusion / total)
                f.write("s\t%s\t%s\t%s\t%s\t0\t1\n" % (ME, inclusion, exclusion, psi))

    def test_script_writes_whippet_style_microexon_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            files = []
            for name, inclusion in (("a1", 90), ("a2", 85), ("b1", 10), ("b2", 12)):
                path = os.path.join(tmp, name + ".corrected.PSI.gz")
                self.write_psi(path, [("chr1_+_1000_1015", inclusion, 100 - inclusion), ("chr1_-_5000_5003", 0, 0)])
                files.append(path)
            microexons = os.path.join(tmp, "out.robustly_detected.txt")
            with open(microexons, "w") as f:
                f.write("ME\tTranscript\nchr1_+_1000_1015\tENST1.2\nchr1_-_5000_5003\tENST2.1\n")
            gtf = os.path.join(tmp, "annotation.gtf")
            with open(gtf, "w") as f:
                f.write('chr1\tx\texon\t1\t2\t.\t+\t.\tgene_id "ENSG1.4"; transcript_id "ENST1.2";\n')

            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                me_delta.main(["-a", ",".join(files[:2]), "-b", ",".join(files[2:]),
                               "--microexons", microexons, "--gtf", gtf])

        lines = out.getvalue().splitlines()
        self.assertEqual(lines[0].split("\t"), me_delta.HEADER)
        self.assertEqual(len(lines), 2)   # the microexon without reads is not reported
        row = dict(zip(me_delta.HEADER, lines[1].split("\t")))
        self.assertEqual((row["exon_ID"], row["Gene"], row["Coord"], row["Strand"]), ("chr1_+_1000_1015", "ENSG1.4", "chr1:1001-1015", "+"))
        self.assertAlmostEqual(float(row["DeltaPsi"]), 0.76, delta=0.05)
        self.assertEqual(row["Probability"], "1.000")


if __name__ == "__main__":
    unittest.main()
