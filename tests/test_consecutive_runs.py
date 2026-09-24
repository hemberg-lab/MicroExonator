"""Tests for reporting runs of adjacent microexons."""

import pathlib
import sys
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

import consecutive_runs  # noqa: E402


class ConsecutiveRunTests(unittest.TestCase):
    def test_introns_touching_another_microexon_join_them_into_one_run(self):
        # Minus strand, transcript order m1 (high) -> m2 -> m3 (low); m4 is on its own.
        m1, m2, m3, m4 = "chr17_-_900_918", "chr17_-_600_630", "chr17_-_300_312", "chr17_-_5000_5009"
        introns = {
            m1: {"chr17:100-2000", "chr17:630-2000"},   # E_L..E_R and m2..E_L
            m2: {"chr17:312-900", "chr17:100-900"},     # m3..m1 and E_R..m1
            m3: {"chr17:100-2000", "chr17:100-600"},    # E_R..E_L and E_R..m2
            m4: {"chr17:4000-6000"},
        }
        runs = consecutive_runs.consecutive_runs(introns)
        self.assertEqual(set(runs), {m1, m2, m3})
        self.assertEqual(runs[m2][0], [m3, m2, m1])
        self.assertEqual(runs[m1][1], [m2])
        self.assertEqual(runs[m2][1], sorted([m1, m3]))


if __name__ == "__main__":
    unittest.main()
