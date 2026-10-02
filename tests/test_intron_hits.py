"""The single bedtools call in intron_hits.py matches the old per-microexon intersect."""

import importlib.util
import pathlib
import shutil
import sys
import unittest

REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
PYBEDTOOLS = importlib.util.find_spec("pybedtools") is not None and shutil.which("bedtools")


@unittest.skipUnless(PYBEDTOOLS, "requires pybedtools and bedtools")
class ContainingIntronsTests(unittest.TestCase):
    def test_matches_one_intersect_per_microexon(self):
        sys.path.insert(0, str(REPOSITORY / "src"))
        from pybedtools import BedTool
        from intron_hits import containing_introns

        introns = BedTool("\n".join([
            "chr1 100 500 SJ 0 +", "chr1 100 300 SJ 0 +", "chr1 150 500 SJ 0 +",
            "chr1 100 500 SJ 0 -", "chr2 100 500 SJ 0 +",
        ]), from_string=True).sort()
        microexons = [
            ("chr1", 200, 215, "+", 15),   # inside three + introns
            ("chr1", 290, 310, "+", 20),   # crosses the end of one of them
            ("chr1", 200, 201, "+", 1),    # 1 nt: zero-length query interval
            ("chr1", 200, 215, "-", 15),   # other strand
            ("chr1", 600, 615, "+", 15),   # no containing intron
            ("chr3", 200, 215, "+", 15),   # chromosome without introns
        ]
        expected = {}
        for key in microexons:
            chrom, start, end, strand, _ = key
            query = BedTool(" ".join([chrom, str(start), str(end - 1), "ME", "0", strand]), from_string=True)
            hits = [str(i).strip("\n") for i in introns.intersect(query, wa=True, s=True, F=1, nonamecheck=True)]
            if hits:
                expected[key] = hits
        self.assertEqual(containing_introns(introns, microexons), expected)
        self.assertEqual(len(expected[microexons[0]]), 3)
        self.assertNotIn(microexons[4], expected)


if __name__ == "__main__":
    unittest.main()
