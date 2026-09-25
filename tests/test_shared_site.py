"""Tests for microexons that share a splice site with a longer exon."""

import pathlib
import sys
import unittest
from collections import defaultdict


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

import shared_site  # noqa: E402


class SharedSiteTests(unittest.TestCase):
    def test_shared_acceptor_uses_down_side_and_shared_donor_the_up_side(self):
        self.assertEqual(shared_site.unique_side_half_reads(40, 12, shared_site.shared_sides(True, False)), 12)
        self.assertEqual(shared_site.unique_side_half_reads(9, 30, shared_site.shared_sides(False, True)), 9)
        self.assertEqual(shared_site.unique_side_half_reads(9, 30, shared_site.shared_sides(True, True)), 9)

    def introns(self):
        # plus strand: U (100-200) -> L (1000-1100) -> D (2000-2100);
        # the microexon is 1000-1015 and shares L's acceptor at 1000.
        estart = defaultdict(set)
        eend = defaultdict(set)
        for intron, low, high in (("chr1:200+1000", 200, 1000), ("chr1:1100+2000", 1100, 2000)):
            estart["chr1_%d" % high].add(intron)
            eend["chr1_%d" % low].add(intron)
        return estart, eend

    def test_far_side_is_the_one_the_microexon_does_not_share(self):
        estart, eend = self.introns()
        far, shared = shared_site.far_and_shared_introns("chr1:1000+1100", "+", True, estart, eend)
        self.assertEqual((far, shared), ({"chr1:1100+2000"}, {"chr1:200+1000"}))
        far, shared = shared_site.far_and_shared_introns("chr1:1000+1100", "+", False, estart, eend)
        self.assertEqual((far, shared), ({"chr1:200+1000"}, {"chr1:1100+2000"}))

    def test_minus_strand_acceptor_is_the_high_coordinate(self):
        estart, eend = self.introns()
        far, shared = shared_site.far_and_shared_introns("chr1:1000-1100", "-", True, estart, eend)
        self.assertEqual((far, shared), ({"chr1:200+1000"}, {"chr1:1100+2000"}))

    def test_each_intron_is_counted_once_and_untagged_far_side_falls_back(self):
        estart, eend = self.introns()
        coverage = defaultdict(int, {"chr1:1100+2000": 50, "chr1:200+1000": 70})
        credited = set()
        counts = shared_site.long_exon_units(
            ["chr1:1000+1100", "chr1:1000+1100"], "+", True, estart, eend,
            coverage, {"chr1:1100+2000", "chr1:200+1000"}, credited)
        self.assertEqual(counts, [50, 0])
        counts = shared_site.long_exon_units(
            ["chr1:1000+1100"], "+", True, estart, eend, coverage, {"chr1:200+1000"}, set())
        self.assertEqual(counts, [70])

    def test_no_shared_side_fallback_when_another_exon_has_a_far_junction(self):
        estart, eend = self.introns()
        eend["chr1_1300"].add("chr1:1300+1900")
        coverage = defaultdict(int, {"chr1:1100+2000": 50, "chr1:200+1000": 70, "chr1:1300+1900": 0})
        tagged = {"chr1:1100+2000", "chr1:200+1000"}   # the second exon's far intron has no tag
        counts = shared_site.long_exon_units(
            ["chr1:1000+1100", "chr1:1000+1300"], "+", True, estart, eend, coverage, tagged, set())
        self.assertEqual(counts, [50, 0])


if __name__ == "__main__":
    unittest.main()
