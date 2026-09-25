"""Tests for keeping one tag per distinct annotated flank path."""

import pathlib
import sys
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

import tag_paths  # noqa: E402


class MaximalFlankTests(unittest.TestCase):
    def test_identical_and_contained_upstream_flanks_are_merged(self):
        # Two transcripts share the last 6 nt before the junction and then diverge.
        flanks = ["AAAACCGGTT", "GGGGCCGGTT", "CCGGTT", "AAAACCGGTT"]
        self.assertEqual(tag_paths.maximal_flanks(flanks, "up"), ["AAAACCGGTT", "GGGGCCGGTT"])

    def test_downstream_flanks_are_compared_from_the_junction(self):
        flanks = ["TTGGCCAAAA", "TTGGCC", "TTGGCCGGGG", "ATGGCC"]
        self.assertEqual(tag_paths.maximal_flanks(flanks, "down"), ["TTGGCCAAAA", "TTGGCCGGGG", "ATGGCC"])

    def test_contained_tags_are_dropped_only_when_junctions_align(self):
        long_tag = (4, "AAAACCCC", 4, "tx1")
        short_same = (2, "AACC", 2, "tx2")      # the middle of long_tag, same junction
        short_shifted = (1, "AAAC", 3, "tx3")   # same letters, junction elsewhere
        other_path = (4, "GGGGCCCC", 4, "tx4")
        kept = tag_paths.maximal_tags([short_same, long_tag, short_shifted, other_path])
        self.assertEqual([t[3] for t in kept], ["tx1", "tx4", "tx3"])

    def test_first_path_keeps_its_transcript_name(self):
        self.assertEqual(tag_paths.path_label("ENSMUST1", 0), "ENSMUST1")
        self.assertEqual(tag_paths.path_label("ENSMUST1", 2), "ENSMUST1#p2")


if __name__ == "__main__":
    unittest.main()
