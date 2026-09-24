"""Tests for choosing one tag alignment per read."""

import pathlib
import sys
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

import select_tag_alignments as sta  # noqa: E402

READ = "A" * 100


def sam(read, tag, pos, nm=0):
    return "\t".join([read, "0", tag, str(pos), "255", "100M", "*", "0", "0", READ, "I" * 100,
                      "XA:i:0", "MD:Z:100", "NM:i:{}".format(nm)]) + "\n"


# The read crosses the microexon junctions of its own inclusion tag (150_ME_150)
# but also fits inside the 150 nt flank of a neighbouring tag, where it crosses nothing.
INSIDE_FLANK = sam("r1", "chr7:100-200|tx|150_150", 21)          # covers 20..120, junction at 150
ON_JUNCTION = sam("r1", "chr7:300-900|tx|150_GAGCTGACGGGACTGAAG_150", 101)


class SelectTagAlignmentTests(unittest.TestCase):
    def test_junction_crossing_alignment_is_preferred_over_flank_only(self):
        for order in ([INSIDE_FLANK, ON_JUNCTION], [ON_JUNCTION, INSIDE_FLANK]):
            self.assertEqual(list(sta.select(order)), [ON_JUNCTION])

    def test_fewer_mismatches_win_among_crossing_alignments(self):
        worse = sam("r2", "chr1:10+500|tx|100_100", 41, nm=2)
        better = sam("r2", "chr1:10+600|tx|100_100", 41, nm=0)
        self.assertEqual(list(sta.select([worse, better])), [better])

    def test_ties_are_broken_the_same_way_every_time(self):
        a = sam("r3", "chr1:10+500|tx|100_100", 41)
        b = sam("r3", "chr1:10+600|tx|100_100", 41)
        first = list(sta.select([a, b]))
        self.assertEqual(len(first), 1)
        self.assertEqual(list(sta.select([a, b])), first)

    def test_reads_without_crossing_alignment_keep_their_first_and_headers_pass(self):
        lines = ["@HD\tVN:1.0\n", INSIDE_FLANK, sam("r4", "chr1:10+500|tx|100_100", 1)]
        self.assertEqual(list(sta.select(lines)), lines)

    def test_one_line_per_read(self):
        lines = [INSIDE_FLANK, ON_JUNCTION, sam("r2", "chr1:10+500|tx|100_100", 41)]
        self.assertEqual([l.split("\t")[0] for l in sta.select(lines)], ["r1", "r2"])


if __name__ == "__main__":
    unittest.main()
