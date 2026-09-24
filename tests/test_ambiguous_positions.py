"""Tests for listing microexons whose position reads cannot resolve."""

import pathlib
import sys
import tempfile
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

import ambiguous_positions  # noqa: E402


def centric_row(candidates, seq="CAG"):
    total = ",".join("{}|{}|0".format(me, score) for me, score in candidates)
    return "\t".join(["x", "tx", "0", "chr17:100-900", "0", str(len(seq)), seq, "1", "0", "0", "0", total]) + "\n"


class AmbiguousPositionTests(unittest.TestCase):
    def test_only_events_with_several_positions_are_listed_with_the_best_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "TOTAL.ME_centric.txt"
            path.write_text(
                centric_row([("chr17_-_200_203", "73.0"), ("chr17_-_600_603", "86.1"), ("chr17_-_400_403", "61.0")])
                + centric_row([("chr1_+_500_512", "80.0")], seq="ACGTACGTACGT")
                + "\t".join(["y", "tx", "0", "chr1:1+2", "0", "0", "", "0", "0", "0", "0", ""]) + "\n"
            )
            events = ambiguous_positions.ambiguous_events(path)
        self.assertEqual(events, {"chr17_-_600_603": ("CAG", ["chr17_-_200_203", "chr17_-_400_403"])})


if __name__ == "__main__":
    unittest.main()
